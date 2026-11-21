#!/usr/bin/env python3
"""Publish this preprocessing project to GitHub without git.

Creates/updates two repositories through the GitHub REST API:

    <owner>/plant-genome-preprocess        complete runnable project
    <owner>/plant-genome-preprocess-skill  skill package (SKILL.md + full pipeline)

Usage:
    setx GH_TOKEN ghp_xxx            (or: $env:GH_TOKEN="ghp_xxx")
    python tools/publish_github.py --token $env:GH_TOKEN --owner <login>
    python tools/publish_github.py --dry-run          # 只打印将上传的文件

The token needs `repo` scope (classic) or Contents+Administration read/write
(fine-grained).  Nothing is uploaded without a valid token.
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import time
import urllib.error
import urllib.request
from pathlib import Path

API = "https://api.github.com"

PROJECT_ROOT = Path(__file__).resolve().parents[2]      # wgdi_pipeline/
SKILL_ROOT = Path(os.environ.get("PGCI_SKILL_ROOT", r"D:\lab-task\07_REPORTS\skill_src\plant-genome-preprocess"))

SKIP_DIRS = {"__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache", ".git"}
SKIP_SUFFIX = {".pyc", ".pyo", ".zip"}


def collect(root: Path, prefix: str = "") -> list[tuple[str, Path]]:
    entries: list[tuple[str, Path]] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if any(part in SKIP_DIRS for part in path.parts):
            continue
        if path.suffix in SKIP_SUFFIX:
            continue
        rel = path.relative_to(root).as_posix()
        entries.append((f"{prefix}{rel}", path))
    return entries


def project_files() -> list[tuple[str, Path]]:
    files = collect(PROJECT_ROOT / "code", "code/")
    for name in ("QUICKSTART.md", "README.md", "requirements-wgdi-validated.txt"):
        p = PROJECT_ROOT / name
        if p.is_file():
            files.append((name, p))
    for p in PROJECT_ROOT.glob("*.sh"):
        files.append((p.name, p))
    for p in PROJECT_ROOT.glob("*.ps1"):
        files.append((p.name, p))
    files += collect(PROJECT_ROOT / "bin", "bin/")
    files += collect(PROJECT_ROOT / "vendor" / "wgdi", "vendor/wgdi/")
    site = PROJECT_ROOT / "vendor" / "sitecustomize.py"
    if site.is_file():
        files.append(("vendor/sitecustomize.py", site))
    return files


def skill_files() -> list[tuple[str, Path]]:
    files = [("SKILL.md", SKILL_ROOT / "SKILL.md")]
    files += collect(SKILL_ROOT / "pipeline", "pipeline/")
    return files


def request(token: str, method: str, url: str, payload: dict | None = None) -> dict:
    data = json.dumps(payload).encode() if payload is not None else None
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", f"Bearer {token}")
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("User-Agent", "pgci-publish")
    req.add_header("X-GitHub-Api-Version", "2022-11-28")
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            body = resp.read()
            return json.loads(body) if body else {}
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", "replace")[:400]
        raise SystemExit(f"{method} {url} -> {exc.code}: {detail}") from exc


def has_ref(token: str, owner: str, name: str, branch: str) -> bool:
    try:
        request(token, "GET", f"{API}/repos/{owner}/{name}/git/ref/heads/{branch}")
        return True
    except SystemExit:
        return False


def seed_empty_repo(token: str, owner: str, name: str, branch: str) -> None:
    """An empty repository rejects blob creation with 409; seed one commit."""
    content = base64.b64encode(b"# seed\n").decode()
    request(token, "PUT", f"{API}/repos/{owner}/{name}/contents/README.md", {
        "message": "chore: initialize repository", "content": content, "branch": branch,
    })
    for _ in range(10):
        if has_ref(token, owner, name, branch):
            return
        time.sleep(1)


def ensure_repo(token: str, owner: str, name: str, description: str, private: bool) -> str:
    try:
        info = request(token, "GET", f"{API}/repos/{owner}/{name}")
    except SystemExit:
        info = None
    if info:
        branch = info.get("default_branch") or "main"
        if not has_ref(token, owner, name, branch):
            seed_empty_repo(token, owner, name, branch)
        return branch
    created = request(token, "POST", f"{API}/user/repos", {
        "name": name, "description": description, "private": private, "auto_init": False,
    })
    time.sleep(2)
    branch = created.get("default_branch") or "main"
    seed_empty_repo(token, owner, name, branch)
    return branch


def push_tree(token: str, owner: str, repo: str, branch: str, files: list[tuple[str, Path]], message: str) -> str:
    entries = []
    total = 0
    for arc, path in files:
        content = base64.b64encode(path.read_bytes()).decode()
        blob = request(token, "POST", f"{API}/repos/{owner}/{repo}/git/blobs", {
            "content": content, "encoding": "base64",
        })
        entries.append({"path": arc, "mode": "100644", "type": "blob", "sha": blob["sha"]})
        total += path.stat().st_size
    tree = request(token, "POST", f"{API}/repos/{owner}/{repo}/git/trees", {"tree": entries})
    parents: list[str] = []
    try:
        ref = request(token, "GET", f"{API}/repos/{owner}/{repo}/git/ref/heads/{branch}")
        parents = [ref["object"]["sha"]]
    except SystemExit:
        parents = []
    commit = request(token, "POST", f"{API}/repos/{owner}/{repo}/git/commits", {
        "message": message, "tree": tree["sha"], "parents": parents,
    })
    if parents:
        request(token, "PATCH", f"{API}/repos/{owner}/{repo}/git/refs/heads/{branch}", {"sha": commit["sha"]})
    else:
        request(token, "POST", f"{API}/repos/{owner}/{repo}/git/refs", {
            "ref": f"refs/heads/{branch}", "sha": commit["sha"],
        })
    print(f"  pushed {len(files)} files, {total / 1024:.0f} KiB, commit {commit['sha'][:10]}")
    return commit["sha"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--token", default=os.environ.get("GH_TOKEN", ""))
    ap.add_argument("--owner", default="")
    ap.add_argument("--private", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    project = project_files()
    skill = skill_files()
    print(f"project files: {len(project)}; skill files: {len(skill)}")
    if args.dry_run:
        for arc, _ in project[:5] + skill[-5:]:
            print("  ", arc)
        return 0
    if not args.token:
        raise SystemExit("no GH_TOKEN; run: $env:GH_TOKEN='ghp_...'; python tools/publish_github.py")
    user = request(args.token, "GET", f"{API}/user")
    owner = args.owner or user["login"]
    print("publishing as", owner)
    branch1 = ensure_repo(args.token, owner, "plant-genome-preprocess",
                          "Plant genome annotation preprocessing for WGDI (complete, tested)", args.private)
    push_tree(args.token, owner, "plant-genome-preprocess", branch1, project,
              "feat: complete plant genome preprocessing pipeline")
    branch2 = ensure_repo(args.token, owner, "plant-genome-preprocess-skill",
                          "Skill: plant genome preprocessing (complete project inside)", args.private)
    push_tree(args.token, owner, "plant-genome-preprocess-skill", branch2, skill,
              "feat: plant-genome-preprocess skill with full pipeline")
    print("done: https://github.com/%s/plant-genome-preprocess" % owner)
    print("done: https://github.com/%s/plant-genome-preprocess-skill" % owner)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
