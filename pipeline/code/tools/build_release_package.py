#!/usr/bin/env python3
"""Build the self-contained pgci_preprocess + 0-6 WGDI workflow source release.

Source layout (this repository):

    wgdi_pipeline/                  pipeline entry scripts, bin/, vendor/
    wgdi_pipeline/code/             preprocessing source (src/, scripts/, tools/, tests/)

Release layout (what the other person receives):

    wgdi_pipeline/
        QUICKSTART.md, *.sh, *.ps1, requirements*.txt, bin/, vendor/
        code/
            src/pgci_preprocess/...
            scripts/ tools/ tests/ docs/
            resources/species_abbreviations.tsv
            pyproject.toml
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

CODE_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_ROOT = CODE_ROOT.parent
sys.path.insert(0, str(CODE_ROOT / "src"))
from pgci_preprocess.version import __version__ as VERSION  # noqa: E402

TOP_LEVEL_FILES = (
    "QUICKSTART.md",
    "run_all_pipeline.ps1",
    "bootstrap_wgdi.ps1",
    "bootstrap_wgdi.sh",
    "requirements-wgdi-validated.txt",
    "0check_data.sh",
    "1run_blast.sh",
    "1run_diamond.sh",
    "2run_pre_wgdi.sh",
    "3blast_to_wgdi.sh",
    "4run_wgdi_different.sh",
    "5run_wgdi_same.sh",
    "6run_wgdi.sh",
    "bin/wgdi",
    "bin/wgdi_entry.py",
    "vendor/sitecustomize.py",
)
CODE_FILES = (
    "README.md",
    "pyproject.toml",
    "resources/species_abbreviations.tsv",
    "resources/format_exceptions.tsv",
)
CODE_GLOBS = (
    "src/pgci_preprocess/*.py",
    "pgci_preprocess/*.py",
    "scripts/*.py",
    "tests/*.py",
    "tools/*.py",
    "tools/*.sh",
    "docs/**/*.md",
)
# WGDI 兼容 shim（4–6 步脚本依赖）；大体积 vendor/_Bio_py311_unused 不打包
PIPELINE_GLOBS = (
    "vendor/wgdi/*.py",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def collect_files() -> list[tuple[Path, str]]:
    """Return (source path, archive-relative path) pairs."""
    files: dict[str, Path] = {}
    for name in TOP_LEVEL_FILES:
        path = PIPELINE_ROOT / name
        if not path.is_file():
            raise FileNotFoundError(path)
        files[name] = path
    for name in CODE_FILES:
        path = CODE_ROOT / name
        if not path.is_file():
            raise FileNotFoundError(path)
        files[f"code/{name}"] = path
    for pattern in CODE_GLOBS:
        for path in CODE_ROOT.glob(pattern):
            if path.is_file():
                files[f"code/{path.relative_to(CODE_ROOT).as_posix()}"] = path
    for pattern in PIPELINE_GLOBS:
        for path in PIPELINE_ROOT.glob(pattern):
            if path.is_file():
                files[path.relative_to(PIPELINE_ROOT).as_posix()] = path
    return [(path, arc) for arc, path in sorted(files.items())]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"refusing to overwrite {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    entries = collect_files()
    manifest_lines = ["path\tsha256\tbytes"]
    for path, arc in entries:
        manifest_lines.append(f"{arc}\t{sha256(path)}\t{path.stat().st_size}")
    info = {
        "name": "pgci-preprocess-wgdi-workflow",
        "version": VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": (
            "Canonical preprocessing source, tests, 0-6 execution scripts, wrappers, "
            "compatibility shims and smoke tools; no raw data or previous result directories."
        ),
        "file_count": len(entries),
        "entry_layout": "wgdi_pipeline/code/<source tree>; launchers at wgdi_pipeline/<name>",
    }
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path, arc in entries:
            archive.write(path, "wgdi_pipeline/" + arc)
        archive.writestr("wgdi_pipeline/PACKAGE_SHA256.tsv", "\n".join(manifest_lines) + "\n")
        archive.writestr("wgdi_pipeline/PACKAGE_INFO.json", json.dumps(info, ensure_ascii=False, indent=2) + "\n")
    with zipfile.ZipFile(output) as archive:
        bad = archive.testzip()
        if bad:
            raise ValueError(f"ZIP CRC validation failed: {bad}")
    print(json.dumps({
        "output": str(output),
        "version": VERSION,
        "source_files": len(entries),
        "zip_entries": len(entries) + 2,
        "bytes": output.stat().st_size,
        "sha256": sha256(output),
        "crc_check": "PASS",
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
