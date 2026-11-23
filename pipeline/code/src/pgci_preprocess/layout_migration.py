"""Move an existing project from ``analysis_v2`` to ``wgdi_input``.

The container was renamed so the project root reads as what it is: one folder
of WGDI-ready species bundles.  Nothing inside a bundle changes — the four
WGDI files and the traceability files keep their bytes — only the directory
name, the index paths that point into it, and the now-redundant root-level
compatibility directories are touched.

The migration is deliberately explicit (``migrate-layout``) rather than
automatic: renaming a directory that holds hundreds of megabytes is a user
decision, and a half-migrated project is worse than an old one.
"""

from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass, field
from pathlib import Path

from pgci_preprocess.project_paths import (
    LEGACY_OUTPUT_DIRNAMES,
    OUTPUT_DIRNAME,
    index_file,
    machine_config_file,
)

INDEX_RELATIVE = Path("resources") / "wgdi_v2_index.tsv"
REPORT_NAME = "data-validation-report.tsv"
BUNDLE_SUFFIXES = (".pep", ".cds", ".new.gff", ".lens")


@dataclass(slots=True)
class LayoutMigrationResult:
    project_root: Path
    old_name: str = ""
    new_name: str = OUTPUT_DIRNAME
    status: str = "NOTHING"
    index_rows_changed: int = 0
    removed_links: list[str] = field(default_factory=list)
    kept_directories: list[str] = field(default_factory=list)
    problems: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.status not in {"BLOCKED", "FAIL"}

    def as_dict(self) -> dict[str, object]:
        return {
            "project_root": str(self.project_root),
            "old_name": self.old_name,
            "new_name": self.new_name,
            "status": self.status,
            "index_rows_changed": self.index_rows_changed,
            "removed_links": list(self.removed_links),
            "kept_directories": list(self.kept_directories),
            "problems": list(self.problems),
        }


def _read_index(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def _relative_bundle_paths(rows: list[dict[str, str]], old_name: str) -> list[str]:
    prefix = f"{old_name}/"
    found: list[str] = []
    for row in rows:
        for value in row.values():
            text = str(value or "")
            if text.startswith(prefix):
                found.append(text)
    return found


def _rewrite_index(
    path: Path, old_name: str, new_name: str, *, dry_run: bool
) -> int:
    fields, rows = _read_index(path)
    prefix = f"{old_name}/"
    changed = 0
    for row in rows:
        row_changed = False
        for key, value in row.items():
            text = str(value or "")
            if text.startswith(prefix):
                row[key] = f"{new_name}/{text[len(prefix):]}"
                row_changed = True
        changed += int(row_changed)
    if changed and not dry_run:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
    return changed


def _looks_like_compat_dir(path: Path, bundle_root: Path) -> bool:
    """A root directory that only mirrors the four WGDI files of one bundle."""

    if not path.is_dir():
        return False
    expected = {f"{path.name}{suffix}" for suffix in BUNDLE_SUFFIXES}
    entries = {item.name for item in path.iterdir()}
    if not entries or not entries.issubset(expected):
        return False
    if not expected.issubset(entries):
        return False
    source_dir = bundle_root / path.name
    for name in expected:
        source = source_dir / name
        candidate = path / name
        try:
            if not source.is_file() or not os.path.samefile(source, candidate):
                return False
        except OSError:
            return False
    return True


def _remove_compat_dirs(
    root: Path, bundle_root: Path, *, dry_run: bool
) -> tuple[list[str], list[str]]:
    removed: list[str] = []
    kept: list[str] = []
    for path in sorted(root.iterdir(), key=lambda item: item.name):
        if not path.is_dir() or path.name.startswith("."):
            continue
        if path.name in {OUTPUT_DIRNAME, *LEGACY_OUTPUT_DIRNAMES}:
            continue
        if (bundle_root / path.name).is_dir() and _looks_like_compat_dir(path, bundle_root):
            removed.append(path.name)
            if not dry_run:
                for item in path.iterdir():
                    item.unlink()
                path.rmdir()
        elif (bundle_root / path.name).is_dir():
            kept.append(path.name)
    return removed, kept


def _update_config(project_root: Path) -> None:
    config = machine_config_file(project_root)
    if not config.is_file():
        return
    try:
        payload = json.loads(config.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return
    if not isinstance(payload, dict):
        return
    paths = payload.setdefault("paths", {})
    if isinstance(paths, dict):
        paths["output_dirname"] = OUTPUT_DIRNAME
        config.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )


def _rewrite_report(project_root: Path, old_name: str, *, dry_run: bool) -> None:
    report = index_file(project_root, REPORT_NAME)
    if not report.is_file():
        return
    text = report.read_text(encoding="utf-8")
    updated = text.replace(f"{old_name}{os.sep}", f"{OUTPUT_DIRNAME}{os.sep}")
    updated = updated.replace(f"{old_name}/", f"{OUTPUT_DIRNAME}/")
    if updated != text and not dry_run:
        report.write_text(updated, encoding="utf-8")


def migrate_layout(
    project_root: Path | str, *, dry_run: bool = False
) -> LayoutMigrationResult:
    """Rename the bundle container and repair every path that points into it."""

    root = Path(project_root).expanduser().resolve()
    result = LayoutMigrationResult(project_root=root)
    old = root / LEGACY_OUTPUT_DIRNAMES[0]
    new = root / OUTPUT_DIRNAME
    index = root / INDEX_RELATIVE

    if not old.is_dir():
        if new.is_dir():
            result.status = "ALREADY"
            result.new_name = new.name
            if index.is_file():
                # A previous attempt may have renamed the directory but not the
                # index; finish that half without touching any data.
                _, rows = _read_index(index)
                stale = _relative_bundle_paths(rows, old.name)
                missing = [
                    item
                    for item in stale
                    if not (root / OUTPUT_DIRNAME / item[len(old.name) + 1 :]).is_file()
                ]
                if missing:
                    result.problems.append(
                        f"index still points at {old.name} but {len(missing)} file(s) "
                        "are missing under the new container"
                    )
                    result.status = "BLOCKED"
                else:
                    result.index_rows_changed = _rewrite_index(
                        index, old.name, OUTPUT_DIRNAME, dry_run=dry_run
                    )
            if not dry_run:
                _update_config(root)
            return result
        result.status = "NOTHING"
        return result

    if new.is_dir():
        result.status = "BLOCKED"
        result.problems.append(
            f"both {old.name}/ and {new.name}/ exist; merge them by hand first"
        )
        return result

    result.old_name = old.name
    if not index.is_file():
        result.status = "BLOCKED"
        result.problems.append(f"missing {index}")
        return result

    _, rows = _read_index(index)
    stale = _relative_bundle_paths(rows, old.name)
    result.index_rows_changed = sum(
        1
        for row in rows
        if any(str(value or "").startswith(f"{old.name}/") for value in row.values())
    )
    missing = [
        item
        for item in stale
        if not (root / item).is_file()
    ]
    if missing:
        result.status = "BLOCKED"
        result.problems.append(
            f"{old.name} index references {len(missing)} missing file(s); "
            "repair the project before migrating"
        )
        return result

    if dry_run:
        result.status = "DRY_RUN"
        result.removed_links, result.kept_directories = _remove_compat_dirs(
            root, old, dry_run=True
        )
        return result

    try:
        old.rename(new)
    except OSError as exc:
        result.status = "FAIL"
        result.problems.append(f"cannot rename {old} -> {new}: {exc}")
        return result

    result.index_rows_changed = _rewrite_index(
        index, old.name, OUTPUT_DIRNAME, dry_run=False
    )
    result.removed_links, result.kept_directories = _remove_compat_dirs(
        root, new, dry_run=False
    )
    _rewrite_report(root, old.name, dry_run=False)
    _update_config(root)
    result.status = "RENAMED"
    return result


def render_text(result: LayoutMigrationResult) -> str:
    lines = [f"项目根 : {result.project_root}"]
    if result.status == "NOTHING":
        lines.append(f"状态   : 没有 {LEGACY_OUTPUT_DIRNAMES[0]}/，无需迁移")
    elif result.status == "ALREADY":
        lines.append(f"状态   : 已是 {OUTPUT_DIRNAME}/ 布局")
    elif result.status == "RENAMED":
        lines.append(f"状态   : {result.old_name}/ -> {OUTPUT_DIRNAME}/")
        lines.append(f"索引   : 更新 {result.index_rows_changed} 行")
        if result.removed_links:
            lines.append(
                "清理   : 移除 "
                + f"{len(result.removed_links)} 个根目录兼容文件夹（硬链接，未删数据）"
            )
    elif result.status == "DRY_RUN":
        lines.append(
            f"状态   : 预演（未改动），将重命名 {result.old_name}/ -> {OUTPUT_DIRNAME}/"
        )
        if result.removed_links:
            lines.append(
                "        并将移除 "
                + f"{len(result.removed_links)} 个根目录兼容文件夹"
            )
    else:
        lines.append(f"状态   : {result.status}")
    if result.kept_directories:
        lines.append(
            "保留   : " + ", ".join(result.kept_directories) + "（含额外文件，未删除）"
        )
    for problem in result.problems:
        lines.append(f"[FAIL] {problem}")
    return "\n".join(lines)
