"""Remove optional artifacts from an existing project, keeping it valid.

Some artifacts are useful while building a dataset and pure weight afterwards:

* ``<Abbr>.all.*`` — the pre-filtering snapshot of every candidate gene.
  Nothing downstream reads it; on the 41-species pilot it was 51% of
  ``analysis_v2`` (2.4 GB). New runs skip it unless ``run --keep-all``.
* root ``<Abbr>/`` compatibility directories — hard links used only by the
  0-6 shell layer. They cost no extra bytes but add one directory per species.

Deleting files alone would leave the index pointing at absent paths, so
pruning drops the files and blanks the affected index cells in one step.
"""

from __future__ import annotations

from pathlib import Path

from pgci_preprocess.artifact_naming import OPTIONAL_ANALYSIS_SUFFIXES
from pgci_preprocess.reports import (
    INDEX_FIELDS,
    read_tsv,
    write_tsv,
)

ALL_INDEX_KEYS = {
    ".all.pep.fa": "all_pep",
    ".all.cds.fa": "all_cds",
    ".all.gff": "all_gff",
    ".all.lens": "all_lens",
}


def _index_path(project_root: Path) -> Path:
    return project_root / "resources" / "wgdi_v2_index.tsv"


def prune_project(
    project_root: Path | str,
    *,
    all_snapshots: bool = False,
    compat_links: bool = False,
    dry_run: bool = False,
) -> dict[str, object]:
    """Drop the requested optional artifacts and re-derive index + manifest."""

    root = Path(project_root).resolve()
    index_path = _index_path(root)
    if not index_path.is_file():
        raise FileNotFoundError(f"no index at {index_path}")

    rows = read_tsv(index_path)
    codes = [row.get("code", "").strip() for row in rows if row.get("code", "").strip()]

    removed_files = 0
    removed_dirs = 0
    removed_bytes = 0
    blanked_rows = 0

    if all_snapshots:
        for row in rows:
            code = (row.get("code") or "").strip()
            if not code:
                continue
            changed = False
            for suffix in OPTIONAL_ANALYSIS_SUFFIXES:
                for layer in ("analysis_v2", "clean_v2"):
                    path = root / layer / code / f"{code}{suffix}"
                    if not path.is_file():
                        continue
                    removed_bytes += path.stat().st_size
                    removed_files += 1
                    if not dry_run:
                        path.unlink()
                key = ALL_INDEX_KEYS.get(suffix)
                if key and row.get(key, ""):
                    # Only blank the cell when the file is really gone.
                    if dry_run or not (root / str(row[key])).is_file():
                        row[key] = ""
                        changed = True
            if changed:
                blanked_rows += 1

    if compat_links:
        for code in codes:
            link_dir = root / code
            if not link_dir.is_dir():
                continue
            entries = [p for p in link_dir.rglob("*") if p.is_file()]
            removed_files += len(entries)
            removed_dirs += 1
            if not dry_run:
                for path in sorted(link_dir.rglob("*"), reverse=True):
                    if path.is_file() or path.is_symlink():
                        path.unlink()
                for path in sorted(link_dir.rglob("*"), reverse=True):
                    if path.is_dir():
                        path.rmdir()
                link_dir.rmdir()

    if not dry_run and all_snapshots:
        write_tsv(index_path, INDEX_FIELDS, rows)

    return {
        "project_root": str(root),
        "dry_run": dry_run,
        "removed_files": removed_files,
        "removed_dirs": removed_dirs,
        "removed_bytes": removed_bytes,
        "index_rows_updated": blanked_rows,
        "note": (
            "compat directories are hard links: deleting them frees no extra "
            "bytes, they only remove clutter"
            if compat_links
            else ""
        ),
    }
