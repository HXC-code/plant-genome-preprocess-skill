"""Project layout resolution.

Raw species folders live directly under ``data/``:

    <project>/data/{show_id}_{Order}_{Family}({Species})/...

Earlier releases nested them under ``data/downloads/``. Projects that still
have that layout keep working: the legacy directory stays authoritative while
it contains species folders, so nothing has to be moved by hand.
"""

from __future__ import annotations

from pathlib import Path

DATA_RELATIVE = Path("data")
LEGACY_DOWNLOADS_RELATIVE = DATA_RELATIVE / "downloads"

# Directories that may legitimately sit under data/ but never hold species.
NON_SPECIOUS_NAMES = frozenset({"downloads", "derived", "normalized"})


def is_species_directory(path: Path) -> bool:
    """A species folder is ``<show_id>_<Order>_...``, keyed by its numeric id."""

    return (
        path.is_dir()
        and path.name not in NON_SPECIOUS_NAMES
        and "_" in path.name
        and path.name.split("_", 1)[0].isdigit()
    )


def downloads_root(project_root: Path | str) -> Path:
    """Return the directory that holds raw species folders.

    ``data/`` wins whenever it contains species folders, so a partially
    migrated project resolves to the new layout. ``data/downloads`` is used
    only while it holds species folders and ``data/`` holds none; an empty
    legacy directory never captures resolution.
    """

    root = Path(project_root)
    data = root / DATA_RELATIVE
    if data.is_dir() and any(is_species_directory(p) for p in data.iterdir()):
        return data
    legacy = root / LEGACY_DOWNLOADS_RELATIVE
    if legacy.is_dir() and any(is_species_directory(p) for p in legacy.iterdir()):
        return legacy
    return data


def species_directories(project_root: Path | str) -> list[Path]:
    """Every species folder under the resolved raw-data root, sorted by show_id."""

    base = downloads_root(project_root)
    if not base.is_dir():
        return []
    return sorted(
        (path for path in base.iterdir() if is_species_directory(path)),
        key=lambda path: path.name,
    )


def species_by_show_id(project_root: Path | str) -> dict[str, Path]:
    return {path.name.split("_", 1)[0]: path for path in species_directories(project_root)}


def relative_to_data(path: Path | str, project_root: Path | str) -> str:
    """Project-relative POSIX path, tolerating both layouts.

    The deepest matching base wins, so a legacy ``data/downloads/<species>``
    file is reported without the ``downloads/`` prefix even when ``data/``
    itself is the active root.
    """

    resolved = Path(path).resolve()
    candidates = (
        Path(project_root) / LEGACY_DOWNLOADS_RELATIVE,
        downloads_root(project_root),
        Path(project_root) / DATA_RELATIVE,
    )
    best: tuple[int, str] | None = None
    for base in candidates:
        try:
            relative = resolved.relative_to(base.resolve())
        except ValueError:
            continue
        depth = len(base.resolve().parts)
        if best is None or depth > best[0]:
            best = (depth, relative.as_posix())
    return best[1] if best else resolved.as_posix()
