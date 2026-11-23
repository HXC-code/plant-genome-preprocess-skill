"""Project layout resolution.

Raw species folders live directly under ``data/``:

    <project>/data/{show_id}_{Order}_{Family}({Species})/...

Earlier releases nested them under ``data/downloads/``. Projects that still
have that layout keep working: the legacy directory stays authoritative while
it contains species folders, so nothing has to be moved by hand.
"""

from __future__ import annotations

import json
from pathlib import Path

DATA_RELATIVE = Path("data")
LEGACY_DOWNLOADS_RELATIVE = DATA_RELATIVE / "downloads"

# One container for every preprocessing product.  Each species gets a folder
# named after its abbreviation, so WGDI can be pointed straight at
# ``wgdi_input/<code>`` and the project root stays readable.
OUTPUT_DIRNAME = "wgdi_input"
# Scan metadata and the small pipeline indexes are named by content, not by
# the historic "resources" label.  A compatibility alias keeps old projects
# working; new projects write into metadata/.
METADATA_DIRNAME = "metadata"
LEGACY_METADATA_DIRNAME = "resources"
MACHINE_CONFIG_RELATIVE = Path(METADATA_DIRNAME) / "machine_paths.json"


def index_file(project_root: Path | str, name: str) -> Path:
    """Resolve a small pipeline index with a lossless legacy fallback.

    Existing metadata wins.  If a project was created while a compatibility
    alias could not be made, the same file may still live in ``resources/`` or
    at the project root; neither location should make the pipeline disagree
    with the file that is actually present.
    """

    root = Path(project_root)
    candidates = (
        root / METADATA_DIRNAME / name,
        root / LEGACY_METADATA_DIRNAME / name,
        root / name,
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    for candidate in candidates:
        if candidate.parent.is_dir():
            return candidate
    return candidates[-1]


def machine_config_file(project_root: Path | str) -> Path:
    """Locate the scan path map in the canonical or legacy metadata folder."""

    return index_file(project_root, "machine_paths.json")


# Earlier releases wrote the same bundles to ``analysis_v2`` (and, before that,
# a second ``clean_v2`` evidence layer).  Those names stay readable so an old
# project keeps working until it is migrated with ``migrate-layout``.
LEGACY_OUTPUT_DIRNAMES = ("analysis_v2",)

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


def configured_data_root(project_root: Path | str) -> Path | None:
    """Raw-data directory recorded by ``discover``, when it lives elsewhere.

    A new machine may keep genomes on another drive while the project (and its
    ``resources/``) sits wherever the user works.  ``machine_paths.json``
    records that location, so the pipeline does not have to copy terabytes or
    rely on a symbolic link.
    """

    config = machine_config_file(project_root)
    if not config.is_file():
        return None
    try:
        payload = json.loads(config.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(payload, dict):
        return None
    value = str((payload.get("paths") or {}).get("data_root") or "").strip()
    if not value:
        return None
    candidate = Path(value).expanduser()
    if not candidate.is_absolute():
        candidate = Path(project_root) / candidate
    return candidate if candidate.is_dir() else None


def downloads_root(project_root: Path | str) -> Path:
    """Return the directory that holds raw species folders.

    A path map written by ``discover`` wins (data may live on another drive),
    then ``data/``.  ``data/`` wins whenever it contains species folders, so a
    partially migrated project resolves to the new layout.
    ``data/downloads`` is used only while it holds species folders and
    ``data/`` holds none; an empty legacy directory never captures resolution.
    """

    root = Path(project_root)
    configured = configured_data_root(root)
    if configured is not None and any(is_species_directory(p) for p in configured.iterdir()):
        return configured
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


def configured_output_dirname(project_root: Path | str) -> str | None:
    """Output directory recorded by ``discover`` (``paths.output_dirname``).

    Only a single path component is accepted: the bundle root always lives
    directly under the project root, next to ``resources/``.
    """

    config = machine_config_file(project_root)
    if not config.is_file():
        return None
    try:
        payload = json.loads(config.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(payload, dict):
        return None
    value = str((payload.get("paths") or {}).get("output_dirname") or "").strip()
    if not value or value in {".", ".."}:
        return None
    if "/" in value or "\\" in value or Path(value).is_absolute():
        return None
    return value


def output_dirname(project_root: Path | str) -> str:
    """Name of the active bundle container.

    An existing directory wins over configuration, so a project is never
    silently redirected while its data sits somewhere else.  ``wgdi_input``
    wins over the legacy ``analysis_v2`` when both exist.
    """

    root = Path(project_root)
    if (root / OUTPUT_DIRNAME).is_dir():
        return OUTPUT_DIRNAME
    for legacy in LEGACY_OUTPUT_DIRNAMES:
        if (root / legacy).is_dir():
            return legacy
    configured = configured_output_dirname(root)
    return configured or OUTPUT_DIRNAME


def output_root(project_root: Path | str) -> Path:
    """Directory holding ``<code>/<code>.*`` bundles for the active layout."""

    return Path(project_root) / output_dirname(project_root)


def is_legacy_output(project_root: Path | str) -> bool:
    """True while the project still uses a pre-``wgdi_input`` container."""

    return output_dirname(project_root) != OUTPUT_DIRNAME


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
