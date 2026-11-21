"""Create a canonical project and optionally register one explicit species bundle."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import tempfile
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Any, Iterable, Mapping

from pgci_preprocess.batch_prepare import ASSEMBLY_MAP_FIELDS, INVENTORY_FIELDS
from pgci_preprocess.species_names import (
    load_species_names,
    resolve_species_code,
    species_name_path,
)


class ProjectInitError(RuntimeError):
    """Raised when a project cannot be initialized without guessing or overwriting."""


# 发布包自带的默认资源：新建项目时自动复制，用户无需再手动找表
BUNDLED_RESOURCES = Path(__file__).resolve().parents[2] / "resources"


@dataclass(frozen=True)
class ProjectInitResult:
    project_root: Path
    inventory_path: Path
    assembly_map_path: Path
    show_id: str
    species_code: str
    species_directory: Path | None
    files: dict[str, str]
    created_paths: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "project_root": str(self.project_root),
            "inventory": str(self.inventory_path),
            "assembly_map": str(self.assembly_map_path),
            "show_id": self.show_id,
            "species_code": self.species_code,
            "species_directory": str(self.species_directory) if self.species_directory else "",
            "files": self.files,
            "created_paths": list(self.created_paths),
        }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _clean_field(value: str, label: str) -> str:
    value = value.strip()
    if not value or any(character in value for character in "\t\r\n\x00"):
        raise ProjectInitError(f"{label} must be non-empty and cannot contain tabs or newlines")
    return value


def _safe_component(value: str, label: str) -> str:
    value = _clean_field(value, label)
    if value in {".", ".."} or "/" in value or "\\" in value:
        raise ProjectInitError(f"{label} cannot contain path separators: {value!r}")
    return value


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ProjectInitError(f"TSV has no header: {path}")
        return list(reader.fieldnames), list(reader)


def _serialize_tsv(fields: Iterable[str], rows: Iterable[Mapping[str, Any]]) -> str:
    output = StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
    return output.getvalue()


def _atomic_write(path: Path, content: str | bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    binary = isinstance(content, bytes)
    mode = "wb" if binary else "w"
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        if binary:
            with os.fdopen(descriptor, mode) as handle:
                handle.write(content)
        else:
            with os.fdopen(descriptor, mode, encoding="utf-8", newline="") as handle:
                handle.write(content)
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def _write_tsv_if_missing(path: Path, fields: list[str], created: list[str]) -> None:
    if path.exists():
        if not path.is_file():
            raise ProjectInitError(f"expected a file: {path}")
        return
    _atomic_write(path, _serialize_tsv(fields, []))
    created.append(str(path))


def _mirror(source: Path, destination: Path) -> str:
    try:
        os.link(source, destination)
        transfer = "hardlink"
    except OSError:
        shutil.copy2(source, destination)
        transfer = "copy"
    if _sha256(source) != _sha256(destination):
        raise ProjectInitError(f"hash mismatch after mirroring: {source}")
    return transfer


_ASSEMBLY_EXTENSIONS = (
    ".genome.fasta",
    ".genome.fa",
    ".genome.fna",
    ".genome.fas",
    ".fasta",
    ".fa",
    ".fna",
    ".fas",
)


def _infer_assembly_id(genome: Path | str) -> str:
    """Derive an assembly id from the genome file name when the user omits one.

    ``GWHBHAK00000000.genome.fasta.gz`` -> ``GWHBHAK00000000`` and
    ``GCF_000001.1_ASM123v1_genomic.fna.gz`` keeps its versioned accession.
    Returns an empty string when nothing sensible can be derived, so the caller
    still reports the explicit error.
    """

    name = Path(str(genome)).name
    for compression in (".gz", ".bz2", ".xz"):
        if name.lower().endswith(compression):
            name = name[: -len(compression)]
    lowered = name.lower()
    for extension in _ASSEMBLY_EXTENSIONS:
        if lowered.endswith(extension):
            name = name[: -len(extension)]
            break
    return name.strip()


def init_project(
    project_root: Path | str,
    *,
    show_id: str = "",
    scientific_name: str = "",
    order: str = "",
    family: str = "",
    assembly_id: str = "",
    sources: dict[str, Path | str] | None = None,
    species_names_path: Path | str | None = None,
) -> ProjectInitResult:
    """Initialize directories or register one species without modifying source files.

    With no species arguments this creates only the standard directory and TSV
    headers. Registration requires explicit GFF/GFF3, CDS and protein paths;
    genome FASTA and assembly identity are optional. No biological role is
    inferred from a suffix.
    """

    root = Path(project_root).expanduser().resolve()
    sources = sources or {}
    unknown_roles = sorted(set(sources) - {"genome", "gff", "cds", "pep"})
    if unknown_roles:
        raise ProjectInitError(f"unknown input roles: {', '.join(unknown_roles)}")
    registration_requested = bool(
        show_id
        or scientific_name
        or order
        or family
        or assembly_id
        or sources
    )
    if root.exists() and not root.is_dir():
        raise ProjectInitError(f"project root is a file: {root}")

    downloads = root / "data" / "downloads"
    resources = root / "resources"
    created: list[str] = []
    for directory in (root, downloads, resources):
        if not directory.exists():
            directory.mkdir(parents=True, exist_ok=True)
            created.append(str(directory))
    inventory_path = resources / "download_inventory.tsv"
    assembly_map_path = resources / "assembly_accessions.tsv"
    _write_tsv_if_missing(inventory_path, INVENTORY_FIELDS, created)
    _write_tsv_if_missing(assembly_map_path, ASSEMBLY_MAP_FIELDS, created)
    configured_names = species_name_path(root, explicit_path=species_names_path)
    if configured_names is None:
        bundled_names = BUNDLED_RESOURCES / "species_abbreviations.tsv"
        if bundled_names.is_file():
            # 无论建空项目还是注册物种都带上缩写表：表里的 207+ 个物种直接沿用
            # 既有缩写，表里没有的新物种按命名规则现算（见 species_codegen）。
            configured_names = bundled_names
    if configured_names is not None:
        source_names = configured_names
        if not source_names.is_file():
            # 用户把 --species-names 指向“项目内尚未生成的路径”时，回落到包内自带的表，
            # 否则第一次建项目会直接报“找不到表”。
            bundled_fallback = BUNDLED_RESOURCES / "species_abbreviations.tsv"
            if bundled_fallback.is_file():
                source_names = bundled_fallback
            else:
                raise ProjectInitError(f"species abbreviation table not found: {source_names}")
        destination_names = resources / "species_abbreviations.tsv"
        if destination_names.exists():
            if _sha256(destination_names) != _sha256(source_names):
                raise ProjectInitError(
                    f"project already has a different species abbreviation table: "
                    f"{destination_names}"
                )
        elif source_names != destination_names:
            shutil.copy2(source_names, destination_names)
            created.append(str(destination_names))

    bundled_exceptions = BUNDLED_RESOURCES / "format_exceptions.tsv"
    destination_exceptions = resources / "format_exceptions.tsv"
    if bundled_exceptions.is_file() and not destination_exceptions.exists():
        shutil.copy2(bundled_exceptions, destination_exceptions)
        created.append(str(destination_exceptions))

    if not registration_requested:
        return ProjectInitResult(
            root, inventory_path, assembly_map_path, "", "", None, {}, tuple(created)
        )

    show_id = _clean_field(show_id, "show_id")
    if not re.fullmatch(r"[0-9]+", show_id):
        raise ProjectInitError("show_id must contain digits only")
    scientific_name = _safe_component(scientific_name, "scientific_name")
    order = _safe_component(order, "order")
    family = _safe_component(family, "family")
    if not str(assembly_id).strip() and "genome" in sources:
        assembly_id = _infer_assembly_id(sources["genome"])
    assembly_id = (
        _clean_field(assembly_id, "assembly_id") if str(assembly_id).strip() else ""
    )
    if assembly_id and not re.fullmatch(r"[A-Za-z0-9_.-]+", assembly_id):
        raise ProjectInitError("assembly_id may contain only letters, digits, dot, underscore and hyphen")
    missing_roles = sorted({"gff", "cds", "pep"} - set(sources))
    if missing_roles:
        raise ProjectInitError(f"missing explicit input roles: {', '.join(missing_roles)}")

    source_paths: dict[str, Path] = {}
    for role in ("genome", "gff", "cds", "pep"):
        if role not in sources:
            continue
        source = Path(sources[role]).expanduser().resolve()
        if not source.is_file() or source.stat().st_size == 0:
            raise ProjectInitError(f"missing or empty {role} input: {source}")
        source_paths[role] = source
    if len(set(source_paths.values())) != len(source_paths):
        raise ProjectInitError("registered input roles must reference distinct source files")

    inventory_fields, inventory_rows = _read_tsv(inventory_path)
    if "species_code" not in inventory_fields:
        insert_at = inventory_fields.index("showid") + 1 if "showid" in inventory_fields else 0
        inventory_fields.insert(insert_at, "species_code")
    catalog = load_species_names(root, explicit_path=species_names_path)
    taken_codes = set(catalog.by_code) if catalog is not None else set()
    for inventory_entry in inventory_rows:
        inventory_show_id = str(inventory_entry.get("showid", "")).strip()
        if inventory_show_id:
            code = resolve_species_code(
                root,
                inventory_show_id,
                latin_name=str(inventory_entry.get("Species", "")),
                explicit_path=species_names_path,
                occupied=taken_codes,
            )
            inventory_entry["species_code"] = code
            taken_codes.add(code.casefold())
    species_code = resolve_species_code(
        root,
        show_id,
        latin_name=scientific_name,
        explicit_path=species_names_path,
        occupied=taken_codes,
    )
    required_inventory = {"showid", "Order", "Family", "Species", "gff", "genome", "cds", "pep"}
    missing_columns = sorted(required_inventory - set(inventory_fields))
    if missing_columns:
        raise ProjectInitError(f"inventory missing columns: {', '.join(missing_columns)}")
    if any(str(row.get("showid", "")).strip() == show_id for row in inventory_rows):
        raise ProjectInitError(f"show_id already registered: {show_id}")
    assembly_fields, assembly_rows = _read_tsv(assembly_map_path)
    if set(ASSEMBLY_MAP_FIELDS) - set(assembly_fields):
        raise ProjectInitError("assembly map has an incompatible header")
    if any(str(row.get("show_id", "")).strip() == show_id for row in assembly_rows):
        raise ProjectInitError(f"assembly mapping already registered: {show_id}")

    species_name = f"{show_id}_{order}_{family}({scientific_name})"
    species_directory = downloads / species_name
    if species_directory.exists():
        raise ProjectInitError(f"species directory already exists: {species_directory}")

    original_inventory = inventory_path.read_bytes()
    original_assembly = assembly_map_path.read_bytes()
    staging = Path(tempfile.mkdtemp(prefix=".init-project-", dir=downloads))
    committed = False
    transfers: dict[str, str] = {}
    try:
        for role, source in source_paths.items():
            destination = staging / source.name
            transfers[role] = _mirror(source, destination)
        os.replace(staging, species_directory)
        committed = True
        inventory_row = {field: "" for field in inventory_fields}
        inventory_row.update(
            {
                "showid": show_id,
                "species_code": species_code,
                "Order": order,
                "Family": family,
                "Species": scientific_name,
                "genome": source_paths["genome"].name if "genome" in source_paths else "",
                "gff": source_paths["gff"].name,
                "cds": source_paths["cds"].name,
                "pep": source_paths["pep"].name,
            }
        )
        inventory_rows.append(inventory_row)
        inventory_rows.sort(key=lambda row: (0, int(row["showid"])) if str(row.get("showid", "")).isdigit() else (1, str(row.get("showid", ""))))
        _atomic_write(inventory_path, _serialize_tsv(inventory_fields, inventory_rows))
        if assembly_id:
            assembly_row = {field: "" for field in assembly_fields}
            assembly_row.update(
                {
                    "show_id": show_id,
                    "assembly_id": assembly_id,
                    "source": "init-project",
                    "evidence_path": str(source_paths["gff"]),
                    "evidence_sha256": _sha256(source_paths["gff"]),
                }
            )
            assembly_rows.append(assembly_row)
            assembly_rows.sort(
                key=lambda row: (
                    (0, int(row["show_id"]))
                    if str(row.get("show_id", "")).isdigit()
                    else (1, str(row.get("show_id", "")))
                )
            )
            _atomic_write(
                assembly_map_path,
                _serialize_tsv(assembly_fields, assembly_rows),
            )
        manifest = {
            "show_id": show_id,
            "species_code": species_code,
            "assembly_id": assembly_id,
            "species_directory": str(species_directory),
            "sources": {
                role: {
                    "source": str(source_paths[role]),
                    "filename": source_paths[role].name,
                    "sha256": _sha256(source_paths[role]),
                    "transfer": transfers[role],
                }
                for role in ("genome", "gff", "cds", "pep")
                if role in source_paths
            },
        }
        _atomic_write(
            resources / "project_init" / f"{species_code}.json",
            json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        )
    except Exception:
        _atomic_write(inventory_path, original_inventory)
        _atomic_write(assembly_map_path, original_assembly)
        if committed and species_directory.exists():
            shutil.rmtree(species_directory)
        raise
    finally:
        if staging.exists():
            shutil.rmtree(staging)

    created.append(str(species_directory))
    created.append(str(inventory_path))
    created.append(str(assembly_map_path))
    return ProjectInitResult(
        root,
        inventory_path,
        assembly_map_path,
        show_id,
        species_code,
        species_directory,
        {role: source_paths[role].name for role in source_paths},
        tuple(created),
    )
