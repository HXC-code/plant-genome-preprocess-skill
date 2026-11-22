"""Rename legacy ``s<show_id>`` bundles to the canonical abbreviation.

The 2026-11-19 naming migration only rewrote metadata (index, QC ``code``
fields, compatibility links).  Directories, file names and, most importantly,
the canonical gene IDs inside ``*.gff`` / ``*.cds`` / ``*.pep`` stayed legacy.
This module performs the complete migration:

1. move ``<layer>/s<show_id>`` to a staging directory;
2. rename every file ``s<show_id>*`` -> ``<code>*``;
3. replace every gene ID ``s<show_id>g<order>`` -> ``<code>g<order>``
   (byte-level streaming rewrite; the numeric suffix is preserved);
4. verify file names and that no legacy token survives;
5. publish the staging directory under ``<layer>/<code>``;
6. refresh the index and compatibility links;
7. optionally move the stale ``<workspace>/s<show_id>`` copies to quarantine.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import uuid
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from pgci_preprocess.artifact_naming import (
    ANALYSIS_REQUIRED_SUFFIXES,
    CLEAN_REQUIRED_SUFFIXES,
    NamingViolation,
    bundle_name_problems,
    first_gene_id,
    legacy_gene_id_regex,
    legacy_prefix,
)
from pgci_preprocess.project_paths import output_dirname, output_root
from pgci_preprocess.reports import (
    read_tsv,
    read_tsv_with_fields,
    write_tsv,
)
from pgci_preprocess.species_names import (
    SpeciesNamingError,
    load_species_names,
)
from pgci_preprocess.validate import (
    BUNDLE_KEYS,
    ROLE_SUFFIX,
    resolve_layout,
    sync_compat_links,
)

LEGACY_LAYERS = ("clean_v2", "analysis_v2")


def _layers(canonical_root: Path) -> tuple[str, ...]:
    """Every directory that may hold a bundle, active container last."""

    return tuple(
        dict.fromkeys((*LEGACY_LAYERS, output_dirname(canonical_root)))
    )
TEXT_SUFFIXES = (
    ".gff",
    ".cds",
    ".pep",
    ".lens",
    ".fa",
    ".tsv",
    ".csv",
)
QC_SUFFIX = ".qc.json"
INDEX_PATH = Path("resources") / "wgdi_v2_index.tsv"
INDEX_PATH_FIELDS = (
    "pep",
    "cds",
    "gff",
    "lens",
    "all_pep",
    "all_cds",
    "all_gff",
    "all_lens",
    "gene_map",
    "chr_map",
    "rejected",
    "qc",
)


@dataclass
class SpeciesMigration:
    show_id: str
    code: str
    status: str
    files_renamed: int = 0
    ids_replaced: int = 0
    bytes_rewritten: int = 0
    problems: list[str] = field(default_factory=list)
    pruned_legacy_dir: str = ""

    def as_dict(self) -> dict[str, object]:
        return {
            "show_id": self.show_id,
            "code": self.code,
            "status": self.status,
            "files_renamed": self.files_renamed,
            "ids_replaced": self.ids_replaced,
            "bytes_rewritten": self.bytes_rewritten,
            "problems": list(self.problems),
            "pruned_legacy_dir": self.pruned_legacy_dir,
        }


@dataclass
class ArtifactMigrationResult:
    project_root: Path
    canonical_root: Path
    report_path: Path
    quarantine_dir: Path | None
    species: list[SpeciesMigration]
    links: list[str]
    index_rows_changed: int
    dry_run: bool

    def as_dict(self) -> dict[str, object]:
        counts: dict[str, int] = {}
        for item in self.species:
            counts[item.status] = counts.get(item.status, 0) + 1
        return {
            "project_root": str(self.project_root),
            "canonical_root": str(self.canonical_root),
            "report_path": str(self.report_path),
            "dry_run": self.dry_run,
            "species_counts": counts,
            "species": [item.as_dict() for item in self.species],
            "links": self.links[:200],
            "index_rows_changed": self.index_rows_changed,
            "quarantine_dir": str(self.quarantine_dir) if self.quarantine_dir else "",
        }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_write_bytes(path: Path, payload: bytes) -> None:
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        temporary.write_bytes(payload)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _id_pattern(path: Path, show_id: str) -> re.Pattern[bytes] | None:
    """Regex that matches a legacy gene ID in its canonical position only.

    * GFF: second column (``<seqid>\t<gene_id>\t``, group 1 = prefix + tab);
    * FASTA: first header token (``>``);
    * TSV: first column (``<canonical_gene_id>\t``).

    Restricting the match to the canonical field avoids rewriting unrelated
    provenance columns (source gene IDs) that happen to look similar.
    """

    sid = str(show_id).encode()
    suffix = path.suffix.lower()
    if suffix == ".gff":
        return re.compile(rb"(?m)^([^\t]*\t)s" + sid + rb"(\d+g\d{4,})\t")
    if suffix in {".fa", ".cds", ".pep"}:
        return re.compile(rb"(?m)^>s" + sid + rb"(\d+g\d{4,})(?=[ \t\r\n])")
    if suffix == ".tsv":
        return re.compile(rb"(?m)^s" + sid + rb"(\d+g\d{4,})\t")
    return None


def _rewrite_file_ids(path: Path, show_id: str, code: str) -> tuple[int, int]:
    """Replace legacy gene IDs and return ``(ids_replaced, bytes_before)``."""

    pattern = _id_pattern(path, show_id)
    if pattern is None:
        return 0, 0
    payload = path.read_bytes()
    if b"s" + str(show_id).encode() not in payload:
        return 0, len(payload)
    if path.suffix.lower() == ".gff":
        # Group 1 ends with the tab that separates seqid and gene ID.
        replacement = rb"\1" + code.encode() + rb"\2\t"
    elif path.suffix.lower() == ".tsv":
        replacement = code.encode() + rb"\1\t"
    else:
        # The FASTA pattern consumes the leading ">"; keep it in the output.
        replacement = rb">" + code.encode() + rb"\1"
    replaced, count = pattern.subn(replacement, payload)
    if count:
        _atomic_write_bytes(path, replaced)
    if pattern.search(replaced):
        raise NamingViolation(
            f"{path}: legacy gene-ID tokens survived the rewrite"
        )
    return count, len(payload)


def _repair_fasta_headers(path: Path, code: str) -> int:
    """Restore the ``>`` marker dropped by an older migration attempt."""

    payload = path.read_bytes()
    header = re.compile(rb"(?m)^" + re.escape(code.encode()) + rb"\d+g\d{4,}(?=[ \t])")
    lines = payload.splitlines(keepends=True)
    changed = 0
    for index, raw in enumerate(lines):
        body = raw.rstrip(b"\r\n")
        if header.match(body) and not raw.startswith(b">"):
            ending = b"\r\n" if raw.endswith(b"\r\n") else b"\n"
            lines[index] = b">" + body + ending
            changed += 1
    if changed:
        _atomic_write_bytes(path, b"".join(lines))
    return changed


def _repair_gff_columns(path: Path, code: str) -> int:
    """Rebuild the gene-ID column of GFF rows damaged by an older migration.

    A previously interrupted rename could leave rows shaped as
    ``<code><seq>\\t<start>\\t<end>\\t<strand>\\t<order>\\t<source>``: the
    original gene ID was dropped.  The gene ID is fully determined by
    ``code + sequence + order``, so the row can be rebuilt deterministically.
    """

    payload = path.read_bytes()
    lines = payload.splitlines(keepends=True)
    changed = 0
    bogus = re.compile(re.escape(code.encode()) + rb"b'(\d+)'g(\d+)")
    for index, raw in enumerate(lines):
        if not raw.strip() or raw.startswith(b"#"):
            continue
        fields = raw.rstrip(b"\r\n").split(b"\t")
        if len(fields) == 7:
            if bogus.fullmatch(fields[1]):
                sequence = fields[0]
                if not sequence.isdigit() or not fields[5].isdigit():
                    raise NamingViolation(f"{path}: cannot rebuild gene ID from {raw[:60]!r}")
                gene_id = f"{code}{sequence.decode()}g{int(fields[5]):04d}".encode()
                fields[1] = gene_id
                ending = b"\r\n" if raw.endswith(b"\r\n") else b"\n"
                lines[index] = b"\t".join(fields) + ending
                changed += 1
            continue
        if len(fields) != 6:
            raise NamingViolation(f"{path}: unexpected GFF column count {len(fields)}")
        sequence_field, order_field = fields[0], fields[4]
        if not sequence_field.startswith(code.encode()):
            raise NamingViolation(
                f"{path}: cannot rebuild gene ID from {sequence_field!r}"
            )
        sequence = sequence_field[len(code):]
        if not sequence.isdigit() or not order_field.isdigit():
            raise NamingViolation(
                f"{path}: cannot rebuild gene ID from {raw[:60]!r}"
            )
        gene_id = f"{code}{sequence.decode('ascii')}g{int(order_field):04d}".encode()
        ending = b"\r\n" if raw.endswith(b"\r\n") else b"\n"
        lines[index] = b"\t".join([sequence, gene_id, *fields[1:]]) + ending
        changed += 1
    if changed:
        _atomic_write_bytes(path, b"".join(lines))
    return changed


def _migrate_qc(path: Path, show_id: str, code: str) -> None:
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError):
        return
    if not isinstance(record, dict):
        return
    record["code"] = code
    record["show_id"] = show_id
    record["artifact_code"] = code
    record["legacy_code"] = legacy_prefix(show_id)
    path.write_text(
        json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def _precheck_legacy_ids(directory: Path, show_id: str) -> None:
    pattern = re.compile(legacy_gene_id_regex(show_id))
    for suffix, kind in (
        (".new.gff", "gff"),
        (".all.gff", "gff"),
        (".cds", "fasta"),
        (".pep", "fasta"),
        (".all.cds.fa", "fasta"),
        (".all.pep.fa", "fasta"),
    ):
        path = directory / f"{legacy_prefix(show_id)}{suffix}"
        if not path.is_file():
            continue
        gene_id = first_gene_id(path, kind)
        if gene_id and not pattern.fullmatch(gene_id):
            raise NamingViolation(
                f"{path}: first gene ID {gene_id!r} is not "
                f"'s{show_id}g<order>'; refusing to rename a mixed bundle"
            )


def _rename_tree(staging: Path, show_id: str, code: str) -> tuple[int, int, int]:
    """Rename files and rewrite IDs inside an already-moved staging bundle."""

    files_renamed = 0
    ids_replaced = 0
    bytes_rewritten = 0
    for path in sorted(staging.rglob("*")):
        if not path.is_file():
            continue
        old_prefix = legacy_prefix(show_id)
        if path.name.startswith(old_prefix):
            target = path.with_name(code + path.name[len(old_prefix):])
            if target.exists():
                raise NamingViolation(f"rename target already exists: {target}")
            path.rename(target)
            path = target
            files_renamed += 1
        suffix = path.suffix.lower()
        if suffix == ".json" and path.name.endswith(QC_SUFFIX):
            _migrate_qc(path, show_id, code)
            continue
        if suffix == ".gff":
            _repair_gff_columns(path, code)
        elif suffix in {".fa", ".cds", ".pep"}:
            _repair_fasta_headers(path, code)
        if suffix in TEXT_SUFFIXES:
            count, size = _rewrite_file_ids(path, show_id, code)
            if count:
                ids_replaced += count
                bytes_rewritten += size
    return files_renamed, ids_replaced, bytes_rewritten


def _migrate_layer(
    canonical_root: Path,
    layer: str,
    show_id: str,
    code: str,
    *,
    dry_run: bool,
    item: SpeciesMigration,
    required: tuple[str, ...],
) -> None:
    legacy_dir = canonical_root / layer / legacy_prefix(show_id)
    canonical_dir = canonical_root / layer / code
    if not legacy_dir.is_dir():
        if canonical_dir.is_dir():
            problems = bundle_name_problems(
                canonical_dir,
                code,
                required_suffixes=required,
                kind_map={".new.gff": "gff", ".cds": "fasta", ".pep": "fasta"},
            )
            if problems:
                item.problems.extend(f"{layer}:{problem}" for problem in problems)
                item.status = "FAIL"
        return
    if canonical_dir.exists():
        item.problems.append(f"{layer}: both legacy and canonical directories exist")
        item.status = "FAIL"
        return
    _precheck_legacy_ids(legacy_dir, show_id)
    if dry_run:
        item.files_renamed += sum(
            1
            for path in legacy_dir.rglob("*")
            if path.is_file() and path.name.startswith(legacy_prefix(show_id))
        )
        return

    staging = canonical_root / layer / f".{code}.migrate-{uuid.uuid4().hex}"
    os.replace(legacy_dir, staging)
    try:
        files_renamed, ids_replaced, bytes_rewritten = _rename_tree(
            staging, show_id, code
        )
        item.files_renamed += files_renamed
        item.ids_replaced += ids_replaced
        item.bytes_rewritten += bytes_rewritten
        problems = bundle_name_problems(
            staging,
            code,
            required_suffixes=required,
            kind_map={".new.gff": "gff", ".cds": "fasta", ".pep": "fasta"},
            check_dir_name=False,
        )
        if problems:
            raise NamingViolation(";".join(problems))
        os.replace(staging, canonical_dir)
        if legacy_dir.is_dir() and not any(legacy_dir.iterdir()):
            legacy_dir.rmdir()
    except Exception:
        if staging.exists() and not canonical_dir.exists():
            os.replace(staging, legacy_dir)
        raise


def _migrate_index(
    canonical_root: Path,
    migrated: dict[str, str],
    legacy_codes: dict[str, str],
) -> int:
    path = canonical_root / INDEX_PATH
    if not path.is_file():
        return 0
    fields, rows = read_tsv_with_fields(path)
    changed = 0
    for row in rows:
        show_id = str(row.get("show_id", "")).strip()
        if not show_id:
            code_field = str(row.get("code", "")).strip()
            if code_field.casefold().startswith("s") and code_field[1:].isdigit():
                show_id = code_field[1:]
        code = migrated.get(show_id)
        if not code:
            continue
        legacy_code = legacy_codes.get(show_id, legacy_prefix(show_id))
        before = dict(row)
        if "show_id" in fields:
            row["show_id"] = show_id
        if "code" in fields:
            row["code"] = code
        if "artifact_code" in fields:
            row["artifact_code"] = code
        for field_name in INDEX_PATH_FIELDS:
            value = str(row.get(field_name, ""))
            if not value:
                continue
            row[field_name] = value.replace(
                f"/{legacy_code}/", f"/{code}/"
            ).replace(f"{legacy_code}.", f"{code}.")
        if row != before:
            changed += 1
    if changed:
        write_tsv(path, fields, rows)
    return changed


def _migrated_codes_from_layers(canonical_root: Path) -> dict[str, str]:
    codes: dict[str, str] = {}
    for layer in _layers(canonical_root):
        root = canonical_root / layer
        if not root.is_dir():
            continue
        for path in root.iterdir():
            if path.is_dir() and re.fullmatch(r"s\d+", path.name):
                codes.setdefault(path.name[1:], "")
    return codes


def _legacy_root_values(canonical_root: Path) -> dict[str, str]:
    """Read ``legacy_code`` / ``artifact_code`` per show_id from the index."""

    values: dict[str, str] = {}
    path = canonical_root / INDEX_PATH
    if not path.is_file():
        return values
    for row in read_tsv(path):
        show_id = str(row.get("show_id", "")).strip()
        if not show_id:
            code = str(row.get("code", "")).strip()
            if code.casefold().startswith("s") and code[1:].isdigit():
                show_id = code[1:]
        if show_id:
            values[show_id] = str(
                row.get("artifact_code") or row.get("code") or ""
            ).strip()
    return values


def _prune_legacy_compat(
    workspace: Path,
    canonical_root: Path,
    code_by_show_id: dict[str, str],
    *,
    dry_run: bool,
) -> tuple[int, Path | None, list[str]]:
    """Move byte-identical stale ``<workspace>/s<show_id>`` copies to quarantine."""

    candidates: list[tuple[Path, Path]] = []
    for show_id, code in code_by_show_id.items():
        legacy_dir = workspace / legacy_prefix(show_id)
        canonical_dir = output_root(canonical_root) / code
        if not legacy_dir.is_dir() or not canonical_dir.is_dir():
            continue
        identical = True
        for key in BUNDLE_KEYS:
            old = legacy_dir / f"{legacy_prefix(show_id)}{ROLE_SUFFIX[key]}"
            new = canonical_dir / f"{code}{ROLE_SUFFIX[key]}"
            if not old.is_file() or not new.is_file() or _sha256(old) != _sha256(new):
                identical = False
                break
        if identical:
            candidates.append((legacy_dir, canonical_dir))
    if not candidates:
        return 0, None, []
    stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    quarantine = workspace / "_quarantine" / f"legacy_compat_dirs_{stamp}"
    messages = [f"{'would move' if dry_run else 'moved'}: {path}" for path, _ in candidates]
    if dry_run:
        return len(candidates), quarantine, messages
    quarantine.mkdir(parents=True, exist_ok=False)
    for path, _ in candidates:
        shutil.move(str(path), str(quarantine / path.name))
    return len(candidates), quarantine, messages


def migrate_artifacts(
    project_root: Path | str,
    *,
    dry_run: bool = False,
    refresh_links: bool = True,
    prune_legacy_dirs: bool = False,
    show_ids: set[str] | None = None,
) -> ArtifactMigrationResult:
    workspace, canonical_root = resolve_layout(project_root)
    names = load_species_names(canonical_root)
    if names is None:
        raise SpeciesNamingError(
            "resources/species_abbreviations.tsv is required for the migration"
        )
    legacy_codes = _legacy_root_values(canonical_root)
    found = sorted(
        set(_migrated_codes_from_layers(canonical_root))
        | {
            show_id
            for show_id, artifact_code in legacy_codes.items()
            if artifact_code.casefold().startswith("s")
        }
    )
    selected = sorted(
        show_id for show_id in found if show_ids is None or show_id in show_ids
    )

    items: list[SpeciesMigration] = []
    migrated: dict[str, str] = {}
    layers = _layers(canonical_root)
    for show_id in selected:
        code = names.code_for(show_id)
        item = SpeciesMigration(show_id=show_id, code=code, status="MIGRATED")
        has_legacy = any(
            (canonical_root / layer / legacy_prefix(show_id)).is_dir()
            for layer in layers
        )
        canonical_complete = (
            canonical_root / output_dirname(canonical_root) / code
        ).is_dir()
        if not has_legacy and not canonical_complete:
            continue
        if not has_legacy and canonical_complete:
            problems: list[str] = []
            for layer in layers:
                required = (
                    CLEAN_REQUIRED_SUFFIXES
                    if layer == "clean_v2"
                    else ANALYSIS_REQUIRED_SUFFIXES
                )
                if not (canonical_root / layer / code).is_dir():
                    continue
                problems.extend(
                    f"{layer}:{problem}"
                    for problem in bundle_name_problems(
                        canonical_root / layer / code,
                        code,
                        required_suffixes=required,
                        kind_map={
                            ".new.gff": "gff",
                            ".cds": "fasta",
                            ".pep": "fasta",
                        },
                    )
                )
            if problems:
                item.status = "FAIL"
                item.problems.extend(problems)
            else:
                item.status = "ALREADY_CANONICAL"
            items.append(item)
            migrated[show_id] = code
            continue
        try:
            for layer in layers:
                _migrate_layer(
                    canonical_root,
                    layer,
                    show_id,
                    code,
                    dry_run=dry_run,
                    item=item,
                    required=(
                        CLEAN_REQUIRED_SUFFIXES
                        if layer == "clean_v2"
                        else ANALYSIS_REQUIRED_SUFFIXES
                    ),
                )
        except (NamingViolation, OSError) as exc:
            item.status = "FAIL"
            item.problems.append(str(exc))
        if item.status != "FAIL":
            item.status = "DRY_RUN" if dry_run else "MIGRATED"
            migrated[show_id] = code
        items.append(item)

    index_changed = 0
    links: list[str] = []
    quarantine_dir: Path | None = None
    pruned = 0
    if not dry_run:
        index_changed = _migrate_index(canonical_root, migrated, legacy_codes)
        if refresh_links:
            links = sync_compat_links(canonical_root, workspace)
    if prune_legacy_dirs and migrated and not dry_run:
        pruned, quarantine_dir, prune_messages = _prune_legacy_compat(
            workspace, canonical_root, migrated, dry_run=dry_run
        )
        links.extend(prune_messages)

    report_path = canonical_root / "resources" / "artifact_naming_migration.json"
    result = ArtifactMigrationResult(
        project_root=workspace,
        canonical_root=canonical_root,
        report_path=report_path,
        quarantine_dir=quarantine_dir,
        species=items,
        links=links,
        index_rows_changed=index_changed,
        dry_run=dry_run,
    )
    if not dry_run:
        report_path.write_text(
            json.dumps(result.as_dict(), ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
    return result
