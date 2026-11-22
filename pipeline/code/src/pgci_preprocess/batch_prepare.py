from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

from pgci_preprocess.annotation_convert import (
    ConversionError,
    convert_annotation,
)
from pgci_preprocess.conversion_common import default_output_prefix
from pgci_preprocess.input_router import (
    InputRoutingError,
    route_annotation,
    select_role_source,
)
from pgci_preprocess.project_paths import downloads_root
from pgci_preprocess.sequence_extract import (
    SequenceExtractionError,
    extract_cds_pep_from_gff,
)
from pgci_preprocess.species_names import (
    legacy_species_code,
    load_species_names,
    resolve_species_code,
)


ROLE_FIELDS = ("genome", "gff", "cds", "pep")
SOURCE_COLUMNS = {
    "genome": "current_genome",
    "gff": "current_gff_or_gtf",
    "cds": "current_cds",
    "pep": "current_pep",
}
INVENTORY_FIELDS = [
    "showid",
    "species_code",
    "Order",
    "Family",
    "Species",
    "gff",
    "genome",
    "cds",
    "pep",
    "rna",
    "gbff",
    "gtf",
]
ASSEMBLY_MAP_FIELDS = [
    "show_id",
    "assembly_id",
    "source",
    "evidence_path",
    "evidence_sha256",
]
STANDARD_ASSEMBLY_RE = re.compile(r"^(?:GC[AF]_\d+\.\d+|GWH[A-Z0-9]+\.\d+)$")
FIGSHARE_ID_RE = re.compile(r"(?:figshare[./])(\d+)", re.IGNORECASE)


class BatchPreparationError(RuntimeError):
    pass


@dataclass(frozen=True)
class PreparedFile:
    show_id: str
    species: str
    bundle: str
    role: str
    source_path: Path
    target_path: Path
    transfer: str
    bytes: int
    sha256: str
    compound_suffix: str
    action: str = "mirror"
    content_format: str = ""
    derived_from: str = ""
    candidate_count: int = 1
    route_note: str = ""


@dataclass(frozen=True)
class BatchPreparationResult:
    candidate_ids: tuple[str, ...]
    prepared_files: tuple[PreparedFile, ...]
    manifest_path: Path | None
    suffix_summary_path: Path | None
    routing_path: Path | None = None
    converted: int = 0
    extracted: int = 0


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise BatchPreparationError(f"TSV has no header: {path}")
        return list(reader.fieldnames), list(reader)


def _write_tsv(
    path: Path,
    fields: Iterable[str],
    rows: Iterable[Mapping[str, Any]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def compound_suffix(filename: str | Path) -> str:
    """Return all trailing suffixes, e.g. .cds.fa.gz rather than only .gz."""
    return "".join(Path(filename).suffixes).lower()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_source(folder: Path, value: str, role: str, show_id: str) -> Path:
    if not value.strip():
        raise BatchPreparationError(f"s{show_id}: missing {role} path")
    try:
        source, _candidates, _detail = select_role_source(
            folder,
            value,
            role,
            show_id,
        )
    except InputRoutingError as exc:
        raise BatchPreparationError(str(exc)) from exc
    return source


def _target_for(
    project_root: Path,
    species_directory: str,
    source: Path,
    source_folder: Path,
) -> tuple[Path, str]:
    relative = source.relative_to(source_folder).as_posix()
    downloads = downloads_root(project_root).resolve()
    target = (downloads / species_directory / relative).resolve()
    try:
        target.relative_to(downloads)
    except ValueError as exc:
        raise BatchPreparationError(f"target escapes downloads directory: {target}") from exc
    return target, relative


def _mirror(source: Path, target: Path, dry_run: bool) -> tuple[str, str]:
    source_hash = _sha256(source)
    if target.exists():
        if not target.is_file() or target.stat().st_size != source.stat().st_size:
            raise BatchPreparationError(f"target already exists with different size: {target}")
        if _sha256(target) != source_hash:
            raise BatchPreparationError(f"target already exists with different hash: {target}")
        return "existing", source_hash
    if dry_run:
        return "dry-run", source_hash
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(source, target)
        transfer = "hardlink"
    except OSError:
        shutil.copy2(source, target)
        transfer = "copy"
    return transfer, _sha256(target)


def _relative_in_species(species_dir: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(species_dir.resolve()).as_posix()
    except ValueError as exc:
        raise BatchPreparationError(
            f"derived output escapes species directory: {path}"
        ) from exc


def _conversion_names(fmt: str, prefix: str, include_genome: bool) -> dict[str, str]:
    if fmt == "gtf":
        return {"gff": f"{prefix}.gff3.gz"}
    if fmt == "gbff":
        names = {
            "gff": f"{prefix}.gff3.gz",
            "cds": f"{prefix}.cds.fa.gz",
            "pep": f"{prefix}.pep.fa.gz",
        }
        if include_genome:
            names["genome"] = f"{prefix}.genome.fa.gz"
        return names
    raise BatchPreparationError(f"unsupported annotation conversion route: {fmt}")


def _reuse_conversion(
    source: Path,
    output_dir: Path,
    prefix: str,
    fmt: str,
    *,
    include_genome: bool,
) -> tuple[dict[str, Path], dict[str, int]]:
    names = _conversion_names(fmt, prefix, include_genome)
    manifest_path = output_dir / f"{prefix}.conversion.json"
    expected_paths = {role: output_dir / name for role, name in names.items()}
    existing = [path for path in expected_paths.values() if path.exists()]
    if not manifest_path.is_file():
        if existing:
            raise BatchPreparationError(
                "derived conversion outputs exist without a verified manifest: "
                + ", ".join(str(path) for path in existing)
            )
        result = convert_annotation(
            source,
            output_dir,
            fmt=fmt,
            prefix=prefix,
            include_genome=include_genome,
        )
        return dict(result.outputs), dict(result.counts)

    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        raise BatchPreparationError(f"invalid conversion manifest: {manifest_path}") from exc
    if manifest.get("format") != fmt or manifest.get("prefix") != prefix:
        raise BatchPreparationError(
            f"conversion manifest does not match requested route: {manifest_path}"
        )
    if manifest.get("complete") is not True:
        raise BatchPreparationError(
            f"incomplete conversion output cannot be registered: {manifest_path}"
        )
    source_record = manifest.get("source") or {}
    if (
        Path(str(source_record.get("path", ""))).resolve() != source.resolve()
        or int(source_record.get("bytes", -1)) != source.stat().st_size
        or str(source_record.get("sha256", "")).lower() != _sha256(source)
    ):
        raise BatchPreparationError(
            f"conversion manifest source changed; remove/reconvert explicitly: {manifest_path}"
        )
    output_records = manifest.get("outputs") or {}
    for role, path in expected_paths.items():
        record = output_records.get(role) or {}
        if not path.is_file() or path.stat().st_size == 0:
            raise BatchPreparationError(f"conversion output missing: {path}")
        if (
            int(record.get("bytes", -1)) != path.stat().st_size
            or str(record.get("sha256", "")).lower() != _sha256(path)
        ):
            raise BatchPreparationError(
                f"conversion output does not match its manifest: {path}"
            )
    return expected_paths, dict(manifest.get("counts") or {})


def _reuse_or_extract(
    genome: Path,
    gff: Path,
    cds_path: Path,
    pep_path: Path,
    manifest_path: Path,
) -> dict[str, int]:
    if cds_path.exists() or pep_path.exists() or manifest_path.exists():
        if not manifest_path.is_file():
            raise BatchPreparationError(
                "derived CDS/PEP exists without an extraction manifest: "
                f"{cds_path}, {pep_path}"
            )
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            raise BatchPreparationError(f"invalid extraction manifest: {manifest_path}") from exc
        if not cds_path.is_file() or not pep_path.is_file():
            raise BatchPreparationError(
                f"incomplete genome+GFF extraction bundle: {manifest_path}"
            )
        for role, source in (("genome", genome), ("gff", gff)):
            record = (manifest.get("sources") or {}).get(role) or {}
            if (
                int(record.get("bytes", -1)) != source.stat().st_size
                or str(record.get("sha256", "")).lower() != _sha256(source)
            ):
                raise BatchPreparationError(
                    f"extraction source changed; do not silently reuse {manifest_path}"
                )
        for role, path in (("cds", cds_path), ("pep", pep_path)):
            record = (manifest.get("outputs") or {}).get(role) or {}
            if (
                int(record.get("bytes", -1)) != path.stat().st_size
                or str(record.get("sha256", "")).lower() != _sha256(path)
            ):
                raise BatchPreparationError(
                    f"extraction output does not match its manifest: {path}"
                )
        return dict(manifest.get("counts") or {})

    try:
        extracted = extract_cds_pep_from_gff(genome, gff, cds_path, pep_path)
    except SequenceExtractionError as exc:
        raise BatchPreparationError(str(exc)) from exc
    manifest = {
        "extractor": "pgci_preprocess.sequence_extract",
        "sources": {
            "genome": {
                "path": str(genome.resolve()),
                "bytes": genome.stat().st_size,
                "sha256": _sha256(genome),
            },
            "gff": {
                "path": str(gff.resolve()),
                "bytes": gff.stat().st_size,
                "sha256": _sha256(gff),
            },
        },
        "outputs": {
            role: {
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for role, path in (("cds", cds_path), ("pep", pep_path))
        },
        "counts": extracted.counts,
        "warnings": list(extracted.warnings),
    }
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return dict(extracted.counts)


def _prepared_file(
    *,
    show_id: str,
    species: str,
    bundle: str,
    role: str,
    source: Path,
    target: Path,
    transfer: str,
    digest: str,
    action: str,
    content_format: str = "",
    derived_from: str = "",
    candidate_count: int = 1,
    route_note: str = "",
) -> PreparedFile:
    return PreparedFile(
        show_id=show_id,
        species=species,
        bundle=bundle,
        role=role,
        source_path=source,
        target_path=target,
        transfer=transfer,
        bytes=target.stat().st_size if target.exists() else source.stat().st_size,
        sha256=digest,
        compound_suffix=compound_suffix(target.name),
        action=action,
        content_format=content_format,
        derived_from=derived_from,
        candidate_count=candidate_count,
        route_note=route_note,
    )


def _discover_processed(project_root: Path) -> set[str]:
    analysis = project_root / "analysis_v2"
    if not analysis.is_dir():
        return set()
    return {
        item.name
        for item in analysis.iterdir()
        if item.is_dir() and not item.name.startswith(".")
    }


def _infer_assembly_id(row: dict[str, str], show_id: str) -> str:
    accession = str(row.get("candidate_accession", "")).strip()
    if STANDARD_ASSEMBLY_RE.fullmatch(accession):
        return accession
    if accession and not accession.isdigit():
        match = FIGSHARE_ID_RE.search(accession)
        if match:
            return f"Figshare_{match.group(1)}_s{show_id}"
        return accession

    public_ids = (
        accession,
        str(row.get("candidate_version", "")).strip(),
        str(row.get("download_url", "")).strip(),
    )
    for value in public_ids:
        match = FIGSHARE_ID_RE.search(value)
        if match:
            return f"Figshare_{match.group(1)}_s{show_id}"
        if value.isdigit():
            return f"Figshare_{value}_s{show_id}"
    raise BatchPreparationError(
        f"s{show_id}: no evidence-backed assembly identifier in candidate metadata"
    )


def prepare_batch(
    project_root: Path | str,
    source_status: Path | str,
    *,
    expected_count: int | None = None,
    dry_run: bool = False,
    batch_tag: str = "current",
    species_names_path: Path | str | None = None,
) -> BatchPreparationResult:
    project_root = Path(project_root).resolve()
    source_status = Path(source_status).resolve()
    if not project_root.is_dir():
        raise BatchPreparationError(f"project root not found: {project_root}")
    if not source_status.is_file():
        raise BatchPreparationError(f"source status table not found: {source_status}")

    fields, status_rows = _read_tsv(source_status)
    required = {"show_id", "Order", "Family", "Species", "status", "folder", *SOURCE_COLUMNS.values()}
    missing_columns = sorted(required - set(fields))
    if missing_columns:
        raise BatchPreparationError(f"status table missing columns: {', '.join(missing_columns)}")

    processed = _discover_processed(project_root)
    names_catalog = load_species_names(project_root, explicit_path=species_names_path)
    taken_codes = set(names_catalog.by_code) if names_catalog is not None else set()
    candidates: list[dict[str, str]] = []
    for row in status_rows:
        show_id = str(row.get("show_id", "")).strip()
        if row.get("status") != "complete" or not show_id:
            continue
        species_code = resolve_species_code(
            project_root,
            show_id,
            latin_name=str(row.get("Species", "")),
            explicit_path=species_names_path,
            occupied=taken_codes,
        )
        taken_codes.add(species_code.casefold())
        if species_code in processed or legacy_species_code(show_id) in processed:
            continue
        candidates.append(row)
    candidates.sort(key=lambda row: int(row["show_id"]))
    if expected_count is not None and len(candidates) != expected_count:
        raise BatchPreparationError(
            f"expected {expected_count} candidates, found {len(candidates)}: "
            + ",".join(
                f"{resolve_species_code(project_root, row['show_id'])}"
                f"(show_id={row['show_id']})"
                for row in candidates
            )
        )
    if not candidates:
        raise BatchPreparationError("no unprocessed complete species found")

    inventory_path = project_root / "resources" / "download_inventory.tsv"
    inventory_fields, inventory_rows = (
        _read_tsv(inventory_path) if inventory_path.is_file() else (INVENTORY_FIELDS, [])
    )
    if "species_code" not in inventory_fields:
        insert_at = inventory_fields.index("showid") + 1 if "showid" in inventory_fields else 0
        inventory_fields.insert(insert_at, "species_code")
    inventory_by_id = {row.get("showid", ""): row for row in inventory_rows}
    for inventory_row in inventory_rows:
        inventory_show_id = str(inventory_row.get("showid", "")).strip()
        if not inventory_show_id:
            continue
        existing_code = str(inventory_row.get("species_code", "")).strip()
        if existing_code:
            taken_codes.add(existing_code.casefold())
            continue
        code = resolve_species_code(
            project_root,
            inventory_show_id,
            latin_name=str(inventory_row.get("Species", "")),
            explicit_path=species_names_path,
            occupied=taken_codes,
        )
        inventory_row["species_code"] = code
        taken_codes.add(code.casefold())
    assembly_map_path = project_root / "resources" / "assembly_accessions.tsv"
    assembly_fields, assembly_rows = (
        _read_tsv(assembly_map_path)
        if assembly_map_path.is_file()
        else (ASSEMBLY_MAP_FIELDS, [])
    )
    assembly_by_id = {row.get("show_id", ""): row for row in assembly_rows}

    prepared: list[PreparedFile] = []
    routing_rows: list[dict[str, object]] = []
    converted_count = 0
    extracted_count = 0
    downloads_root_path = downloads_root(project_root).resolve()
    downloads_root_path.mkdir(parents=True, exist_ok=True)
    for row in candidates:
        show_id = str(row["show_id"]).strip()
        species = str(row["Species"]).strip()
        species_code = resolve_species_code(
            project_root,
            show_id,
            latin_name=species,
            explicit_path=species_names_path,
            occupied=taken_codes,
        )
        taken_codes.add(species_code.casefold())
        source_folder = Path(row["folder"]).resolve()
        if not source_folder.is_dir():
            raise BatchPreparationError(f"s{show_id}: species folder missing: {source_folder}")

        inventory = inventory_by_id.get(show_id)
        if inventory is None:
            inventory = {
                "showid": show_id,
                "species_code": species_code,
                "Order": row["Order"],
                "Family": row["Family"],
                "Species": species,
                "gff": "",
                "genome": "",
                "cds": "",
                "pep": "",
                "rna": "",
                "gbff": "",
                "gtf": "",
            }
            inventory_by_id[show_id] = inventory
            inventory_rows.append(inventory)
        else:
            inventory["species_code"] = species_code

        selected_sources: dict[str, Path] = {}
        selected_candidates: dict[str, int] = {}
        roles: dict[str, str] = {}
        bundle = str(row.get("current_bundle", ""))
        for role in ROLE_FIELDS:
            raw_value = str(row.get(SOURCE_COLUMNS[role], "")).strip()
            if not raw_value:
                if role in {"genome", "cds", "pep"}:
                    continue
                for fallback in (
                    f"current_{role}",
                    role,
                    f"{role}_file",
                ):
                    raw_value = str(row.get(fallback, "")).strip()
                    if raw_value:
                        break
            if not raw_value:
                continue
            try:
                source, candidate_count, detail = select_role_source(
                    source_folder,
                    raw_value,
                    role,
                    show_id,
                )
            except InputRoutingError as exc:
                raise BatchPreparationError(str(exc)) from exc
            target, relative = _target_for(
                project_root,
                source_folder.name,
                source,
                source_folder,
            )
            transfer, digest = _mirror(source, target, dry_run)
            selected_sources[role] = source
            selected_candidates[role] = candidate_count
            roles[role] = relative
            prepared.append(
                _prepared_file(
                    show_id=show_id,
                    species=species,
                    bundle=bundle,
                    role=role,
                    source=source,
                    target=target,
                    transfer=transfer,
                    digest=digest,
                    action="mirror",
                    content_format=detail,
                    candidate_count=candidate_count,
                )
            )

        direct_roles = dict(roles)
        if "gff" not in selected_sources:
            raise BatchPreparationError(
                f"s{show_id}: no annotation selected in current_gff_or_gtf "
                "(or a gtf/gbff fallback column)"
            )

        annotation_source = selected_sources["gff"]
        try:
            annotation_route = route_annotation(annotation_source)
        except InputRoutingError as exc:
            raise BatchPreparationError(f"s{show_id}: {exc}") from exc
        runnable_gff_source = annotation_source
        species_dir = downloads_root_path / source_folder.name
        derived_dir = species_dir / "derived"

        if annotation_route.route == "direct_gff":
            pass
        elif annotation_route.route in {"convert_gtf", "convert_gbff"}:
            fmt = "gtf" if annotation_route.route == "convert_gtf" else "gbff"
            prefix = default_output_prefix(annotation_source.name) + f".from_{fmt}"
            include_genome = fmt == "gbff" and "genome" not in selected_sources
            names = _conversion_names(fmt, prefix, include_genome)
            if dry_run:
                conversion_outputs = {role: derived_dir / name for role, name in names.items()}
                conversion_counts: dict[str, int] = {}
            else:
                try:
                    conversion_outputs, conversion_counts = _reuse_conversion(
                        annotation_source,
                        derived_dir,
                        prefix,
                        fmt,
                        include_genome=include_genome,
                    )
                except (BatchPreparationError, ConversionError) as exc:
                    raise BatchPreparationError(
                        f"s{show_id}: {fmt.upper()} conversion failed: {exc}"
                    ) from exc
                converted_count += 1
            for role, output in conversion_outputs.items():
                relative = _relative_in_species(species_dir, output)
                roles[role] = relative
                if not dry_run:
                    prepared.append(
                        _prepared_file(
                            show_id=show_id,
                            species=species,
                            bundle=bundle,
                            role=role,
                            source=annotation_source,
                            target=output,
                            transfer=f"converted_{fmt}",
                            digest=_sha256(output),
                            action=f"convert_{fmt}",
                            content_format=annotation_route.content_format,
                            derived_from=selected_sources["gff"].name,
                            candidate_count=selected_candidates["gff"],
                            route_note=annotation_route.note,
                        )
                    )
                elif role == "gff":
                    prepared.append(
                        _prepared_file(
                            show_id=show_id,
                            species=species,
                            bundle=bundle,
                            role=role,
                            source=annotation_source,
                            target=output,
                            transfer=f"dry_run_convert_{fmt}",
                            digest="",
                            action=f"convert_{fmt}",
                            content_format=annotation_route.content_format,
                            derived_from=selected_sources["gff"].name,
                            candidate_count=selected_candidates["gff"],
                            route_note=annotation_route.note,
                        )
                    )
            runnable_gff_source = conversion_outputs["gff"]
            inventory["gff"] = roles["gff"]
            inventory[fmt] = direct_roles.get("gff", "")
            if fmt == "gbff":
                for role in ("cds", "pep", "genome"):
                    if role in conversion_outputs:
                        inventory[role] = roles[role]
        else:
            raise BatchPreparationError(
                f"s{show_id}: {annotation_route.content_format} is metadata-only and "
                "cannot be registered as WGDI coordinates"
            )

        # A GFF with complete CDS rows plus a genome is enough to derive the
        # missing FASTA pair.  If only one FASTA role was supplied, derive both
        # so the registered bundle is guaranteed to share one protein-key set.
        if annotation_route.route != "convert_gbff" and (
            "cds" not in roles or "pep" not in roles
        ):
            if "genome" not in selected_sources:
                raise BatchPreparationError(
                    f"s{show_id}: CDS/PEP are missing and no genome FASTA is "
                    "available for genome+GFF extraction"
                )
            extraction_prefix = default_output_prefix(runnable_gff_source.name) + ".derived"
            cds_target = derived_dir / f"{extraction_prefix}.cds.fa"
            pep_target = derived_dir / f"{extraction_prefix}.pep.fa"
            extraction_manifest = derived_dir / f"{extraction_prefix}.extraction.json"
            if not dry_run:
                try:
                    _reuse_or_extract(
                        selected_sources["genome"],
                        runnable_gff_source,
                        cds_target,
                        pep_target,
                        extraction_manifest,
                    )
                except BatchPreparationError as exc:
                    raise BatchPreparationError(
                        f"s{show_id}: genome+GFF extraction failed: {exc}"
                    ) from exc
                extracted_count += 1
            for role, output in (("cds", cds_target), ("pep", pep_target)):
                roles[role] = _relative_in_species(species_dir, output)
                if dry_run:
                    prepared.append(
                        _prepared_file(
                            show_id=show_id,
                            species=species,
                            bundle=bundle,
                            role=role,
                            source=selected_sources["gff"],
                            target=output,
                            transfer="dry_run_extract",
                            digest="",
                            action="extract_genome_gff",
                            content_format="genome_gff",
                            derived_from=selected_sources["gff"].name,
                        )
                    )
                else:
                    prepared.append(
                        _prepared_file(
                            show_id=show_id,
                            species=species,
                            bundle=bundle,
                            role=role,
                            source=selected_sources["gff"],
                            target=output,
                            transfer="extracted",
                            digest=_sha256(output),
                            action="extract_genome_gff",
                            content_format="genome_gff",
                            derived_from=selected_sources["gff"].name,
                        )
                    )

        for role in ("gff", "genome", "cds", "pep"):
            if role in roles:
                inventory[role] = roles[role]
        routing_rows.append(
            {
                "show_id": show_id,
                "species_code": species_code,
                "species": species,
                "source_folder": str(source_folder),
                "annotation_source": str(annotation_source),
                "suffix_hint": annotation_route.suffix_hint,
                "content_format": annotation_route.content_format,
                "route": annotation_route.route,
                "registered_gff": inventory.get("gff", ""),
                "registered_cds": inventory.get("cds", ""),
                "registered_pep": inventory.get("pep", ""),
                "registered_genome": inventory.get("genome", ""),
                "selected_candidates": selected_candidates.get("gff", 1),
                "note": annotation_route.note,
            }
        )

        assembly_id = _infer_assembly_id(row, show_id)
        gff_path = annotation_source
        existing_assembly = assembly_by_id.get(show_id)
        if existing_assembly and existing_assembly.get("assembly_id") != assembly_id:
            raise BatchPreparationError(
                f"s{show_id}: assembly map conflict "
                f"{existing_assembly.get('assembly_id')} != {assembly_id}"
            )
        if existing_assembly is None:
            existing_assembly = {
                "show_id": show_id,
                "assembly_id": assembly_id,
                "source": str(row.get("candidate_source", "")).strip(),
                "evidence_path": str(gff_path),
                "evidence_sha256": _sha256(gff_path),
            }
            assembly_by_id[show_id] = existing_assembly
            assembly_rows.append(existing_assembly)

    if not dry_run:
        def numeric(row: Mapping[str, Any]):
            value = str(row.get("showid", ""))
            return (0, int(value)) if value.isdigit() else (1, value)

        inventory_rows.sort(key=numeric)
        _write_tsv(inventory_path, inventory_fields, inventory_rows)
        assembly_rows.sort(
            key=lambda row: (0, int(row["show_id"])) if str(row.get("show_id", "")).isdigit() else (1, str(row.get("show_id", "")))
        )
        _write_tsv(assembly_map_path, assembly_fields, assembly_rows)

    manifest_path: Path | None = None
    suffix_summary_path: Path | None = None
    routing_path: Path | None = None
    if not dry_run:
        batch_dir = project_root / "resources" / "preprocess_batches"
        manifest_path = batch_dir / f"{batch_tag}_manifest.tsv"
        suffix_summary_path = batch_dir / f"{batch_tag}_suffix_summary.tsv"
        routing_path = batch_dir / f"{batch_tag}_routing.tsv"
        manifest_fields = [
            "show_id",
            "species_code",
            "species",
            "bundle",
            "role",
            "source_path",
            "target_path",
            "transfer",
            "bytes",
            "sha256",
            "compound_suffix",
            "action",
            "content_format",
            "derived_from",
            "candidate_count",
            "route_note",
        ]
        _write_tsv(
            manifest_path,
            manifest_fields,
            (
                {
                    "show_id": item.show_id,
                    "species_code": str(
                        inventory_by_id.get(item.show_id, {}).get(
                            "species_code",
                            "",
                        )
                    ),
                    "species": item.species,
                    "bundle": item.bundle,
                    "role": item.role,
                    "source_path": str(item.source_path),
                    "target_path": str(item.target_path),
                    "transfer": item.transfer,
                    "bytes": item.bytes,
                    "sha256": item.sha256,
                    "compound_suffix": item.compound_suffix,
                    "action": item.action,
                    "content_format": item.content_format,
                    "derived_from": item.derived_from,
                    "candidate_count": item.candidate_count,
                    "route_note": item.route_note,
                }
                for item in prepared
            ),
        )
        _write_tsv(
            routing_path,
            [
                "show_id",
                "species_code",
                "species",
                "source_folder",
                "annotation_source",
                "suffix_hint",
                "content_format",
                "route",
                "registered_gff",
                "registered_cds",
                "registered_pep",
                "registered_genome",
                "selected_candidates",
                "note",
            ],
            routing_rows,
        )
        suffix_counts = Counter((item.role, item.compound_suffix) for item in prepared)
        _write_tsv(
            suffix_summary_path,
            ["role", "compound_suffix", "files"],
            (
                {"role": role, "compound_suffix": suffix, "files": count}
                for (role, suffix), count in sorted(suffix_counts.items())
            ),
        )

    return BatchPreparationResult(
        candidate_ids=tuple(str(row["show_id"]) for row in candidates),
        prepared_files=tuple(prepared),
        manifest_path=manifest_path,
        suffix_summary_path=suffix_summary_path,
        routing_path=routing_path,
        converted=converted_count,
        extracted=extracted_count,
    )
