from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import uuid
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

from pgci_preprocess.artifact_naming import (
    ANALYSIS_REQUIRED_SUFFIXES,
    CLEAN_REQUIRED_SUFFIXES,
    bundle_name_problems,
)
from pgci_preprocess.artifact_migration import migrate_artifacts
from pgci_preprocess.annotation_convert import ConversionError, detect_conversion_format
from pgci_preprocess.fasta import (
    normalize_evm_protein_sequences,
    open_text,
    read_cds_fasta,
    read_pep_fasta,
)
from pgci_preprocess.gff_formats import (
    GffFormatError,
    apply_gwh_accessions,
    backfill_regions_from_transcripts,
    mrna_accession_map,
    normalize_cds_protein_ids,
    parse_with_adapter,
)
from pgci_preprocess.input_router import InputRoutingError, select_role_source
from pgci_preprocess.models import Anomaly, Candidate, Region
from pgci_preprocess.reports import (
    finalize_project,
    read_tsv as read_report_tsv,
)
from pgci_preprocess.selection import (
    choose_representative,
    region_sort_key,
    select_analysis_sequences,
)
from pgci_preprocess.sequence import validate_cds_pep
from pgci_preprocess.species_names import (
    legacy_species_code,
    load_species_names,
    resolve_show_ids,
)
from pgci_preprocess.version import __version__

ACCESSION_RE = re.compile(
    r"^((?:GC[AF]_\d+\.\d+|GWH[A-Z0-9]+\.\d+))(?:[_\.])"
)
GENOME_LENGTH_RE = re.compile(r"(?:^|\s)Len=(\d+)", re.IGNORECASE)
CHROMOSOME_LABEL_RE = re.compile(r"(?:^|\s)Chromosome\s+([^\s]+)", re.IGNORECASE)
# Derived from the package version so that every ``*.qc.json`` records the same
# release the CLI reports.  See ``version.py``.
PIPELINE_VERSION = __version__


@dataclass(slots=True)
class SpeciesInput:
    show_id: str
    code: str
    species: str
    order: str
    family: str
    directory: Path
    gff_name: str
    cds_name: str
    pep_name: str
    genome_name: str
    accession: str
    assembly_id: str = ""

    @property
    def legacy_code(self) -> str:
        return legacy_species_code(self.show_id)


@dataclass(slots=True)
class PipelineOptions:
    project_root: Path
    clean_dirname: str = "clean_v2"
    analysis_dirname: str = "analysis_v2"
    min_genes_per_sequence: int = 50
    target_gene_coverage: float = 0.90
    max_sequences: int = 200
    # Default must match the CLI: never replace an existing validated bundle
    # unless the caller opts in explicitly with overwrite=True.
    overwrite: bool = False
    gff_format: str = "auto"
    protein_normalization: str = "auto"
    accession_policy: str = "strict"
    backfill_regions_from_transcripts: bool = True
    format_exceptions: dict[str, dict[str, bool]] | None = None
    species_names_path: Path | None = None
    # Legacy ``s<show_id>`` bundles are renamed to the canonical abbreviation
    # before they are reused, so preprocessing never reports name-compliant
    # output while the content still carries legacy gene IDs.
    migrate_legacy_artifacts: bool = True


class PipelineFailure(RuntimeError):
    def __init__(self, reason: str, detail: str = "") -> None:
        super().__init__(detail or reason)
        self.reason = reason
        self.detail = detail


_SHOW_ID_RE = re.compile(r"^(\d+)_")


def _accession_from_filename(name: str) -> str:
    match = ACCESSION_RE.match(name)
    if match:
        return match.group(1)
    # 无标准 accession 的物种（如 JGI 的 3282_/4193_）：用文件名前缀的 show-id 兜底
    sid = _SHOW_ID_RE.match(name)
    return f"JGI_{sid.group(1)}" if sid else ""


ASSEMBLY_MAP_RELATIVE = Path("resources") / "assembly_accessions.tsv"
FORMAT_EXCEPTIONS_RELATIVE = Path("resources") / "format_exceptions.tsv"


def load_format_exceptions(project_root: Path | str) -> dict[str, dict[str, bool]]:
    """Load evidence-backed per-species format switches.

    The canonical package contains no species IDs.  Structural exceptions that
    cannot be inferred safely (for example NCBI segmented gene rows that must be
    merged) are stored in ``resources/format_exceptions.tsv`` with a reason
    column, and recorded in each species' QC record.
    """

    path = Path(project_root) / FORMAT_EXCEPTIONS_RELATIVE
    mapping: dict[str, dict[str, bool]] = {}
    if not path.is_file():
        return mapping
    for row in read_report_tsv(path):
        show_id = str(row.get("show_id") or "").strip()
        if not show_id:
            continue
        mapping[show_id] = {
            "merge_split_genes": _truthy(row.get("merge_split_genes")),
            "gene_span_only": _truthy(row.get("gene_span_only")),
            "gene_cds_direct": _truthy(row.get("gene_cds_direct")),
        }
    return mapping


def _truthy(value: object) -> bool:
    return str(value or "").strip().casefold() in {"1", "true", "yes", "y"}


def load_assembly_map(project_root: Path | str) -> dict[str, str]:
    """Load official assembly identity overrides.

    Returns {show_id: assembly_id}. The table is an evidence-backed provenance
    record (official GFF header / source manifest), not guessed accessions.
    Missing or empty table is allowed: accession validation falls back to the
    filename rule, preserving fail-closed behaviour.
    """
    project_root = Path(project_root)
    path = project_root / ASSEMBLY_MAP_RELATIVE
    mapping: dict[str, str] = {}
    if not path.is_file():
        return mapping
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            show_id = (row.get("show_id") or "").strip()
            assembly_id = (row.get("assembly_id") or "").strip()
            if show_id and assembly_id:
                mapping[show_id] = assembly_id
    return mapping


def discover_species(
    project_root: Path | str,
    *,
    species_names_path: Path | str | None = None,
) -> list[SpeciesInput]:
    project_root = Path(project_root)
    assembly_map = load_assembly_map(project_root)
    names = load_species_names(
        project_root,
        explicit_path=species_names_path,
    )
    inventory = project_root / "resources" / "download_inventory.tsv"
    if not inventory.is_file():
        raise FileNotFoundError(f"download inventory not found: {inventory}")

    downloads = project_root / "data" / "downloads"
    directories = {
        path.name.split("_", 1)[0]: path
        for path in downloads.iterdir()
        if path.is_dir() and "_" in path.name
    }
    species: list[SpeciesInput] = []
    taken_codes = set(names.by_code) if names is not None else set()
    with inventory.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            show_id = row.get("showid", "").strip()
            if not show_id or show_id not in directories:
                continue
            selected_names: dict[str, str] = {}
            for role in ("gff", "genome", "cds", "pep"):
                value = row.get(role, "").strip()
                if not value:
                    selected_names[role] = ""
                    continue
                if "|" not in value:
                    selected_names[role] = value
                    continue
                try:
                    selected, _count, _detail = select_role_source(
                        directories[show_id],
                        value,
                        role,
                        show_id,
                    )
                except InputRoutingError as exc:
                    raise FileNotFoundError(str(exc)) from exc
                selected_names[role] = selected.resolve().relative_to(
                    directories[show_id].resolve()
                ).as_posix()
            files = [
                selected_names["gff"],
                selected_names["cds"],
                selected_names["pep"],
            ]
            accessions = {_accession_from_filename(name) for name in files if name}
            accessions.discard("")
            accession = next(iter(accessions)) if len(accessions) == 1 else "NA"
            species_name = row.get("Species", "").strip()
            code = row.get("species_code", "").strip()
            if not code:
                if names is not None:
                    from pgci_preprocess.species_codegen import derive_code

                    if show_id in names.by_show_id:
                        code = names.code_for(show_id, latin_name=species_name)
                    elif species_name:
                        code = derive_code(species_name, taken_codes)
                    else:
                        code = f"s{show_id}"
                else:
                    code = f"s{show_id}"
            taken_codes.add(code.casefold())
            species.append(
                SpeciesInput(
                    show_id=show_id,
                    code=code,
                    species=species_name,
                    order=row.get("Order", "").strip(),
                    family=row.get("Family", "").strip(),
                    directory=directories[show_id],
                    gff_name=selected_names["gff"],
                    cds_name=selected_names["cds"],
                    pep_name=selected_names["pep"],
                    genome_name=selected_names["genome"],
                    accession=accession,
                    assembly_id=assembly_map.get(show_id, ""),
                )
            )
    species.sort(key=lambda item: int(item.show_id) if item.show_id.isdigit() else 10**12)
    return species


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _remove_staging(path: Path, parent: Path) -> None:
    path = path.resolve()
    parent = parent.resolve()
    if path.parent != parent or not path.name.startswith("."):
        raise ValueError(f"refusing to remove unexpected staging path: {path}")
    if path.exists():
        shutil.rmtree(path)


def _replace_directory(staging: Path, destination: Path, overwrite: bool) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not overwrite:
        raise FileExistsError(f"output already exists: {destination}")
    backup: Path | None = None
    if destination.exists():
        backup = destination.with_name(f".{destination.name}.previous-{uuid.uuid4().hex}")
        os.replace(destination, backup)
    try:
        os.replace(staging, destination)
    except Exception:
        if backup is not None and backup.exists() and not destination.exists():
            os.replace(backup, destination)
        raise
    if backup is not None and backup.exists():
        shutil.rmtree(backup)


def _write_fasta(path: Path, records: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            canonical_id = str(record["canonical_gene_id"])
            source_gene = str(record["source_gene_id"])
            source_protein = str(record["source_protein_id"])
            handle.write(f">{canonical_id} source_gene={source_gene} source_protein={source_protein}\n")
            sequence = str(record["sequence"])
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index : index + 80] + "\n")


def _write_gff(path: Path, entries: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for entry in entries:
            handle.write(
                "\t".join(
                    [
                        str(entry["clean_sequence_id"]),
                        str(entry["canonical_gene_id"]),
                        str(entry["start"]),
                        str(entry["end"]),
                        str(entry["strand"]),
                        str(entry["order"]),
                        str(entry["source_gene_id"]),
                    ]
                )
                + "\n"
            )


def _write_lens(path: Path, region_rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in region_rows:
            handle.write(
                f"{row['clean_sequence_id']}\t{row['length']}\t{row['gene_count']}\n"
            )


def _read_fasta_ids(path: Path) -> list[str]:
    with path.open(encoding="utf-8") as handle:
        return [
            line[1:].split(None, 1)[0]
            for line in handle
            if line.startswith(">") and len(line) > 1
        ]


def _read_lens(path: Path) -> dict[str, tuple[int, int]]:
    result: dict[str, tuple[int, int]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 3:
                raise PipelineFailure("invalid_lens", f"{path.name}:{line_number}")
            sequence_id, length, gene_count = fields
            result[sequence_id] = (int(length), int(gene_count))
    return result


def _validate_bundle(
    gff_path: Path,
    lens_path: Path,
    cds_path: Path,
    pep_path: Path,
) -> dict[str, object]:
    pep_ids = _read_fasta_ids(pep_path)
    cds_ids = _read_fasta_ids(cds_path)
    if len(pep_ids) != len(set(pep_ids)):
        raise PipelineFailure("duplicate_pep_id", str(pep_path))
    if len(cds_ids) != len(set(cds_ids)):
        raise PipelineFailure("duplicate_cds_id", str(cds_path))

    lens = _read_lens(lens_path)
    gff_ids: list[str] = []
    orders: dict[str, list[int]] = defaultdict(list)
    per_sequence_counts: Counter[str] = Counter()
    with gff_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 7:
                raise PipelineFailure("invalid_gff_column_count", f"{gff_path.name}:{line_number}")
            sequence_id, gene_id, start, end, strand, order_text, _source_gene = fields
            if sequence_id not in lens:
                raise PipelineFailure("gff_sequence_missing_from_lens", sequence_id)
            start_int, end_int, order = int(start), int(end), int(order_text)
            if start_int < 1 or end_int < start_int or end_int > lens[sequence_id][0]:
                raise PipelineFailure("invalid_gff_coordinates", gene_id)
            if strand not in {"+", "-"}:
                raise PipelineFailure("invalid_strand", gene_id)
            orders[sequence_id].append(order)
            per_sequence_counts[sequence_id] += 1
            gff_ids.append(gene_id)

    if len(gff_ids) != len(set(gff_ids)):
        raise PipelineFailure("duplicate_gff_id", gff_path.name)
    if set(gff_ids) != set(cds_ids) or set(gff_ids) != set(pep_ids):
        raise PipelineFailure(
            "bundle_id_mismatch",
            f"gff={len(gff_ids)} cds={len(cds_ids)} pep={len(pep_ids)}",
        )
    for sequence_id, sequence_orders in orders.items():
        if sequence_orders != list(range(1, len(sequence_orders) + 1)):
            raise PipelineFailure("non_contiguous_order", sequence_id)
    for sequence_id, (_length, expected_count) in lens.items():
        if per_sequence_counts[sequence_id] != expected_count:
            raise PipelineFailure("lens_gene_count_mismatch", sequence_id)
    return {
        "n_genes": len(gff_ids),
        "n_sequences": len(lens),
    }


def _add_rejection(
    rows: list[dict[str, str]],
    scope: str,
    identifier: str,
    reason: str,
    detail: str = "",
) -> None:
    rows.append(
        {
            "scope": scope,
            "id": identifier,
            "reason": reason,
            "detail": detail,
        }
    )


def _discover_accession(species: SpeciesInput, policy: str = "strict") -> str:
    """Resolve the assembly identity, fail-closed unless legacy policy is set.

    Strict policy (default):
      1) explicit assembly map entry (evidence-backed official label);
      2) a single accession agreed by all three input filenames.

    Legacy policy reproduces the historical conversion script: an assembly map
    entry wins, otherwise the status-table accession is used, and finally
    ``NA``.  The chosen policy is written into QC so relaxed provenance is
    never silent.
    """

    if policy not in {"strict", "legacy"}:
        raise PipelineFailure("invalid_accession_policy", policy)
    filename_matches = {
        match.group(1)
        for match in (
            ACCESSION_RE.match(species.gff_name),
            ACCESSION_RE.match(species.cds_name),
            ACCESSION_RE.match(species.pep_name),
        )
        if match
    }

    mapped = (species.assembly_id or "").strip()
    if mapped:
        if policy == "strict" and filename_matches and filename_matches != {mapped}:
            raise PipelineFailure(
                "accession_mismatch",
                f"map={mapped} conflicts with filenames={sorted(filename_matches)}",
            )
        return mapped

    if policy == "legacy":
        return (species.accession or "").strip() or "NA"

    if len(filename_matches) == 1:
        return next(iter(filename_matches))
    if not filename_matches:
        # 手写注释/GWH 等文件名不含 GCA/GCF accession：没有证据说明冲突，
        # 记为 NA，由 assembly_accessions.tsv 提供权威值。
        return "NA"
    raise PipelineFailure(
        "accession_mismatch",
        f"no assembly map and filenames do not agree: {sorted(filename_matches)}",
    )


def _read_genome_regions(path: Path) -> dict[str, Region]:
    regions: dict[str, Region] = {}
    current_id = ""
    current_length = 0
    current_description = ""
    counted_bases = 0
    source_order = 0

    def flush() -> None:
        nonlocal current_id, current_length, counted_bases, source_order
        if not current_id:
            return
        length = current_length or counted_bases
        if length <= 0 or current_id in regions:
            return
        description = current_description
        lowered = description.casefold()
        if "chloroplast" in lowered:
            genome_type = "chloroplast"
        elif "mitochondrion" in lowered or "mitochondrial" in lowered:
            genome_type = "mitochondrion"
        elif "plastid" in lowered:
            genome_type = "plastid"
        else:
            genome_type = "genomic"
        chromosome_match = CHROMOSOME_LABEL_RE.search(description)
        chromosome = chromosome_match.group(1) if chromosome_match else ""
        is_primary = genome_type == "genomic" and (
            bool(chromosome_match) or "OriSeqID=Chr" in description
        )
        regions[current_id] = Region(
            seqid=current_id,
            length=length,
            genome_type=genome_type,
            chromosome=chromosome,
            taxon_id="",
            source_order=source_order,
            is_primary=is_primary,
        )
        source_order += 1

    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                current_id = header.split(None, 1)[0] if header else ""
                current_description = header
                length_match = GENOME_LENGTH_RE.search(header)
                current_length = int(length_match.group(1)) if length_match else 0
                counted_bases = 0
            elif current_id:
                counted_bases += len("".join(line.split()))
    flush()
    return regions


def _prepare_staging(
    options: PipelineOptions,
    code: str,
) -> tuple[Path, Path, Path, Path]:
    clean_root = options.project_root / options.clean_dirname
    analysis_root = options.project_root / options.analysis_dirname
    clean_root.mkdir(parents=True, exist_ok=True)
    analysis_root.mkdir(parents=True, exist_ok=True)
    suffix = uuid.uuid4().hex
    clean_staging = clean_root / f".{code}.tmp-{suffix}"
    analysis_staging = analysis_root / f".{code}.tmp-{suffix}"
    clean_staging.mkdir()
    analysis_staging.mkdir()
    return clean_root, analysis_root, clean_staging, analysis_staging


def _existing_qc(
    species: SpeciesInput,
    options: PipelineOptions,
) -> tuple[dict[str, object], str] | None:
    """Return the prior QC record when a complete canonical bundle exists.

    Only the canonical abbreviation directory is accepted: a legacy
    ``s<show_id>`` bundle is migrated by ``run_project`` before this check, and
    an un-migratable legacy bundle must be re-processed instead of reported as
    a name-compliant SKIPPED result.
    """
    artifact_code = species.code
    analysis_final = options.project_root / options.analysis_dirname / artifact_code
    clean_final = options.project_root / options.clean_dirname / artifact_code
    qc_path = analysis_final / f"{artifact_code}.qc.json"
    if not qc_path.is_file():
        return None
    try:
        qc = json.loads(qc_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError):
        return None
    if qc.get("status") != "OK" or str(qc.get("show_id", "")) != species.show_id:
        return None
    required = [
        analysis_final / f"{artifact_code}.pep",
        analysis_final / f"{artifact_code}.cds",
        analysis_final / f"{artifact_code}.new.gff",
        analysis_final / f"{artifact_code}.lens",
        clean_final / f"{artifact_code}.gene_map.tsv",
        clean_final / f"{artifact_code}.chr_map.tsv",
    ]
    if any(not path.is_file() or path.stat().st_size == 0 for path in required):
        return None
    analysis_problems = bundle_name_problems(
        analysis_final,
        artifact_code,
        required_suffixes=ANALYSIS_REQUIRED_SUFFIXES,
        kind_map={".new.gff": "gff", ".cds": "fasta", ".pep": "fasta"},
    )
    clean_problems = bundle_name_problems(
        clean_final,
        artifact_code,
        required_suffixes=CLEAN_REQUIRED_SUFFIXES,
    )
    if analysis_problems or clean_problems:
        return None
    return qc, artifact_code


def _existing_outputs_present(species: SpeciesInput, options: PipelineOptions) -> bool:
    return any(
        (
            options.project_root / options.analysis_dirname / code
        ).exists()
        or (options.project_root / options.clean_dirname / code).exists()
        for code in (species.code, species.legacy_code)
    )


def _skipped_result(
    species: SpeciesInput,
    qc: dict[str, object],
    artifact_code: str,
) -> dict[str, object]:
    selection = qc.get("sequence_selection") or {}
    all_bundle = qc.get("all_bundle") or {}
    wgdi_bundle = qc.get("wgdi_bundle") or {}
    selection = selection if isinstance(selection, dict) else {}
    all_bundle = all_bundle if isinstance(all_bundle, dict) else {}
    wgdi_bundle = wgdi_bundle if isinstance(wgdi_bundle, dict) else {}
    nuclear_taxa = qc.get("nuclear_taxa")
    anomalies = qc.get("anomalies")
    nuclear_taxa = nuclear_taxa if isinstance(nuclear_taxa, list) else []
    anomalies = anomalies if isinstance(anomalies, list) else []
    return {
        "show_id": species.show_id,
        "code": species.code,
        "artifact_code": artifact_code,
        "species": species.species,
        "accession": qc.get("accession", species.accession),
        "status": "SKIPPED",
        "assembly_build": qc.get("assembly_build", ""),
        "n_regions": qc.get("n_regions", ""),
        "n_representative_genes": qc.get("n_representative_genes", ""),
        "n_all_genes": all_bundle.get("n_genes", ""),
        "n_wgdi_genes": wgdi_bundle.get("n_genes", ""),
        "n_wgdi_sequences": wgdi_bundle.get("n_sequences", ""),
        "sequence_mode": selection.get("mode", ""),
        "retain_ratio": selection.get("ratio", ""),
        "wgdi_ready": selection.get("wgdi_ready", ""),
        "rejection_counts": qc.get("rejection_counts") or {},
        "nuclear_taxa": ",".join(str(item) for item in nuclear_taxa),
        "n_anomalies": len(anomalies),
        "reason": "",
        "detail": "existing_valid_output_no_overwrite",
    }


_WINDOWS_PATH_LIMIT = 250


def _guard_windows_path_budget(species: SpeciesInput, options: PipelineOptions) -> None:
    """Fail fast with an actionable error instead of a cryptic OSError.

    Windows still enforces the legacy 260-character path limit unless long
    paths are enabled machine-wide.  Species folders plus per-species outputs
    stack enough segments that a deep project root silently exceeds it: the
    directory gets created but files cannot be written, which surfaces much
    later as ``FileNotFoundError``.
    """

    if os.name != "nt":
        return
    candidates = [
        species.directory,
        species.directory / species.gff_name,
        species.directory / species.cds_name,
        species.directory / species.pep_name,
        options.project_root / options.analysis_dirname / species.code,
        options.project_root / options.clean_dirname / species.code,
    ]
    for path in candidates:
        length = len(str(path))
        if length <= _WINDOWS_PATH_LIMIT:
            continue
        raise PipelineFailure(
            "path_too_long_for_windows",
            f"{path} is {length} characters; the Windows limit is about 260. "
            "Move the project to a shorter root, shorten the species folder name, "
            "or enable long-path support for the machine.",
        )


def process_species(species: SpeciesInput, options: PipelineOptions) -> dict[str, object]:
    if not options.overwrite:
        existing_qc = _existing_qc(species, options)
        if existing_qc is not None:
            qc, artifact_code = existing_qc
            return _skipped_result(species, qc, artifact_code)
        if _existing_outputs_present(species, options):
            return {
                "show_id": species.show_id,
                "code": species.code,
                "species": species.species,
                "accession": species.accession,
                "status": "FAIL",
                "reason": "existing_output_invalid_no_overwrite",
                "detail": str(
                    options.project_root / options.analysis_dirname / species.code
                ),
            }
    try:
        _guard_windows_path_budget(species, options)
    except PipelineFailure as exc:
        return {
            "show_id": species.show_id,
            "code": species.code,
            "species": species.species,
            "accession": species.accession,
            "status": "FAIL",
            "reason": exc.reason,
            "detail": exc.detail,
        }
    clean_root, analysis_root, clean_staging, analysis_staging = _prepare_staging(
        options, species.code
    )
    clean_final = clean_root / species.code
    analysis_final = analysis_root / species.code
    try:
        result = _process_species_into(
            species,
            options,
            clean_staging,
            analysis_staging,
        )
        _replace_directory(clean_staging, clean_final, options.overwrite)
        _replace_directory(analysis_staging, analysis_final, options.overwrite)
        return result
    except PipelineFailure as exc:
        _remove_staging(clean_staging, clean_root)
        _remove_staging(analysis_staging, analysis_root)
        return {
            "show_id": species.show_id,
            "code": species.code,
            "species": species.species,
            "accession": species.accession,
            "status": "FAIL",
            "reason": exc.reason,
            "detail": exc.detail,
        }
    except GffFormatError as exc:
        _remove_staging(clean_staging, clean_root)
        _remove_staging(analysis_staging, analysis_root)
        return {
            "show_id": species.show_id,
            "code": species.code,
            "species": species.species,
            "accession": species.accession,
            "status": "FAIL",
            "reason": "gff_format_error",
            "detail": str(exc),
        }
    except Exception:
        _remove_staging(clean_staging, clean_root)
        _remove_staging(analysis_staging, analysis_root)
        raise


def _process_species_into(
    species: SpeciesInput,
    options: PipelineOptions,
    clean_dir: Path,
    analysis_dir: Path,
) -> dict[str, object]:
    accession = _discover_accession(species, options.accession_policy)
    gff_path = species.directory / species.gff_name
    cds_path = species.directory / species.cds_name
    pep_path = species.directory / species.pep_name
    genome_path = (
        species.directory / species.genome_name
        if species.genome_name
        else None
    )
    for path in (gff_path, cds_path, pep_path):
        if not path.is_file() or path.stat().st_size == 0:
            if not species.cds_name or not species.pep_name:
                raise PipelineFailure(
                    "sequence_extraction_required",
                    "CDS/PEP are absent; run `pgci_preprocess prepare "
                    "--source-status <table>` with genome+GFF to derive them",
                )
            raise PipelineFailure("missing_or_empty_input", str(path))
    try:
        annotation_kind = detect_conversion_format(gff_path)
    except ConversionError:
        annotation_kind = ""
    if annotation_kind in {"gtf", "gbff"}:
        fmt = annotation_kind
        raise PipelineFailure(
            "conversion_required",
            f"raw {fmt.upper()} cannot be run directly; use "
            f"`pgci_preprocess convert-annotation --input \"{gff_path}\" "
            f"--output-dir <derived> --format {fmt}` and then register the "
            "converted GFF/bundle, or let `prepare`/`batch` route it automatically",
        )

    if options.format_exceptions is None:
        options.format_exceptions = load_format_exceptions(options.project_root)
    format_flags = options.format_exceptions.get(species.show_id, {})

    parsed_result = parse_with_adapter(
        gff_path,
        requested_format=options.gff_format,
        merge_split_genes=bool(format_flags.get("merge_split_genes")),
    )
    parsed = parsed_result.data
    gwh_mapped = 0
    if parsed_result.detection.has_transcript:
        gwh_mapped = apply_gwh_accessions(parsed, mrna_accession_map(gff_path))
    cds_prefix_changed = normalize_cds_protein_ids(parsed)
    fatal_gff = [item for item in parsed.anomalies if item.level == "fatal"]
    if fatal_gff:
        raise PipelineFailure(fatal_gff[0].code, fatal_gff[0].detail)

    region_backfill: dict[str, list[str]] = {
        "genome": [],
        "transcript": list(parsed_result.region_backfill),
    }
    missing_region_seqids = {
        transcript.seqid
        for transcript in parsed.transcripts.values()
        if transcript.seqid not in parsed.regions
    }
    synthetic_seqids = {
        seqid for seqid, region in parsed.regions.items() if region.synthetic
    }
    if genome_path is not None and genome_path.is_file() and genome_path.stat().st_size > 0:
        genome_regions = _read_genome_regions(genome_path)
        if genome_regions:
            # Authoritative genome FASTA lengths supersede synthetic max-end
            # regions (transcript-only / simplified adapters).  Real GFF
            # region rows and ##sequence-region directives are never touched.
            for seqid in sorted(synthetic_seqids):
                if seqid in genome_regions:
                    parsed.regions[seqid] = genome_regions[seqid]
                    region_backfill["genome"].append(seqid)
            for seqid in sorted(missing_region_seqids):
                if seqid in genome_regions and seqid not in parsed.regions:
                    parsed.regions[seqid] = genome_regions[seqid]
                    region_backfill["genome"].append(seqid)
            for seqid, gff_region in sorted(parsed.regions.items()):
                genome_region = genome_regions.get(seqid)
                if genome_region is None or gff_region.synthetic:
                    continue
                if gff_region.length > genome_region.length:
                    raise PipelineFailure(
                        "gff_genome_length_mismatch",
                        f"{seqid}: gff={gff_region.length} genome={genome_region.length}",
                    )
                if gff_region.length != genome_region.length:
                    parsed.anomalies.append(
                        Anomaly(
                            "warning",
                            "gff_genome_length_mismatch",
                            f"{seqid}: gff={gff_region.length} genome={genome_region.length}",
                        )
                    )
            if region_backfill["genome"]:
                parsed.anomalies.append(
                    Anomaly(
                        "warning",
                        "region_backfilled_from_genome",
                        f"seqids={','.join(sorted(region_backfill['genome'])[:5])}",
                    )
                )
    unresolved = {
        seqid for seqid in missing_region_seqids if seqid not in parsed.regions
    }
    if unresolved and options.backfill_regions_from_transcripts:
        region_backfill["transcript"].extend(
            backfill_regions_from_transcripts(parsed, reason="gff_missing_regions")
        )
        unresolved = {
            seqid for seqid in missing_region_seqids if seqid not in parsed.regions
        }
    if unresolved:
        raise PipelineFailure(
            "region_unresolved",
            ",".join(sorted(unresolved)[:5]),
        )

    nuclear_taxa = {
        region.taxon_id
        for region in parsed.regions.values()
        if region.taxon_id and not region.is_organelle
    }
    if len(nuclear_taxa) > 1:
        raise PipelineFailure(
            "multiple_nuclear_taxa",
            f"{species.species}: taxon IDs {sorted(nuclear_taxa)}",
        )

    cds = read_cds_fasta(cds_path)
    pep = read_pep_fasta(pep_path)
    protein_normalization = {
        "mode": options.protein_normalization,
        "applied": "none",
        "records": 0,
        "trailing_stop": 0,
        "internal_wildcard": 0,
    }
    if options.protein_normalization in {"auto", "evm-dot"}:
        has_dots = any("." in sequence for sequence in pep.sequences.values())
        if options.protein_normalization == "evm-dot" or has_dots:
            report = normalize_evm_protein_sequences(pep)
            protein_normalization.update(
                {
                    "applied": "evm_dot" if report["records"] else "none",
                    "records": report["records"],
                    "trailing_stop": report["trailing_stop"],
                    "internal_wildcard": report["internal_wildcard"],
                }
            )
    fatal_fasta = [
        item for item in [*cds.anomalies, *pep.anomalies] if item.level == "fatal"
    ]
    if fatal_fasta:
        raise PipelineFailure(fatal_fasta[0].code, fatal_fasta[0].detail)

    # Annotations whose CDS rows carry no protein_id (e.g. Phytozome GFF3) are
    # resolved against the FASTA keys themselves: the transcript id or its Name
    # must match the SAME key in both CDS and PEP. JGI polypeptide files append a
    # ".p" suffix to the transcript key; that alias is recorded, never guessed.
    protein_id_from_fasta_key = 0
    pep_alias: dict[str, str] = {}
    for transcript in parsed.transcripts.values():
        if transcript.protein_ids:
            continue
        matched: dict[str, str] = {}
        for value in (transcript.source_id, transcript.name):
            if not value or value not in cds.sequences:
                continue
            if value in pep.sequences:
                matched[value] = value
            elif f"{value}.p" in pep.sequences:
                matched[value] = f"{value}.p"
        if len(matched) == 1:
            protein_id, pep_key = next(iter(matched.items()))
            transcript.protein_ids = {protein_id}
            if pep_key != protein_id:
                pep_alias[protein_id] = pep_key
            protein_id_from_fasta_key += 1

    rejected: list[dict[str, str]] = []
    candidates: list[Candidate] = []
    reason_counts: Counter[str] = Counter()
    for transcript_id, transcript in parsed.transcripts.items():
        gene = parsed.genes.get(transcript.gene_id)
        if gene is None:
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                "gene_missing",
                transcript.gene_id,
            )
            reason_counts["gene_missing"] += 1
            continue
        region = parsed.regions.get(transcript.seqid)
        if region is None:
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                "region_missing",
                transcript.seqid,
            )
            reason_counts["region_missing"] += 1
            continue
        if region.is_organelle:
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                "organelle_sequence",
                f"{transcript.seqid}:{region.genome_type}",
            )
            reason_counts["organelle_sequence"] += 1
            continue
        if transcript.strand not in {"+", "-"}:
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                "invalid_strand",
                transcript.strand,
            )
            reason_counts["invalid_strand"] += 1
            continue
        if len(transcript.protein_ids) != 1:
            reason = "protein_id_missing" if not transcript.protein_ids else "multiple_protein_ids"
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                reason,
                ",".join(sorted(transcript.protein_ids)),
            )
            reason_counts[reason] += 1
            continue
        protein_id = next(iter(transcript.protein_ids))
        pep_key = pep_alias.get(protein_id, protein_id)
        cds_sequence = cds.sequences.get(protein_id)
        pep_sequence = pep.sequences.get(pep_key)
        if cds_sequence is None or pep_sequence is None:
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                "sequence_missing",
                f"cds={protein_id in cds.sequences} pep={pep_key in pep.sequences}",
            )
            reason_counts["sequence_missing"] += 1
            continue
        check = validate_cds_pep(
            cds_sequence,
            pep_sequence,
            translation_table=transcript.translation_table,
        )
        if not check.valid:
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                check.reason,
                protein_id,
            )
            reason_counts[check.reason] += 1
            continue
        if (
            transcript.start < 1
            or transcript.end < transcript.start
            or transcript.end > region.length
        ):
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                "coordinate_out_of_region",
                f"{transcript.start}-{transcript.end}/{region.length}",
            )
            reason_counts["coordinate_out_of_region"] += 1
            continue
        candidates.append(
            Candidate(
                source_gene_id=gene.source_id,
                source_gene_name=gene.name,
                gene_biotype=gene.biotype,
                source_transcript_id=transcript_id,
                source_protein_id=protein_id,
                seqid=transcript.seqid,
                start=transcript.start,
                end=transcript.end,
                strand=transcript.strand,
                cds_length=len(cds_sequence),
                protein_length=len(pep_sequence.rstrip("*")),
                translation_exact=check.exact,
                translation_mismatches=check.mismatches,
            )
        )

    protein_to_genes: dict[str, set[str]] = defaultdict(set)
    for candidate in candidates:
        protein_to_genes[candidate.source_protein_id].add(candidate.source_gene_id)
    ambiguous_proteins = {
        protein_id for protein_id, gene_ids in protein_to_genes.items() if len(gene_ids) > 1
    }
    if ambiguous_proteins:
        retained: list[Candidate] = []
        for candidate in candidates:
            if candidate.source_protein_id in ambiguous_proteins:
                _add_rejection(
                    rejected,
                    "transcript",
                    candidate.source_transcript_id,
                    "protein_id_reused_across_genes",
                    candidate.source_protein_id,
                )
                reason_counts["protein_id_reused_across_genes"] += 1
            else:
                retained.append(candidate)
        candidates = retained

    by_gene: dict[str, list[Candidate]] = defaultdict(list)
    for candidate in candidates:
        by_gene[candidate.source_gene_id].append(candidate)
    representatives: list[Candidate] = []
    for gene_id in sorted(by_gene):
        group = by_gene[gene_id]
        representative = choose_representative(group)
        representatives.append(representative)
        for candidate in group:
            if candidate is not representative:
                _add_rejection(
                    rejected,
                    "transcript",
                    candidate.source_transcript_id,
                    "alternative_isoform",
                    f"representative={representative.source_transcript_id}",
                )
                reason_counts["alternative_isoform"] += 1
    translation_tables = {
        transcript.source_id: transcript.translation_table
        for transcript in parsed.transcripts.values()
    }

    if not representatives:
        # An annotation whose every transcript was rejected as organelle-only
        # (chloroplast/mitochondrion submissions such as the GWH GeSeq stubs)
        # is a data-availability problem, not a processing bug.  Report it as
        # such so batch triage never mistakes it for a converter failure.
        if reason_counts and set(reason_counts) == {"organelle_sequence"}:
            raise PipelineFailure(
                "organelle_only_annotation",
                f"transcripts={reason_counts['organelle_sequence']}",
            )
        raise PipelineFailure("no_valid_representative_genes")

    gene_counts: Counter[str] = Counter(candidate.seqid for candidate in representatives)
    sequence_selection = select_analysis_sequences(
        parsed.regions,
        gene_counts,
        options.min_genes_per_sequence,
        options.target_gene_coverage,
        options.max_sequences,
    )
    if not sequence_selection.keep:
        raise PipelineFailure("empty_wgdi_sequence_set")

    nuclear_with_genes = [
        parsed.regions[seqid]
        for seqid, count in gene_counts.items()
        if count > 0 and seqid in parsed.regions and not parsed.regions[seqid].is_organelle
    ]
    nuclear_with_genes.sort(key=region_sort_key)
    clean_sequence_id = {
        region.seqid: index for index, region in enumerate(nuclear_with_genes, 1)
    }

    by_sequence: dict[str, list[Candidate]] = defaultdict(list)
    for candidate in representatives:
        by_sequence[candidate.seqid].append(candidate)

    entries: list[dict[str, object]] = []
    for region in nuclear_with_genes:
        rows = sorted(
            by_sequence[region.seqid],
            key=lambda candidate: (
                candidate.start,
                candidate.end,
                candidate.source_gene_id,
            ),
        )
        for order, candidate in enumerate(rows, 1):
            canonical_id = (
                f"{species.code}{clean_sequence_id[candidate.seqid]}"
                f"g{order:04d}"
            )
            entries.append(
                {
                    **asdict(candidate),
                    "canonical_gene_id": canonical_id,
                    "clean_sequence_id": clean_sequence_id[candidate.seqid],
                    "order": order,
                    "selected_for_wgdi": candidate.seqid in sequence_selection.keep,
                    "cds_sequence": cds.sequences[candidate.source_protein_id],
                    "pep_sequence": pep.sequences[pep_alias.get(candidate.source_protein_id, candidate.source_protein_id)].rstrip("*"),
                }
            )

    selected_entries = [
        entry for entry in entries if bool(entry["selected_for_wgdi"])
    ]
    if not selected_entries:
        raise PipelineFailure("empty_wgdi_gene_set")

    all_region_rows: list[dict[str, object]] = []
    selected_region_rows: list[dict[str, object]] = []
    for region in nuclear_with_genes:
        row = {
            "clean_sequence_id": clean_sequence_id[region.seqid],
            "source_seqid": region.seqid,
            "length": region.length,
            "gene_count": gene_counts[region.seqid],
            "selected": region.seqid in sequence_selection.keep,
        }
        all_region_rows.append(row)
        if row["selected"]:
            selected_region_rows.append(row)

    all_cds_records = [
        {
            "canonical_gene_id": entry["canonical_gene_id"],
            "source_gene_id": entry["source_gene_id"],
            "source_protein_id": entry["source_protein_id"],
            "sequence": entry["cds_sequence"],
        }
        for entry in entries
    ]
    all_pep_records = [
        {
            "canonical_gene_id": entry["canonical_gene_id"],
            "source_gene_id": entry["source_gene_id"],
            "source_protein_id": entry["source_protein_id"],
            "sequence": entry["pep_sequence"],
        }
        for entry in entries
    ]
    selected_cds_records = [
        record
        for entry, record in zip(entries, all_cds_records)
        if bool(entry["selected_for_wgdi"])
    ]
    selected_pep_records = [
        record
        for entry, record in zip(entries, all_pep_records)
        if bool(entry["selected_for_wgdi"])
    ]

    code = species.code
    all_cds_path = analysis_dir / f"{code}.all.cds.fa"
    all_pep_path = analysis_dir / f"{code}.all.pep.fa"
    all_gff_path = analysis_dir / f"{code}.all.gff"
    all_lens_path = analysis_dir / f"{code}.all.lens"
    cds_out = analysis_dir / f"{code}.cds"
    pep_out = analysis_dir / f"{code}.pep"
    gff_out = analysis_dir / f"{code}.new.gff"
    lens_out = analysis_dir / f"{code}.lens"
    _write_fasta(all_cds_path, all_cds_records)
    _write_fasta(all_pep_path, all_pep_records)
    _write_gff(all_gff_path, entries)
    _write_lens(all_lens_path, all_region_rows)
    _write_fasta(cds_out, selected_cds_records)
    _write_fasta(pep_out, selected_pep_records)
    _write_gff(gff_out, selected_entries)
    _write_lens(lens_out, selected_region_rows)

    all_validation = _validate_bundle(
        all_gff_path, all_lens_path, all_cds_path, all_pep_path
    )
    selected_validation = _validate_bundle(gff_out, lens_out, cds_out, pep_out)

    gene_map_path = clean_dir / f"{code}.gene_map.tsv"
    gene_map_fields = [
        "canonical_gene_id",
        "clean_sequence_id",
        "canonical_order",
        "source_gene_id",
        "source_gene_name",
        "source_transcript_id",
        "source_protein_id",
        "source_seqid",
        "start",
        "end",
        "strand",
        "gene_biotype",
        "cds_length",
        "protein_length",
        "translation_exact",
        "translation_mismatches",
        "translation_table",
        "selected_for_wgdi",
    ]
    with gene_map_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=gene_map_fields, delimiter="\t")
        writer.writeheader()
        for entry in entries:
            writer.writerow(
                {
                    "canonical_gene_id": entry["canonical_gene_id"],
                    "clean_sequence_id": entry["clean_sequence_id"],
                    "canonical_order": entry["order"],
                    "source_gene_id": entry["source_gene_id"],
                    "source_gene_name": entry["source_gene_name"],
                    "source_transcript_id": entry["source_transcript_id"],
                    "source_protein_id": entry["source_protein_id"],
                    "source_seqid": entry["seqid"],
                    "start": entry["start"],
                    "end": entry["end"],
                    "strand": entry["strand"],
                    "gene_biotype": entry["gene_biotype"],
                    "cds_length": entry["cds_length"],
                    "protein_length": entry["protein_length"],
                    "translation_exact": "Y" if entry["translation_exact"] else "N",
                    "translation_mismatches": entry["translation_mismatches"],
                    "translation_table": translation_tables[
                        str(entry["source_transcript_id"])
                    ],
                    "selected_for_wgdi": "Y" if entry["selected_for_wgdi"] else "N",
                }
            )

    chr_map_path = clean_dir / f"{code}.chr_map.tsv"
    with chr_map_path.open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "source_seqid",
            "clean_sequence_id",
            "length",
            "genome_type",
            "taxon_id",
            "is_primary",
            "is_organelle",
            "gene_count",
            "in_wgdi_set",
            "exclusion_reason",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for region in sorted(
            parsed.regions.values(),
            key=lambda item: (
                item.is_organelle,
                *region_sort_key(item),
            ),
        ):
            count = gene_counts.get(region.seqid, 0)
            clean_id = clean_sequence_id.get(region.seqid, "")
            reason = ""
            if region.is_organelle:
                reason = "organelle_sequence"
            elif count == 0:
                reason = "no_valid_protein_gene"
            elif region.seqid not in sequence_selection.keep:
                reason = "outside_wgdi_sequence_set"
            writer.writerow(
                {
                    "source_seqid": region.seqid,
                    "clean_sequence_id": clean_id,
                    "length": region.length,
                    "genome_type": region.genome_type,
                    "taxon_id": region.taxon_id,
                    "is_primary": "Y" if region.is_primary else "N",
                    "is_organelle": "Y" if region.is_organelle else "N",
                    "gene_count": count,
                    "in_wgdi_set": "Y" if region.seqid in sequence_selection.keep else "N",
                    "exclusion_reason": reason,
                }
            )

    rejected_path = clean_dir / f"{code}.rejected.tsv"
    with rejected_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["scope", "id", "reason", "detail"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rejected)

    input_info = {}
    input_paths = [("gff", gff_path), ("cds", cds_path), ("pep", pep_path)]
    if genome_path is not None and genome_path.is_file():
        input_paths.append(("genome", genome_path))
    for label, path in input_paths:
        input_info[label] = {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": _sha256(path),
        }
    qc = {
        "show_id": species.show_id,
        "code": species.code,
        "species": species.species,
        "order": species.order,
        "family": species.family,
        "accession": accession,
        "assembly_id_source": (
            "assembly_accessions.tsv"
            if (species.assembly_id or "").strip()
            else "legacy_status_table"
            if options.accession_policy == "legacy"
            else "filename"
        ),
        "assembly_build": parsed.assembly_build,
        "annotation_source": parsed.annotation_source,
        "annotation_level": parsed.annotation_level,
        "split_genes_merged": parsed.split_genes_merged,
        "protein_id_from_fasta_key": protein_id_from_fasta_key,
        "pep_key_alias": len(pep_alias),
        "format_detection": parsed_result.detection.as_dict(),
        "format_policy": {
            "requested_gff_format": options.gff_format,
            "selected_gff_format": parsed_result.detection.format,
            "protein_normalization": options.protein_normalization,
            "accession_policy": options.accession_policy,
            "transcript_region_backfill": options.backfill_regions_from_transcripts,
        },
        "format_exceptions": {
            "merge_split_genes": bool(format_flags.get("merge_split_genes")),
            "gene_span_only": bool(format_flags.get("gene_span_only"))
            or parsed_result.detection.format == "gene_span",
            "gene_cds_direct": bool(format_flags.get("gene_cds_direct"))
            or parsed_result.detection.format == "gene_cds_direct",
            "transcript_only_gff": parsed_result.detection.format == "transcript_only",
            "simplified5_gff": parsed_result.detection.format == "simplified5",
            "gwh_accessions": gwh_mapped,
            "cds_protein_id_prefix_stripped": cds_prefix_changed,
            "region_backfilled_from_genome": len(region_backfill["genome"]),
            "region_backfilled_from_transcript": len(region_backfill["transcript"]),
            "protein_dot_normalized": protein_normalization["records"],
        },
        "region_backfill": region_backfill,
        "protein_normalization": protein_normalization,
        "pipeline_version": PIPELINE_VERSION,
        "inputs": input_info,
        "n_regions": len(parsed.regions),
        "n_genes": len(parsed.genes),
        "n_transcripts": len(parsed.transcripts),
        "n_candidates": len(candidates),
        "n_representative_genes": len(representatives),
        "n_organelle_genes": len(
            {
                transcript.gene_id
                for transcript in parsed.transcripts.values()
                if transcript.seqid in parsed.regions
                and parsed.regions[transcript.seqid].is_organelle
            }
        ),
        "n_cds_records": cds.counts["records"],
        "n_pep_records": pep.counts["records"],
        "n_cds_duplicate_conflicts": cds.counts["duplicate_conflict"],
        "n_pep_duplicate_conflicts": pep.counts["duplicate_conflict"],
        "rejection_counts": dict(sorted(reason_counts.items())),
        "nuclear_taxa": sorted(nuclear_taxa),
        "sequence_selection": {
            "mode": sequence_selection.mode,
            "selected_genes": sequence_selection.selected_genes,
            "total_genes": sequence_selection.total_genes,
            "ratio": sequence_selection.ratio,
            "wgdi_ready": sequence_selection.ready,
            "min_genes_per_sequence": options.min_genes_per_sequence,
            "target_gene_coverage": options.target_gene_coverage,
            "max_sequences": options.max_sequences,
        },
        "all_bundle": all_validation,
        "wgdi_bundle": selected_validation,
        "gff_feature_counts": dict(sorted(parsed.feature_counts.items())),
        "anomalies": [asdict(item) for item in [*parsed.anomalies, *cds.anomalies, *pep.anomalies]],
        "artifact_naming": {
            "code": code,
            "table": "resources/species_abbreviations.tsv",
            "rule": "<code>/<code>.<role> and gene IDs <code><seq>g<order>",
            "verified": True,
        },
        "status": "OK",
    }
    with (analysis_dir / f"{code}.qc.json").open("w", encoding="utf-8") as handle:
        json.dump(qc, handle, ensure_ascii=False, indent=2)
    shutil.copy2(gene_map_path, analysis_dir / f"{code}.gene_map.tsv")
    shutil.copy2(chr_map_path, analysis_dir / f"{code}.chr_map.tsv")

    naming_problems = [
        f"analysis:{problem}"
        for problem in bundle_name_problems(
            analysis_dir,
            code,
            required_suffixes=ANALYSIS_REQUIRED_SUFFIXES,
            kind_map={".new.gff": "gff", ".cds": "fasta", ".pep": "fasta"},
            check_dir_name=False,
        )
    ] + [
        f"clean:{problem}"
        for problem in bundle_name_problems(
            clean_dir,
            code,
            required_suffixes=CLEAN_REQUIRED_SUFFIXES,
            check_dir_name=False,
        )
    ]
    if naming_problems:
        raise PipelineFailure("naming_mismatch", ";".join(naming_problems))

    return {
        "show_id": species.show_id,
        "code": species.code,
        "species": species.species,
        "accession": accession,
        "status": "OK",
        "assembly_build": parsed.assembly_build,
        "n_regions": len(parsed.regions),
        "n_representative_genes": len(representatives),
        "n_all_genes": all_validation["n_genes"],
        "n_wgdi_genes": selected_validation["n_genes"],
        "n_wgdi_sequences": selected_validation["n_sequences"],
        "sequence_mode": sequence_selection.mode,
        "retain_ratio": sequence_selection.ratio,
        "wgdi_ready": sequence_selection.ready,
        "rejection_counts": dict(sorted(reason_counts.items())),
        "nuclear_taxa": ",".join(sorted(nuclear_taxa)),
        "n_anomalies": len(parsed.anomalies) + len(cds.anomalies) + len(pep.anomalies),
    }


def _promote_legacy_bundles(
    options: PipelineOptions,
    species: list[SpeciesInput],
) -> None:
    """Rename any legacy ``s<show_id>`` bundle of the selected species once.

    The rename is a pure naming migration (same sequences, same gene order)
    and is verified inside :mod:`pgci_preprocess.artifact_migration`.  It runs
    before the per-species loop so the index/manifest and compatibility links
    are refreshed once instead of once per species.
    """

    selected = [
        item
        for item in species
        if any(
            (options.project_root / layer / item.legacy_code).is_dir()
            for layer in (options.clean_dirname, options.analysis_dirname)
        )
    ]
    if not selected:
        return
    try:
        result = migrate_artifacts(
            options.project_root,
            show_ids={item.show_id for item in selected},
            refresh_links=True,
        )
    except Exception as exc:  # migration must never hide the real pipeline error
        print(f"[WARN] legacy artifact rename skipped: {exc}", flush=True)
        return
    migrated = [item for item in result.species if item.status == "MIGRATED"]
    failed = [item for item in result.species if item.status == "FAIL"]
    if migrated:
        print(
            f"[NAMING] renamed {len(migrated)} legacy bundle(s) to the "
            "canonical abbreviation table",
            flush=True,
        )
    for item in failed:
        print(
            f"[NAMING][FAIL] {item.show_id} -> {item.code}: "
            + "; ".join(item.problems),
            flush=True,
        )


def run_project(
    project_root: Path | str,
    *,
    show_ids: set[str] | None = None,
    options: PipelineOptions | None = None,
) -> list[dict[str, object]]:
    project_root = Path(project_root)
    options = options or PipelineOptions(project_root=project_root)
    all_species = discover_species(
        project_root,
        species_names_path=options.species_names_path,
    )
    species = all_species
    if show_ids is not None:
        selected_show_ids = resolve_show_ids(
            project_root,
            show_ids,
            explicit_path=options.species_names_path,
        )
        species = [item for item in species if item.show_id in selected_show_ids]
    selected_codes = {
        code
        for item in species
        for code in (item.code, item.legacy_code)
    }

    results: list[dict[str, object]] = []
    if options.migrate_legacy_artifacts and not options.overwrite:
        _promote_legacy_bundles(options, species)
    for item in species:
        try:
            result = process_species(item, options)
        except Exception as exc:
            result = {
                "show_id": item.show_id,
                "code": item.code,
                "species": item.species,
                "accession": item.accession,
                "status": "ERROR",
                "reason": type(exc).__name__,
                "detail": str(exc),
            }
        results.append(result)
        print(
            f"[{result['status']}] {item.code} show_id={item.show_id} {item.species} "
            f"genes={result.get('n_wgdi_genes', '')} "
            f"seqs={result.get('n_wgdi_sequences', '')} "
            f"ready={result.get('wgdi_ready', '')} "
            f"reason={result.get('reason', '')}",
            flush=True,
        )

    species_rows = [
        {
            "show_id": item.show_id,
            "code": item.code,
            "species": item.species,
            "order": item.order,
            "family": item.family,
            "accession": item.accession,
            "gff": item.gff_name,
            "cds": item.cds_name,
            "pep": item.pep_name,
            "genome": item.genome_name,
            "directory": item.directory.name,
        }
        for item in all_species
    ]
    if species:
        finalize_project(
            project_root,
            species_rows=species_rows,
            results=results,
            processed_codes=selected_codes,
            analysis_dirname=options.analysis_dirname,
            clean_dirname=options.clean_dirname,
        )
    return results

