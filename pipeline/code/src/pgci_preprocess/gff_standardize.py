from __future__ import annotations

import gzip
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

from pgci_preprocess.fasta import open_text
from pgci_preprocess.gff_formats import GffFormatError, detect_gff_format
from pgci_preprocess.reports import read_tsv_with_fields, write_tsv
from pgci_preprocess.species_names import resolve_show_ids, resolve_species_code


TRANSCRIPT_TYPES = {
    "mRNA",
    "transcript",
    "pseudogenic_transcript",
    "lnc_RNA",
    "tRNA",
    "rRNA",
    "snRNA",
    "snoRNA",
    "ncRNA",
    "miRNA",
}
GENE_TYPES = {"gene", "pseudogene"}


class GffStandardizationError(RuntimeError):
    pass


@dataclass(frozen=True)
class StandardizationResult:
    converted: bool
    source_path: Path
    output_path: Path | None
    feature_counts: dict[str, int]
    transcript_count: int
    gene_count: int


def _attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        if not item:
            continue
        key, separator, value = item.partition("=")
        if separator:
            result[key.strip()] = unquote(value.strip())
    return result


def _attribute_text(values: dict[str, str]) -> str:
    return ";".join(f"{key}={value}" for key, value in values.items() if value)


def _strip_rna_prefix(value: str) -> str:
    return value[4:] if value.startswith("rna-") else value


def _gff_rows(path: Path):
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if raw_line.startswith("#") or not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise GffStandardizationError(
                    f"{path.name}:{line_number}: expected 9 columns, found {len(fields)}"
                )
            yield line_number, raw_line.rstrip("\r\n"), fields


def _inspect(path: Path) -> tuple[Counter[str], dict[str, tuple[str, str, str, str]]]:
    counts: Counter[str] = Counter()
    transcripts: dict[str, tuple[str, str, str, str]] = {}
    for line_number, _raw, fields in _gff_rows(path):
        seqid, _source, feature, start, end, _score, strand, _phase, attr_text = fields
        counts[feature] += 1
        if feature not in TRANSCRIPT_TYPES:
            continue
        attributes = _attributes(attr_text)
        transcript_id = _strip_rna_prefix(attributes.get("ID", "").strip())
        if not transcript_id:
            raise GffStandardizationError(
                f"{path.name}:{line_number}: transcript has no ID"
            )
        if transcript_id in transcripts:
            raise GffStandardizationError(
                f"{path.name}:{line_number}: duplicate transcript ID {transcript_id}"
            )
        transcripts[transcript_id] = (seqid, start, end, strand)
    return counts, transcripts


def standardize_transcript_only_gff(
    source_path: Path | str,
    output_path: Path | str,
) -> StandardizationResult:
    source = Path(source_path)
    output = Path(output_path)
    counts, transcripts = _inspect(source)
    gene_count = sum(counts[item] for item in GENE_TYPES)
    if gene_count:
        return StandardizationResult(
            converted=False,
            source_path=source,
            output_path=None,
            feature_counts=dict(sorted(counts.items())),
            transcript_count=len(transcripts),
            gene_count=gene_count,
        )
    if not transcripts:
        raise GffStandardizationError(
            f"{source.name}: no gene rows and no supported transcript rows"
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wt", encoding="utf-8", newline="\n") as destination:
        destination.write("##gff-version 3\n")
        destination.write(
            "# normalized: transcript-only EVM GFF expanded to gene/mRNA/CDS hierarchy\n"
        )
        for line_number, raw_line, fields in _gff_rows(source):
            seqid, source_name, feature, start, end, score, strand, phase, attr_text = fields
            attributes = _attributes(attr_text)
            if feature in TRANSCRIPT_TYPES:
                transcript_id = _strip_rna_prefix(attributes.get("ID", "").strip())
                gene_id = f"gene-{transcript_id}"
                gene = [
                    seqid,
                    source_name,
                    "gene",
                    start,
                    end,
                    score,
                    strand,
                    ".",
                    _attribute_text({"ID": gene_id, "Name": transcript_id}),
                ]
                destination.write("\t".join(gene) + "\n")
                attributes = dict(attributes)
                attributes["ID"] = f"rna-{transcript_id}"
                attributes["Parent"] = gene_id
                normalized = [
                    seqid,
                    source_name,
                    feature,
                    start,
                    end,
                    score,
                    strand,
                    phase,
                    _attribute_text(attributes),
                ]
                destination.write("\t".join(normalized) + "\n")
                continue

            parent_text = attributes.get("Parent", "")
            if parent_text:
                parents = [
                    _strip_rna_prefix(item.strip())
                    for item in parent_text.split(",")
                    if item.strip()
                ]
                if feature == "CDS" and len(parents) != 1:
                    raise GffStandardizationError(
                        f"{source.name}:{line_number}: CDS must have exactly one transcript parent"
                    )
                if parents and all(item in transcripts for item in parents):
                    attributes["Parent"] = ",".join(f"rna-{item}" for item in parents)
                    if feature == "CDS" and not (
                        attributes.get("protein_id")
                        or attributes.get("Protein_Accession")
                    ):
                        attributes["protein_id"] = parents[0]
                    fields[8] = _attribute_text(attributes)
                    destination.write("\t".join(fields) + "\n")
                    continue
            destination.write(raw_line + "\n")

    return StandardizationResult(
        converted=True,
        source_path=source,
        output_path=output,
        feature_counts=dict(sorted(counts.items())),
        transcript_count=len(transcripts),
        gene_count=gene_count,
    )


@dataclass(frozen=True, slots=True)
class SpeciesNormalizationResult:
    show_id: str
    species_code: str
    status: str
    source_path: Path
    output_path: Path | None
    detail: str = ""


def normalize_gff_species(
    project_root: Path | str,
    show_ids: list[str] | tuple[str, ...],
    *,
    update_inventory: bool = True,
    overwrite: bool = False,
    species_names_path: Path | str | None = None,
) -> list[SpeciesNormalizationResult]:
    """Normalize transcript-only GFF files for selected inventory species.

    This is the shared implementation behind the ``normalize-gff`` subcommand
    and the legacy ``standardize_transcript_gff.py`` wrapper.  Inventory rows
    are only rewritten when at least one file was converted and
    ``update_inventory`` is true.
    """

    project = Path(project_root).resolve()
    inventory_path = project / "resources" / "download_inventory.tsv"
    if not inventory_path.is_file():
        raise GffStandardizationError(f"download inventory not found: {inventory_path}")
    fields, inventory = read_tsv_with_fields(inventory_path)
    if "species_code" not in fields:
        insert_at = fields.index("showid") + 1 if "showid" in fields else 0
        fields.insert(insert_at, "species_code")
    by_id = {row.get("showid", ""): row for row in inventory}
    requested_ids = resolve_show_ids(
        project,
        show_ids,
        explicit_path=species_names_path,
    )
    downloads = project / "data" / "downloads"
    directories = {
        path.name.split("_", 1)[0]: path
        for path in downloads.iterdir()
        if path.is_dir() and "_" in path.name
    }

    results: list[SpeciesNormalizationResult] = []
    converted = 0
    inventory_changed = False
    for show_id in sorted(requested_ids, key=lambda value: int(value) if value.isdigit() else 10**12):
        row = by_id.get(show_id)
        directory = directories.get(show_id)
        if row is None or directory is None:
            results.append(
                SpeciesNormalizationResult(
                    show_id=show_id,
                    species_code="",
                    status="ERROR",
                    source_path=Path(""),
                    output_path=None,
                    detail="missing inventory row or species directory",
                )
            )
            continue
        species_code = resolve_species_code(
            project,
            show_id,
            latin_name=str(row.get("Species", "")),
            explicit_path=species_names_path,
        )
        row["species_code"] = species_code
        inventory_changed = True
        relative = Path(row.get("gff", ""))
        source = (directory / relative).resolve()
        if not source.is_file():
            results.append(
                SpeciesNormalizationResult(
                    show_id=show_id,
                    species_code=species_code,
                    status="ERROR",
                    source_path=source,
                    output_path=None,
                    detail="GFF file not found",
                )
            )
            continue
        try:
            detection = detect_gff_format(source)
        except GffFormatError as exc:
            results.append(
                SpeciesNormalizationResult(
                    show_id=show_id,
                    species_code=species_code,
                    status="ERROR",
                    source_path=source,
                    output_path=None,
                    detail=str(exc),
                )
            )
            continue
        if detection.format != "transcript_only":
            results.append(
                SpeciesNormalizationResult(
                    show_id=show_id,
                    species_code=species_code,
                    status="UNCHANGED",
                    source_path=source,
                    output_path=None,
                    detail=f"format={detection.format}",
                )
            )
            continue

        output_relative = Path("normalized") / f"{species_code}.transcript_only.normalized.gff.gz"
        output = directory / output_relative
        if output.exists() and not overwrite:
            if Path(row.get("gff", "")) != output_relative:
                row["gff"] = output_relative.as_posix()
                inventory_changed = True
            results.append(
                SpeciesNormalizationResult(
                    show_id=show_id,
                    species_code=species_code,
                    status="EXISTING",
                    source_path=source,
                    output_path=output,
                    detail="normalized GFF already exists",
                )
            )
            continue
        try:
            result = standardize_transcript_only_gff(source, output)
        except GffStandardizationError as exc:
            results.append(
                SpeciesNormalizationResult(
                    show_id=show_id,
                    species_code=species_code,
                    status="ERROR",
                    source_path=source,
                    output_path=None,
                    detail=str(exc),
                )
            )
            continue
        if result.converted:
            row["gff"] = output_relative.as_posix()
            converted += 1
            inventory_changed = True
        results.append(
            SpeciesNormalizationResult(
                show_id=show_id,
                species_code=species_code,
                status="CONVERTED" if result.converted else "UNCHANGED",
                source_path=source,
                output_path=result.output_path,
                detail=f"genes={result.gene_count} transcripts={result.transcript_count}",
            )
        )

    if inventory_changed and update_inventory:
        write_tsv(inventory_path, fields, inventory)
    return results
