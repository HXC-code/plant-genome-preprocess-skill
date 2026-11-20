"""Non-destructive GTF, GBFF and functional-annotation conversion.

The converters in this module standardize *content*, not merely filename
suffixes.  Raw downloads are never modified.  Outputs are staged in a temporary
sibling directory, hashed, and committed only after the conversion succeeds.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import shutil
import tempfile
from collections.abc import Iterable
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO
from urllib.parse import unquote

from pgci_preprocess.conversion_common import (
    ConversionError,
    default_output_prefix,
    gff_escape,
    sha256_file,
)
from pgci_preprocess.fasta import open_text
from pgci_preprocess.gff import parse_gff
from pgci_preprocess.gff_formats import detect_gff_format
from pgci_preprocess.version import __version__

# Derived from the package version; see ``version.py``.
CONVERTER_VERSION = __version__
FORMATS = {"auto", "gtf", "gbff", "functional-go"}
TRANSCRIPT_FEATURES = {
    "mRNA",
    "transcript",
    "pseudogenic_transcript",
    "lnc_RNA",
    "ncRNA",
    "tRNA",
    "rRNA",
    "snRNA",
    "snoRNA",
    "miRNA",
}
CHILD_FEATURES = {
    "exon",
    "CDS",
    "start_codon",
    "stop_codon",
    "UTR",
    "five_prime_UTR",
    "three_prime_UTR",
    "Selenocysteine",
}
GO_RE = re.compile(r"^GO:\d{7}$")
FUNCTIONAL_SOURCES = frozenset({"GO", "NR", "SP", "TR", "SWISSPROT", "TREMBL", "KEGG", "KOG", "PFAM"})


@dataclass(frozen=True, slots=True)
class ConversionResult:
    source_path: Path
    format: str
    output_dir: Path
    prefix: str
    outputs: dict[str, Path]
    manifest_path: Path
    counts: dict[str, int]
    warnings: tuple[str, ...]
    complete: bool = True


@dataclass(slots=True)
class _GeneSpan:
    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    source: str
    name: str = ""
    biotype: str = ""
    locus_tag: str = ""
    dbxref: str = ""


@dataclass(slots=True)
class _TranscriptSpan:
    transcript_id: str
    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    source: str
    feature_type: str = "transcript"
    biotype: str = ""
    has_cds: bool = False
    protein_ids: set[str] = field(default_factory=set)
    translation_tables: set[int] = field(default_factory=set)


def parse_gtf_attributes(text: str) -> dict[str, str]:
    """Parse GTF ``key \"value\"`` and GFF3 ``key=value`` attributes."""

    attributes: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item and (" " not in item.split("=", 1)[0]):
            key, value = item.split("=", 1)
            attributes[key.strip()] = unquote(value.strip().strip('"'))
            continue
        key, _separator, value = item.partition(" ")
        if key:
            attributes[key] = value.strip().strip('"')
    return attributes


def _gff_attributes(values: Iterable[tuple[str, object]]) -> str:
    return ";".join(
        f"{key}={gff_escape(value)}"
        for key, value in values
        if value is not None and str(value) != ""
    )


def _row(
    seqid: str,
    source: str,
    feature: str,
    start: int | str,
    end: int | str,
    score: str,
    strand: str,
    phase: str,
    attributes: str,
) -> str:
    return "\t".join(
        (seqid, source or "pgci", feature, str(start), str(end), score or ".", strand, phase or ".", attributes)
    )


def _update_gene(
    genes: dict[str, _GeneSpan],
    gene_id: str,
    *,
    seqid: str,
    start: int,
    end: int,
    strand: str,
    source: str,
    attributes: dict[str, str],
) -> _GeneSpan:
    current = genes.get(gene_id)
    if current is None:
        current = _GeneSpan(
            gene_id=gene_id,
            seqid=seqid,
            start=start,
            end=end,
            strand=strand,
            source=source,
            name=attributes.get("gene") or attributes.get("Name") or gene_id,
            biotype=attributes.get("gene_biotype") or attributes.get("gbkey") or "",
            locus_tag=attributes.get("locus_tag", ""),
            dbxref=attributes.get("db_xref") or attributes.get("Dbxref") or "",
        )
        genes[gene_id] = current
        return current
    if (current.seqid, current.strand) != (seqid, strand):
        raise ConversionError(
            f"GTF gene {gene_id!r} occurs on incompatible locations: "
            f"{current.seqid}/{current.strand} and {seqid}/{strand}"
        )
    current.start = min(current.start, start)
    current.end = max(current.end, end)
    return current


def _update_transcript(
    transcripts: dict[str, _TranscriptSpan],
    transcript_id: str,
    gene_id: str,
    *,
    seqid: str,
    start: int,
    end: int,
    strand: str,
    source: str,
    feature_type: str,
    attributes: dict[str, str],
) -> _TranscriptSpan:
    current = transcripts.get(transcript_id)
    if current is None:
        current = _TranscriptSpan(
            transcript_id=transcript_id,
            gene_id=gene_id,
            seqid=seqid,
            start=start,
            end=end,
            strand=strand,
            source=source,
            feature_type=feature_type,
            biotype=attributes.get("transcript_biotype") or attributes.get("gbkey") or "",
        )
        transcripts[transcript_id] = current
        return current
    if (current.gene_id, current.seqid, current.strand) != (gene_id, seqid, strand):
        raise ConversionError(
            f"GTF transcript {transcript_id!r} has conflicting parents/locations"
        )
    current.start = min(current.start, start)
    current.end = max(current.end, end)
    if feature_type in TRANSCRIPT_FEATURES:
        current.feature_type = feature_type
    return current


def _scan_gtf(source: Path) -> tuple[dict[str, _GeneSpan], dict[str, _TranscriptSpan], Counter[str]]:
    genes: dict[str, _GeneSpan] = {}
    transcripts: dict[str, _TranscriptSpan] = {}
    counts: Counter[str] = Counter()
    with open_text(source) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if raw_line.startswith("#") or not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                counts["malformed_rows"] += 1
                continue
            seqid, source_name, feature, start_text, end_text, _score, strand, _phase, attr_text = fields
            try:
                start, end = int(start_text), int(end_text)
            except ValueError:
                counts["malformed_rows"] += 1
                continue
            if start < 1 or end < start:
                counts["malformed_rows"] += 1
                continue
            attributes = parse_gtf_attributes(attr_text)
            gene_id = (attributes.get("gene_id") or attributes.get("locus_tag") or attributes.get("gene") or "").strip()
            transcript_id = (attributes.get("transcript_id") or "").strip()
            counts["input_rows"] += 1
            counts[f"input_{feature}"] += 1
            if feature == "gene" and not gene_id:
                counts["gene_missing_id"] += 1
                continue
            if feature in CHILD_FEATURES and (not gene_id or not transcript_id):
                counts["child_missing_parent_id"] += 1
                continue
            if not gene_id:
                counts["rows_without_gene_id"] += 1
                continue
            gene = _update_gene(
                genes,
                gene_id,
                seqid=seqid,
                start=start,
                end=end,
                strand=strand,
                source=source_name,
                attributes=attributes,
            )
            if feature == "gene":
                gene.biotype = attributes.get("gene_biotype") or gene.biotype
            if not transcript_id:
                continue
            transcript = _update_transcript(
                transcripts,
                transcript_id,
                gene_id,
                seqid=seqid,
                start=start,
                end=end,
                strand=strand,
                source=source_name,
                feature_type=feature if feature in TRANSCRIPT_FEATURES else "transcript",
                attributes=attributes,
            )
            if feature == "CDS":
                transcript.has_cds = True
                protein_id = attributes.get("protein_id", "").strip()
                if protein_id:
                    transcript.protein_ids.add(protein_id)
                else:
                    counts["cds_without_protein_id"] += 1
                raw_table = attributes.get("transl_table", "").strip()
                if raw_table:
                    try:
                        transcript.translation_tables.add(int(raw_table))
                    except ValueError:
                        counts["invalid_translation_table"] += 1
    return genes, transcripts, counts


def convert_gtf_to_gff3(
    source_path: Path | str,
    output_path: Path | str,
) -> tuple[dict[str, int], list[str]]:
    """Convert NCBI/EMBL-style GTF into a gene/mRNA/exon/CDS GFF3 hierarchy.

    The input is read twice so exon/CDS rows are not retained in memory.  GTF
    has no authoritative sequence lengths, so no fake ``region`` lengths are
    emitted; the main pipeline backfills them from the matching genome FASTA.
    """

    source, output = Path(source_path), Path(output_path)
    genes, transcripts, counts = _scan_gtf(source)
    if not genes or not transcripts or not counts["input_CDS"]:
        raise ConversionError(
            f"{source.name}: GTF requires gene/transcript structure and CDS rows; "
            f"genes={len(genes)} transcripts={len(transcripts)} cds={counts['input_CDS']}"
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    exon_counter = 0
    child_counter = 0
    with gzip.open(output, "wt", encoding="utf-8", newline="\n") as destination:
        destination.write("##gff-version 3\n")
        destination.write(f"# pgci-converter: GTF to GFF3; source={source.name}\n")
        destination.write("# sequence lengths are intentionally obtained from the genome FASTA\n")
        ordered_genes = sorted(genes.values(), key=lambda item: (item.seqid, item.start, item.end, item.gene_id))
        for gene in ordered_genes:
            destination.write(
                _row(
                    gene.seqid,
                    gene.source,
                    "gene",
                    gene.start,
                    gene.end,
                    ".",
                    gene.strand,
                    ".",
                    _gff_attributes(
                        [
                            ("ID", f"gene-{gene.gene_id}"),
                            ("Name", gene.name),
                            ("gene_biotype", gene.biotype),
                            ("locus_tag", gene.locus_tag),
                            ("Dbxref", gene.dbxref),
                        ]
                    ),
                )
                + "\n"
            )
        ordered_transcripts = sorted(
            transcripts.values(),
            key=lambda item: (item.seqid, item.start, item.end, item.transcript_id),
        )
        for transcript in ordered_transcripts:
            feature_type = "mRNA" if transcript.has_cds else transcript.feature_type
            destination.write(
                _row(
                    transcript.seqid,
                    transcript.source,
                    feature_type,
                    transcript.start,
                    transcript.end,
                    ".",
                    transcript.strand,
                    ".",
                    _gff_attributes(
                        [
                            ("ID", f"rna-{transcript.transcript_id}"),
                            ("Parent", f"gene-{transcript.gene_id}"),
                            ("Name", transcript.transcript_id),
                            ("transcript_biotype", transcript.biotype),
                        ]
                    ),
                )
                + "\n"
            )

        with open_text(source) as handle:
            for raw_line in handle:
                if raw_line.startswith("#") or not raw_line.strip():
                    continue
                fields = raw_line.rstrip("\r\n").split("\t")
                if len(fields) != 9:
                    continue
                seqid, source_name, feature, start, end, score, strand, phase, attr_text = fields
                if feature not in CHILD_FEATURES:
                    continue
                attributes = parse_gtf_attributes(attr_text)
                gene_id = (attributes.get("gene_id") or attributes.get("locus_tag") or attributes.get("gene") or "").strip()
                transcript_id = (attributes.get("transcript_id") or "").strip()
                if not gene_id or transcript_id not in transcripts:
                    continue
                parent = f"rna-{transcript_id}"
                if feature == "exon":
                    exon_counter += 1
                    number = attributes.get("exon_number") or str(exon_counter)
                    feature_id = f"exon-{transcript_id}-{number}-{exon_counter}"
                    counts["exons"] += 1
                elif feature == "CDS":
                    protein_id = attributes.get("protein_id", "").strip()
                    feature_id = f"cds-{protein_id or transcript_id}"
                    counts["cds_segments"] += 1
                else:
                    child_counter += 1
                    feature_id = f"{feature}-{transcript_id}-{child_counter}"
                    counts["other_child_features"] += 1
                output_attributes = [
                    ("ID", feature_id),
                    ("Parent", parent),
                    ("protein_id", attributes.get("protein_id", "") if feature == "CDS" else ""),
                    ("transl_table", attributes.get("transl_table", "") if feature == "CDS" else ""),
                    ("product", attributes.get("product", "")),
                    ("Dbxref", attributes.get("db_xref") or attributes.get("Dbxref") or ""),
                    ("exon_number", attributes.get("exon_number", "")),
                ]
                destination.write(
                    _row(
                        seqid,
                        source_name,
                        feature,
                        start,
                        end,
                        score,
                        strand,
                        phase,
                        _gff_attributes(output_attributes),
                    )
                    + "\n"
                )

    if counts["input_CDS"] and not counts["cds_segments"]:
        raise ConversionError(
            f"{source.name}: GTF contained {counts['input_CDS']} CDS rows but none could "
            "be attached to a transcript; refusing to emit a CDS-free GFF3"
        )

    detection = detect_gff_format(output)
    parsed = parse_gff(output)
    fatal = [item for item in parsed.anomalies if item.level == "fatal"]
    if fatal:
        raise ConversionError(f"converted GFF3 failed validation: {fatal[0].code}: {fatal[0].detail}")
    if detection.format != "standard" or len(parsed.transcripts) != len(transcripts):
        raise ConversionError("converted GFF3 hierarchy could not be read back consistently")

    counts["genes"] = len(genes)
    counts["transcripts"] = len(transcripts)
    counts["transcripts_with_protein_id"] = sum(bool(item.protein_ids) for item in transcripts.values())
    counts["sequences"] = len({item.seqid for item in transcripts.values()})
    warnings: list[str] = []
    for key, label in (
        ("malformed_rows", "malformed GTF rows skipped"),
        ("child_missing_parent_id", "child rows missing gene_id/transcript_id skipped"),
        ("cds_without_protein_id", "CDS rows without protein_id; FASTA-key matching will be required"),
        ("invalid_translation_table", "CDS rows with a non-integer transl_table; table 1 was used"),
    ):
        if counts[key]:
            warnings.append(f"{label}: {counts[key]}")
    warnings.append("GTF contains no authoritative sequence lengths; genome FASTA is required for region backfill")
    return dict(sorted(counts.items())), warnings


def _functional_records(handle: TextIO):
    """Join physical lines belonging to one unquoted seven-column record.

    The source tables contain occasional literal newlines inside descriptions.
    They are not CSV-quoted, so a record is complete only when its accumulated
    tab fields reach the seven-column schema.
    """

    buffer: list[str] | None = None
    start_line = 0
    physical_lines = 0
    for line_number, raw_line in enumerate(handle, 1):
        if buffer is None and not raw_line.strip():
            continue
        parts = raw_line.rstrip("\r\n").split("\t")
        if buffer is None:
            buffer = parts
            start_line = line_number
            physical_lines = 0
        else:
            buffer[-1] += "\n" + parts[0]
            buffer.extend(parts[1:])
        physical_lines += 1
        if len(buffer) == 7:
            yield start_line, buffer, physical_lines
            buffer = None
        elif len(buffer) > 7:
            raise ConversionError(
                f"functional annotation record starting at line {start_line} has more than 7 tab-separated columns"
            )
    if buffer is not None:
        raise ConversionError(
            f"functional annotation record starting at line {start_line} is truncated; expected 7 columns, found {len(buffer)}"
        )


def normalize_functional_go_table(
    source_path: Path | str,
    output_path: Path | str,
) -> tuple[dict[str, int], list[str]]:
    """Normalize the project's headerless seven-column GO table.

    These files are named ``*.annotations.csv.gz`` but are actually tab
    separated.  They contain functional annotations, not genomic coordinates,
    so conversion to GFF would be scientifically invalid.
    """

    source, output = Path(source_path), Path(output_path)
    counts: Counter[str] = Counter()
    features: set[str] = set()
    terms: set[str] = set()
    species: set[str] = set()
    sources_seen: Counter[str] = Counter()
    output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(source) as handle, gzip.open(output, "wt", encoding="utf-8", newline="\n") as destination:
        destination.write("feature_id\tannotation_id\tevalue\tscore\tdescription\tspecies\tsource\n")
        for line_number, fields, physical_lines in _functional_records(handle):
            feature_id, annotation_id, evalue, score, description, species_name, source_name = fields
            try:
                float(evalue)
                float(score)
            except ValueError as exc:
                raise ConversionError(
                    f"{source.name}:{line_number}: evalue and score must be numeric"
                ) from exc
            if (
                not feature_id.strip()
                or not annotation_id.strip()
                or not species_name.strip()
                or not source_name.strip()
            ):
                raise ConversionError(
                    f"{source.name}:{line_number}: unrecognized functional annotation schema"
                )
            escaped_fields = [
                value.replace("\\", "\\\\").replace("\t", "\\t").replace("\r", "\\r").replace("\n", "\\n")
                for value in fields
            ]
            destination.write("\t".join(escaped_fields) + "\n")
            counts["rows"] += 1
            counts["physical_lines"] += physical_lines
            if physical_lines > 1:
                counts["continuation_rows"] += 1
            features.add(feature_id)
            terms.add(annotation_id)
            species.add(species_name)
            sources_seen[source_name.strip()] += 1
    if not counts["rows"]:
        raise ConversionError(f"{source.name}: no functional annotation rows")
    counts["features"] = len(features)
    counts["terms"] = len(terms)
    counts["species"] = len(species)
    counts["sources"] = len(sources_seen)
    for source_name, number in sorted(sources_seen.items()):
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", source_name)
        counts[f"source_{safe_name}"] = number
    warnings = [
        "functional annotation table normalized; it is metadata, not WGDI coordinate input",
        "annotation sources: " + ", ".join(f"{key}={value}" for key, value in sorted(sources_seen.items())),
    ]
    if counts["continuation_rows"]:
        warnings.append(
            f"literal newlines inside {counts['continuation_rows']} descriptions were preserved as escaped \\n text"
        )
    return dict(sorted(counts.items())), warnings


def _looks_functional_row(fields: list[str]) -> bool:
    if len(fields) < 5 or len(fields) > 7:
        return False
    if not fields[0].strip() or not fields[1].strip():
        return False
    try:
        float(fields[2])
        float(fields[3])
    except (ValueError, IndexError):
        return False
    annotation_id = fields[1].strip()
    source_name = fields[6].strip().upper() if len(fields) == 7 else ""
    if len(fields) == 7 and source_name in FUNCTIONAL_SOURCES:
        return True
    return bool(
        GO_RE.fullmatch(annotation_id)
        or re.match(r"^(?:[NXW]P_|YP_|AP_|A0A|Q[A-Z0-9])", annotation_id)
    )


def detect_conversion_format(source_path: Path | str) -> str:
    source = Path(source_path)
    if not source.is_file() or source.stat().st_size == 0:
        raise ConversionError(f"input file not found or empty: {source}")
    sampled_data_rows = 0
    try:
        with open_text(source) as handle:
            for raw_line in handle:
                if not raw_line.strip():
                    continue
                if raw_line.lstrip().startswith("LOCUS"):
                    return "gbff"
                if raw_line.startswith("#"):
                    continue
                fields = raw_line.rstrip("\r\n").split("\t")
                if _looks_functional_row(fields):
                    return "functional-go"
                if len(fields) == 9:
                    sampled_data_rows += 1
                    if re.search(
                        r"(?:^|;)\s*(?:gene_id|transcript_id|locus_tag|gene|gene_name|"
                        r"transcript_name|exon_number)\s+\"",
                        fields[8],
                    ):
                        return "gtf"
                    if sampled_data_rows >= 64:
                        break
    except (OSError, EOFError, ValueError) as exc:
        raise ConversionError(
            f"{source.name}: compressed input cannot be read: {exc}"
        ) from exc
    raise ConversionError(
        f"{source.name}: unrecognized input; supported content types are GTF, GenBank/GBFF, "
        "and the seven-column functional GO table"
    )


def _expected_names(fmt: str, prefix: str, include_genome: bool) -> dict[str, str]:
    if fmt == "gtf":
        return {"gff": f"{prefix}.gff3.gz"}
    if fmt == "gbff":
        outputs = {
            "gff": f"{prefix}.gff3.gz",
            "cds": f"{prefix}.cds.fa.gz",
            "pep": f"{prefix}.pep.fa.gz",
        }
        if include_genome:
            outputs["genome"] = f"{prefix}.genome.fa.gz"
        return outputs
    if fmt == "functional-go":
        return {"functional_annotations": f"{prefix}.functional.tsv.gz"}
    raise ConversionError(f"unsupported conversion format: {fmt}")


def convert_annotation(
    source_path: Path | str,
    output_dir: Path | str,
    *,
    fmt: str = "auto",
    prefix: str = "",
    overwrite: bool = False,
    include_genome: bool = True,
    max_records: int | None = None,
) -> ConversionResult:
    """Convert one annotation input and write a checksummed manifest.

    ``max_records`` is a smoke-test option for GBFF only.  Its outputs are
    explicitly marked ``complete=false`` and must not be released as a complete
    species bundle.
    """

    source = Path(source_path).resolve()
    destination_dir = Path(output_dir).resolve()
    if fmt not in FORMATS:
        raise ConversionError(f"unsupported format {fmt!r}; choose from {sorted(FORMATS)}")
    detected = detect_conversion_format(source) if fmt == "auto" else fmt
    if detected == fmt == "auto":  # defensive; auto is resolved above
        raise ConversionError("could not resolve input format")
    if max_records is not None and (detected != "gbff" or max_records < 1):
        raise ConversionError("--max-records is only valid for GBFF and must be positive")
    output_prefix = prefix.strip() or default_output_prefix(source)
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", output_prefix):
        raise ConversionError("output prefix may contain only letters, digits, dot, underscore and hyphen")
    expected = _expected_names(detected, output_prefix, include_genome)
    manifest_name = f"{output_prefix}.conversion.json"
    final_paths = {role: destination_dir / name for role, name in expected.items()}
    manifest_path = destination_dir / manifest_name
    existing = [path for path in [*final_paths.values(), manifest_path] if path.exists()]
    if existing and not overwrite:
        raise ConversionError(
            "output already exists (default is no overwrite): " + ", ".join(str(path) for path in existing)
        )
    if destination_dir.exists() and not destination_dir.is_dir():
        raise ConversionError(f"output directory is a file: {destination_dir}")
    destination_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".pgci-convert-", dir=destination_dir.parent))
    committed = False
    try:
        if detected == "gtf":
            counts, warnings = convert_gtf_to_gff3(source, staging / expected["gff"])
        elif detected == "functional-go":
            counts, warnings = normalize_functional_go_table(source, staging / expected["functional_annotations"])
        else:
            from pgci_preprocess.gbff_convert import convert_gbff_to_bundle

            counts, warnings = convert_gbff_to_bundle(
                source,
                staging,
                output_prefix,
                include_genome=include_genome,
                max_records=max_records,
            )
        for name in expected.values():
            path = staging / name
            if not path.is_file() or path.stat().st_size == 0:
                raise ConversionError(f"converter did not produce a non-empty output: {name}")
        output_metadata = {
            role: {
                "file": name,
                "bytes": (staging / name).stat().st_size,
                "sha256": sha256_file(staging / name),
            }
            for role, name in sorted(expected.items())
        }
        manifest = {
            "converter": "pgci_preprocess.convert-annotation",
            "converter_version": CONVERTER_VERSION,
            "source": {
                "path": str(source),
                "bytes": source.stat().st_size,
                "sha256": sha256_file(source),
            },
            "format": detected,
            "prefix": output_prefix,
            "complete": max_records is None,
            "max_records": max_records,
            "counts": counts,
            "warnings": warnings,
            "outputs": output_metadata,
        }
        staging_manifest = staging / manifest_name
        staging_manifest.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        destination_dir.mkdir(parents=True, exist_ok=True)
        for name in [*expected.values(), manifest_name]:
            os.replace(staging / name, destination_dir / name)
        committed = True
    finally:
        if not committed:
            shutil.rmtree(staging, ignore_errors=True)
        elif staging.exists():
            shutil.rmtree(staging, ignore_errors=True)
    return ConversionResult(
        source_path=source,
        format=detected,
        output_dir=destination_dir,
        prefix=output_prefix,
        outputs=final_paths,
        manifest_path=manifest_path,
        counts=counts,
        warnings=tuple(warnings),
        complete=max_records is None,
    )
