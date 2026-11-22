"""Explicit GFF structure detection and adapters.

The canonical GFF3 parser in :mod:`pgci_preprocess.gff` only understands
nine-column gene/mRNA/CDS hierarchies.  Historical inputs used by this project
also appear as:

* simplified five-column GFF (``seqid start end strand gene_id``);
* transcript-only nine-column GFF (mRNA rows, no gene rows);
* gene-only nine-column GFF (one gene row per representative transcript);
* gene + CDS nine-column GFF (CDS parented directly to the gene).

This module turns those structural variants into the same ``GffData`` model
without import-time monkeypatching.  Every transformation is reported back to
the pipeline so it can be written into the per-species QC record.
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from pgci_preprocess.fasta import open_text
from pgci_preprocess.gff import (
    GENE_TYPES,
    ORGANELLE_TYPES,
    TRANSCRIPT_TYPES,
    parse_attributes,
    parse_gff,
    strip_prefix,
)
from pgci_preprocess.models import Anomaly, Gene, GffData, Region, Transcript

FORMATS = {
    "standard",
    "transcript_only",
    "simplified5",
    "gene_span",
    "gene_cds_direct",
}
CHROM_RE = re.compile(r"^[Cc]hr(?:omosome)?[_]?0*(\d+)\b")
#: Chromosome labels embedded in a prefixed seqid, e.g. ``Ajul.Chr01``,
#: ``scaffold_Chr3`` or ``sample|chromosome12``.  Real assemblies often ship
#: FASTA records whose seqid carries the chromosome label while the header
#: description is empty; strict chromosome mode must still recognise them.
CHROM_SEGMENT_RE = re.compile(r"(?:^|[._\s|])[Cc]hr(?:omosome)?[_]?0*(\d+)\b")
CDS_PREFIX_RE = re.compile(r"^cds[.\-]")


class GffFormatError(RuntimeError):
    pass


@dataclass(frozen=True, slots=True)
class GffDetection:
    format: str
    n_columns: int
    has_gene: bool
    has_transcript: bool
    has_cds: bool
    rows: int
    note: str = ""

    def as_dict(self) -> dict[str, object]:
        return {
            "format": self.format,
            "n_columns": self.n_columns,
            "has_gene": self.has_gene,
            "has_transcript": self.has_transcript,
            "has_cds": self.has_cds,
            "rows": self.rows,
            "note": self.note,
        }


@dataclass(slots=True)
class ParsedGff:
    data: GffData
    detection: GffDetection
    requested_format: str
    region_backfill: list[str] = field(default_factory=list)
    gwh_accessions: int = 0
    cds_prefixes_stripped: int = 0


def chromosome_from_seqid(seqid: str) -> str:
    """Return a chromosome label from common seqid spellings, else empty.

    Both bare (``Chr1``, ``chromosome_02``) and prefixed (``Ajul.Chr01``,
    ``sample|Chr3``) spellings are recognised.  Accessions such as
    ``NW_026555332.1`` and ``GWHEUBM00000024.1`` deliberately return empty:
    a scaffold accession is not chromosome-mounted evidence.
    """

    match = CHROM_RE.match(seqid) or CHROM_SEGMENT_RE.search(seqid)
    return match.group(1) if match else ""


def _rows(path: Path):
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if raw_line.startswith("#") or not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            yield line_number, raw_line.rstrip("\r\n"), fields


def detect_gff_format(path: Path | str) -> GffDetection:
    """Inspect the file once and classify its structure.

    Column count is authoritative: five columns is simplified GFF, nine
    columns is GFF3.  For nine-column input the presence of gene, transcript and
    CDS rows decides which adapter is required.  Filename suffixes are never
    used for classification.
    """

    path = Path(path)
    columns: Counter[int] = Counter()
    features: Counter[str] = Counter()
    rows = 0
    for _line_number, _raw, fields in _rows(path):
        rows += 1
        columns[len(fields)] += 1
        if len(fields) == 9:
            features[fields[2]] += 1
        elif len(fields) == 5:
            features["gene"] += 1

    if not rows:
        raise GffFormatError(f"{path.name}: no GFF data rows")
    if set(columns) == {5}:
        return GffDetection("simplified5", 5, True, False, False, rows, "five-column simplified GFF")
    if set(columns) == {9}:
        has_gene = any(feature in GENE_TYPES for feature in features)
        has_transcript = any(feature in TRANSCRIPT_TYPES for feature in features)
        has_cds = bool(features["CDS"])
        if has_gene and has_transcript:
            return GffDetection("standard", 9, True, True, has_cds, rows, "gene/mRNA hierarchy")
        if has_gene and has_cds and not has_transcript:
            return GffDetection(
                "gene_cds_direct",
                9,
                True,
                False,
                True,
                rows,
                "CDS parented directly to gene rows",
            )
        if has_gene and not has_transcript and not has_cds:
            return GffDetection("gene_span", 9, True, False, False, rows, "gene rows only")
        if has_transcript and not has_gene:
            return GffDetection(
                "transcript_only",
                9,
                False,
                True,
                has_cds,
                rows,
                "mRNA rows without a gene layer",
            )
        return GffDetection("standard", 9, has_gene, has_transcript, has_cds, rows, "unclassified nine-column GFF")
    raise GffFormatError(f"{path.name}: mixed or unsupported column counts {dict(columns)}")


def _region_from_length(seqid: str, length: int, source_order: int) -> Region:
    chromosome = chromosome_from_seqid(seqid)
    return Region(
        seqid=seqid,
        length=max(1, length),
        genome_type="genomic",
        chromosome=chromosome,
        taxon_id="",
        source_order=source_order,
        is_primary=bool(chromosome),
        synthetic=True,
    )


def parse_simplified5_gff(path: Path | str) -> GffData:
    """Parse ``seqid start end strand gene_id`` rows into gene/transcript pairs."""

    path = Path(path)
    genes: dict[str, Gene] = {}
    transcripts: dict[str, Transcript] = {}
    counts: Counter[str] = Counter()
    anomalies: list[Anomaly] = []
    max_end: dict[str, int] = {}
    skipped = 0

    for line_number, _raw, fields in _rows(path):
        if len(fields) != 5:
            skipped += 1
            continue
        seqid, start_text, end_text, strand, gene_id = fields
        gene_id = gene_id.strip()
        try:
            start, end = int(start_text), int(end_text)
        except ValueError:
            anomalies.append(
                Anomaly("warning", "invalid_coordinate", f"line={line_number} {start_text}-{end_text}")
            )
            continue
        if not gene_id or start < 1 or end < start:
            anomalies.append(Anomaly("warning", "invalid_simplified_row", f"line={line_number}"))
            continue
        if strand not in {"+", "-"}:
            anomalies.append(Anomaly("warning", "invalid_strand", f"{gene_id}:{strand}"))
            continue
        counts["gene"] += 1
        counts["mRNA"] += 1
        max_end[seqid] = max(max_end.get(seqid, 0), end)
        if gene_id in genes:
            anomalies.append(Anomaly("warning", "duplicate_gene", f"{gene_id} line={line_number}"))
            continue
        genes[gene_id] = Gene(
            source_id=gene_id,
            seqid=seqid,
            start=start,
            end=end,
            strand=strand,
            name=gene_id,
            biotype="protein_coding",
            feature_type="gene",
        )
        transcripts[gene_id] = Transcript(
            source_id=gene_id,
            gene_id=gene_id,
            seqid=seqid,
            start=start,
            end=end,
            strand=strand,
            feature_type="mRNA",
            transcript_id=gene_id,
            name=gene_id,
            protein_ids={gene_id},
            cds_count=1,
        )

    if skipped:
        anomalies.append(Anomaly("warning", "simplified_column_count", f"skipped={skipped}"))
    regions = {
        seqid: _region_from_length(seqid, length, order)
        for order, (seqid, length) in enumerate(max_end.items())
    }
    anomalies.append(
        Anomaly(
            "warning",
            "region_backfilled_from_simplified_gff",
            f"regions={len(regions)}",
        )
    )
    return GffData(
        regions=regions,
        genes=genes,
        transcripts=transcripts,
        feature_counts=counts,
        anomalies=anomalies,
        annotation_level="simplified5",
    )


def parse_transcript_only_gff(path: Path | str) -> GffData:
    """Parse transcript-only GFF3 into a synthetic one-transcript gene layer.

    The transcript ID doubles as the synthetic gene ID.  CDS ``protein_id`` /
    ``Protein_Accession`` values are authoritative; the transcript ID is only a
    fallback when the CDS rows provide no protein key.
    """

    path = Path(path)
    regions: dict[str, Region] = {}
    genes: dict[str, Gene] = {}
    transcripts: dict[str, Transcript] = {}
    counts: Counter[str] = Counter()
    anomalies: list[Anomaly] = []
    max_end: dict[str, int] = {}
    cds_counts: Counter[str] = Counter()
    exon_counts: Counter[str] = Counter()
    cds_protein_ids: dict[str, set[str]] = defaultdict(set)
    cds_translation_tables: dict[str, set[int]] = defaultdict(set)
    source_order = 0

    for line_number, _raw, fields in _rows(path):
        if len(fields) != 9:
            continue
        seqid, _source, feature, start_text, end_text, _score, strand, _phase, attr_text = fields
        counts[feature] += 1
        attributes = parse_attributes(attr_text)
        if feature == "region":
            try:
                start, end = int(start_text), int(end_text)
            except ValueError:
                continue
            length = end - start + 1
            if seqid not in regions:
                genome_type = attributes.get("genome", "genomic").strip().lower() or "genomic"
                chromosome = attributes.get("chromosome", "").strip()
                regions[seqid] = Region(
                    seqid=seqid,
                    length=length,
                    genome_type=genome_type,
                    chromosome=chromosome,
                    taxon_id="",
                    source_order=source_order,
                    is_primary=genome_type not in ORGANELLE_TYPES and bool(chromosome),
                )
                source_order += 1
            continue
        if feature in TRANSCRIPT_TYPES:
            transcript_id = strip_prefix(attributes.get("ID", "").strip(), "rna-")
            if not transcript_id:
                anomalies.append(
                    Anomaly("warning", "transcript_without_id", f"line={line_number}")
                )
                continue
            try:
                start, end = int(start_text), int(end_text)
            except ValueError:
                continue
            if start < 1 or end < start:
                continue
            if transcript_id in transcripts:
                anomalies.append(
                    Anomaly("warning", "duplicate_transcript", f"{transcript_id} line={line_number}")
                )
                continue
            gene = Gene(
                source_id=transcript_id,
                seqid=seqid,
                start=start,
                end=end,
                strand=strand,
                name=attributes.get("Name") or transcript_id,
                biotype=attributes.get("gene_biotype") or attributes.get("gbkey") or "protein_coding",
                feature_type="gene",
            )
            genes[transcript_id] = gene
            transcripts[transcript_id] = Transcript(
                source_id=transcript_id,
                gene_id=transcript_id,
                seqid=seqid,
                start=start,
                end=end,
                strand=strand,
                feature_type=feature,
                transcript_id=attributes.get("transcript_id", ""),
                name=attributes.get("Name", ""),
                protein_ids=set(),
            )
            max_end[seqid] = max(max_end.get(seqid, 0), end)
            counts["synthetic_gene"] += 1
        elif feature == "CDS":
            protein_id = (
                attributes.get("protein_id", "")
                or attributes.get("Protein_Accession", "")
            ).strip()
            raw_id = attributes.get("ID", "").strip()
            if not protein_id and raw_id.startswith("cds-"):
                protein_id = raw_id[len("cds-") :]
            raw_table = attributes.get("transl_table", "").strip()
            for parent in attributes.get("Parent", "").split(","):
                parent = strip_prefix(parent.strip(), "rna-")
                if parent:
                    cds_counts[parent] += 1
                    if protein_id:
                        cds_protein_ids[parent].add(protein_id)
                    if raw_table:
                        try:
                            cds_translation_tables[parent].add(int(raw_table))
                        except ValueError:
                            anomalies.append(
                                Anomaly(
                                    "warning",
                                    "invalid_translation_table",
                                    f"{parent}:{raw_table}",
                                )
                            )
        elif feature == "exon":
            for parent in attributes.get("Parent", "").split(","):
                parent = strip_prefix(parent.strip(), "rna-")
                if parent:
                    exon_counts[parent] += 1

    for transcript_id, transcript in transcripts.items():
        transcript.cds_count = cds_counts.get(transcript_id, 0)
        transcript.exon_count = exon_counts.get(transcript_id, 0)
        transcript.protein_ids = cds_protein_ids.get(transcript_id, set()) or {
            transcript_id
        }
        tables = cds_translation_tables.get(transcript_id, set())
        if len(tables) > 1:
            anomalies.append(
                Anomaly(
                    "fatal",
                    "conflicting_translation_tables",
                    f"{transcript_id}:{sorted(tables)}",
                )
            )
        elif tables:
            transcript.translation_table = next(iter(tables))

    for seqid, length in max_end.items():
        if seqid not in regions:
            regions[seqid] = _region_from_length(seqid, length, len(regions))

    if not genes:
        raise GffFormatError(f"{path.name}: transcript-only GFF contains no supported transcript rows")
    anomalies.append(
        Anomaly(
            "warning",
            "synthetic_gene_layer",
            f"genes={len(genes)} source=transcript-only",
        )
    )
    return GffData(
        regions=regions,
        genes=genes,
        transcripts=transcripts,
        feature_counts=counts,
        anomalies=anomalies,
        annotation_level="transcript_only",
    )


def backfill_regions_from_transcripts(data: GffData, *, reason: str) -> list[str]:
    """Create missing regions from transcript coordinates.

    This is a structural fallback, not a QC relaxation: it is recorded as a
    warning on the GFF data and returned to the pipeline for the QC record.
    Existing region rows are never overwritten.
    """

    missing: list[str] = []
    max_end: dict[str, int] = defaultdict(int)
    for transcript in data.transcripts.values():
        max_end[transcript.seqid] = max(max_end[transcript.seqid], transcript.end)
        if transcript.seqid not in data.regions:
            missing.append(transcript.seqid)
    backfilled: list[str] = []
    for seqid in dict.fromkeys(missing):
        if seqid in data.regions:
            continue
        data.regions[seqid] = _region_from_length(
            seqid,
            max_end[seqid],
            len(data.regions),
        )
        backfilled.append(seqid)
    if backfilled:
        data.anomalies.append(
            Anomaly(
                "warning",
                "region_backfilled_from_transcript",
                f"reason={reason} seqids={','.join(sorted(backfilled)[:5])}",
            )
        )
    return backfilled


def mrna_accession_map(path: Path | str) -> dict[str, str]:
    """Map GWH-style mRNA IDs to the ``Accession`` used in CDS/PEP FASTA."""

    path = Path(path)
    mapping: dict[str, str] = {}
    for _line_number, _raw, fields in _rows(path):
        if len(fields) != 9 or fields[2] not in TRANSCRIPT_TYPES:
            continue
        attributes = parse_attributes(fields[8])
        transcript_id = strip_prefix(attributes.get("ID", "").strip(), "rna-")
        accession = attributes.get("Accession", "").strip()
        if transcript_id and accession and transcript_id != accession:
            mapping[transcript_id] = accession
    return mapping


def apply_gwh_accessions(data: GffData, mapping: dict[str, str]) -> int:
    """Fall back to the mRNA ``Accession`` only when no protein id is known.

    GWH GFF3 rows usually carry an authoritative ``Protein_Accession``
    (``GWHP...``) that matches the ``Protein=`` key extracted from the CDS and
    PEP FASTA headers.  Overwriting it with the mRNA accession (``GWHT...``)
    made every transcript fail ``sequence_missing`` (regression found by the
    Golden-109 isolated rebuild on s2224/s4290).  The mapping is therefore a
    fallback for annotations with no protein id at all, never a replacement.
    """

    changed = 0
    for transcript in data.transcripts.values():
        if transcript.protein_ids:
            continue
        accession = mapping.get(transcript.source_id)
        if accession:
            transcript.protein_ids = {accession}
            changed += 1
    if changed:
        data.anomalies.append(
            Anomaly("warning", "gwh_accession_mapped", f"transcripts={changed}")
        )
    return changed


def normalize_cds_protein_ids(data: GffData) -> int:
    """Strip legacy ``cds.``/``cds-`` prefixes from CDS-derived protein IDs."""

    changed = 0
    for transcript in data.transcripts.values():
        normalized = {CDS_PREFIX_RE.sub("", protein_id) for protein_id in transcript.protein_ids}
        normalized.discard("")
        if normalized != transcript.protein_ids:
            transcript.protein_ids = normalized
            changed += 1
    if changed:
        data.anomalies.append(
            Anomaly("warning", "cds_protein_id_prefix_stripped", f"transcripts={changed}")
        )
    return changed


def parse_with_adapter(
    path: Path | str,
    *,
    requested_format: str = "auto",
    merge_split_genes: bool = False,
) -> ParsedGff:
    """Parse any supported GFF structure according to an explicit policy."""

    detection = detect_gff_format(path)
    actual = detection.format
    if requested_format not in {"auto", *FORMATS}:
        raise GffFormatError(f"unknown GFF format policy: {requested_format}")
    selected = actual if requested_format == "auto" else requested_format
    if selected == "simplified5":
        data = parse_simplified5_gff(path)
    elif selected == "transcript_only":
        data = parse_transcript_only_gff(path)
    elif selected == "gene_span":
        data = parse_gff(
            path,
            merge_split_genes=merge_split_genes,
            gene_span_only=True,
        )
    elif selected == "gene_cds_direct":
        data = parse_gff(
            path,
            merge_split_genes=merge_split_genes,
            gene_cds_direct=True,
        )
    else:
        data = parse_gff(path, merge_split_genes=merge_split_genes)
    return ParsedGff(
        data=data,
        detection=detection,
        requested_format=requested_format,
    )
