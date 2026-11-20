from __future__ import annotations

import re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

from pgci_preprocess.fasta import open_text
from pgci_preprocess.models import Anomaly, Gene, GffData, Region, Transcript

REGION_TYPES = {"region"}
GENE_TYPES = {"gene", "pseudogene"}
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
ORGANELLE_TYPES = {"chloroplast", "mitochondrion", "plastid", "apicoplast"}
UNKNOWN_VALUES = {"", ".", "na", "n/a", "nan", "none", "null", "unknown"}
TAXON_RE = re.compile(r"(?:^|,)taxon:(\d+)")


def parse_attributes(field: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in field.split(";"):
        if not item:
            continue
        key, separator, value = item.partition("=")
        key = key.strip()
        if not key:
            continue
        attributes[key] = unquote(value.strip()) if separator else ""
    return attributes


def extract_taxon(attributes: dict[str, str]) -> str:
    match = TAXON_RE.search(attributes.get("Dbxref", ""))
    return match.group(1) if match else ""


def strip_prefix(value: str, prefix: str) -> str:
    return value[len(prefix) :] if value.startswith(prefix) else value


def _signed_int(value: str, line_number: int, anomalies: list[Anomaly]) -> int | None:
    try:
        return int(value)
    except ValueError:
        anomalies.append(
            Anomaly("warning", "invalid_coordinate", f"line={line_number} value={value!r}")
        )
        return None


def _register_region(
    region: Region,
    regions: dict[str, Region],
    anomalies: list[Anomaly],
) -> None:
    previous = regions.get(region.seqid)
    if previous and previous.length != region.length:
        anomalies.append(
            Anomaly(
                "fatal",
                "conflicting_region_length",
                f"{region.seqid}: {previous.length} != {region.length}",
            )
        )
        return
    if previous:
        if (
            previous.genome_type == "genomic"
            and not previous.chromosome
            and not previous.taxon_id
            and not previous.is_primary
        ):
            previous.genome_type = region.genome_type
            previous.chromosome = region.chromosome
            previous.taxon_id = region.taxon_id
            previous.is_primary = region.is_primary
            return
        anomalies.append(
            Anomaly("warning", "duplicate_region", f"duplicate region row for {region.seqid}")
        )
        return
    regions[region.seqid] = region


def parse_gff(
    path: Path | str,
    *,
    merge_split_genes: bool = False,
    gene_span_only: bool = False,
    gene_cds_direct: bool = False,
) -> GffData:
    """Parse a GFF into genes/transcripts.

    merge_split_genes: when the same gene id occurs on the same seqid and strand
      more than once (NCBI segmented gene rows), merge into one gene spanning
      min(start)..max(end) and record the merge count. True conflicts (different
      seqid or strand) still fail.
    gene_span_only: for annotations with gene rows only (no transcript/exon
      structure), synthesize one transcript per gene whose protein id equals the
      gene id, so CDS/PEP keyed by gene id can be validated.
    gene_cds_direct: for annotations that attach CDS rows straight to a gene
      (Parent=gene-...) without any mRNA layer, synthesize one transcript per
      gene. The protein id is taken from the CDS rows themselves, never invented.
    """
    path = Path(path)
    regions: dict[str, Region] = {}
    genes: dict[str, Gene] = {}
    transcripts: dict[str, Transcript] = {}
    transcript_links: dict[str, set[str]] = defaultdict(set)
    transcript_cds_counts: Counter[str] = Counter()
    transcript_exon_counts: Counter[str] = Counter()
    transcript_translation_tables: dict[str, set[int]] = defaultdict(set)
    counts: Counter[str] = Counter()
    anomalies: list[Anomaly] = []
    assembly_build = ""
    annotation_source = ""
    source_order = 0

    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if raw_line.startswith("#"):
                if raw_line.startswith("##sequence-region"):
                    parts = raw_line.rstrip("\n").split()
                    if len(parts) >= 4:
                        try:
                            length = int(parts[3])
                        except ValueError:
                            anomalies.append(
                                Anomaly(
                                    "warning",
                                    "invalid_sequence_region",
                                    f"line={line_number} {raw_line.strip()}",
                                )
                            )
                            continue
                        _register_region(
                            Region(
                                seqid=parts[1],
                                length=length,
                                genome_type="genomic",
                                chromosome="",
                                taxon_id="",
                                source_order=source_order,
                                is_primary=False,
                            ),
                            regions,
                            anomalies,
                        )
                elif raw_line.startswith("#!genome-build "):
                    assembly_build = raw_line.split(None, 1)[1].strip()
                elif raw_line.startswith("#!annotation-source "):
                    annotation_source = raw_line.split(None, 1)[1].strip()
                continue

            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                counts["malformed_line"] += 1
                anomalies.append(
                    Anomaly("warning", "gff_column_count", f"line={line_number} columns={len(fields)}")
                )
                continue

            seqid, _source, feature_type, start_text, end_text, _score, strand, _phase, attr_text = fields
            counts[feature_type] += 1
            attributes = parse_attributes(attr_text)

            if feature_type in REGION_TYPES:
                start = _signed_int(start_text, line_number, anomalies)
                end = _signed_int(end_text, line_number, anomalies)
                if start is None or end is None or start < 1 or end < start:
                    continue
                genome_type = attributes.get("genome", "genomic").strip().lower() or "genomic"
                chromosome = attributes.get("chromosome", "").strip()
                chromosome_known = chromosome.casefold() not in UNKNOWN_VALUES
                is_primary = (
                    genome_type not in ORGANELLE_TYPES
                    and (genome_type == "chromosome" or chromosome_known)
                )
                region = Region(
                    seqid=seqid,
                    length=end - start + 1,
                    genome_type=genome_type,
                    chromosome=chromosome,
                    taxon_id=extract_taxon(attributes),
                    source_order=source_order,
                    is_primary=is_primary,
                )
                _register_region(region, regions, anomalies)
                source_order += 1
                continue

            if feature_type in GENE_TYPES:
                gene_id = strip_prefix(attributes.get("ID", ""), "gene-")
                if not gene_id:
                    anomalies.append(
                        Anomaly("warning", "gene_without_id", f"{seqid}:{start_text}-{end_text}")
                    )
                    continue
                start = _signed_int(start_text, line_number, anomalies)
                end = _signed_int(end_text, line_number, anomalies)
                if start is None or end is None or start < 1 or end < start:
                    continue
                biotype = (
                    attributes.get("gene_biotype")
                    or attributes.get("gbkey")
                    or ("pseudogene" if feature_type == "pseudogene" else "unknown")
                )
                name = attributes.get("Name") or attributes.get("gene") or gene_id
                gene = Gene(
                    source_id=gene_id,
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand=strand,
                    name=name,
                    biotype=biotype,
                    feature_type=feature_type,
                )
                previous_gene = genes.get(gene_id)
                if previous_gene:
                    previous_region = regions.get(previous_gene.seqid)
                    current_region = regions.get(gene.seqid)
                    both_organelle = bool(
                        previous_region
                        and current_region
                        and previous_region.is_organelle
                        and current_region.is_organelle
                    )
                    same_span = (
                        previous_gene.seqid, previous_gene.start, previous_gene.end
                    ) == (gene.seqid, gene.start, gene.end)
                    same_location = (
                        previous_gene.seqid == gene.seqid
                        and previous_gene.strand == gene.strand
                    )
                    if merge_split_genes and same_location and not same_span:
                        # NCBI splits one gene into multiple same-strand intervals.
                        merged_start = min(previous_gene.start, gene.start)
                        merged_end = max(previous_gene.end, gene.end)
                        genes[gene_id] = Gene(
                            source_id=previous_gene.source_id,
                            seqid=previous_gene.seqid,
                            start=merged_start,
                            end=merged_end,
                            strand=previous_gene.strand,
                            name=previous_gene.name,
                            biotype=previous_gene.biotype,
                            feature_type=previous_gene.feature_type,
                        )
                        counts["split_gene_merged"] += 1
                        continue
                    level = "warning" if same_span or both_organelle else "fatal"
                    anomalies.append(
                        Anomaly(
                            level,
                            "duplicate_organelle_gene" if both_organelle else "duplicate_gene",
                            f"{gene_id}:{gene.seqid}:{gene.start}-{gene.end} line={line_number}",
                        )
                    )
                    continue
                genes[gene_id] = gene
                continue

            if feature_type in TRANSCRIPT_TYPES:
                transcript_id = strip_prefix(attributes.get("ID", ""), "rna-")
                parent = attributes.get("Parent", "").split(",", 1)[0].strip()
                gene_id = strip_prefix(parent, "gene-")
                if not transcript_id or not gene_id:
                    anomalies.append(
                        Anomaly(
                            "warning",
                            "transcript_without_id_or_parent",
                            f"{seqid}:{start_text}-{end_text}",
                        )
                    )
                    continue
                start = _signed_int(start_text, line_number, anomalies)
                end = _signed_int(end_text, line_number, anomalies)
                if start is None or end is None or start < 1 or end < start:
                    continue
                transcript_record = Transcript(
                    source_id=transcript_id,
                    gene_id=gene_id,
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand=strand,
                    feature_type=feature_type,
                    transcript_id=attributes.get("transcript_id", ""),
                    name=attributes.get("Name", ""),
                )
                previous_transcript = transcripts.get(transcript_id)
                if previous_transcript:
                    anomalies.append(
                        Anomaly("warning", "duplicate_transcript", f"{transcript_id} line={line_number}")
                    )
                    continue
                transcripts[transcript_id] = transcript_record
                continue

            if feature_type == "CDS":
                parents: list[str] = []
                for item in attributes.get("Parent", "").split(","):
                    item = item.strip()
                    if not item:
                        continue
                    item = strip_prefix(item, "rna-")
                    if gene_cds_direct:
                        item = strip_prefix(item, "gene-")
                    parents.append(item)
                protein_id = (
                    attributes.get("protein_id", "")
                    or attributes.get("Protein_Accession", "")
                )
                translation_table = 1
                raw_table = attributes.get("transl_table", "").strip()
                if raw_table:
                    try:
                        translation_table = int(raw_table)
                    except ValueError:
                        anomalies.append(
                            Anomaly(
                                "warning",
                                "invalid_translation_table",
                                f"line={line_number} value={raw_table!r}",
                            )
                        )
                if not protein_id:
                    raw_id = attributes.get("ID", "")
                    # only a real "cds-" style ID is a protein key; Phytozome CDS
                    # rows carry "<mRNA>.CDS.n" which must NOT be mistaken for one
                    if raw_id.startswith("cds-"):
                        protein_id = raw_id[len("cds-"):]
                for transcript_id in parents:
                    if protein_id:
                        transcript_links[transcript_id].add(protein_id)
                    else:
                        transcript_links[transcript_id]
                    transcript_cds_counts[transcript_id] += 1
                    transcript_translation_tables[transcript_id].add(translation_table)
                if not parents:
                    anomalies.append(
                        Anomaly("warning", "cds_without_parent", f"{seqid}:{start_text}-{end_text}")
                    )
                continue

            if feature_type == "exon":
                for parent in attributes.get("Parent", "").split(","):
                    parent = strip_prefix(parent.strip(), "rna-")
                    if parent:
                        transcript_exon_counts[parent] += 1

    if gene_span_only and not transcripts and genes:
        # gene-only annotation: one representative per gene, protein id == gene id
        for gene_id, gene in genes.items():
            transcripts[gene_id] = Transcript(
                source_id=gene_id,
                gene_id=gene_id,
                seqid=gene.seqid,
                start=gene.start,
                end=gene.end,
                strand=gene.strand,
                feature_type="gene_span",
                transcript_id="",
                protein_ids={gene_id},
            )
        counts["synthesized_gene_span_transcripts"] = len(transcripts)

    if gene_cds_direct and genes:
        # CDS rows hang directly off the gene: build one transcript per uncovered
        # gene from the CDS rows' own protein ids (no invented identifiers).
        covered = {transcript.gene_id for transcript in transcripts.values()}
        synthesized = 0
        for gene_id, gene in genes.items():
            if gene_id in covered:
                continue
            links = transcript_links.get(gene_id)
            if not links:
                continue
            transcripts[gene_id] = Transcript(
                source_id=gene_id,
                gene_id=gene_id,
                seqid=gene.seqid,
                start=gene.start,
                end=gene.end,
                strand=gene.strand,
                feature_type="gene_cds",
                transcript_id="",
                protein_ids=set(links),
                cds_count=transcript_cds_counts.get(gene_id, 0),
            )
            synthesized += 1
        counts["synthesized_gene_cds_transcripts"] = synthesized

    for transcript_id, protein_ids in transcript_links.items():
        linked_transcript = transcripts.get(transcript_id)
        if linked_transcript is None:
            anomalies.append(
                Anomaly(
                    "warning",
                    "cds_parent_missing_transcript",
                    f"transcript={transcript_id}",
                )
            )
            continue
        linked_transcript.protein_ids = {
            protein_id for protein_id in protein_ids if protein_id
        }
        linked_transcript.cds_count = transcript_cds_counts[transcript_id]
        linked_transcript.exon_count = transcript_exon_counts[transcript_id]
        tables = transcript_translation_tables.get(transcript_id, set())
        if len(tables) > 1:
            anomalies.append(
                Anomaly(
                    "fatal",
                    "conflicting_translation_tables",
                    f"{transcript_id}:{sorted(tables)}",
                )
            )
        elif tables:
            linked_transcript.translation_table = next(iter(tables))

    return GffData(
        regions=regions,
        genes=genes,
        transcripts=transcripts,
        feature_counts=counts,
        anomalies=anomalies,
        assembly_build=assembly_build,
        annotation_source=annotation_source,
        split_genes_merged=counts.get("split_gene_merged", 0),
        annotation_level=(
            "gene_cds_direct"
            if gene_cds_direct
            else "gene_span_only"
            if gene_span_only
            else "exon_structured"
        ),
    )

