"""GenBank/GBFF to WGDI-oriented coding annotation and FASTA bundle.

Biopython is intentionally imported lazily: GTF and functional-table conversion
remain standard-library only.  The output is a coding-sequence bundle, not a
lossless representation of every GenBank feature.
"""

from __future__ import annotations

import gzip
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from pgci_preprocess.conversion_common import (
    ConversionError,
    fasta_record,
    gff_escape,
)
from pgci_preprocess.fasta import open_text, read_cds_fasta, read_pep_fasta
from pgci_preprocess.gff import parse_gff
from pgci_preprocess.gff_formats import detect_gff_format
from pgci_preprocess.sequence import validate_cds_pep

GENE_FEATURES = {"gene", "pseudogene"}
TRANSCRIPT_FEATURES = {"mRNA", "transcript"}


@dataclass(slots=True)
class _Gene:
    key: str
    identity: str
    seqid: str
    start: int
    end: int
    strand: str
    name: str
    locus_tag: str = ""
    dbxref: str = ""
    merged_features: int = 1


@dataclass(slots=True)
class _Coding:
    protein_id: str
    transcript_id: str
    gene: _Gene
    mrna: Any | None
    cds_parts: list[tuple[int, int, int]]
    exon_parts: list[tuple[int, int, int]]
    nucleotide: str
    protein: str
    translation_table: int
    synthesized_transcript: bool = False


def _values(feature: Any, key: str) -> list[str]:
    return [str(value).strip() for value in feature.qualifiers.get(key, []) if str(value).strip()]


def _first(feature: Any, key: str) -> str:
    values = _values(feature, key)
    return values[0] if values else ""


def _gene_ids(feature: Any) -> list[str]:
    return [value.split(":", 1)[1] for value in _values(feature, "db_xref") if value.startswith("GeneID:")]


def _identity(feature: Any, fallback: str = "") -> str:
    locus = _first(feature, "locus_tag")
    if locus:
        return locus
    gene_ids = _gene_ids(feature)
    if gene_ids:
        return f"GeneID:{gene_ids[0]}"
    gene = _first(feature, "gene")
    return gene or fallback


def _strand(feature: Any) -> int:
    return int(feature.location.strand or 1)


def _span(feature: Any) -> tuple[int, int]:
    return int(feature.location.start), int(feature.location.end)


def _contains(outer: Any, inner: Any) -> bool:
    outer_start, outer_end = _span(outer)
    inner_start, inner_end = _span(inner)
    return outer_start <= inner_start and inner_end <= outer_end


def _parts(feature: Any) -> list[tuple[int, int, int]]:
    location = feature.location
    raw_parts = getattr(location, "parts", [location])
    return [(int(part.start), int(part.end), int(part.strand or _strand(feature))) for part in raw_parts]


def _compatible_parts(cds: Any, mrna: Any) -> bool:
    mrna_parts = _parts(mrna)
    for start, end, strand in _parts(cds):
        if not any(
            outer_start <= start and end <= outer_end and outer_strand == strand
            for outer_start, outer_end, outer_strand in mrna_parts
        ):
            return False
    return True


def _build_genes(record: Any, counts: Counter[str]) -> dict[str, _Gene]:
    genes: dict[str, _Gene] = {}
    for feature in record.features:
        if feature.type not in GENE_FEATURES:
            continue
        identity = _identity(feature)
        if not identity:
            counts["gene_without_identity"] += 1
            continue
        start, end = _span(feature)
        strand = _strand(feature)
        key = f"{identity}|{strand}"
        current = genes.get(key)
        if current is None:
            genes[key] = _Gene(
                key=key,
                identity=identity,
                seqid=record.id,
                start=start,
                end=end,
                strand="+" if strand > 0 else "-",
                name=_first(feature, "gene") or identity,
                locus_tag=_first(feature, "locus_tag"),
                dbxref=",".join(_values(feature, "db_xref")),
            )
        else:
            current.start = min(current.start, start)
            current.end = max(current.end, end)
            current.merged_features += 1
            counts["gene_features_merged"] += 1
    return genes


def _product_variant(feature: Any) -> str:
    product = _first(feature, "product")
    match = re.search(r"\b(?:isoform|transcript variant)\s+([A-Za-z0-9_.-]+)", product, re.IGNORECASE)
    return match.group(1).casefold() if match else ""


def _match_mrna(
    cds: Any,
    gene: _Gene,
    mrnas: list[tuple[str, Any]],
    excluded_mrna_ids: set[int] | None = None,
) -> Any | None:
    excluded = excluded_mrna_ids or set()
    explicit = _first(cds, "transcript_id")
    candidates = [
        feature
        for key, feature in mrnas
        if key == gene.key
        and id(feature) not in excluded
        and _strand(feature) == _strand(cds)
        and _contains(feature, cds)
    ]
    if explicit:
        explicit_matches = [
            feature
            for _key, feature in mrnas
            if _first(feature, "transcript_id") == explicit
            and id(feature) not in excluded
            and _strand(feature) == _strand(cds)
            and _contains(feature, cds)
        ]
        if len(explicit_matches) == 1:
            return explicit_matches[0]
        explicit_matches = [feature for feature in candidates if _first(feature, "transcript_id") == explicit]
        if len(explicit_matches) == 1:
            return explicit_matches[0]
    compatible = [feature for feature in candidates if _compatible_parts(cds, feature)]
    variant = _product_variant(cds)
    if variant:
        variant_matches = [feature for feature in compatible if _product_variant(feature) == variant]
        if len(variant_matches) == 1:
            return variant_matches[0]
        if len(variant_matches) > 1:
            return None
        if any(_product_variant(feature) for feature in compatible):
            return None
    if len(compatible) == 1:
        return compatible[0]
    return None


def _adjusted_cds_parts(feature: Any, offset: int) -> list[tuple[int, int, int]]:
    parts = _parts(feature)
    if not offset:
        return parts
    start, end, strand = parts[0]
    if strand > 0:
        start += offset
    else:
        end -= offset
    if end <= start:
        return []
    parts[0] = (start, end, strand)
    return parts


def _gff_row(
    seqid: str,
    feature: str,
    start: int,
    end: int,
    strand: str,
    attributes: list[tuple[str, object]],
    phase: str = ".",
) -> str:
    attr = ";".join(
        f"{key}={gff_escape(value)}"
        for key, value in attributes
        if value is not None and str(value) != ""
    )
    return f"{seqid}\tGenBank\t{feature}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attr}\n"


def _unique_gene_id(
    gene: _Gene,
    record_id: str,
    used: dict[str, tuple[str, str, int]],
    counts: Counter[str],
) -> str:
    internal = (record_id, gene.key, gene.start)
    previous = used.get(gene.identity)
    if previous is None:
        used[gene.identity] = internal
        return gene.identity
    if previous[0] == record_id and previous[1] == gene.key:
        return gene.identity
    candidate = f"{gene.identity}__{record_id}"
    if candidate != gene.identity:
        counts["gene_ids_disambiguated"] += 1
    used[candidate] = internal
    return candidate


def convert_gbff_to_bundle(
    source_path: Path | str,
    output_dir: Path | str,
    prefix: str,
    *,
    include_genome: bool = True,
    max_records: int | None = None,
    strict_translation: bool = True,
) -> tuple[dict[str, int], list[str]]:
    """Extract coding GFF3, CDS/PEP FASTA and optionally genome FASTA.

    ``strict_translation=False`` skips a CDS whose extracted sequence and
    ``/translation`` qualifier disagree (real NCBI files contain a few such
    records) instead of aborting the whole file; skipped records are counted
    and reported as warnings.
    """

    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise ConversionError(
            "GBFF conversion requires Biopython; install with "
            "`python -m pip install \"pgci-preprocess[gbff]\"` or `python -m pip install biopython>=1.83`"
        ) from exc

    source = Path(source_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    gff_path = output_dir / f"{prefix}.gff3.gz"
    cds_path = output_dir / f"{prefix}.cds.fa.gz"
    pep_path = output_dir / f"{prefix}.pep.fa.gz"
    genome_path = output_dir / f"{prefix}.genome.fa.gz"
    counts: Counter[str] = Counter()
    used_proteins: dict[str, tuple[str, str]] = {}
    used_gene_ids: dict[str, tuple[str, str, int]] = {}
    transcript_proteins: dict[str, str] = {}
    assigned_mrnas: dict[int, str] = {}

    with gzip.open(gff_path, "wt", encoding="utf-8", newline="\n") as gff_handle, gzip.open(
        cds_path, "wt", encoding="utf-8", newline="\n"
    ) as cds_handle, gzip.open(pep_path, "wt", encoding="utf-8", newline="\n") as pep_handle:
        genome_handle = (
            gzip.open(genome_path, "wt", encoding="utf-8", newline="\n") if include_genome else None
        )
        try:
            with open_text(source) as input_handle:
                for record in SeqIO.parse(input_handle, "genbank"):
                    if max_records is not None and counts["records"] >= max_records:
                        break
                    counts["records"] += 1
                    if not record.id or not record.seq:
                        raise ConversionError(f"GBFF record {counts['records']} has no ID or sequence")
                    if genome_handle is not None:
                        fasta_record(genome_handle, record.id, str(record.seq).upper())
                        counts["genome_records"] += 1
                    source_feature = next((item for item in record.features if item.type == "source"), None)
                    chromosome = _first(source_feature, "chromosome") if source_feature else ""
                    if chromosome.casefold() in {"", "unknown", "na", "n/a"}:
                        chromosome = ""
                    organelle = (_first(source_feature, "organelle") if source_feature else "").casefold()
                    genome_type = "genomic"
                    if "mitochond" in organelle:
                        genome_type = "mitochondrion"
                    elif "plastid" in organelle or "chloroplast" in organelle:
                        genome_type = "plastid"
                    elif chromosome:
                        genome_type = "chromosome"
                    taxon = next(
                        (
                            value.split(":", 1)[1]
                            for value in (_values(source_feature, "db_xref") if source_feature else [])
                            if value.startswith("taxon:")
                        ),
                        "",
                    )
                    gff_handle.write(f"##sequence-region {record.id} 1 {len(record.seq)}\n")
                    gff_handle.write(
                        _gff_row(
                            record.id,
                            "region",
                            1,
                            len(record.seq),
                            "+",
                            [
                                ("ID", record.id),
                                ("genome", genome_type),
                                ("chromosome", chromosome),
                                ("Dbxref", f"taxon:{taxon}" if taxon else ""),
                            ],
                        )
                    )
                    counts["regions"] += 1

                    genes = _build_genes(record, counts)
                    mrnas: list[tuple[str, Any]] = []
                    for feature in record.features:
                        if feature.type not in TRANSCRIPT_FEATURES:
                            continue
                        identity = _identity(feature)
                        key = f"{identity}|{_strand(feature)}" if identity else ""
                        mrnas.append((key, feature))

                    coding: list[_Coding] = []
                    for feature in record.features:
                        if feature.type != "CDS":
                            continue
                        counts["cds_features"] += 1
                        protein_id = _first(feature, "protein_id")
                        if not protein_id:
                            counts["cds_missing_protein_id_skipped"] += 1
                            continue
                        if protein_id in used_proteins:
                            previous = used_proteins[protein_id]
                            raise ConversionError(
                                f"duplicate protein_id {protein_id!r}: {previous} and {record.id}"
                            )
                        identity = _identity(feature, fallback=protein_id)
                        gene_key = f"{identity}|{_strand(feature)}"
                        gene = genes.get(gene_key)
                        if gene is None:
                            start, end = _span(feature)
                            gene = _Gene(
                                key=gene_key,
                                identity=identity,
                                seqid=record.id,
                                start=start,
                                end=end,
                                strand="+" if _strand(feature) > 0 else "-",
                                name=_first(feature, "gene") or identity,
                                locus_tag=_first(feature, "locus_tag"),
                                dbxref=",".join(_values(feature, "db_xref")),
                            )
                            genes[gene_key] = gene
                            counts["genes_synthesized_from_cds"] += 1
                        mrna = _match_mrna(feature, gene, mrnas, set(assigned_mrnas))
                        transcript_id = _first(feature, "transcript_id") or (
                            _first(mrna, "transcript_id") if mrna is not None else ""
                        ) or protein_id
                        previous_protein = transcript_proteins.get(transcript_id)
                        if previous_protein is not None and previous_protein != protein_id:
                            raise ConversionError(
                                f"transcript {transcript_id!r} maps to multiple proteins: "
                                f"{previous_protein}, {protein_id}"
                            )
                        try:
                            offset = max(0, int(_first(feature, "codon_start") or "1") - 1)
                        except ValueError:
                            raise ConversionError(f"{protein_id}: invalid codon_start")
                        nucleotide = str(feature.extract(record.seq)).upper()
                        if offset:
                            nucleotide = nucleotide[offset:]
                            counts["codon_start_adjusted"] += 1
                        cds_parts = _adjusted_cds_parts(feature, offset)
                        if not nucleotide or not cds_parts:
                            counts["empty_cds_skipped"] += 1
                            continue
                        raw_table = _first(feature, "transl_table") or "1"
                        try:
                            translation_table = int(raw_table)
                        except ValueError as exc:
                            raise ConversionError(
                                f"{protein_id}: invalid transl_table {raw_table!r}"
                            ) from exc
                        translation = _first(feature, "translation")
                        if translation:
                            protein = "".join(translation.split()).upper()
                            check = validate_cds_pep(
                                nucleotide,
                                protein,
                                translation_table=translation_table,
                            )
                            if not check.valid:
                                if not strict_translation:
                                    counts["translation_inconsistent_skipped"] += 1
                                    continue
                                raise ConversionError(
                                    f"{protein_id}: extracted CDS and translation qualifier "
                                    f"are inconsistent ({check.reason})"
                                )
                            counts["translation_from_qualifier"] += 1
                            counts["translation_qualifier_verified"] += 1
                        else:
                            protein = str(
                                feature.extract(record.seq)[offset:].translate(
                                    table=translation_table
                                )
                            ).upper()
                            counts["translation_derived"] += 1
                        protein = protein.rstrip("*")
                        if not protein:
                            counts["empty_translation_skipped"] += 1
                            continue
                        exon_parts = _parts(mrna) if mrna is not None else list(cds_parts)
                        synthesized = mrna is None
                        if synthesized:
                            counts["transcripts_synthesized_from_cds"] += 1
                        item = _Coding(
                            protein_id=protein_id,
                            transcript_id=transcript_id,
                            gene=gene,
                            mrna=mrna,
                            cds_parts=cds_parts,
                            exon_parts=exon_parts,
                            nucleotide=nucleotide,
                            protein=protein,
                            translation_table=translation_table,
                            synthesized_transcript=synthesized,
                        )
                        coding.append(item)
                        if mrna is not None:
                            assigned_mrnas[id(mrna)] = protein_id
                        used_proteins[protein_id] = (record.id, transcript_id)
                        transcript_proteins[transcript_id] = protein_id
                        fasta_record(cds_handle, protein_id, nucleotide)
                        fasta_record(pep_handle, protein_id, protein)
                        counts["proteins"] += 1
                        counts["cds_segments"] += len(cds_parts)
                        counts["exon_segments"] += len(exon_parts)

                    coding.sort(key=lambda item: (min(part[0] for part in item.cds_parts), item.protein_id))
                    emitted_genes: dict[str, str] = {}
                    for item in coding:
                        gene = item.gene
                        gene_id = _unique_gene_id(gene, record.id, used_gene_ids, counts)
                        if gene.key not in emitted_genes:
                            tx_start = min(part[0] for part in item.exon_parts)
                            tx_end = max(part[1] for part in item.exon_parts)
                            gene.start = min(gene.start, tx_start)
                            gene.end = max(gene.end, tx_end)
                            gff_handle.write(
                                _gff_row(
                                    record.id,
                                    "gene",
                                    gene.start + 1,
                                    gene.end,
                                    gene.strand,
                                    [
                                        ("ID", f"gene-{gene_id}"),
                                        ("Name", gene.name),
                                        ("gene_biotype", "protein_coding"),
                                        ("locus_tag", gene.locus_tag),
                                        ("Dbxref", gene.dbxref),
                                    ],
                                )
                            )
                            emitted_genes[gene.key] = gene_id
                            counts["genes"] += 1
                        gene_id = emitted_genes[gene.key]
                        tx_start = min(part[0] for part in item.exon_parts)
                        tx_end = max(part[1] for part in item.exon_parts)
                        tx_strand = "+" if item.exon_parts[0][2] > 0 else "-"
                        gff_handle.write(
                            _gff_row(
                                record.id,
                                "mRNA",
                                tx_start + 1,
                                tx_end,
                                tx_strand,
                                [
                                    ("ID", f"rna-{item.transcript_id}"),
                                    ("Parent", f"gene-{gene_id}"),
                                    ("Name", item.transcript_id),
                                    ("transcript_biotype", "mRNA"),
                                ],
                            )
                        )
                        counts["transcripts"] += 1
                        for index, (start, end, strand) in enumerate(
                            sorted(item.exon_parts, key=lambda part: part[0]), 1
                        ):
                            gff_handle.write(
                                _gff_row(
                                    record.id,
                                    "exon",
                                    start + 1,
                                    end,
                                    "+" if strand > 0 else "-",
                                    [
                                        ("ID", f"exon-{item.transcript_id}-{index}"),
                                        ("Parent", f"rna-{item.transcript_id}"),
                                        ("exon_number", index),
                                    ],
                                )
                            )
                        for start, end, strand in sorted(item.cds_parts, key=lambda part: part[0]):
                            gff_handle.write(
                                _gff_row(
                                    record.id,
                                    "CDS",
                                    start + 1,
                                    end,
                                    "+" if strand > 0 else "-",
                                    [
                                        ("ID", f"cds-{item.protein_id}"),
                                        ("Parent", f"rna-{item.transcript_id}"),
                                        ("protein_id", item.protein_id),
                                        ("transl_table", item.translation_table),
                                    ],
                                )
                            )
        finally:
            if genome_handle is not None:
                genome_handle.close()

    if not counts["proteins"] and max_records is None:
        raise ConversionError(f"{source.name}: no usable coding sequences with protein_id and translation")
    if not counts["proteins"]:
        counts["partial_without_coding_sequences"] = 1
    if max_records is not None:
        counts["partial_input"] = 1

    detection = detect_gff_format(gff_path)
    parsed = parse_gff(gff_path)
    fatal = [item for item in parsed.anomalies if item.level == "fatal"]
    if fatal:
        raise ConversionError(f"converted GBFF GFF3 failed validation: {fatal[0].code}: {fatal[0].detail}")
    cds_data = read_cds_fasta(cds_path)
    pep_data = read_pep_fasta(pep_path)
    gff_proteins = {
        protein_id
        for transcript in parsed.transcripts.values()
        for protein_id in transcript.protein_ids
    }
    if not (
        gff_proteins == set(cds_data.sequences) == set(pep_data.sequences)
    ):
        raise ConversionError(
            "GBFF output ID mismatch across GFF/CDS/PEP: "
            f"gff={len(gff_proteins)} cds={len(cds_data.sequences)} pep={len(pep_data.sequences)}"
        )
    if detection.format != "standard" or len(parsed.transcripts) != counts["transcripts"]:
        raise ConversionError("converted GBFF hierarchy could not be read back consistently")
    for key in ("records", "regions", "genes", "transcripts", "proteins", "cds_features", "genome_records"):
        counts[key] += 0

    warnings: list[str] = []
    for key, label in (
        ("partial_without_coding_sequences", "selected GBFF smoke-test records contain no coding sequences; this is not a full-file result"),
        ("cds_missing_protein_id_skipped", "CDS features skipped because protein_id was absent"),
        ("empty_cds_skipped", "empty CDS features skipped"),
        ("empty_translation_skipped", "empty translations skipped"),
        ("transcripts_synthesized_from_cds", "transcripts explicitly synthesized from CDS because no unique mRNA matched"),
        ("genes_synthesized_from_cds", "gene rows explicitly synthesized from CDS"),
        ("gene_ids_disambiguated", "gene IDs disambiguated using record IDs"),
        ("gene_features_merged", "same-identity split gene features merged"),
        ("translation_derived", "proteins translated from CDS because the translation qualifier was absent"),
        ("translation_inconsistent_skipped", "CDS skipped because extracted sequence and /translation disagreed"),
    ):
        if counts[key]:
            warnings.append(f"{label}: {counts[key]}")
    if max_records is not None:
        warnings.append(f"partial smoke-test conversion limited to {max_records} GBFF records; not release-ready")
    warnings.append("GBFF output is WGDI-oriented coding annotation; non-CDS features are not a lossless GenBank export")
    return dict(sorted(counts.items())), warnings
