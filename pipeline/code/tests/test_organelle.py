from __future__ import annotations

from collections import Counter

import pytest

from pgci_preprocess.models import GffData, Region, Transcript
from pgci_preprocess.organelle import infer_synthetic_organelles


def _data(
    *,
    seqid: str,
    gene_count: int,
    gff_span: int,
    cds_length: int = 900,
    synthetic: bool = True,
) -> tuple[GffData, dict[str, int]]:
    region = Region(
        seqid=seqid,
        length=gene_count * 2_000,
        genome_type="genomic",
        chromosome="",
        taxon_id="",
        source_order=0,
        is_primary=False,
        synthetic=synthetic,
    )
    transcripts: dict[str, Transcript] = {}
    cds_lengths: dict[str, int] = {}
    for index in range(gene_count):
        transcript_id = f"rna-generic-{index:04d}"
        start = index * 2_000 + 1
        transcripts[transcript_id] = Transcript(
            source_id=transcript_id,
            gene_id=transcript_id,
            seqid=seqid,
            start=start,
            end=start + gff_span - 1,
            strand="+",
            feature_type="mRNA",
            protein_ids={transcript_id},
        )
        cds_lengths[transcript_id] = cds_length
    return (
        GffData(
            regions={seqid: region},
            genes={},
            transcripts=transcripts,
            feature_counts=Counter({"gene": gene_count}),
            anomalies=[],
            annotation_level="simplified5",
        ),
        cds_lengths,
    )


def test_compact_intronless_synthetic_region_is_inferred_as_organelle() -> None:
    data, cds_lengths = _data(seqid="NC_999999.1", gene_count=60, gff_span=900)
    report = infer_synthetic_organelles(data, cds_lengths)
    assert report["inferred"][0]["seqid"] == "NC_999999.1"
    assert report["inferred"][0]["reason"] == "dense_compact_synthetic_region"
    assert data.regions["NC_999999.1"].is_organelle
    assert any(
        item.code == "organelle_sequence_inferred" for item in data.anomalies
    )


def test_dense_multiexon_region_is_reported_for_review_not_excluded() -> None:
    data, cds_lengths = _data(seqid="NC_888888.1", gene_count=60, gff_span=4_000)
    report = infer_synthetic_organelles(data, cds_lengths)
    assert report["inferred"] == []
    assert report["review_candidates"][0]["seqid"] == "NC_888888.1"
    assert not data.regions["NC_888888.1"].is_organelle


def test_authoritative_region_is_never_reclassified_by_heuristic() -> None:
    data, cds_lengths = _data(
        seqid="NC_777777.1",
        gene_count=60,
        gff_span=900,
        synthetic=False,
    )
    report = infer_synthetic_organelles(data, cds_lengths)
    assert report["inferred"] == []
    assert not data.regions["NC_777777.1"].is_organelle


def test_custom_pattern_marks_synthetic_sequence() -> None:
    data, cds_lengths = _data(seqid="PTG_SCAFFOLD_99", gene_count=3, gff_span=5_000)
    report = infer_synthetic_organelles(
        data,
        cds_lengths,
        custom_pattern=r"^PTG_",
    )
    assert report["inferred"][0]["reason"] == "custom_organelle_pattern"
    assert data.regions["PTG_SCAFFOLD_99"].genome_type == "plastid"


def test_invalid_custom_pattern_is_rejected() -> None:
    data, cds_lengths = _data(seqid="scaffold1", gene_count=2, gff_span=900)
    with pytest.raises(ValueError, match="invalid organelle pattern"):
        infer_synthetic_organelles(
            data,
            cds_lengths,
            custom_pattern="[",
        )


def _named_data(
    *,
    seqid: str,
    gene_names: list[str],
    length: int,
    gff_span: int = 4_000,
    cds_length: int = 900,
) -> tuple[GffData, dict[str, int]]:
    region = Region(
        seqid=seqid,
        length=length,
        genome_type="genomic",
        chromosome="",
        taxon_id="",
        source_order=0,
        is_primary=False,
        synthetic=True,
    )
    transcripts: dict[str, Transcript] = {}
    cds_lengths: dict[str, int] = {}
    for index, name in enumerate(gene_names):
        start = index * 5_000 + 1
        transcripts[name] = Transcript(
            source_id=name,
            gene_id=name,
            seqid=seqid,
            start=start,
            end=start + gff_span - 1,
            strand="+",
            feature_type="mRNA",
            protein_ids={name},
        )
        cds_lengths[name] = cds_length
    return (
        GffData(
            regions={seqid: region},
            genes={},
            transcripts=transcripts,
            feature_counts=Counter({"gene": len(gene_names)}),
            anomalies=[],
            annotation_level="simplified5",
        ),
        cds_lengths,
    )


def test_mitochondrial_core_gene_cluster_is_inferred() -> None:
    names = ["nad5", "nad1", "cob", "cox1", "atp6"] + [
        f"hypo-{i:03d}" for i in range(35)
    ]
    data, cds_lengths = _named_data(
        seqid="scaffold_12", gene_names=names, length=400_000
    )
    report = infer_synthetic_organelles(data, cds_lengths)
    assert report["inferred"][0]["seqid"] == "scaffold_12"
    assert report["inferred"][0]["genome_type"] == "mitochondrion"
    assert report["inferred"][0]["reason"] == "mitochondrial_gene_content"
    assert {"nad5", "cob", "cox1", "atp6"} <= set(
        report["inferred"][0]["mt_core_genes"]
    )
    assert data.regions["scaffold_12"].is_organelle


def test_single_numt_copy_is_reviewed_not_excluded() -> None:
    names = ["cox1"] + [f"nuclear-{i:03d}" for i in range(60)]
    data, cds_lengths = _named_data(
        seqid="scaffold_3", gene_names=names, length=3_000_000
    )
    report = infer_synthetic_organelles(data, cds_lengths)
    assert report["inferred"] == []
    assert not data.regions["scaffold_3"].is_organelle
    assert report["review_candidates"][0]["genome_type"] == "mitochondrion"
    assert (
        report["review_candidates"][0]["reason"]
        == "mitochondrial_gene_content_requires_review"
    )


def test_opaque_large_sparse_scaffold_is_not_guessed() -> None:
    data, cds_lengths = _data(
        seqid="scaffold_99", gene_count=40, gff_span=4_000
    )
    data.regions["scaffold_99"].length = 1_500_000
    report = infer_synthetic_organelles(data, cds_lengths)
    assert report["inferred"] == []
    assert report["review_candidates"] == []
    assert not data.regions["scaffold_99"].is_organelle
    assert "limits" in report


def test_mt_core_gene_threshold_is_configurable() -> None:
    names = ["nad5", "cob", "cox1"] + [f"hypo-{i:03d}" for i in range(30)]
    data, cds_lengths = _named_data(
        seqid="scaffold_5", gene_names=names, length=300_000
    )
    default = infer_synthetic_organelles(data, cds_lengths)
    assert default["inferred"][0]["genome_type"] == "mitochondrion"

    data2, cds2 = _named_data(
        seqid="scaffold_5", gene_names=names, length=300_000
    )
    strict = infer_synthetic_organelles(data2, cds2, mt_min_core_genes=5)
    assert strict["inferred"] == []
    assert strict["review_candidates"][0]["genome_type"] == "mitochondrion"


def test_chrm_identifier_is_mitochondrion() -> None:
    data, cds_lengths = _data(seqid="chrM", gene_count=3, gff_span=900)
    report = infer_synthetic_organelles(data, cds_lengths)
    assert report["inferred"][0]["genome_type"] == "mitochondrion"
    assert report["inferred"][0]["reason"] == "sequence_or_gene_identifier"
