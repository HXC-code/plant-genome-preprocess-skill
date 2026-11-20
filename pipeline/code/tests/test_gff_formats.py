from __future__ import annotations

from collections import Counter
from pathlib import Path

from pgci_preprocess.fasta import normalize_evm_protein_sequences
from pgci_preprocess.gff import parse_gff
from pgci_preprocess.gff_formats import (
    apply_gwh_accessions,
    backfill_regions_from_transcripts,
    detect_gff_format,
    mrna_accession_map,
    normalize_cds_protein_ids,
    parse_transcript_only_gff,
    parse_with_adapter,
)
from pgci_preprocess.models import FastaData, GffData, Transcript
from tests.support import build_multiformat_project


def _gff_path(project: Path, show_id: str, name: str) -> Path:
    return project / "data" / "downloads" / f"{show_id}_Testales_Testaceae(Testus {['primus','secundus','tertius','quartus'][int(show_id)-1]})" / name


def test_detect_all_supported_formats(tmp_path: Path) -> None:
    project = build_multiformat_project(tmp_path)["project"]
    cases = {
        "1": ("standard", "GCF_000000001.1_test1_genomic.gff"),
        "2": ("transcript_only", "GCF_000000002.1_test2_genomic.gff"),
        "3": ("simplified5", "GCF_000000003.1_test3.gff3"),
        "4": ("standard", "GCF_000000004.1_test4_genomic.gff"),
    }
    for show_id, (expected, filename) in cases.items():
        detection = detect_gff_format(_gff_path(project, show_id, filename))
        assert detection.format == expected
        assert detection.rows > 0


def test_transcript_only_and_simplified_parsers(tmp_path: Path) -> None:
    project = build_multiformat_project(tmp_path)["project"]
    transcript_only = parse_transcript_only_gff(
        _gff_path(project, "2", "GCF_000000002.1_test2_genomic.gff")
    )
    assert len(transcript_only.genes) == 2
    assert len(transcript_only.transcripts) == 2
    assert transcript_only.regions["1"].length == 308
    assert transcript_only.transcripts["txA"].protein_ids == {"txA"}
    assert transcript_only.annotation_level == "transcript_only"

    simplified = parse_with_adapter(
        _gff_path(project, "3", "GCF_000000003.1_test3.gff3")
    ).data
    assert len(simplified.genes) == 3
    assert len(simplified.transcripts) == 3
    assert simplified.regions["1"].length == 308
    assert simplified.regions["2"].length == 108
    assert simplified.annotation_level == "simplified5"


def test_transcript_only_preserves_explicit_cds_protein_id(tmp_path: Path) -> None:
    gff = tmp_path / "explicit_protein.gff"
    gff.write_text(
        "1\tsource\tmRNA\t10\t30\t.\t+\t.\tID=tx1;\n"
        "1\tsource\tCDS\t10\t18\t.\t+\t0\t"
        "ID=cds-P1;Parent=tx1;protein_id=P1;transl_table=11\n",
        encoding="utf-8",
    )
    data = parse_transcript_only_gff(gff)
    assert data.transcripts["tx1"].protein_ids == {"P1"}
    assert data.transcripts["tx1"].translation_table == 11


def test_gwh_accession_mapping(tmp_path: Path) -> None:
    gff = tmp_path / "gwh.gff"
    gff.write_text(
        "\n".join(
            [
                "1\tsource\tgene\t10\t30\t.\t+\t.\tID=gene1",
                "1\tsource\tmRNA\t10\t30\t.\t+\t.\tID=mrna1;Parent=gene1;Accession=protA",
                "1\tsource\tCDS\t10\t18\t.\t+\t0\tID=cds1;Parent=mrna1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    mapping = mrna_accession_map(gff)
    assert mapping == {"mrna1": "protA"}
    data = parse_gff(gff)
    changed = apply_gwh_accessions(data, mapping)
    assert changed == 1
    assert data.transcripts["mrna1"].protein_ids == {"protA"}


def test_gwh_accession_mapping_never_overwrites_protein_accession(tmp_path: Path) -> None:
    """Regression: Golden-109 s2224/s4290 lost every gene to sequence_missing.

    GWH GFF3 CDS rows carry ``Protein_Accession=GWHP...`` which matches the
    ``Protein=`` key in both CDS and PEP FASTA.  The mRNA ``Accession``
    (``GWHT...``) must only be a fallback when no protein id exists.
    """

    gff = tmp_path / "gwh_protein_accession.gff"
    gff.write_text(
        "\n".join(
            [
                "1\tsource\tgene\t10\t30\t.\t+\t.\tID=GGL01G000001;Accession=GWHGTEST000001.1",
                "1\tsource\tmRNA\t10\t30\t.\t+\t.\tID=GGL01G000001.mRNA1;Parent=GGL01G000001;Accession=GWHTTEST000001.1",
                "1\tsource\tCDS\t10\t30\t.\t+\t0\tID=cds1;Parent=GGL01G000001.mRNA1;Protein_Accession=GWHPTEST000001.1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    data = parse_gff(gff)
    assert data.transcripts["GGL01G000001.mRNA1"].protein_ids == {"GWHPTEST000001.1"}
    changed = apply_gwh_accessions(data, mrna_accession_map(gff))
    assert changed == 0
    assert data.transcripts["GGL01G000001.mRNA1"].protein_ids == {"GWHPTEST000001.1"}


def test_cds_prefix_normalization() -> None:
    data = GffData(
        regions={},
        genes={},
        transcripts={
            "t1": Transcript(
                source_id="t1",
                gene_id="g1",
                seqid="1",
                start=1,
                end=9,
                strand="+",
                feature_type="mRNA",
                protein_ids={"cds.p1", "cds-p2", "plain"},
            )
        },
        feature_counts=Counter(),
        anomalies=[],
    )
    changed = normalize_cds_protein_ids(data)
    assert changed == 1
    assert data.transcripts["t1"].protein_ids == {"p1", "p2", "plain"}


def test_transcript_region_backfill() -> None:
    data = GffData(
        regions={},
        genes={},
        transcripts={
            "t1": Transcript(
                source_id="t1",
                gene_id="g1",
                seqid="chr1",
                start=10,
                end=120,
                strand="+",
                feature_type="mRNA",
                protein_ids={"p1"},
            )
        },
        feature_counts=Counter(),
        anomalies=[],
    )
    backfilled = backfill_regions_from_transcripts(data, reason="test")
    assert backfilled == ["chr1"]
    assert data.regions["chr1"].length == 120
    assert data.regions["chr1"].is_primary is True
    assert any(item.code == "region_backfilled_from_transcript" for item in data.anomalies)


def test_gene_span_and_gene_cds_direct_detection(tmp_path: Path) -> None:
    gene_span = tmp_path / "gene_span.gff"
    gene_span.write_text(
        "1\tsource\tgene\t10\t20\t.\t+\t.\tID=g1\n",
        encoding="utf-8",
    )
    assert detect_gff_format(gene_span).format == "gene_span"

    gene_cds = tmp_path / "gene_cds.gff"
    gene_cds.write_text(
        "\n".join(
            [
                "1\tsource\tgene\t10\t30\t.\t+\t.\tID=g1",
                "1\tsource\tCDS\t10\t18\t.\t+\t0\tID=cds1;Parent=g1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    assert detect_gff_format(gene_cds).format == "gene_cds_direct"


def test_gene_cds_direct_merge_split_flag_is_forwarded(tmp_path: Path) -> None:
    gff = tmp_path / "gene_cds_split.gff"
    gff.write_text(
        "\n".join(
            [
                "1\tsource\tgene\t100\t120\t.\t+\t.\tID=g1",
                "1\tsource\tgene\t200\t220\t.\t+\t.\tID=g1",
                "1\tsource\tCDS\t100\t117\t.\t+\t0\tID=cds1;Parent=g1",
                "1\tsource\tCDS\t200\t217\t.\t+\t0\tID=cds2;Parent=g1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    parsed = parse_with_adapter(gff, merge_split_genes=True)
    assert parsed.detection.format == "gene_cds_direct"
    assert parsed.data.split_genes_merged == 1
    assert parsed.data.genes["g1"].start == 100
    assert parsed.data.genes["g1"].end == 220


def test_evm_protein_dot_normalization() -> None:
    data = FastaData(
        sequences={"p1": "MAA.", "p2": "MA.A.", "p3": "MAA"},
        counts=Counter(),
        anomalies=[],
    )
    report = normalize_evm_protein_sequences(data)
    assert report == {"records": 2, "trailing_stop": 2, "internal_wildcard": 1}
    assert data.sequences == {"p1": "MAA*", "p2": "MAXA*", "p3": "MAA"}
    assert data.counts["dot_normalized_records"] == 2
