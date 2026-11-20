from __future__ import annotations

from pathlib import Path

import pytest

from pgci_preprocess.fasta import read_cds_fasta, read_pep_fasta
from pgci_preprocess.sequence_extract import (
    SequenceExtractionError,
    extract_cds_pep_from_gff,
)


def test_extracts_plus_and_minus_cds_with_shared_ids(tmp_path: Path) -> None:
    genome = tmp_path / "genome.fa"
    genome.write_text(
        ">chr1\n"
        + "A" * 9
        + "ATGGCTTAA"
        + "A" * 82
        + "\n"
        + ">chr2\n"
        + "A" * 49
        + "TTAAGCCAT"
        + "A" * 42
        + "\n",
        encoding="utf-8",
    )
    gff = tmp_path / "annotation.gff"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tsrc\tregion\t1\t100\t.\t+\t.\tID=chr1;chromosome=1\n"
        "chr1\tsrc\tgene\t10\t18\t.\t+\t.\tID=g1\n"
        "chr1\tsrc\tmRNA\t10\t18\t.\t+\t.\tID=rna-t1;Parent=g1\n"
        "chr1\tsrc\tCDS\t10\t18\t.\t+\t0\tID=cds-P1;Parent=rna-t1;protein_id=P1\n"
        "chr2\tsrc\tregion\t1\t100\t.\t+\t.\tID=chr2;chromosome=2\n"
        "chr2\tsrc\tgene\t50\t58\t.\t-\t.\tID=g2\n"
        "chr2\tsrc\tmRNA\t50\t58\t.\t-\t.\tID=rna-t2;Parent=g2\n"
        "chr2\tsrc\tCDS\t50\t58\t.\t-\t0\tID=cds-P2;Parent=rna-t2;protein_id=P2\n",
        encoding="utf-8",
    )
    result = extract_cds_pep_from_gff(
        genome,
        gff,
        tmp_path / "derived.cds.fa",
        tmp_path / "derived.pep.fa",
    )
    cds = read_cds_fasta(result.cds_path)
    pep = read_pep_fasta(result.pep_path)
    assert cds.sequences == {"P1": "ATGGCTTAA", "P2": "ATGGCTTAA"}
    assert pep.sequences == {"P1": "MA", "P2": "MA"}
    assert result.counts["proteins"] == 2


def test_extraction_rejects_cds_coordinates_outside_genome(tmp_path: Path) -> None:
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "A" * 20 + "\n", encoding="utf-8")
    gff = tmp_path / "annotation.gff"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tsrc\tregion\t1\t20\t.\t+\t.\tID=chr1;chromosome=1\n"
        "chr1\tsrc\tgene\t18\t30\t.\t+\t.\tID=g1\n"
        "chr1\tsrc\tmRNA\t18\t30\t.\t+\t.\tID=rna-t1;Parent=g1\n"
        "chr1\tsrc\tCDS\t18\t30\t.\t+\t0\tParent=rna-t1;protein_id=P1\n",
        encoding="utf-8",
    )
    with pytest.raises(SequenceExtractionError, match="CDS exceeds genome"):
        extract_cds_pep_from_gff(
            genome,
            gff,
            tmp_path / "derived.cds.fa",
            tmp_path / "derived.pep.fa",
        )
