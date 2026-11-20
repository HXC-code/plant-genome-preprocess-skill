from __future__ import annotations

import gzip
from pathlib import Path

from pgci_preprocess.gff_standardize import standardize_transcript_only_gff


def test_expands_transcript_only_gff(tmp_path: Path) -> None:
    source = tmp_path / "input.gff.gz"
    output = tmp_path / "normalized.gff.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        handle.write("chr1\tEVM\tCDS\t10\t20\t.\t+\t0\tParent=tx1;\n")
        handle.write("chr1\tEVM\tmRNA\t10\t50\t.\t+\t.\tID=tx1;\n")

    result = standardize_transcript_only_gff(source, output)
    assert result.converted
    with gzip.open(output, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert "ID=gene-tx1" in text
    assert "ID=rna-tx1;Parent=gene-tx1" in text
    assert "Parent=rna-tx1;protein_id=tx1" in text


def test_transcript_only_normalization_never_overwrites_protein_id(tmp_path: Path) -> None:
    source = tmp_path / "input.gff.gz"
    output = tmp_path / "normalized.gff.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        handle.write("chr1\tEVM\tCDS\t10\t18\t.\t+\t0\tParent=tx1;protein_id=P1\n")
        handle.write("chr1\tEVM\tmRNA\t10\t50\t.\t+\t.\tID=tx1;Name=alpha\n")

    standardize_transcript_only_gff(source, output)
    with gzip.open(output, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert "Parent=rna-tx1;protein_id=P1" in text
    assert "ID=rna-tx1" in text
    assert "Parent=gene-tx1" in text
    assert "Name=alpha" in text
