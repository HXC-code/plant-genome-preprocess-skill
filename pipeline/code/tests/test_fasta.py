from __future__ import annotations

import bz2
import gzip
import lzma
import zipfile
from pathlib import Path

from pgci_preprocess.fasta import (
    extract_cds_id,
    open_text,
    read_cds_fasta,
    read_pep_fasta,
)


def test_read_ncbi_cds_uses_protein_id(tmp_path: Path) -> None:
    path = tmp_path / "sample.cds.fna.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(
            ">lcl|chr1_cds_XP_1.1_1 [gene=G1] [protein_id=XP_1.1] [location=1..9]\n"
            "ATGGCTTAA\n"
            ">lcl|chr1_cds_XP_2.1_1 [gene=G2] [protein_id=XP_2.1]\n"
            "ATGAAATAA\n"
        )

    result = read_cds_fasta(path)
    assert result.sequences == {"XP_1.1": "ATGGCTTAA", "XP_2.1": "ATGAAATAA"}
    assert result.counts["records"] == 2


def test_read_pep_detects_duplicate_conflict_and_empty_record(tmp_path: Path) -> None:
    path = tmp_path / "sample.pep.faa"
    path.write_text(
        ">P1 description one\nMA\n>P1 description two\nMW\n>P2 empty\n>P3\nMAA\n",
        encoding="utf-8",
    )

    result = read_pep_fasta(path)
    assert result.sequences["P1"] == "MA"
    assert result.sequences["P3"] == "MAA"
    assert result.counts["duplicate_conflict"] == 1
    assert result.counts["empty_sequence"] == 1


def test_cds_fallback_extracts_protein_id_from_cds_marker() -> None:
    """The cds marker fallback must work without a protein_id qualifier.

    The released 2.3.0 pattern used ``\\s`` and ``\\d`` inside a raw string, so
    it matched literal backslashes and never fired as documented; headers fell
    through to the raw first token and every pairing then failed as
    ``sequence_missing``.
    """

    assert extract_cds_id("lcl|NC_001.1_cds_XP_001792555.1_1") == "XP_001792555.1"
    assert extract_cds_id("lcl|scaffold_1_cds_ABC123_2") == "ABC123"
    assert extract_cds_id("cds_XP_001792555.1") == "XP_001792555.1"
    assert extract_cds_id("cds_ABC123_1 extra description") == "ABC123"
    assert extract_cds_id("cds_ABC123") == "ABC123"
    assert extract_cds_id("cds.ENSOP_00123") == "ENSOP_00123"
    assert extract_cds_id("cds-GA277554.1") == "GA277554.1"
    # A long numeric tail is part of the identifier, not a CDS index.
    assert extract_cds_id("cds_ENSOP_00123") == "ENSOP_00123"


def test_cds_fallback_keeps_explicit_qualifiers_and_plain_ids() -> None:
    assert (
        extract_cds_id("lcl|chr1_cds_XP_1.1_1 [gene=G1] [protein_id=XP_1.1]")
        == "XP_1.1"
    )
    assert extract_cds_id("sp|P12345.1|G1 Protein=G1") == "G1"
    # No cds marker: the identifier is the first token, unchanged behaviour.
    assert extract_cds_id("P1 description one") == "P1"
    assert extract_cds_id("XP_001792555.1") == "XP_001792555.1"


def test_read_cds_fasta_pairs_with_plain_pep_ids(tmp_path: Path) -> None:
    cds_path = tmp_path / "x.cds.fna"
    cds_path.write_text(
        ">lcl|scaffold_1_cds_XP_1.1_1\nATGGCTTAA\n"
        ">lcl|scaffold_1_cds_XP_2.1_1\nATGAAATAA\n",
        encoding="utf-8",
    )
    result = read_cds_fasta(cds_path)
    assert set(result.sequences) == {"XP_1.1", "XP_2.1"}


def test_open_text_detects_compression_by_magic(tmp_path: Path) -> None:
    text = ">P1\nMAA\n"
    mislabeled_gzip = tmp_path / "mislabeled.fa"
    with gzip.open(mislabeled_gzip, "wt", encoding="utf-8") as handle:
        handle.write(text)
    bzip = tmp_path / "sample.fa.bz2"
    with bz2.open(bzip, "wt", encoding="utf-8") as handle:
        handle.write(text)
    xz = tmp_path / "sample.fa.xz"
    with lzma.open(xz, "wt", encoding="utf-8") as handle:
        handle.write(text)
    zipped = tmp_path / "sample.fa.zip"
    with zipfile.ZipFile(zipped, "w") as archive:
        archive.writestr("sample.fa", text)
    for path in (mislabeled_gzip, bzip, xz, zipped):
        with open_text(path) as handle:
            assert handle.read().splitlines() == text.splitlines()

