from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest

from pgci_preprocess.batch_prepare import (
    BatchPreparationError,
    compound_suffix,
    prepare_batch,
)


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _write_gzip(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def _write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_compound_suffix() -> None:
    assert compound_suffix("demo.cds.fa.gz") == ".cds.fa.gz"
    assert compound_suffix("demo.cds") == ".cds"
    assert compound_suffix("demo.genome.fa.gz") == ".genome.fa.gz"


def test_prepare_batch_mirrors_and_updates_inventory(tmp_path: Path) -> None:
    project = tmp_path / "project"
    source = tmp_path / "source" / "1_Testales_Testaceae(Testus plantus)"
    status = tmp_path / "current_gap_jobs.tsv"
    (project / "data" / "downloads").mkdir(parents=True)
    (project / "resources").mkdir(parents=True)
    source.mkdir(parents=True)

    _write_gzip(source / "demo.genome.fa.gz", ">chr1\n" + "A" * 100 + "\n")
    _write_gzip(
        source / "demo.gff3.gz",
        "##gff-version 3\n"
        "chr1\tsrc\tregion\t1\t100\t.\t+\t.\tID=chr1;chromosome=1\n"
        "chr1\tsrc\tgene\t10\t27\t.\t+\t.\tID=gene1\n"
        "chr1\tsrc\tmRNA\t10\t27\t.\t+\t.\tID=rna-tx1;Parent=gene1\n"
        "chr1\tsrc\tCDS\t10\t27\t.\t+\t0\tID=cds-P1;Parent=rna-tx1;protein_id=P1\n",
    )
    _write_gzip(source / "demo.cds.fa.gz", ">P1\n" + "ATG" * 6 + "\n")
    _write_gzip(source / "demo.pep.fa.gz", ">P1\nMMMMMM\n")

    fields = [
        "show_id",
        "Order",
        "Family",
        "Species",
        "status",
        "folder",
        "current_bundle",
        "current_genome",
        "current_gff_or_gtf",
        "current_cds",
        "current_pep",
            "candidate_accession",
            "candidate_source",
            "candidate_version",
            "download_url",
    ]
    _write_tsv(
        status,
        fields,
        [
            {
                "show_id": "1",
                "Order": "Testales",
                "Family": "Testaceae",
                "Species": "Testus plantus",
                "status": "complete",
                "folder": str(source),
                "current_bundle": "alt_assembly",
                "current_genome": "demo.genome.fa.gz",
                "current_gff_or_gtf": "demo.gff3.gz",
                "current_cds": "demo.cds.fa.gz",
                "current_pep": "demo.pep.fa.gz",
                "candidate_accession": "TESTASSEMBLY1",
                "candidate_source": "Test source",
                "candidate_version": "1",
                "download_url": "https://example.invalid/test",
            }
        ],
    )
    _write_tsv(
        project / "resources" / "download_inventory.tsv",
        ["showid", "Order", "Family", "Species", "gff", "genome", "cds", "pep", "rna", "gbff", "gtf"],
        [],
    )

    result = prepare_batch(project, status, expected_count=1, batch_tag="test")
    assert result.candidate_ids == ("1",)
    assert len(result.prepared_files) == 4
    target = project / "data" / "downloads" / source.name
    assert (target / "demo.gff3.gz").is_file()
    assert result.manifest_path and result.manifest_path.is_file()
    assert result.suffix_summary_path and result.suffix_summary_path.is_file()

    with (project / "resources" / "download_inventory.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["showid"] == "1"
    assert row["genome"] == "demo.genome.fa.gz"
    assert row["gff"] == "demo.gff3.gz"
    assert row["cds"] == "demo.cds.fa.gz"
    assert row["pep"] == "demo.pep.fa.gz"
    with (project / "resources" / "assembly_accessions.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        assembly = next(csv.DictReader(handle, delimiter="\t"))
    assert assembly["show_id"] == "1"
    assert assembly["assembly_id"] == "TESTASSEMBLY1"
    assert assembly["evidence_path"].endswith("demo.gff3.gz")


def test_prepare_batch_rejects_count_mismatch(tmp_path: Path) -> None:
    project = tmp_path / "project"
    (project / "data" / "downloads").mkdir(parents=True)
    (project / "resources").mkdir(parents=True)
    status = tmp_path / "empty.tsv"
    _write_tsv(
        status,
        [
            "show_id",
            "Order",
            "Family",
            "Species",
            "status",
            "folder",
            "current_genome",
            "current_gff_or_gtf",
            "current_cds",
            "current_pep",
        ],
        [],
    )
    with pytest.raises(BatchPreparationError, match="expected 14 candidates, found 0"):
        prepare_batch(project, status, expected_count=14)


def test_prepare_routes_gtf_converts_registers_and_extracts_missing_fastas(
    tmp_path: Path,
) -> None:
    project = tmp_path / "project"
    source = tmp_path / "source" / "2_Testales_Testaceae(Testus secundus)"
    status = tmp_path / "status.tsv"
    (project / "data" / "downloads").mkdir(parents=True)
    (project / "resources").mkdir(parents=True)
    source.mkdir(parents=True)

    genome = "A" * 9 + "ATGGCTTAA" + "A" * 82
    _write_gzip(source / "demo.genome.fa.gz", ">chr1\n" + genome + "\n")
    _write_gzip(source / "demo.pro.fa.gz", ">spurious\nMMMM\n")
    _write_gzip(
        source / "demo.gtf.gz",
        "chr1\tsource\tgene\t10\t18\t.\t+\t.\tgene_id \"G1\";\n"
        "chr1\tsource\ttranscript\t10\t18\t.\t+\t.\t"
        "gene_id \"G1\"; transcript_id \"TX1\";\n"
        "chr1\tsource\tCDS\t10\t18\t.\t+\t0\t"
        "gene_id \"G1\"; transcript_id \"TX1\"; protein_id \"P1\";\n",
    )
    _write_tsv(
        status,
        [
            "show_id",
            "Order",
            "Family",
            "Species",
            "status",
            "folder",
            "current_bundle",
            "current_genome",
            "current_gff_or_gtf",
            "current_cds",
            "current_pep",
            "candidate_accession",
        ],
        [
            {
                "show_id": "2",
                "Order": "Testales",
                "Family": "Testaceae",
                "Species": "Testus secundus",
                "status": "complete",
                "folder": str(source),
                "current_bundle": "primary",
                "current_genome": "demo.pro.fa.gz|demo.genome.fa.gz",
                "current_gff_or_gtf": "missing.gff|demo.gtf.gz",
                "current_cds": "",
                "current_pep": "",
                "candidate_accession": "GCA_000000002.1",
            }
        ],
    )
    _write_tsv(
        project / "resources" / "download_inventory.tsv",
        ["showid", "Order", "Family", "Species", "gff", "genome", "cds", "pep", "rna", "gbff", "gtf"],
        [],
    )

    result = prepare_batch(project, status, expected_count=1, batch_tag="route")
    assert result.converted == 1
    assert result.extracted == 1
    assert result.routing_path and result.routing_path.is_file()
    with (project / "resources" / "download_inventory.tsv").open(
        encoding="utf-8",
        newline="",
    ) as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["gff"].startswith("derived/") and row["gff"].endswith(".gff3.gz")
    assert row["cds"].startswith("derived/") and row["cds"].endswith(".cds.fa")
    assert row["pep"].startswith("derived/") and row["pep"].endswith(".pep.fa")
    assert row["genome"] == "demo.genome.fa.gz"
    assert row["gtf"] == "demo.gtf.gz"


def test_prepare_rejects_ambiguous_pipe_candidates(tmp_path: Path) -> None:
    project = tmp_path / "project"
    source = tmp_path / "source" / "3_Testales_Testaceae(Testus tertius)"
    source.mkdir(parents=True)
    (project / "data" / "downloads").mkdir(parents=True)
    (project / "resources").mkdir(parents=True)
    for name in ("one.genome.fa", "two.genome.fa"):
        _write(source / name, ">chr1\n" + "A" * 100 + "\n")
    _write(
        source / "demo.gff",
        "##gff-version 3\n"
        "chr1\tsrc\tregion\t1\t100\t.\t+\t.\tID=chr1;chromosome=1\n",
    )
    _write(source / "demo.cds.fa", ">P1\nATGGCTTAA\n")
    _write(source / "demo.pep.fa", ">P1\nMA\n")
    status = tmp_path / "status.tsv"
    _write_tsv(
        status,
        [
            "show_id",
            "Order",
            "Family",
            "Species",
            "status",
            "folder",
            "current_genome",
            "current_gff_or_gtf",
            "current_cds",
            "current_pep",
            "candidate_accession",
        ],
        [
            {
                "show_id": "3",
                "Order": "Testales",
                "Family": "Testaceae",
                "Species": "Testus tertius",
                "status": "complete",
                "folder": str(source),
                "current_genome": "one.genome.fa|two.genome.fa",
                "current_gff_or_gtf": "demo.gff",
                "current_cds": "demo.cds.fa",
                "current_pep": "demo.pep.fa",
                "candidate_accession": "GCA_000000003.1",
            }
        ],
    )
    _write_tsv(
        project / "resources" / "download_inventory.tsv",
        ["showid", "Order", "Family", "Species", "gff", "genome", "cds", "pep", "rna", "gbff", "gtf"],
        [],
    )
    with pytest.raises(BatchPreparationError, match="ambiguous genome candidates"):
        prepare_batch(project, status, expected_count=1)
