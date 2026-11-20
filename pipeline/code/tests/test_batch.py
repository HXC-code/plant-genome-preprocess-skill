from __future__ import annotations

import csv
from pathlib import Path

from pgci_preprocess.batch import run_batch
from pgci_preprocess.legacy_convert import run_legacy_conversion
from pgci_preprocess.pipeline import PipelineOptions
from pgci_preprocess.reports import read_tsv


def _write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_bundle(directory: Path, show_id: str) -> dict[str, str]:
    directory.mkdir(parents=True)
    prefix = f"GCF_00000000{show_id}.1_test{show_id}"
    gff = (
        "##gff-version 3\n"
        "1\tRefSeq\tregion\t1\t1000\t.\t+\t.\tID=1;genome=genomic;chromosome=1\n"
        "1\tRefSeq\tgene\t100\t120\t.\t+\t.\tID=gene1;Name=gene1;gene_biotype=protein_coding\n"
        "1\tRefSeq\tmRNA\t100\t120\t.\t+\t.\tID=rna-tx1;Parent=gene1;transcript_id=tx1\n"
        "1\tRefSeq\tCDS\t100\t108\t.\t+\t0\tID=cds1;Parent=rna-tx1;protein_id=prot1\n"
    )
    (directory / f"{prefix}_genomic.gff").write_text(gff, encoding="utf-8")
    (directory / f"{prefix}_genome.fa").write_text(">1\n" + "N" * 1000 + "\n", encoding="utf-8")
    (directory / f"{prefix}_cds.fa").write_text(">prot1\nATGGCTTAA\n", encoding="utf-8")
    (directory / f"{prefix}_protein.faa").write_text(">prot1\nMA\n", encoding="utf-8")
    return {
        "genome": f"{prefix}_genome.fa",
        "gff": f"{prefix}_genomic.gff",
        "cds": f"{prefix}_cds.fa",
        "pep": f"{prefix}_protein.faa",
    }


def test_batch_prepare_run_validate(tmp_path: Path) -> None:
    project = tmp_path / "project"
    source = tmp_path / "downloads" / "7_Testales_Testaceae(Testus septimus)"
    b = _write_bundle(source, "7")
    (project / "data" / "downloads").mkdir(parents=True)
    (project / "resources").mkdir(parents=True)
    (project / "resources" / "download_inventory.tsv").write_text(
        "showid\tOrder\tFamily\tSpecies\tgff\tgenome\tcds\tpep\trna\tgbff\tgtf\n",
        encoding="utf-8",
    )
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
            "current_bundle",
            "current_genome",
            "current_gff_or_gtf",
            "current_cds",
            "current_pep",
            "candidate_accession",
        ],
        [
            {
                "show_id": "7",
                "Order": "Testales",
                "Family": "Testaceae",
                "Species": "Testus septimus",
                "status": "complete",
                "folder": str(source),
                "current_bundle": "primary",
                "current_genome": b["genome"],
                "current_gff_or_gtf": b["gff"],
                "current_cds": b["cds"],
                "current_pep": b["pep"],
                "candidate_accession": "GCF_000000007.1",
            }
        ],
    )

    result = run_batch(
        project,
        status,
        expected_count=1,
        batch_tag="test",
        pipeline_options=None,
    )
    assert result.candidates == ("7",)
    assert result.prepared_files >= 4
    assert result.run_results and result.run_results[0]["status"] == "OK"
    assert result.validation is not None and result.validation.ok
    assert result.ok
    assert {row["code"] for row in read_tsv(project / "resources" / "wgdi_v2_index.tsv")} == {"s7"}


def test_batch_count_mismatch_is_a_failure(tmp_path: Path) -> None:
    project = tmp_path / "project"
    (project / "data" / "downloads").mkdir(parents=True)
    (project / "resources").mkdir(parents=True)
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
        ],
        [],
    )
    result = run_batch(project, status, expected_count=1)
    assert result.prepare_failed is True
    assert result.ok is False
    assert "expected 1 candidates" in result.skipped_prepare_reason


def test_batch_auto_converts_gtf_extracts_fastas_and_runs(tmp_path: Path) -> None:
    project = tmp_path / "project"
    source = tmp_path / "downloads" / "8_Testales_Testaceae(Testus octavus)"
    source.mkdir(parents=True)
    (project / "data" / "downloads").mkdir(parents=True)
    (project / "resources").mkdir(parents=True)
    (project / "resources" / "download_inventory.tsv").write_text(
        "showid\tOrder\tFamily\tSpecies\tgff\tgenome\tcds\tpep\trna\tgbff\tgtf\n",
        encoding="utf-8",
    )
    genome = "A" * 9 + "ATGGCTTAA" + "A" * 82
    (source / "GCA_000000008.1_genomic.fna").write_text(
        ">chr1 Chromosome 1\n" + genome + "\n",
        encoding="utf-8",
    )
    (source / "GCA_000000008.1_genomic.gtf").write_text(
        "chr1\tsource\tgene\t10\t18\t.\t+\t.\tgene_id \"G1\";\n"
        "chr1\tsource\ttranscript\t10\t18\t.\t+\t.\t"
        "gene_id \"G1\"; transcript_id \"TX1\";\n"
        "chr1\tsource\tCDS\t10\t18\t.\t+\t0\t"
        "gene_id \"G1\"; transcript_id \"TX1\"; protein_id \"P1\";\n",
        encoding="utf-8",
    )
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
            "current_bundle",
            "current_genome",
            "current_gff_or_gtf",
            "current_cds",
            "current_pep",
            "candidate_accession",
        ],
        [
            {
                "show_id": "8",
                "Order": "Testales",
                "Family": "Testaceae",
                "Species": "Testus octavus",
                "status": "complete",
                "folder": str(source),
                "current_bundle": "primary",
                "current_genome": "GCA_000000008.1_genomic.fna",
                "current_gff_or_gtf": "GCA_000000008.1_genomic.gtf",
                "current_cds": "",
                "current_pep": "",
                "candidate_accession": "GCA_000000008.1",
            }
        ],
    )
    result = run_batch(
        project,
        status,
        expected_count=1,
        batch_tag="gtf-auto",
        pipeline_options=PipelineOptions(
            project_root=project,
            min_genes_per_sequence=1,
        ),
        require_compat_links=False,
    )
    assert result.converted == 1
    assert result.extracted == 1
    assert result.run_results and result.run_results[0]["status"] == "OK", result.run_results
    assert result.validation is not None and result.validation.ok
    assert result.ok


def test_legacy_conversion_folder_scan_and_finalize(tmp_path: Path) -> None:
    project = tmp_path / "project"
    resources = project / "resources"
    resources.mkdir(parents=True)
    downloads = tmp_path / "downloads"
    source = downloads / "99_GAP_Testaceae(Testus nonaginta)"
    _write_bundle(source, "99")
    status = tmp_path / "legacy.tsv"
    _write_tsv(
        status,
        ["show_id", "Species", "Order", "Family", "folder", "candidate_accession"],
        [
            {
                "show_id": "99",
                "Species": "Testus nonaginta",
                "Order": "Testales",
                "Family": "Testaceae",
                "folder": str(source),
                "candidate_accession": "GCF_000000099.1",
            }
        ],
    )

    results = run_legacy_conversion(
        project,
        status,
        downloads_root=downloads,
        overwrite=False,
    )
    assert len(results) == 1
    assert results[0]["status"] == "OK"
    assert {row["code"] for row in read_tsv(resources / "wgdi_v2_index.tsv")} == {"s99"}
    assert {row["code"] for row in read_tsv(resources / "stage_preprocess_v2_summary.tsv")} == {"s99"}
    assert (project / "analysis_v2" / "s99" / "s99.qc.json").is_file()
