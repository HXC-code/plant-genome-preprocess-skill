from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from pgci_preprocess.pipeline import (
    PipelineOptions,
    _read_genome_regions,
    discover_species,
    process_species,
    run_project,
)
from pgci_preprocess.reports import read_tsv
from pgci_preprocess.validate import ValidationOptions, validate_project


def test_genome_regions_recognise_chromosome_seqids_without_description(
    tmp_path: Path,
) -> None:
    genome = tmp_path / "genome.fa"
    genome.write_text(
        ">Chr1\nACGT\n"
        ">Ajul.Chr01\nACGTACGT\n"
        ">NW_026555332.1 unplaced genomic scaffold\nACGT\n",
        encoding="utf-8",
    )
    regions = _read_genome_regions(genome)
    assert regions["Chr1"].is_primary is True
    assert regions["Chr1"].chromosome == "1"
    assert regions["Ajul.Chr01"].is_primary is True
    assert regions["Ajul.Chr01"].chromosome == "1"
    assert regions["NW_026555332.1"].is_primary is False


def _write_gz(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def _make_project(tmp_path: Path, *, nuclear_taxon_2: bool = False) -> tuple[Path, Path]:
    project = tmp_path / "project"
    species_dir = project / "data" / "downloads" / "1_Testales_Testaceae(Testus plantus)"
    resources = project / "resources"
    species_dir.mkdir(parents=True)
    resources.mkdir(parents=True)
    prefix = "GCF_000000001.1_Demo"
    extra_region = ""
    if nuclear_taxon_2:
        extra_region = (
            "chr3\tsrc\tregion\t1\t300\t.\t+\t.\t"
            "ID=chr3:1..300;Dbxref=taxon:999;chromosome=3;genome=chromosome\n"
            "chr3\tsrc\tgene\t10\t100\t.\t+\t.\t"
            "ID=gene-X;Name=X;gene_biotype=protein_coding\n"
            "chr3\tsrc\tmRNA\t10\t100\t.\t+\t.\tID=rna-X1;Parent=gene-X\n"
            "chr3\tsrc\tCDS\t10\t90\t.\t+\t0\tID=cds-X;Parent=rna-X1;protein_id=X1\n"
        )
    gff = (
        "##sequence-region chr1 1 1000\n"
        "chr1\tsrc\tregion\t1\t1000\t.\t+\t.\t"
        "ID=chr1:1..1000;Dbxref=taxon:123;chromosome=1;genome=chromosome\n"
        "chr1\tsrc\tgene\t10\t200\t.\t+\t.\tID=gene-A;Name=A;gene_biotype=protein_coding\n"
        "chr1\tsrc\tmRNA\t10\t200\t.\t+\t.\tID=rna-A1;Parent=gene-A\n"
        "chr1\tsrc\tCDS\t10\t18\t.\t+\t0\tID=cds-P1;Parent=rna-A1;protein_id=P1\n"
        "chr1\tsrc\tgene\t300\t500\t.\t-\t.\tID=gene-B;Name=B;gene_biotype=protein_coding\n"
        "chr1\tsrc\tmRNA\t300\t500\t.\t-\t.\tID=rna-B1;Parent=gene-B\n"
        "chr1\tsrc\tCDS\t300\t309\t.\t-\t0\tID=cds-P2;Parent=rna-B1;protein_id=P2\n"
        "chr1\tsrc\tmRNA\t300\t520\t.\t-\t.\tID=rna-B2;Parent=gene-B\n"
        "chr1\tsrc\tCDS\t300\t315\t.\t-\t0\tID=cds-P3;Parent=rna-B2;protein_id=P3\n"
        "chr1\tsrc\tgene\t600\t700\t.\t+\t.\tID=gene-NC;Name=NC;gene_biotype=lncRNA\n"
        "chr1\tsrc\tlnc_RNA\t600\t700\t.\t+\t.\tID=rna-NC;Parent=gene-NC\n"
        "chr2\tsrc\tregion\t1\t500\t.\t+\t.\t"
        "ID=chr2:1..500;Dbxref=taxon:123;chromosome=2;genome=chromosome\n"
        "chr2\tsrc\tgene\t5\t90\t.\t+\t.\tID=gene-C;Name=C;gene_biotype=protein_coding\n"
        "chr2\tsrc\tmRNA\t5\t90\t.\t+\t.\tID=rna-C1;Parent=gene-C\n"
        "chr2\tsrc\tCDS\t5\t13\t.\t+\t0\tID=cds-P4;Parent=rna-C1;protein_id=P4\n"
        "organelle\tsrc\tregion\t1\t200\t.\t+\t.\t"
        "ID=organelle:1..200;Dbxref=taxon:123;genome=mitochondrion\n"
        "organelle\tsrc\tgene\t1\t30\t.\t+\t.\t"
        "ID=gene-O;Name=O;gene_biotype=protein_coding\n"
        "organelle\tsrc\tmRNA\t1\t30\t.\t+\t.\tID=rna-O1;Parent=gene-O\n"
        "organelle\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds-P5;Parent=rna-O1;protein_id=P5\n"
        + extra_region
    )
    cds = (
        ">lcl|chr1_cds_P1_1 [gene=A] [protein_id=P1]\nATGGCTTAA\n"
        ">lcl|chr1_cds_P2_1 [gene=B] [protein_id=P2]\nATGATGTAA\n"
        ">lcl|chr1_cds_P3_1 [gene=B] [protein_id=P3]\nATGATGATGTAA\n"
        ">lcl|chr2_cds_P4_1 [gene=C] [protein_id=P4]\nATGAAATAA\n"
        ">lcl|organelle_cds_P5_1 [gene=O] [protein_id=P5]\nATGAAATAA\n"
        + (">lcl|chr3_cds_X1_1 [gene=X] [protein_id=X1]\nATGAAATAA\n" if nuclear_taxon_2 else "")
    )
    pep = (
        ">P1 alpha\nMA\n"
        ">P2 beta\nMM\n"
        ">P3 beta long\nMMMM\n"
        ">P4 gamma\nMK\n"
        ">P5 organelle\nMK\n"
        + (">X1 foreign\nMK\n" if nuclear_taxon_2 else "")
    )
    _write_gz(species_dir / f"{prefix}_genomic.gff.gz", gff)
    _write_gz(species_dir / f"{prefix}_cds_from_genomic.fna.gz", cds)
    _write_gz(species_dir / f"{prefix}_protein.faa.gz", pep)
    inventory = resources / "download_inventory.tsv"
    with inventory.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["showid", "Order", "Family", "Species", "gff", "genome", "cds", "pep", "rna", "gbff", "gtf"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "showid": "1",
                "Order": "Testales",
                "Family": "Testaceae",
                "Species": "Testus plantus",
                "gff": f"{prefix}_genomic.gff.gz",
                "genome": "",
                "cds": f"{prefix}_cds_from_genomic.fna.gz",
                "pep": f"{prefix}_protein.faa.gz",
                "rna": "",
                "gbff": "",
                "gtf": "",
            }
        )
    return project, inventory


def test_end_to_end_builds_consistent_wgdi_bundle(tmp_path: Path) -> None:
    project, _ = _make_project(tmp_path)
    species = discover_species(project)[0]
    options = PipelineOptions(
        project_root=project,
        min_genes_per_sequence=0,
        target_gene_coverage=0.9,
        max_sequences=200,
    )

    result = process_species(species, options)

    assert result["status"] == "OK"
    analysis = project / "wgdi_input" / "s1"
    assert (analysis / "s1.new.gff").is_file()
    assert (analysis / "s1.lens").is_file()
    assert (analysis / "s1.cds").is_file()
    assert (analysis / "s1.pep").is_file()

    gff_lines = (analysis / "s1.new.gff").read_text(encoding="utf-8").splitlines()
    assert len(gff_lines) == 3
    assert gff_lines[0].split("\t") == ["1", "s11g0001", "10", "200", "+", "1", "A"]
    assert any(line.split("\t")[1] == "s11g0002" for line in gff_lines)
    assert all("P2" not in line for line in gff_lines)
    assert not any("O" == line.split("\t")[-1] for line in gff_lines)

    lens = (analysis / "s1.lens").read_text(encoding="utf-8").splitlines()
    assert lens == ["1\t1000\t2", "2\t500\t1"]

    with (analysis / "s1.pep").open(encoding="utf-8") as handle:
        pep_ids = [line[1:].split()[0] for line in handle if line.startswith(">")]
    with (analysis / "s1.cds").open(encoding="utf-8") as handle:
        cds_ids = [line[1:].split()[0] for line in handle if line.startswith(">")]
    gff_ids = [line.split("\t")[1] for line in gff_lines]
    assert pep_ids == cds_ids == gff_ids


def test_validate_rejects_not_ready_bundle_unless_explicitly_allowed(
    tmp_path: Path,
) -> None:
    project, _ = _make_project(tmp_path)
    run_project(
        project,
        options=PipelineOptions(
            project_root=project,
            min_genes_per_sequence=0,
            target_gene_coverage=0.9,
            max_sequences=200,
        ),
    )
    qc_path = project / "wgdi_input" / "s1" / "s1.qc.json"
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    qc["sequence_selection"]["wgdi_ready"] = "no"
    qc["sequence_selection"]["ratio"] = 0.25
    qc_path.write_text(json.dumps(qc), encoding="utf-8")

    strict = validate_project(
        ValidationOptions(
            project_root=project,
            require_compat_links=False,
        )
    )
    assert not strict.ok
    assert "wgdi_not_ready_low_coverage" in strict.rows[0]["problems"]

    archival = validate_project(
        ValidationOptions(
            project_root=project,
            require_compat_links=False,
            allow_not_ready=True,
        )
    )
    assert archival.ok


def test_all_snapshot_is_opt_in_and_dedupes_gene_map(tmp_path: Path) -> None:
    """The pre-filtering .all.* snapshot is opt-in; one output layer only."""

    project, _ = _make_project(tmp_path)
    results = run_project(
        project,
        options=PipelineOptions(project_root=project, min_genes_per_sequence=0),
    )
    code = str(results[0]["code"])
    analysis = project / "wgdi_input" / code

    # four-piece set is present, the pre-filtering snapshot is not written
    for suffix in (".pep", ".cds", ".new.gff", ".lens"):
        assert (analysis / f"{code}{suffix}").is_file(), suffix
    for suffix in (".all.pep.fa", ".all.cds.fa", ".all.gff", ".all.lens"):
        assert not (analysis / f"{code}{suffix}").exists(), suffix

    # single-layer layout: traceability files sit next to the four-piece set,
    # each exactly once; no clean_v2 evidence layer is created
    for suffix in (".gene_map.tsv", ".chr_map.tsv", ".removed_list.tsv", ".chromosome_map.tsv"):
        assert (analysis / f"{code}{suffix}").is_file(), suffix
    assert not (project / "clean_v2").exists()

    # the index must not point at files that were deliberately not written
    row = next(
        item
        for item in read_tsv(project / "resources" / "wgdi_v2_index.tsv")
        if item["code"] == code
    )
    assert row["all_pep"] == "" and row["all_lens"] == ""
    assert row["pep"].endswith(f"{code}.pep")

    assert results[0]["n_all_genes"] >= results[0]["n_wgdi_genes"]


def test_keep_all_writes_the_snapshot_and_indexes_it(tmp_path: Path) -> None:
    project, _ = _make_project(tmp_path)
    results = run_project(
        project,
        options=PipelineOptions(
            project_root=project,
            min_genes_per_sequence=0,
            write_all_sequences=True,
        ),
    )
    code = str(results[0]["code"])
    analysis = project / "wgdi_input" / code
    for suffix in (".all.pep.fa", ".all.cds.fa", ".all.gff", ".all.lens"):
        assert (analysis / f"{code}{suffix}").is_file(), suffix
    row = next(
        item
        for item in read_tsv(project / "resources" / "wgdi_v2_index.tsv")
        if item["code"] == code
    )
    assert row["all_pep"].endswith(f"{code}.all.pep.fa")


def test_run_shows_chromosome_inventory_to_the_caller(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    project, _ = _make_project(tmp_path)
    run_project(
        project,
        options=PipelineOptions(project_root=project, min_genes_per_sequence=0),
    )
    out = capsys.readouterr().out
    assert "chromosomes " in out
    assert ".chromosome_map.tsv" in out
    # per-chromosome gene counts must be visible, not only the sequence count
    assert "1:" in out
    # the inventory line must stay ASCII so a cp1252 console cannot crash run
    inventory = [
        line for line in out.splitlines() if "(id:genes)" in line
    ]
    assert inventory
    assert all(line.isascii() for line in inventory)


def test_run_survives_a_non_utf8_windows_console(tmp_path: Path) -> None:
    project, _ = _make_project(tmp_path)
    src = Path(__file__).resolve().parents[1] / "src"
    env = {**os.environ, "PYTHONIOENCODING": "cp1252", "PYTHONPATH": str(src)}
    env.pop("PYTHONUTF8", None)
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pgci_preprocess",
            "run",
            "--project-root",
            str(project),
            "--min-genes-per-sequence",
            "0",
            "--min-genes-per-chromosome",
            "0",
        ],
        env=env,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_pipeline_quarantines_multiple_nuclear_taxa(tmp_path: Path) -> None:
    project, _ = _make_project(tmp_path, nuclear_taxon_2=True)
    species = discover_species(project)[0]
    options = PipelineOptions(
        project_root=project,
        min_genes_per_sequence=0,
        target_gene_coverage=0.9,
        max_sequences=200,
    )

    result = process_species(species, options)

    assert result["status"] == "FAIL"
    assert result["reason"] == "multiple_nuclear_taxa"
    assert not (project / "wgdi_input" / "s1" / "s1.new.gff").exists()


def test_run_project_writes_summary_index_and_checksums(tmp_path: Path) -> None:
    project, _ = _make_project(tmp_path)
    options = PipelineOptions(
        project_root=project,
        min_genes_per_sequence=0,
        target_gene_coverage=0.9,
        max_sequences=200,
    )

    results = run_project(project, show_ids={"1"}, options=options)

    assert results[0]["status"] == "OK"
    assert (project / "resources" / "stage_preprocess_v2_summary.tsv").is_file()
    assert (project / "resources" / "wgdi_v2_index.tsv").is_file()
    # no SHA-256 manifest is written any more
    assert not (project / "resources" / "wgdi_v2_manifest_sha256.tsv").exists()


def test_partial_run_preserves_unselected_manifest_and_summary_rows(tmp_path: Path) -> None:
    project, _ = _make_project(tmp_path)
    resources = project / "resources"
    index_fields = [
        "code",
        "species",
        "accession",
        "wgdi_ready",
        "pep",
        "cds",
        "gff",
        "lens",
        "all_pep",
        "all_cds",
        "all_gff",
        "all_lens",
        "gene_map",
        "chr_map",
        "rejected",
        "qc",
    ]
    with (resources / "wgdi_v2_index.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=index_fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({"code": "s999", "species": "Old species"})
    summary_fields = [
        "show_id",
        "code",
        "species",
        "accession",
        "status",
        "assembly_build",
        "n_regions",
        "n_representative_genes",
        "n_all_genes",
        "n_wgdi_genes",
        "n_wgdi_sequences",
        "sequence_mode",
        "retain_ratio",
        "wgdi_ready",
        "nuclear_taxa",
        "n_anomalies",
        "rejection_counts",
        "reason",
        "detail",
    ]
    with (resources / "stage_preprocess_v2_summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({"show_id": "999", "code": "s999", "status": "OK"})

    options = PipelineOptions(
        project_root=project,
        min_genes_per_sequence=0,
        target_gene_coverage=0.9,
        max_sequences=200,
    )
    run_project(project, show_ids={"1"}, options=options)

    index_rows = list(
        csv.DictReader(
            (resources / "wgdi_v2_index.tsv").open(encoding="utf-8"),
            delimiter="\t",
        )
    )
    summary_rows = list(
        csv.DictReader(
            (resources / "stage_preprocess_v2_summary.tsv").open(encoding="utf-8"),
            delimiter="\t",
        )
    )
    assert {row["code"] for row in index_rows} == {"s1", "s999"}
    assert {row["code"] for row in summary_rows} == {"s1", "s999"}


def test_pipeline_options_default_to_no_overwrite(tmp_path: Path) -> None:
    """Library callers must get the same protection the CLI provides by default.

    2.3.0 left ``overwrite=True`` in the dataclass while the CLI defaulted to no
    overwrite, so any script calling ``run_project``/``process_species`` directly
    silently replaced a validated bundle.
    """

    assert PipelineOptions(project_root=tmp_path).overwrite is False

    project, _ = _make_project(tmp_path)
    species = discover_species(project)[0]
    options = PipelineOptions(
        project_root=project,
        min_genes_per_sequence=0,
        target_gene_coverage=0.9,
        max_sequences=200,
    )
    assert options.overwrite is False

    first = process_species(species, options)
    assert first["status"] == "OK"
    gff = project / "wgdi_input" / "s1" / "s1.new.gff"
    before = gff.read_bytes()

    second = process_species(species, options)
    assert second["status"] == "SKIPPED"
    assert second["detail"] == "existing_valid_output_no_overwrite"
    assert gff.read_bytes() == before

    options.overwrite = True
    third = process_species(species, options)
    assert third["status"] == "OK"
    assert gff.read_bytes() == before


def test_run_project_default_options_preserve_existing_output(tmp_path: Path) -> None:
    project, _ = _make_project(tmp_path)
    assert run_project(project, show_ids={"1"})[0]["status"] == "OK"
    assert run_project(project, show_ids={"1"})[0]["status"] == "SKIPPED"


def test_organelle_only_annotation_has_a_distinct_failure_reason(tmp_path: Path) -> None:
    """Chloroplast-only submissions are a data gap, not a converter bug.

    Several raw download folders hold only GeSeq chloroplast stubs or tiny
    mitochondrial records.  Every transcript is rejected as organelle, and the
    batch report must say so explicitly instead of the generic
    ``no_valid_representative_genes``.
    """

    project = tmp_path / "organelle_project"
    species_dir = project / "data" / "downloads" / "1_Testales_Testaceae(Testus plantus)"
    resources = project / "resources"
    species_dir.mkdir(parents=True)
    resources.mkdir(parents=True)
    gff = (
        "cp\tsrc\tregion\t1\t2000\t.\t+\t.\t"
        "ID=cp:1..2000;Dbxref=taxon:123;genome=chloroplast\n"
        "cp\tsrc\tgene\t10\t200\t.\t+\t.\tID=gene-A;Name=psbA;gene_biotype=protein_coding\n"
        "cp\tsrc\tmRNA\t10\t200\t.\t+\t.\tID=rna-A1;Parent=gene-A\n"
        "cp\tsrc\tCDS\t10\t18\t.\t+\t0\tID=cds-P1;Parent=rna-A1;protein_id=P1\n"
    )
    _write_gz(species_dir / "GWHTEST00000000.1.gff.gz", gff)
    _write_gz(
        species_dir / "GWHTEST00000000.1.CDS.fasta.gz",
        ">P1\nATGGCTTAA\n",
    )
    _write_gz(species_dir / "GWHTEST00000000.1.Protein.faa.gz", ">P1\nMA\n")
    (resources / "download_inventory.tsv").write_text(
        "showid\tOrder\tFamily\tSpecies\tgff\tgenome\tcds\tpep\n"
        "1\tTestales\tTestaceae\tTestus plantus\tGWHTEST00000000.1.gff.gz\t\t"
        "GWHTEST00000000.1.CDS.fasta.gz\tGWHTEST00000000.1.Protein.faa.gz\n",
        encoding="utf-8",
    )

    species = discover_species(project)[0]
    result = process_species(species, PipelineOptions(project_root=project))
    assert result["status"] == "FAIL"
    assert result["reason"] == "organelle_only_annotation"
    assert result["detail"] == "transcripts=1"
    assert not (project / "wgdi_input" / "s1").exists()

