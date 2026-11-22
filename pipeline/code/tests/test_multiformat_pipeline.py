from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from pgci_preprocess.gff_formats import detect_gff_format
from pgci_preprocess.gff_standardize import normalize_gff_species
from pgci_preprocess.pipeline import PipelineOptions, discover_species, run_project
from pgci_preprocess.reports import read_tsv, write_tsv
from tests.support import build_multiformat_project

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if not SRC_DIR.is_dir():
    SRC_DIR = PROJECT_ROOT / "scripts"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
LEGACY_WRAPPER = SCRIPTS_DIR / "preprocess_wgdi_v2.py"
if not LEGACY_WRAPPER.is_file():
    LEGACY_WRAPPER = SCRIPTS_DIR / "preprocess.py"


def _run_cli(project: Path, *extra: str, module: bool = True) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SRC_DIR)
    if module:
        command = [sys.executable, "-m", "pgci_preprocess", "run", "--project-root", str(project)]
    else:
        command = [sys.executable, str(LEGACY_WRAPPER), "--project-root", str(project)]
    command.extend(extra)
    return subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        env=env,
        text=True,
        capture_output=True,
        timeout=180,
    )


def _qc(project: Path, code: str) -> dict[str, object]:
    return json.loads((project / "analysis_v2" / code / f"{code}.qc.json").read_text(encoding="utf-8"))


def test_run_multiformat_project(tmp_path: Path) -> None:
    project = build_multiformat_project(tmp_path)["project"]
    results = run_project(
        project,
        options=PipelineOptions(project_root=project, min_genes_per_sequence=1),
    )
    by_code = {str(result["code"]): result for result in results}
    assert set(by_code) == {"s1", "s2", "s3", "s4"}
    assert all(result["status"] == "OK" for result in results)

    assert by_code["s1"]["n_representative_genes"] == 3
    assert by_code["s1"]["n_wgdi_genes"] == 3
    assert by_code["s1"]["n_wgdi_sequences"] == 1
    assert by_code["s1"]["wgdi_ready"] == "yes"
    assert by_code["s1"]["rejection_counts"]["organelle_sequence"] == 1

    assert by_code["s2"]["n_representative_genes"] == 2
    assert by_code["s2"]["n_wgdi_sequences"] == 1
    qc2 = _qc(project, "s2")
    assert qc2["format_detection"]["format"] == "transcript_only"
    assert qc2["format_exceptions"]["transcript_only_gff"] is True
    assert qc2["n_regions"] == 1
    # Authoritative genome length (2000) supersedes the synthetic max-end (308).
    assert qc2["region_backfill"]["genome"] == ["1"]
    assert (project / "analysis_v2" / "s2" / "s2.lens").read_text(encoding="utf-8").splitlines() == ["1\t2000\t2"]

    assert by_code["s3"]["n_representative_genes"] == 3
    assert by_code["s3"]["n_wgdi_sequences"] == 2
    qc3 = _qc(project, "s3")
    assert qc3["format_detection"]["format"] == "simplified5"
    assert qc3["format_exceptions"]["simplified5_gff"] is True

    assert by_code["s4"]["n_representative_genes"] == 2
    qc4 = _qc(project, "s4")
    assert qc4["region_backfill"]["transcript"] == ["1", "2"]
    assert qc4["format_exceptions"]["region_backfilled_from_transcript"] == 2
    assert qc4["format_detection"]["format"] == "standard"

    index = read_tsv(project / "resources" / "wgdi_v2_index.tsv")
    assert [row["code"] for row in index] == ["s1", "s2", "s3", "s4"]
    for row in index:
        for key in ("pep", "cds", "gff", "lens", "qc"):
            assert (project / row[key]).is_file(), (row["code"], key)


def test_inventory_pipe_candidates_are_resolved_by_role(tmp_path: Path) -> None:
    paths = build_multiformat_project(tmp_path)
    inventory = paths["project"] / "resources" / "download_inventory.tsv"
    rows = read_tsv(inventory)
    row = next(item for item in rows if item["showid"] == "1")
    original_gff = row["gff"]
    original_cds = row["cds"]
    original_pep = row["pep"]
    row["gff"] = f"missing.gff|{original_gff}"
    row["cds"] = f"missing.cds|{original_cds}"
    row["pep"] = f"missing.pep|{original_pep}"
    fields = list(rows[0])
    write_tsv(inventory, fields, rows)

    species = {item.show_id: item for item in discover_species(paths["project"])}
    assert species["1"].gff_name == original_gff
    assert species["1"].cds_name == original_cds
    assert species["1"].pep_name == original_pep


def test_normalize_gff_expands_transcript_only(tmp_path: Path) -> None:
    project = build_multiformat_project(tmp_path)["project"]
    results = normalize_gff_species(project, ["2"])
    assert len(results) == 1
    assert results[0].status == "CONVERTED"
    assert results[0].output_path and results[0].output_path.is_file()
    assert detect_gff_format(results[0].output_path).format == "standard"
    inventory = {row["showid"]: row for row in read_tsv(project / "resources" / "download_inventory.tsv")}
    assert inventory["2"]["gff"].startswith("normalized/")

    run = run_project(
        project,
        show_ids={"2"},
        options=PipelineOptions(project_root=project, min_genes_per_sequence=1),
    )
    assert run[0]["status"] == "OK"
    assert run[0]["n_representative_genes"] == 2


def test_cli_and_legacy_wrapper_are_equivalent(tmp_path: Path) -> None:
    fixture = build_multiformat_project(tmp_path)["project"]
    direct = tmp_path / "direct"
    cli = tmp_path / "cli"
    wrapper = tmp_path / "wrapper"
    shutil.copytree(fixture, direct)
    shutil.copytree(fixture, cli)
    shutil.copytree(fixture, wrapper)

    direct_results = run_project(
        direct,
        options=PipelineOptions(project_root=direct, min_genes_per_sequence=1),
    )
    assert all(result["status"] == "OK" for result in direct_results)

    cli_result = _run_cli(cli, "--min-genes-per-sequence", "1", "--min-genes-per-chromosome", "0")
    assert cli_result.returncode == 0, cli_result.stderr
    wrapper_result = _run_cli(wrapper, "--min-genes-per-sequence", "1", "--min-genes-per-chromosome", "0", module=False)
    assert wrapper_result.returncode == 0, wrapper_result.stderr

    for code in ("s1", "s2", "s3", "s4"):
        for suffix in (".new.gff", ".lens", ".cds", ".pep"):
            left = (direct / "analysis_v2" / code / f"{code}{suffix}").read_bytes()
            middle = (cli / "analysis_v2" / code / f"{code}{suffix}").read_bytes()
            right = (wrapper / "analysis_v2" / code / f"{code}{suffix}").read_bytes()
            assert left == middle == right, (code, suffix)
        left_qc = _qc(direct, code)
        middle_qc = _qc(cli, code)
        right_qc = _qc(wrapper, code)
        for field in (
            "status",
            "code",
            "accession",
            "n_regions",
            "n_genes",
            "n_transcripts",
            "n_representative_genes",
            "rejection_counts",
            "sequence_selection",
            "wgdi_bundle",
            "format_detection",
            "format_exceptions",
            "region_backfill",
            "protein_normalization",
        ):
            assert left_qc[field] == middle_qc[field] == right_qc[field], (code, field)

    assert [row["code"] for row in read_tsv(direct / "resources" / "wgdi_v2_index.tsv")] == [
        row["code"] for row in read_tsv(cli / "resources" / "wgdi_v2_index.tsv")
    ] == [row["code"] for row in read_tsv(wrapper / "resources" / "wgdi_v2_index.tsv")]


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _fa(records: dict[str, str]) -> str:
    return "".join(f">{k}\n{v}\n" for k, v in records.items())


def test_synthetic_regions_yield_to_authoritative_genome_lengths(tmp_path: Path) -> None:
    """Regression: transcript-only/simplified max-end regions must not shadow
    authoritative genome FASTA lengths (Golden-109 lens-only DIFF adjudication).

    s1 transcript-only + genome -> lens from genome (2000), not max end (308)
    s2 transcript-only, no genome -> lens from max transcript end (308)
    s3 standard GFF with a real region row (1500) + genome (2000) -> row kept
    s4 transcript-only + genome with a mitochondrion record -> organelle call
    """
    import csv

    project = tmp_path / "plant-genome-resources"
    downloads = project / "data" / "downloads"
    resources = project / "resources"
    resources.mkdir(parents=True)
    tx_gff = "\n".join(
        [
            "##gff-version 3",
            "1\tEVM\tmRNA\t100\t117\t.\t+\t.\tID=txA;Name=txA",
            "1\tEVM\tCDS\t100\t117\t.\t+\t0\tID=cds.txA;Parent=txA",
            "1\tEVM\tmRNA\t300\t308\t.\t+\t.\tID=txB;Name=txB",
            "1\tEVM\tCDS\t300\t308\t.\t+\t0\tID=cds.txB;Parent=txB",
            "",
        ]
    )
    cds = _fa({"txA": "ATGGCTGCTGCTGCTTAA", "txB": "ATGAAATAA"})
    pep = _fa({"txA": "MAAAA", "txB": "MK"})
    rows = []
    for sid, species, genome in (
        ("1", "Testus primus", "GCF_000000001.1_t1_genome.fa"),
        ("2", "Testus secundus", ""),
        ("3", "Testus tertius", "GCF_000000003.1_t3_genome.fa"),
        ("4", "Testus quartus", "GCF_000000004.1_t4_genome.fa"),
    ):
        d = downloads / f"{sid}_Testales_Testaceae({species})"
        d.mkdir(parents=True)
        gff_name = f"GCF_00000000{sid}.1_t{sid}_genomic.gff"
        if sid == "3":
            gff = "\n".join(
                [
                    "##gff-version 3",
                    "1\tsrc\tregion\t1\t1500\t.\t+\t.\tID=1;genome=genomic",
                    "1\tsrc\tgene\t100\t117\t.\t+\t.\tID=gA",
                    "1\tsrc\tmRNA\t100\t117\t.\t+\t.\tID=rna-txA;Parent=gA;transcript_id=txA",
                    "1\tsrc\tCDS\t100\t117\t.\t+\t0\tID=cds-txA;Parent=rna-txA;protein_id=txA",
                    "",
                ]
            )
            _write(d / f"GCF_00000000{sid}.1_t{sid}_cds.fa", _fa({"txA": "ATGGCTGCTGCTGCTTAA"}))
            _write(d / f"GCF_00000000{sid}.1_t{sid}_protein.faa", _fa({"txA": "MAAAA"}))
        elif sid == "4":
            gff = tx_gff + "MT\tEVM\tmRNA\t10\t18\t.\t+\t.\tID=txM;Name=txM\nMT\tEVM\tCDS\t10\t18\t.\t+\t0\tID=cds.txM;Parent=txM\n"
            _write(d / f"GCF_00000000{sid}.1_t{sid}_cds.fa", _fa({"txA": "ATGGCTGCTGCTGCTTAA", "txB": "ATGAAATAA", "txM": "ATGAAATAA"}))
            _write(d / f"GCF_00000000{sid}.1_t{sid}_protein.faa", _fa({"txA": "MAAAA", "txB": "MK", "txM": "MK"}))
        else:
            gff = tx_gff
            _write(d / f"GCF_00000000{sid}.1_t{sid}_cds.fa", cds)
            _write(d / f"GCF_00000000{sid}.1_t{sid}_protein.faa", pep)
        _write(d / gff_name, gff)
        if genome:
            text = ">1\n" + "N" * 2000 + "\n"
            if sid == "4":
                text += ">MT Testus quartus mitochondrion\n" + "N" * 300 + "\n"
            _write(d / genome, text)
        rows.append(
            {
                "showid": sid, "Order": "Testales", "Family": "Testaceae", "Species": species,
                "gff": gff_name, "genome": genome,
                "cds": f"GCF_00000000{sid}.1_t{sid}_cds.fa",
                "pep": f"GCF_00000000{sid}.1_t{sid}_protein.faa",
                "rna": "", "gbff": "", "gtf": "",
            }
        )
    fields = ["showid", "Order", "Family", "Species", "gff", "genome", "cds", "pep", "rna", "gbff", "gtf"]
    with (resources / "download_inventory.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    results = run_project(
        project, options=PipelineOptions(project_root=project, min_genes_per_sequence=1)
    )
    assert all(r["status"] == "OK" for r in results), results

    def lens(code: str) -> list[str]:
        return (project / "analysis_v2" / code / f"{code}.lens").read_text(encoding="utf-8").splitlines()

    assert lens("s1") == ["1\t2000\t2"]
    assert _qc(project, "s1")["region_backfill"]["genome"] == ["1"]
    assert lens("s2") == ["1\t308\t2"]
    assert _qc(project, "s2")["region_backfill"]["genome"] == []
    assert lens("s3") == ["1\t1500\t1"]
    assert _qc(project, "s3")["region_backfill"]["genome"] == []
    assert lens("s4") == ["1\t2000\t2"]
    qc4 = _qc(project, "s4")
    assert qc4["rejection_counts"]["organelle_sequence"] == 1
    assert qc4["n_organelle_genes"] == 1
