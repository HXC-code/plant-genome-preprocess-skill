from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
WORKFLOW_ROOT = PROJECT_ROOT.parent
BUILD_REAL_FIXTURE = PROJECT_ROOT / "tools" / "build_real_wgdi_fixture.py"
COLLECT_REAL = PROJECT_ROOT / "tools" / "collect_real_data_validation.py"


def _bundle(root: Path, code: str) -> Path:
    bundle = root / code
    bundle.mkdir(parents=True)
    rows = []
    pep = []
    cds = []
    for index in range(1, 11):
        gene = f"{code}1g{index:04d}"
        start = 100 + index * 300
        rows.append(
            f"1\t{gene}\t{start}\t{start + 89}\t+\t{index}\t{gene}\n"
        )
        pep.append(f">{gene}\n{'M' * 30}\n")
        cds.append(f">{gene}\n{'ATG' * 30}\n")
    (bundle / f"{code}.new.gff").write_text("".join(rows), encoding="utf-8")
    (bundle / f"{code}.lens").write_text("1\t50000\t10\n", encoding="utf-8")
    (bundle / f"{code}.pep").write_text("".join(pep), encoding="utf-8")
    (bundle / f"{code}.cds").write_text("".join(cds), encoding="utf-8")
    return bundle


def test_real_fixture_builder_extracts_selected_chromosome(tmp_path: Path) -> None:
    left = _bundle(tmp_path / "source", "Aaa")
    right = _bundle(tmp_path / "source", "Bbb")
    pal2nal = tmp_path / "pal2nal.pl"
    pal2nal.write_text("#!/usr/bin/env perl\n", encoding="utf-8")
    output = tmp_path / "fixture"
    env = {**os.environ, "WGDI_PAL2NAL": str(pal2nal), "PYTHONIOENCODING": "utf-8"}
    result = subprocess.run(
        [
            sys.executable,
            str(BUILD_REAL_FIXTURE),
            "--workflow-root",
            str(WORKFLOW_ROOT),
            "--left-bundle",
            str(left),
            "--left-chromosome",
            "1",
            "--right-bundle",
            str(right),
            "--right-chromosome",
            "1",
            "--max-genes",
            "5",
            "--output-root",
            str(output),
        ],
        text=True,
        encoding="utf-8",
        capture_output=True,
        check=False,
        env=env,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert (output / "Aaa" / "Aaa.new.gff").read_text(encoding="utf-8").count("\n") == 5
    assert (output / "external" / "pal2nal.v14" / "pal2nal.pl").is_file()
    manifest = json.loads((output / "fixture_manifest.json").read_text(encoding="utf-8"))
    assert manifest["selection"]["method"] == "real preprocessing output subset"


def test_real_validation_collector_combines_project_and_smoke_receipt(
    tmp_path: Path,
) -> None:
    project = tmp_path / "project"
    resources = project / "resources"
    analysis = project / "analysis_v2" / "Ttt"
    resources.mkdir(parents=True)
    analysis.mkdir(parents=True)
    (resources / "wgdi_v2_index.tsv").write_text(
        "show_id\tcode\tspecies\taccession\twgdi_ready\tqc\n"
        "1\tTtt\tTestus plantus\tGCF_TEST\t"
        "yes\tanalysis_v2/Ttt/Ttt.qc.json\n",
        encoding="utf-8",
    )
    (project / "data-validation-report.tsv").write_text(
        "code\tstatus\tgenes\tsequences\nTtt\tPASS\t10\t1\n",
        encoding="utf-8",
    )
    (resources / "naming_audit.tsv").write_text(
        "show_id\tcode\tstatus\n1\tTtt\tPASS\n",
        encoding="utf-8",
    )
    (analysis / "Ttt.qc.json").write_text(
        json.dumps(
            {
                "n_genes": 12,
                "n_transcripts": 14,
                "n_representative_genes": 11,
                "rejection_counts": {"alternative_isoform": 3},
                "inputs": {
                    role: {"sha256": role * 8}
                    for role in ("gff", "cds", "pep", "genome")
                },
            }
        ),
        encoding="utf-8",
    )
    smoke = tmp_path / "smoke.json"
    smoke.write_text(
        json.dumps(
            {
                "status": "PASS",
                "fixture": {"method": "real preprocessing output subset"},
                "comparisons": [{"status": "PASS"}],
            }
        ),
        encoding="utf-8",
    )
    output = tmp_path / "real-report.json"
    result = subprocess.run(
        [
            sys.executable,
            str(COLLECT_REAL),
            "--project-root",
            str(project),
            "--smoke-report",
            str(smoke),
            "--output",
            str(output),
        ],
        text=True,
        encoding="utf-8",
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads(output.read_text(encoding="utf-8"))
    assert report["status"] == "PASS"
    assert report["species"][0]["wgdi_genes"] == 10
