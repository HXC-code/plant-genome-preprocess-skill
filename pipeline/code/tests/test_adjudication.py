from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from pgci_preprocess.adjudication import (
    gate_status,
    read_adjudication,
    write_adjudication,
)
from pgci_preprocess.pipeline import PipelineOptions, run_project
from pgci_preprocess.validate import ValidationOptions, validate_project
from tests.test_pipeline import _make_project

SRC = Path(__file__).resolve().parents[1] / "src"


def _qc_project(tmp_path: Path, n: int = 7) -> Path:
    project = tmp_path / "proj"
    qc_dir = project / "wgdi_input" / "X"
    qc_dir.mkdir(parents=True)
    (qc_dir / "X.qc.json").write_text(
        json.dumps(
            {
                "code": "X",
                "status": "OK",
                "pipeline_version": "3.6.0",
                "n_wgdi_sequences": n,
            }
        ),
        encoding="utf-8",
    )
    return project


def test_missing_then_ok(tmp_path: Path) -> None:
    project = _qc_project(tmp_path)
    assert gate_status(project, "X") == "missing"
    write_adjudication(project, "X", confirmed_n=7, decision="accept", reviewer="reviewer-A")
    assert gate_status(project, "X") == "ok"
    record = read_adjudication(project, "X")
    assert record["confirmed_n"] == 7
    assert record["reviewer"] == "reviewer-A"
    assert record["qc_sha256"]


def test_count_mismatch_blocks(tmp_path: Path) -> None:
    project = _qc_project(tmp_path)
    write_adjudication(project, "X", confirmed_n=8, decision="accept", reviewer="T")
    assert gate_status(project, "X") == "count_mismatch"


def test_reject_and_adjust_block(tmp_path: Path) -> None:
    project = _qc_project(tmp_path)
    write_adjudication(project, "X", confirmed_n=7, decision="reject", reviewer="T")
    assert gate_status(project, "X") == "rejected"
    write_adjudication(project, "X", confirmed_n=7, decision="adjust", reviewer="T")
    assert gate_status(project, "X") == "pending_adjustment"


def test_rerun_invalidates_stale_signoff(tmp_path: Path) -> None:
    project = _qc_project(tmp_path)
    write_adjudication(project, "X", confirmed_n=7, decision="accept", reviewer="T")
    assert gate_status(project, "X") == "ok"
    qc_path = project / "wgdi_input" / "X" / "X.qc.json"
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    qc["recomputed"] = True
    qc_path.write_text(json.dumps(qc), encoding="utf-8")
    assert gate_status(project, "X") == "stale"


def test_reviewer_required(tmp_path: Path) -> None:
    project = _qc_project(tmp_path)
    with pytest.raises(ValueError, match="reviewer"):
        write_adjudication(
            project, "X", confirmed_n=7, decision="accept", reviewer="  "
        )


def test_validate_gate_requires_adjudication(tmp_path: Path) -> None:
    project, _ = _make_project(tmp_path)
    run_project(
        project,
        options=PipelineOptions(project_root=project, min_genes_per_sequence=0),
    )
    strict = validate_project(
        ValidationOptions(
            project_root=project,
            require_compat_links=False,
            require_adjudication=True,
        )
    )
    assert not strict.ok
    assert any("adjudication_missing" in str(row["problems"]) for row in strict.rows)

    code = next((project / "wgdi_input").iterdir()).name
    qc = json.loads(
        (project / "wgdi_input" / code / f"{code}.qc.json").read_text(encoding="utf-8")
    )
    write_adjudication(
        project,
        code,
        confirmed_n=qc["n_wgdi_sequences"],
        decision="accept",
        reviewer="reviewer-A",
    )
    accepted = validate_project(
        ValidationOptions(
            project_root=project,
            require_compat_links=False,
            require_adjudication=True,
        )
    )
    assert accepted.ok


def test_cli_chromosomes_adjudicate_status(tmp_path: Path) -> None:
    project, _ = _make_project(tmp_path)
    run_project(
        project,
        options=PipelineOptions(project_root=project, min_genes_per_sequence=0),
    )
    code = next((project / "wgdi_input").iterdir()).name
    qc = json.loads(
        (project / "wgdi_input" / code / f"{code}.qc.json").read_text(encoding="utf-8")
    )
    env = {**os.environ, "PYTHONPATH": str(SRC), "PYTHONIOENCODING": "utf-8"}

    def run(*args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, "-m", "pgci_preprocess", *args],
            text=True,
            encoding="utf-8",
            capture_output=True,
            check=False,
            env=env,
        )

    chrom = run("chromosomes", "--project-root", str(project), "--species", code)
    assert chrom.returncode == 0, chrom.stderr
    assert "条染色体" in chrom.stdout

    status = run(
        "adjudication-status", "--project-root", str(project), "--species", code
    )
    assert status.stdout.strip() == f"{code}\tmissing"

    sign = run(
        "adjudicate",
        "--project-root",
        str(project),
        "--species",
        code,
        "--confirmed-n",
        str(qc["n_wgdi_sequences"]),
        "--decision",
        "accept",
        "--reviewer",
        "reviewer-A",
    )
    assert sign.returncode == 0, sign.stderr

    status2 = run(
        "adjudication-status", "--project-root", str(project), "--species", code
    )
    assert status2.stdout.strip() == f"{code}\tok"