from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from pgci_preprocess.readiness import (
    extract_qc_readiness,
    read_readiness_map,
    read_species_readiness,
)


def test_extract_qc_readiness_validates_value() -> None:
    assert extract_qc_readiness({"sequence_selection": {"wgdi_ready": "no"}}) == "no"
    assert (
        extract_qc_readiness({"sequence_selection": {"wgdi_ready": "marginal"}})
        == "marginal"
    )
    assert extract_qc_readiness({"sequence_selection": {"wgdi_ready": "maybe"}}) == ""
    assert extract_qc_readiness({"sequence_selection": {}}) == ""
    assert extract_qc_readiness({"sequence_selection": "no"}) == ""
    assert extract_qc_readiness({}) == ""


def test_stray_top_level_key_is_not_misread(tmp_path: Path) -> None:
    # A text grep would match the first "wgdi_ready" anywhere; the structural
    # reader must use sequence_selection, not a stray top-level/nested key.
    qc = {
        "wgdi_ready": "yes",
        "sequence_selection": {"wgdi_ready": "no", "ratio": 0.1},
    }
    qc_dir = tmp_path / "analysis_v2" / "X"
    qc_dir.mkdir(parents=True)
    (qc_dir / "X.qc.json").write_text(json.dumps(qc), encoding="utf-8")
    assert read_species_readiness(tmp_path, "X") == "no"


def test_read_species_readiness_from_qc(tmp_path: Path) -> None:
    qc_dir = tmp_path / "analysis_v2" / "Y"
    qc_dir.mkdir(parents=True)
    (qc_dir / "Y.qc.json").write_text(
        json.dumps({"sequence_selection": {"wgdi_ready": "marginal"}}),
        encoding="utf-8",
    )
    assert read_species_readiness(tmp_path, "Y") == "marginal"


def test_read_species_readiness_index_fallback(tmp_path: Path) -> None:
    resources = tmp_path / "resources"
    resources.mkdir()
    (resources / "wgdi_v2_index.tsv").write_text(
        "show_id\tcode\twgdi_ready\n1\tZ\tno\n", encoding="utf-8"
    )
    assert read_species_readiness(tmp_path, "Z") == "no"
    assert read_readiness_map(tmp_path, ["Z", "Absent"]) == {
        "Z": "no",
        "Absent": "",
    }


def test_read_species_readiness_missing(tmp_path: Path) -> None:
    assert read_species_readiness(tmp_path, "Nobody") == ""


def test_cli_readiness_emits_tsv(tmp_path: Path) -> None:
    for code, ready in (("A", "yes"), ("B", "no")):
        qc_dir = tmp_path / "analysis_v2" / code
        qc_dir.mkdir(parents=True)
        (qc_dir / f"{code}.qc.json").write_text(
            json.dumps({"sequence_selection": {"wgdi_ready": ready}}),
            encoding="utf-8",
        )
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pgci_preprocess",
            "readiness",
            "--project-root",
            str(tmp_path),
            "--species",
            "A,B,C",
        ],
        text=True,
        encoding="utf-8",
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ["A\tyes", "B\tno", "C\t"]