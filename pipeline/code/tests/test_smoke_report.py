from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
COLLECTOR = PROJECT_ROOT / "tools" / "collect_wgdi_smoke_report.py"
if not COLLECTOR.is_file():
    import pytest
    pytest.skip("WGDI smoke tools are only distributed with the full workflow", allow_module_level=True)


def _build_receipt_fixture(root: Path) -> None:
    comparison = root / "wgdi_output" / "results"
    comparison.mkdir(parents=True)
    state = root / "wgdi_output" / "state"
    state.mkdir(parents=True)
    (root / "fixture_manifest.json").write_text(
        json.dumps({"selection": {"selected_genes": 10, "selected_sequences": 1}}),
        encoding="utf-8",
    )
    settings = root / "wgdi_input" / "blast" / "diamond_run_settings.tsv"
    settings.parent.mkdir(parents=True)
    settings.write_text(
        "key\tvalue\nexecutable\tdiamond\nversion\ttest\nevalue\t1e-5\nsensitivity\tvery-sensitive\n",
        encoding="utf-8",
    )
    # Single authoritative BLAST table under wgdi_input/blast; markers carry
    # the resume fingerprints and live in wgdi_output/state.
    source_blast = root / "wgdi_input" / "blast" / "s1_pep" / "s1_s1.blast"
    source_blast.parent.mkdir(parents=True)
    source_blast.write_text("q\ts\t100\t3\t0\t0\t1\t3\t1\t3\t1e-5\t10\n", encoding="utf-8")
    for stage in (
        "dotplot",
        "collinearity",
        "ks",
        "blockinfo",
        "correspondence",
        "blockks",
        "kspeaks",
        "peaksfit",
    ):
        (state / f"s1.done.{stage}").write_text("core;blast=24-1000\n", encoding="utf-8")
    for name in ("s1.png", "s1.icl", "s1.ks", "s1.bi.csv", "s1.bi.c.csv", "s1.bk.png", "s1.kp.png", "s1.pf.png"):
        (comparison / name).write_text("nonempty\n", encoding="utf-8")
    (root / "smoke_status.tsv").write_text(
        "status\tPASS\ndiamond_evalue\t1e-5\ndiamond_sensitivity\tvery-sensitive\n",
        encoding="utf-8",
    )


def test_smoke_receipt_requires_all_outputs_and_single_blast_copy(tmp_path: Path) -> None:
    root = tmp_path / "smoke"
    root.mkdir()
    _build_receipt_fixture(root)
    passed = subprocess.run(
        [sys.executable, str(COLLECTOR), str(root)],
        text=True,
        capture_output=True,
        timeout=120,
    )
    assert passed.returncode == 0, passed.stdout + passed.stderr
    report = json.loads((root / "smoke_report.json").read_text(encoding="utf-8"))
    assert report["status"] == "PASS"
    assert report["schema_version"] == "wgdi-smoke-v4"
    assert set(report["environment"]["tools"]) == {
        "diamond",
        "blastp",
        "makeblastdb",
        "muscle",
        "yn00",
        "perl",
        "pal2nal",
    }
    assert len(report["comparisons"][0]["outputs"]) == 8
    assert report["comparisons"][0]["auto_pep_blast_present"] is True

    (root / "wgdi_output" / "results" / "s1.pf.png").unlink()
    failed = subprocess.run(
        [sys.executable, str(COLLECTOR), str(root)],
        text=True,
        capture_output=True,
        timeout=120,
    )
    assert failed.returncode == 1
    assert json.loads((root / "smoke_report.json").read_text(encoding="utf-8"))["status"] == "FAIL"


def test_smoke_receipt_core_mode_requires_only_dotplot_and_collinearity(
    tmp_path: Path,
) -> None:
    root = tmp_path / "smoke_core"
    root.mkdir()
    _build_receipt_fixture(root)
    (root / "smoke_status.tsv").write_text(
        "status\tPASS\ndiamond_evalue\t1e-5\ndiamond_sensitivity\tvery-sensitive\nks_enabled\tfalse\n",
        encoding="utf-8",
    )
    comparison = root / "wgdi_output" / "results"
    for name in ("s1.ks", "s1.bi.csv", "s1.bi.c.csv", "s1.bk.png", "s1.kp.png", "s1.pf.png"):
        (comparison / name).unlink()
    result = subprocess.run(
        [sys.executable, str(COLLECTOR), str(root)],
        text=True,
        capture_output=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads((root / "smoke_report.json").read_text(encoding="utf-8"))
    assert report["status"] == "PASS"
    assert report["ks_enabled"] is False
    assert len(report["comparisons"][0]["outputs"]) == 2
