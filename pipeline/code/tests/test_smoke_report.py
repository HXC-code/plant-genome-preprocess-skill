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
    comparison = root / "Auto_wgdi" / "s1_base" / "s1_s1"
    comparison.mkdir(parents=True)
    (root / "fixture_manifest.json").write_text(
        json.dumps({"selection": {"selected_genes": 10, "selected_sequences": 1}}),
        encoding="utf-8",
    )
    settings = root / "Auto_pep" / "diamond_run_settings.tsv"
    settings.parent.mkdir(parents=True)
    settings.write_text(
        "key\tvalue\nexecutable\tdiamond\nversion\ttest\nevalue\t1e-5\nsensitivity\tvery-sensitive\n",
        encoding="utf-8",
    )
    source_blast = root / "Auto_pep" / "s1_pep" / "s1_s1.blast"
    source_blast.parent.mkdir(parents=True)
    source_blast.write_text("q\ts\t100\t3\t0\t0\t1\t3\t1\t3\t1e-5\t10\n", encoding="utf-8")
    (comparison / "s1_s1.blast").write_bytes(source_blast.read_bytes())
    (comparison / ".wgdi.done").write_text("done\n", encoding="utf-8")
    for name in ("s1.png", "s1.icl", "s1.ks", "s1.bi.csv", "s1.bi.c.csv", "s1.bk.png", "s1.kp.png", "s1.pf.png"):
        (comparison / name).write_text("nonempty\n", encoding="utf-8")
    (root / "smoke_status.tsv").write_text(
        "status\tPASS\ndiamond_evalue\t1e-5\ndiamond_sensitivity\tvery-sensitive\n",
        encoding="utf-8",
    )


def test_smoke_receipt_requires_all_outputs_and_preserves_blast(tmp_path: Path) -> None:
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
    assert report["schema_version"] == "wgdi-smoke-v2"
    assert set(report["environment"]["tools"]) == {
        "diamond",
        "blastp",
        "makeblastdb",
        "muscle",
        "yn00",
        "perl",
    }
    assert len(report["comparisons"][0]["outputs"]) == 8
    assert report["comparisons"][0]["auto_pep_blast_preserved_and_identical"] is True

    (root / "Auto_wgdi" / "s1_base" / "s1_s1" / "s1.pf.png").unlink()
    failed = subprocess.run(
        [sys.executable, str(COLLECTOR), str(root)],
        text=True,
        capture_output=True,
        timeout=120,
    )
    assert failed.returncode == 1
    assert json.loads((root / "smoke_report.json").read_text(encoding="utf-8"))["status"] == "FAIL"
