from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from pgci_preprocess import envcheck

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"


def test_pandas3_is_a_blocker_and_foreign_wgdi_a_warning() -> None:
    pins = {"wgdi": "0.75", "pandas": "2.3.3", "numpy": "2.4.6"}
    problems = envcheck.evaluate_conflicts(
        {"pandas": "3.0.6", "wgdi": "0.71", "numpy": "2.4.6"}, pins, "3.8.1"
    )
    by_item = {p["item"]: p for p in problems}
    assert by_item["pandas"]["level"] == "blocker"
    assert by_item["wgdi"]["level"] == "warn"
    assert "numpy" not in by_item  # matches pin -> silent


def test_clean_env_has_no_conflicts() -> None:
    pins = {"wgdi": "0.75", "pandas": "2.3.3"}
    assert envcheck.evaluate_conflicts({"pandas": "2.3.3", "wgdi": "0.75"}, pins) == []


def test_shadowing_source_tree_is_reported() -> None:
    problems = envcheck.evaluate_conflicts(
        {"pgci-preprocess": "3.7.0"}, {"wgdi": "0.75"}, "3.8.1"
    )
    assert any(p["item"] == "pgci-preprocess" and p["level"] == "warn" for p in problems)


def test_parse_requirements_ignores_comments() -> None:
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "req.txt"
        path.write_text("# comment\nwgdi==0.75\npandas==2.3.3\n\n", encoding="utf-8")
        assert envcheck.parse_requirements(path) == {"wgdi": "0.75", "pandas": "2.3.3"}


def _report(plan_inputs: dict) -> dict:
    base = {
        "interpreter": {"reachable": True, "major": 3, "minor": 13},
        "conflicts": [],
        "existing_env": {"present": False},
        "pins": {"wgdi": "0.75", "pandas": "2.3.3"},
        "network": {"reachable": True},
    }
    base.update(plan_inputs)
    return base


def test_choose_plan_reuse_install_offline_blocked() -> None:
    assert envcheck.choose_plan(_report({})) == "install"
    assert envcheck.choose_plan(_report({"network": {"reachable": False}})) == "offline"
    assert envcheck.choose_plan(_report({
        "existing_env": {"present": True, "wgdi": "0.75", "pandas": "2.3.3"},
    })) == "reuse"
    assert envcheck.choose_plan(_report({
        "conflicts": [{"level": "blocker", "item": "pandas", "detail": "x"}],
    })) == "blocked"
    assert envcheck.choose_plan(_report({
        "interpreter": {"reachable": True, "major": 3, "minor": 10},
    })) == "blocked"
    # an existing env with pandas 3 must not be reused
    assert envcheck.choose_plan(_report({
        "existing_env": {"present": True, "wgdi": "0.75", "pandas": "3.0.6"},
    })) != "reuse"


def test_skipped_probe_is_not_claimed_reachable(tmp_path: Path) -> None:
    report = envcheck.scan(tmp_path, probe_index=False)
    assert report["network"]["skipped"] is True
    assert report["plan"] == "offline"
    text = envcheck.render(report)
    assert "未探测" in text
    assert "索引不可达 →" not in text


def test_preflight_cli_reports_and_exits_by_plan(tmp_path: Path) -> None:
    env = {"PYTHONPATH": str(SRC_DIR), "PYTHONIOENCODING": "utf-8"}
    import os

    result = subprocess.run(
        [sys.executable, "-m", "pgci_preprocess", "preflight",
         "--project-root", str(tmp_path), "--no-network", "--json"],
        capture_output=True, text=True, encoding="utf-8", errors="replace",
        env={**os.environ, **env}, timeout=120,
    )
    assert result.returncode in {0, 2, 3}, result.stderr
    report = json.loads(result.stdout)
    assert report["plan"] in {"reuse", "install", "offline", "blocked"}
    assert "interpreter" in report and "conflicts" in report


def test_bootstrap_check_only_does_not_install(tmp_path: Path) -> None:
    import os

    tool = PROJECT_ROOT / "tools" / "bootstrap_wgdi.py"
    result = subprocess.run(
        [sys.executable, str(tool), "--project-root", str(tmp_path),
         "--check-only", "--offline"],
        capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=120,
        env={**os.environ, "PYTHONIOENCODING": "utf-8"},
    )
    assert result.returncode in {0, 2, 3}, result.stderr
    assert "环境预扫" in result.stdout
    assert not (tmp_path / ".wgdi-env").exists()


def test_bootstrap_offline_creates_nothing_and_gives_source_mode(tmp_path: Path) -> None:
    import os

    tool = PROJECT_ROOT / "tools" / "bootstrap_wgdi.py"
    result = subprocess.run(
        [sys.executable, str(tool), "--project-root", str(tmp_path), "--offline"],
        capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=120,
        env={**os.environ, "PYTHONIOENCODING": "utf-8"},
    )
    assert result.returncode == 0, result.stderr
    assert "PYTHONPATH" in result.stdout
    assert "--no-build-isolation" in result.stdout
    assert not (tmp_path / ".wgdi-env").exists()


def test_classify_failure_buckets() -> None:
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location("bw", PROJECT_ROOT / "tools" / "bootstrap_wgdi.py")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.classify_failure("Could not fetch URL https://pypi.org") == "network"
    assert module.classify_failure("Connection timed out") == "network"
    assert module.classify_failure("No matching distribution found for numpy==2.4.6") == "wheel"
    assert module.classify_failure("ResolutionImpossible: conflict") == "conflict"
    assert module.classify_failure("something else") == "unknown"
