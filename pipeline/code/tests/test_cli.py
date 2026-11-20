from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if not SRC_DIR.is_dir():
    SRC_DIR = PROJECT_ROOT / "scripts"


def _expected_version() -> str:
    """Read the version the source tree declares, independent of the runner."""

    text = (SRC_DIR / "pgci_preprocess" / "version.py").read_text(encoding="utf-8")
    match = re.search(r'^__version__\s*=\s*"([^"]+)"', text, re.MULTILINE)
    assert match, "version.py must define __version__"
    return match.group(1)


def _run(args: list[str]) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SRC_DIR)
    return subprocess.run(
        [sys.executable, "-m", "pgci_preprocess", *args],
        cwd=PROJECT_ROOT,
        env=env,
        text=True,
        capture_output=True,
        timeout=120,
    )


@pytest.mark.parametrize(
    "command",
    ["prepare", "init-project", "doctor", "normalize-gff", "convert-annotation", "run", "validate", "batch", "audit", "convert", "source-manifest", "migrate-species-names"],
)
def test_subcommand_help(command: str) -> None:
    result = _run([command, "--help"])
    assert result.returncode == 0
    assert "usage: pgci_preprocess" in result.stdout


def test_top_level_help_and_version() -> None:
    help_result = _run(["--help"])
    assert help_result.returncode == 0
    assert "prepare" in help_result.stdout
    version_result = _run(["--version"])
    assert version_result.returncode == 0
    assert _expected_version() in version_result.stdout


def test_module_entry_works_without_pythonpath() -> None:
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    result = subprocess.run(
        [sys.executable, "-m", "pgci_preprocess", "--version"],
        cwd=SRC_DIR,
        env=env,
        text=True,
        capture_output=True,
        timeout=120,
    )
    assert result.returncode == 0
    assert _expected_version() in result.stdout


def test_unknown_command_is_an_error() -> None:
    result = _run(["not-a-command"])
    assert result.returncode == 2
    assert "unknown command" in result.stderr


def test_no_command_is_an_error() -> None:
    result = _run([])
    assert result.returncode == 2
    assert "usage: pgci_preprocess" in result.stdout


def test_run_with_empty_project_returns_two(tmp_path: Path) -> None:
    project = tmp_path / "empty"
    (project / "data" / "downloads").mkdir(parents=True)
    (project / "resources").mkdir(parents=True)
    (project / "resources" / "download_inventory.tsv").write_text(
        "showid\tOrder\tFamily\tSpecies\tgff\tgenome\tcds\tpep\n",
        encoding="utf-8",
    )
    result = _run(["run", "--project-root", str(project)])
    assert result.returncode == 2
    assert "no species selected" in result.stderr
    assert not (project / "resources" / "wgdi_v2_manifest_sha256.tsv").exists()


def test_normalize_gff_without_species_returns_two() -> None:
    result = _run(["normalize-gff", "--species", ""])
    assert result.returncode == 2
