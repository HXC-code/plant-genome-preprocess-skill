from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


WORKFLOW_ROOT = Path(__file__).resolve().parents[2]
GENERATOR = WORKFLOW_ROOT / "code" / "tools" / "generate_wgdi_configs.py"
SYNTHETIC_FIXTURE = (
    WORKFLOW_ROOT / "code" / "tools" / "build_synthetic_wgdi_smoke.py"
)


def _workspace(tmp_path: Path) -> tuple[Path, Path]:
    project = tmp_path / "workflow"
    project.mkdir()
    (project / "folder_list.txt").write_text("A\nB\n", encoding="utf-8")
    for left, right in (("A", "A"), ("A", "B"), ("B", "A"), ("B", "B")):
        (project / "Auto_wgdi" / f"{left}_base" / f"{left}_{right}").mkdir(
            parents=True
        )
    fake_wgdi = project / "bin" / "wgdi"
    fake_wgdi.parent.mkdir()
    fake_wgdi.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
    return project, fake_wgdi


def _generate(
    project: Path, fake_wgdi: Path, mode: str
) -> subprocess.CompletedProcess[str]:
    env = {
        **os.environ,
        "WGDI": str(fake_wgdi),
        "WGDI_TANDEM": "false",
        "WGDI_EVALUE": "1e-5",
        "PYTHONIOENCODING": "utf-8",
    }
    return subprocess.run(
        [
            sys.executable,
            str(GENERATOR),
            "--project-root",
            str(project),
            "--mode",
            mode,
        ],
        text=True,
        encoding="utf-8",
        capture_output=True,
        check=False,
        env=env,
    )


def test_compatibility_wrappers_share_one_implementation() -> None:
    blast = (WORKFLOW_ROOT / "1run_blast.sh").read_text(encoding="utf-8")
    diamond = (WORKFLOW_ROOT / "1run_diamond.sh").read_text(encoding="utf-8")
    assert "run_homology.sh" in blast
    assert "HOMOLOGY_ENGINE=blast" in blast
    assert "run_homology.sh" in diamond
    assert "HOMOLOGY_ENGINE=diamond" in diamond
    assert len(blast.splitlines()) <= 6
    assert len(diamond.splitlines()) <= 6

    different = (WORKFLOW_ROOT / "4run_wgdi_different.sh").read_text(encoding="utf-8")
    same = (WORKFLOW_ROOT / "5run_wgdi_same.sh").read_text(encoding="utf-8")
    assert "generate_wgdi_configs.py" in different
    assert "--mode different" in different
    assert "generate_wgdi_configs.py" in same
    assert "--mode same" in same
    assert len(different.splitlines()) <= 7
    assert len(same.splitlines()) <= 7


def test_different_config_is_generated_once_and_removes_tandem(
    tmp_path: Path,
) -> None:
    project, fake_wgdi = _workspace(tmp_path)
    result = _generate(project, fake_wgdi, "different")
    assert result.returncode == 0, result.stdout + result.stderr
    config = (
        project / "Auto_wgdi" / "A_base" / "A_B" / "total.conf"
    ).read_text(encoding="utf-8")
    assert config.count("tandem = false") == 3
    assert "ks_col = ks_NG86" in config
    assert "cds_file = all.cds" in config
    assert "gff1 = A.new.gff" in config
    run = (project / "Auto_wgdi" / "A_base" / "A_B" / "run.sh").read_text(
        encoding="utf-8"
    )
    assert "WGDI_TANDEM_ARGS=(--remove-tandem)" in run
    assert "run_stage dotplot" in run
    assert not (project / "Auto_wgdi" / "A_base" / "A_A" / "total.conf").exists()


def test_same_config_uses_single_species_and_is_generated_once(tmp_path: Path) -> None:
    project, fake_wgdi = _workspace(tmp_path)
    result = _generate(project, fake_wgdi, "same")
    assert result.returncode == 0, result.stdout + result.stderr
    config = (
        project / "Auto_wgdi" / "B_base" / "B_B" / "total.conf"
    ).read_text(encoding="utf-8")
    assert config.count("tandem = false") == 3
    assert "cds_file = B.cds" in config
    assert "gff1 = B.new.gff" in config
    assert "gff2 = B.new.gff" in config
    assert "blast = B_B.blast" in config
    assert not (project / "Auto_wgdi" / "A_base" / "A_B" / "total.conf").exists()


def test_wgdi_entry_is_the_only_runtime_compatibility_patch() -> None:
    text = (WORKFLOW_ROOT / "bin" / "wgdi_entry.py").read_text(encoding="utf-8")
    assert "isinstance(value, bool)" in text
    assert "_load_conf_remove_tandem" in text
    assert not (WORKFLOW_ROOT / "vendor").exists()


def test_bootstrap_wrappers_share_one_python_implementation() -> None:
    shell = (WORKFLOW_ROOT / "bootstrap_wgdi.sh").read_text(encoding="utf-8")
    powershell = (WORKFLOW_ROOT / "bootstrap_wgdi.ps1").read_text(encoding="utf-8")
    assert "code/tools/bootstrap_wgdi.py" in shell
    assert "code/tools/bootstrap_wgdi.py" in powershell
    assert len(shell.splitlines()) < 35
    assert len(powershell.splitlines()) < 35


def test_synthetic_wgdi_fixture_is_self_contained(tmp_path: Path) -> None:
    output = tmp_path / "synthetic"
    env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
    result = subprocess.run(
        [
            sys.executable,
            str(SYNTHETIC_FIXTURE),
            "--workflow-root",
            str(WORKFLOW_ROOT),
            "--output-root",
            str(output),
            "--genes-per-species",
            "20",
            "--chromosomes",
            "2",
            "--protein-length",
            "40",
        ],
        text=True,
        encoding="utf-8",
        capture_output=True,
        check=False,
        env=env,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert (output / "fixture_manifest.json").is_file()
    assert (output / "run_homology.sh").is_file()
    assert (output / "code" / "tools" / "pal2nal.v14" / "pal2nal.pl").is_file()
    assert (output / "SynA" / "SynA.pep").is_file()
    assert (output / "SynB" / "SynB.new.gff").is_file()
