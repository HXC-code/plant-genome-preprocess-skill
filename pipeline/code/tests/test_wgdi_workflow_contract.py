from __future__ import annotations

import json
import os
import runpy
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


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
    assert len(different.splitlines()) <= 9
    assert len(same.splitlines()) <= 9


def test_generated_run_script_uses_a_bash_safe_windows_path(tmp_path: Path) -> None:
    namespace = runpy.run_path(str(GENERATOR))
    bash_path = namespace["_bash_path"]
    assert bash_path(Path(r"D:\project\bin\wgdi")) == "/d/project/bin/wgdi"
    project, fake_wgdi = _workspace(tmp_path)
    result = _generate(project, fake_wgdi, "same")
    assert result.returncode == 0, result.stdout + result.stderr
    run_text = (
        project / "Auto_wgdi" / "A_base" / "A_A" / "run.sh"
    ).read_text(encoding="utf-8")
    assert "\\\\" not in run_text.split("WGDI_ENTRY=", 1)[1].splitlines()[0]


def test_homology_skips_not_ready_species_and_infrastructure() -> None:
    text = (WORKFLOW_ROOT / "run_homology.sh").read_text(encoding="utf-8")
    assert "wgdi_ready=no" in text
    assert "analysis_v2|bin|clean_v2|code|data" in text
    assert "HOMOLOGY_ALLOW_NOT_READY" in text
    assert "cleanup_tmp_map" in text
    assert "|| true" in text
    assert "part.${BASHPID:-$$}" in text
    assert "pgci_preprocess readiness" in text
    assert "species_readiness_grep" in text
    six = (WORKFLOW_ROOT / "6run_wgdi.sh").read_text(encoding="utf-8")
    assert "check-wgdi-outputs" in six


def _homology_bash() -> Path | None:
    candidates = [
        Path(os.environ.get("GIT_BASH", "")),
        Path("D:/Git/bin/bash.exe"),
        Path("C:/Program Files/Git/bin/bash.exe"),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    bash = shutil.which("bash")
    return Path(bash) if bash and "System32" not in bash else None


def _homology_workspace(tmp_path: Path) -> tuple[Path, Path]:
    project = tmp_path / "homology"
    project.mkdir()
    for code, ready in (("Rdy", "yes"), ("Not", "no"), ("Mar", "marginal")):
        species = project / code
        species.mkdir()
        for suffix in (".pep", ".cds", ".new.gff", ".lens"):
            (species / f"{code}{suffix}").write_text("x\n", encoding="utf-8")
        qc_dir = project / "analysis_v2" / code
        qc_dir.mkdir(parents=True)
        (qc_dir / f"{code}.qc.json").write_text(
            json.dumps({"sequence_selection": {"wgdi_ready": ready}}),
            encoding="utf-8",
        )
    # A species with the four-piece set but no readiness record at all.
    unknown = project / "Unk"
    unknown.mkdir()
    for suffix in (".pep", ".cds", ".new.gff", ".lens"):
        (unknown / f"Unk{suffix}").write_text("x\n", encoding="utf-8")
    # Infrastructure directories must never be scanned as species.
    (project / "bin").mkdir()
    fake_diamond = project / "fake_diamond"
    fake_diamond.write_bytes(
        b"#!/usr/bin/env bash\n"
        b'case "$1" in\n'
        b'  version) echo "diamond version 2.2.8 (fake)" ;;\n'
        b"  makedb) : ;;\n"
        b"  blastp)\n"
        b'    out=""; prev=""\n'
        b'    for a in "$@"; do\n'
        b'      [[ "$prev" == "--out" ]] && out="$a"\n'
        b'      prev="$a"\n'
        b"    done\n"
        b"    printf 'q1\\ts1\\t100.0\\t50\\t0\\t0\\t1\\t50\\t1\\t50\\t1e-30\\t100\\n' > \"$out\" ;;\n"
        b"esac\n"
        b"exit 0\n"
    )
    fake_diamond.chmod(0o755)
    shutil.copy2(WORKFLOW_ROOT / "run_homology.sh", project / "run_homology.sh")
    return project, fake_diamond


def test_homology_readiness_gate_executes(tmp_path: Path) -> None:
    bash = _homology_bash()
    if bash is None:
        pytest.skip("Git Bash or Linux bash is required to execute run_homology.sh")
    project, fake_diamond = _homology_workspace(tmp_path)
    env = {
        **os.environ,
        "DIAMOND": str(fake_diamond),
        "HOMOLOGY_ENGINE": "diamond",
        "PYTHONIOENCODING": "utf-8",
    }
    result = subprocess.run(
        [str(bash), str(project / "run_homology.sh")],
        cwd=project,
        text=True,
        encoding="utf-8",
        capture_output=True,
        check=False,
        env=env,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    mapped = (project / "species_map.tsv").read_text(encoding="utf-8")
    assert "Rdy\t" in mapped and "Mar\t" in mapped
    assert "Not\t" not in mapped and "Unk\t" not in mapped
    skipped = (project / "skipped_species.txt").read_text(encoding="utf-8")
    assert "Not：wgdi_ready=no" in skipped
    assert "Unk：缺少可读取的 wgdi_ready 记录" in skipped
    assert "bin" not in skipped

    override = subprocess.run(
        [str(bash), str(project / "run_homology.sh")],
        cwd=project,
        text=True,
        encoding="utf-8",
        capture_output=True,
        check=False,
        env={**env, "HOMOLOGY_ALLOW_NOT_READY": "true"},
    )
    assert override.returncode == 0, override.stdout + override.stderr
    mapped = (project / "species_map.tsv").read_text(encoding="utf-8")
    for code in ("Rdy", "Not", "Mar", "Unk"):
        assert f"{code}\t" in mapped


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


def test_wgdi_launcher_normalizes_windows_tmpdir() -> None:
    text = (WORKFLOW_ROOT / "bin" / "wgdi").read_text(encoding="utf-8")
    assert "normalize_shell_path" in text
    assert r"^([A-Za-z]):[\\/](.*)$" in text
    assert 'export TMPDIR="$tmp_root"' in text


def test_bootstrap_wrappers_share_one_python_implementation() -> None:
    shell = (WORKFLOW_ROOT / "bootstrap_wgdi.sh").read_text(encoding="utf-8")
    powershell = (WORKFLOW_ROOT / "bootstrap_wgdi.ps1").read_text(encoding="utf-8")
    assert "code/tools/bootstrap_wgdi.py" in shell
    assert "code/tools/bootstrap_wgdi.py" in powershell
    assert len(shell.splitlines()) < 35
    assert len(powershell.splitlines()) < 35


def test_synthetic_wgdi_fixture_is_self_contained(tmp_path: Path) -> None:
    output = tmp_path / "synthetic"
    fake_pal2nal = tmp_path / "pal2nal.pl"
    fake_pal2nal.write_text("#!/usr/bin/env perl\n", encoding="utf-8")
    env = {
        **os.environ,
        "PYTHONIOENCODING": "utf-8",
        "WGDI_PAL2NAL": str(fake_pal2nal),
    }
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
    assert (output / "external" / "pal2nal.v14" / "pal2nal.pl").is_file()
    assert (output / "SynA" / "SynA.pep").is_file()
    assert (output / "SynB" / "SynB.new.gff").is_file()
    for code in ("SynA", "SynB"):
        qc = json.loads(
            (output / "analysis_v2" / code / f"{code}.qc.json").read_text(
                encoding="utf-8"
            )
        )
        assert qc["sequence_selection"]["wgdi_ready"] == "yes"
