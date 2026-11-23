from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from pgci_preprocess.adjudication import write_adjudication


WORKFLOW_ROOT = Path(__file__).resolve().parents[2]
GENERATOR = WORKFLOW_ROOT / "code" / "tools" / "generate_wgdi_configs.py"
SYNTHETIC_FIXTURE = (
    WORKFLOW_ROOT / "code" / "tools" / "build_synthetic_wgdi_smoke.py"
)


def _workspace(tmp_path: Path) -> Path:
    """A minimal homology-output workspace: scan map + bundles + BLAST tables."""
    project = tmp_path / "workflow"
    project.mkdir()
    meta = project / "metadata"
    meta.mkdir()
    (meta / "folder_list.txt").write_text("A\nB\n", encoding="utf-8")
    (meta / "species_map.tsv").write_text(
        "A\twgdi_input/A\tA.pep\tA.cds\tA.new.gff\tA.lens\n"
        "B\twgdi_input/B\tB.pep\tB.cds\tB.new.gff\tB.lens\n",
        encoding="utf-8",
    )
    for code in ("A", "B"):
        bundle = project / "wgdi_input" / code
        bundle.mkdir(parents=True)
        for suffix in ("pep", "cds", "new.gff", "lens"):
            (bundle / f"{code}.{suffix}").write_text("x\n", encoding="utf-8")
    for left in ("A", "B"):
        pep_dir = project / "wgdi_input" / "blast" / f"{left}_pep"
        pep_dir.mkdir(parents=True)
        for right in ("A", "B"):
            (pep_dir / f"{left}_{right}.blast").write_text(
                "q\ts\t100\t3\t0\t0\t1\t3\t1\t3\t1e-5\t10\n", encoding="utf-8"
            )
    return project


def _generate(
    project: Path, mode: str, ks: str = "false"
) -> subprocess.CompletedProcess[str]:
    env = {
        **os.environ,
        "WGDI_TANDEM": "false",
        "WGDI_EVALUE": "1e-5",
        "WGDI_KS": ks,
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

    gen = (WORKFLOW_ROOT / "2gen_wgdi_conf.sh").read_text(encoding="utf-8")
    assert "generate_wgdi_configs.py" in gen
    assert "--mode" in gen
    # Thin entry point: argument check + one delegation call, nothing else.
    assert len(gen.splitlines()) <= 32
    for legacy in (
        "2run_pre_wgdi.sh",
        "3blast_to_wgdi.sh",
        "4run_wgdi_different.sh",
        "5run_wgdi_same.sh",
        "6run_wgdi.sh",
    ):
        assert not (WORKFLOW_ROOT / legacy).exists(), legacy


def test_runner_contract() -> None:
    runner = (WORKFLOW_ROOT / "3run_wgdi.sh").read_text(encoding="utf-8")
    assert "check-wgdi-outputs" in runner
    assert "species_map.tsv" in runner  # bundle locations from the scan map
    assert "wgdi_input/config" in runner
    assert "wgdi_output/results" in runner
    assert 'run_stage "$prefix" dotplot' in runner
    assert 'run_stage "$prefix" collinearity' in runner
    assert 'run_stage "$prefix" ks ' in runner
    assert 'run_stage "$prefix" peaksfit' in runner
    assert "--remove-tandem" in runner
    assert "$state_dir/all.cds" in runner  # combined Ks inputs rebuilt when stale
    assert 'content_digest "$config"' in runner
    assert 'core;config=' in runner
    assert 'launcher=' in runner


def test_shell_entry_points_are_lf_and_syntax_valid() -> None:
    paths = [
        WORKFLOW_ROOT / name
        for name in (
            "0check_data.sh",
            "_pgci_env.sh",
            "1run_blast.sh",
            "1run_diamond.sh",
            "2gen_wgdi_conf.sh",
            "3run_wgdi.sh",
            "bootstrap_wgdi.sh",
            "run_homology.sh",
        )
    ]
    paths += [WORKFLOW_ROOT / "bin" / "wgdi", Path(__file__).parents[1] / "tools" / "run_wgdi_smoke.sh"]
    for path in paths:
        content = path.read_bytes()
        assert b"\r\n" not in content, f"{path.name} must stay LF-only"

    bash = _homology_bash()
    if bash is None:
        return
    for path in paths:
        result = subprocess.run(
            [str(bash), "-n", str(path)],
            text=True,
            capture_output=True,
            check=False,
        )
        assert result.returncode == 0, f"{path.name}: {result.stderr}"


def test_windows_controller_locates_skill_tools_and_ks_dependencies() -> None:
    controller = (WORKFLOW_ROOT / "run_all_pipeline.ps1").read_text(encoding="utf-8")

    assert 'Join-Path $PSScriptRoot "code\\src"' in controller
    assert 'Join-Path $PSScriptRoot "bin"' in controller
    assert 'if ($WithKs) {' in controller


def test_wgdi_launcher_does_not_hardcode_another_linux_environment() -> None:
    launcher = (WORKFLOW_ROOT / "bin" / "wgdi").read_text(encoding="utf-8")

    assert "WGDI_PYTHON" in launcher
    assert "CONDA_PREFIX" in launcher
    assert "$HOME/miniconda" not in launcher
    assert "/home/" not in launcher
    assert "miniconda3/envs/*" not in launcher


def test_generated_config_is_machine_independent(tmp_path: Path) -> None:
    project = _workspace(tmp_path)
    # Cross-species mode: no self-comparison folder may be created.
    result = _generate(project, "different")
    assert result.returncode == 0, result.stdout + result.stderr
    config = (
        project / "wgdi_input" / "config" / "A_B.conf"
    ).read_text(encoding="utf-8")
    assert not re.search(r"[A-Za-z]:[\\/]", config)
    assert "blast = wgdi_input/blast/A_pep/A_B.blast" in config
    assert "gff1 = wgdi_input/A/A.new.gff" in config
    assert "gff2 = wgdi_input/B/B.new.gff" in config
    assert "lens1 = wgdi_input/A/A.lens" in config
    assert "lens2 = wgdi_input/B/B.lens" in config
    assert "savefile = wgdi_output/results/A_B.icl" in config
    # Lean folder: only the config -- no bundle copies, no per-comparison run.sh.
    assert {entry.name for entry in (project / "wgdi_output").iterdir()} == {
        "results",
        "logs",
        "state",
    }
    assert not (project / "wgdi_input" / "config" / "A.conf").exists()


def test_generator_follows_species_map_for_non_bundle_layout(
    tmp_path: Path,
) -> None:
    # Smoke workspaces keep four-piece sets at the project root; the scan map
    # drives the relative references, so both layouts work unchanged.
    project = tmp_path / "smoke"
    project.mkdir()
    meta = project / "metadata"
    meta.mkdir()
    (meta / "folder_list.txt").write_text("A\nB\n", encoding="utf-8")
    (meta / "species_map.tsv").write_text(
        "A\tA\tA.pep\tA.cds\tA.new.gff\tA.lens\n"
        "B\tB\tB.pep\tB.cds\tB.new.gff\tB.lens\n",
        encoding="utf-8",
    )
    for code in ("A", "B"):
        species_dir = project / code
        species_dir.mkdir()
        for suffix in ("pep", "cds", "new.gff", "lens"):
            (species_dir / f"{code}.{suffix}").write_text("x\n", encoding="utf-8")
        pep_dir = project / "wgdi_input" / "blast" / f"{code}_pep"
        pep_dir.mkdir(parents=True)
        for other in ("A", "B"):
            (pep_dir / f"{code}_{other}.blast").write_text("q\ts\n", encoding="utf-8")
    result = _generate(project, "all")
    assert result.returncode == 0, result.stdout + result.stderr
    config = (
        project / "wgdi_input" / "config" / "A_B.conf"
    ).read_text(encoding="utf-8")
    assert "gff1 = A/A.new.gff" in config
    assert "gff2 = B/B.new.gff" in config
    assert "blast = wgdi_input/blast/A_pep/A_B.blast" in config


def test_generator_fails_fast_on_missing_bundle(tmp_path: Path) -> None:
    project = _workspace(tmp_path)
    (project / "wgdi_input" / "B" / "B.lens").unlink()
    result = _generate(project, "all")
    assert result.returncode == 1
    assert "B.lens" in result.stderr


def test_generator_fails_fast_on_missing_blast(tmp_path: Path) -> None:
    project = _workspace(tmp_path)
    (project / "wgdi_input" / "blast" / "A_pep" / "A_B.blast").unlink()
    result = _generate(project, "all")
    assert result.returncode == 1
    assert "A_B.blast" in result.stderr


def test_dotted_species_codes_keep_full_result_prefix(tmp_path: Path) -> None:
    project = tmp_path / "dotted"
    project.mkdir()
    meta = project / "metadata"
    meta.mkdir()
    codes = ("A.v1", "B.v2")
    (meta / "folder_list.txt").write_text("\n".join(codes) + "\n", encoding="utf-8")
    (meta / "species_map.tsv").write_text(
        "".join(
            f"{code}\twgdi_input/{code}\t{code}.pep\t{code}.cds\t{code}.new.gff\t{code}.lens\n"
            for code in codes
        ),
        encoding="utf-8",
    )
    for left in codes:
        bundle = project / "wgdi_input" / left
        bundle.mkdir(parents=True)
        for suffix in ("pep", "cds", "new.gff", "lens"):
            (bundle / f"{left}.{suffix}").write_text("x\n", encoding="utf-8")
        pep_dir = project / "wgdi_input" / "blast" / f"{left}_pep"
        pep_dir.mkdir(parents=True)
        for right in codes:
            (pep_dir / f"{left}_{right}.blast").write_text("q\ts\n", encoding="utf-8")

    result = _generate(project, "different")
    assert result.returncode == 0, result.stdout + result.stderr
    config = (project / "wgdi_input" / "config" / "A.v1_B.v2.conf").read_text(
        encoding="utf-8"
    )
    assert "savefig = wgdi_output/results/A.v1_B.v2.png" in config
    assert "savefile = wgdi_output/results/A.v1_B.v2.icl" in config


def test_homology_skips_not_ready_species_and_infrastructure() -> None:
    text = (WORKFLOW_ROOT / "run_homology.sh").read_text(encoding="utf-8")
    assert "wgdi_ready=no" in text
    assert all(name in text for name in ("analysis_v2", "bin", "clean_v2", "code", "data"))
    assert "HOMOLOGY_ALLOW_NOT_READY" in text
    assert "cleanup_tmp_map" in text
    assert "|| true" in text
    assert "part.${BASHPID:-$$}" in text
    assert "pgci_preprocess readiness" in text
    assert "species_readiness_grep" in text
    runner = (WORKFLOW_ROOT / "3run_wgdi.sh").read_text(encoding="utf-8")
    assert "check-wgdi-outputs" in runner


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
    for code, ready in (("Rdy", "yes"), ("Not", "no"), ("Mar", "marginal"), ("Unsig", "yes")):
        species = project / code
        species.mkdir()
        for suffix in (".pep", ".cds", ".new.gff", ".lens"):
            (species / f"{code}{suffix}").write_text("x\n", encoding="utf-8")
        qc_dir = project / "analysis_v2" / code
        qc_dir.mkdir(parents=True)
        (qc_dir / f"{code}.qc.json").write_text(
            json.dumps(
                {
                    "sequence_selection": {"wgdi_ready": ready},
                    "pipeline_version": "3.6.0",
                    "n_wgdi_sequences": 1,
                }
            ),
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
        b"      prev=\"$a\"\n"
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


def test_homology_adjudication_gate_executes(tmp_path: Path) -> None:
    bash = _homology_bash()
    if bash is None:
        pytest.skip("Git Bash or Linux bash is required to execute run_homology.sh")
    project, fake_diamond = _homology_workspace(tmp_path)
    for code in ("Rdy", "Mar"):
        write_adjudication(
            project, code, confirmed_n=1, decision="accept", reviewer="reviewer-A"
        )
    env = {
        **os.environ,
        "DIAMOND": str(fake_diamond),
        "HOMOLOGY_ENGINE": "diamond",
        "HOMOLOGY_REQUIRE_ADJUDICATION": "true",
        "PYTHON": sys.executable,
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
    mapped = (project / "resources" / "species_map.tsv").read_text(encoding="utf-8")
    assert "Rdy\t" in mapped and "Mar\t" in mapped
    # ready but un-adjudicated must be blocked by the hard gate
    assert "Unsig\t" not in mapped
    assert "Not\t" not in mapped and "Unk\t" not in mapped
    skipped = (project / "resources" / "skipped_species.txt").read_text(encoding="utf-8")
    assert "Unsig" in skipped and "裁决" in skipped


def test_different_config_is_generated_once(tmp_path: Path) -> None:
    project = _workspace(tmp_path)
    result = _generate(project, "different")
    assert result.returncode == 0, result.stdout + result.stderr
    config = (
        project / "wgdi_input" / "config" / "A_B.conf"
    ).read_text(encoding="utf-8")
    # SOP default: core collinearity only; the Ks branch is opt-in.
    assert "[dotplot]" in config
    assert "[collinearity]" in config
    assert "[ks]" not in config
    assert "[correspondence]" not in config
    assert "ks_col" not in config
    assert not (project / "wgdi_input" / "config" / "A.conf").exists()


def test_ks_branch_is_opt_in_and_adds_all_six_stages(tmp_path: Path) -> None:
    project = _workspace(tmp_path)
    result = _generate(project, "different", ks="true")
    assert result.returncode == 0, result.stdout + result.stderr
    config = (
        project / "wgdi_input" / "config" / "A_B.conf"
    ).read_text(encoding="utf-8")
    assert config.count("tandem = false") == 3
    assert "ks_col = ks_NG86" in config
    assert "cds_file = wgdi_output/state/A_B/all.cds" in config
    assert "pep_file = wgdi_output/state/A_B/all.pep" in config
    for section in ("[ks]", "[blockinfo]", "[correspondence]", "[blockks]", "[kspeaks]", "[peaksfit]"):
        assert section in config
    result_same = _generate(project, "same", ks="true")
    assert result_same.returncode == 0, result_same.stdout + result_same.stderr
    same_config = (
        project / "wgdi_input" / "config" / "A.conf"
    ).read_text(encoding="utf-8")
    assert "cds_file = wgdi_input/A/A.cds" in same_config
    assert "pep_file = wgdi_input/A/A.pep" in same_config


def test_invalid_ks_flag_is_rejected(tmp_path: Path) -> None:
    project = _workspace(tmp_path)
    result = _generate(project, "same", ks="maybe")
    assert result.returncode == 1
    assert "WGDI_KS" in result.stderr


def test_same_config_uses_single_species_and_is_generated_once(tmp_path: Path) -> None:
    project = _workspace(tmp_path)
    result = _generate(project, "same")
    assert result.returncode == 0, result.stdout + result.stderr
    config = (
        project / "wgdi_input" / "config" / "B.conf"
    ).read_text(encoding="utf-8")
    assert "[ks]" not in config
    assert "gff1 = wgdi_input/B/B.new.gff" in config
    assert "gff2 = wgdi_input/B/B.new.gff" in config
    assert "blast = wgdi_input/blast/B_pep/B_B.blast" in config
    assert "savefile = wgdi_output/results/B.icl" in config
    assert not (project / "wgdi_input" / "config" / "A_B.conf").exists()


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
    assert len(shell.splitlines()) < 60
    assert len(powershell.splitlines()) < 60
    # restricted-network / offline / scan-only paths must reach the shared impl
    assert "--index-url" in shell and "--index-url" in powershell
    assert "--offline" in shell and "--offline" in powershell
    assert "--check-only" in shell and "--check-only" in powershell
    assert "--core-only" in shell and "--core-only" in powershell


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
