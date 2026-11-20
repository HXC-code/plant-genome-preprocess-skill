from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from pgci_preprocess.doctor import run_doctor
from pgci_preprocess.pipeline import discover_species
from pgci_preprocess.project_init import ProjectInitError, init_project

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if not SRC_DIR.is_dir():
    SRC_DIR = PROJECT_ROOT / "scripts"


def _sources(tmp_path: Path) -> dict[str, Path]:
    source = tmp_path / "converted"
    source.mkdir(exist_ok=True)
    values = {
        "genome": ("GCF_000000001.1_genomic.fna", ">chr1\nATGC\n"),
        "gff": ("GCF_000000001.1_genomic.gff3", "##gff-version 3\n"),
        "cds": ("GCF_000000001.1_cds.fa", ">gene1\nATG\n"),
        "pep": ("GCF_000000001.1_protein.faa", ">gene1\nM\n"),
    }
    result: dict[str, Path] = {}
    for role, (name, content) in values.items():
        path = source / name
        path.write_text(content, encoding="utf-8")
        result[role] = path
    return result


def _register(tmp_path: Path):
    return init_project(
        tmp_path / "project",
        show_id="42",
        scientific_name="Example species",
        order="Example order",
        family="Example family",
        assembly_id="GCF_000000001.1",
        sources=_sources(tmp_path),
    )


def test_init_empty_scaffold_and_doctor(tmp_path: Path) -> None:
    root = tmp_path / "new" / "nested" / "project"
    result = init_project(root)
    assert (root / "data" / "downloads").is_dir()
    assert result.inventory_path.is_file()
    assert result.assembly_map_path.is_file()
    report = run_doctor(root, allow_empty=True)
    assert report["status"] in {"PASS", "WARN"}
    assert report["species_count"] == 0
    assert report["counts"]["FAIL"] == 0


def test_register_species_and_duplicate_is_fail_closed(tmp_path: Path) -> None:
    result = _register(tmp_path)
    species = discover_species(result.project_root)
    assert len(species) == 1
    assert species[0].show_id == "42"
    assert species[0].assembly_id == "GCF_000000001.1"
    assert (result.species_directory / result.files["pep"]).is_file()
    manifest = result.project_root / "resources" / "project_init" / f"{result.species_code}.json"
    assert json.loads(manifest.read_text(encoding="utf-8"))["show_id"] == "42"
    before = result.inventory_path.read_bytes()
    with pytest.raises(ProjectInitError, match="already registered"):
        init_project(
            result.project_root,
            show_id="42",
            scientific_name="Example species",
            order="Example order",
            family="Example family",
            assembly_id="GCF_000000001.1",
            sources=_sources(tmp_path),
        )
    assert result.inventory_path.read_bytes() == before


def test_registration_uses_species_abbreviation_table(tmp_path: Path) -> None:
    result = init_project(
        tmp_path / "mapped-project",
        show_id="2720",
        scientific_name="Lens culinaris",
        order="Fabales",
        family="Fabaceae",
        assembly_id="GWHEUBR00000000.1",
        sources=_sources(tmp_path),
        species_names_path=PROJECT_ROOT / "resources" / "species_abbreviations.tsv",
    )
    assert result.species_code == "Lcu"
    assert (result.project_root / "resources" / "species_abbreviations.tsv").is_file()
    assert (result.project_root / "resources" / "project_init" / "Lcu.json").is_file()
    roots = discover_species(
        result.project_root,
        species_names_path=PROJECT_ROOT / "resources" / "species_abbreviations.tsv",
    )
    assert roots[0].code == "Lcu"
    assert roots[0].show_id == "2720"


def test_init_can_copy_abbreviation_table_without_registering_species(
    tmp_path: Path,
) -> None:
    root = tmp_path / "mapping-only"
    result = init_project(
        root,
        species_names_path=PROJECT_ROOT / "resources" / "species_abbreviations.tsv",
    )
    assert result.show_id == ""
    assert (root / "resources" / "species_abbreviations.tsv").is_file()


def test_doctor_detects_missing_registered_input(tmp_path: Path) -> None:
    result = _register(tmp_path)
    (result.species_directory / result.files["pep"]).unlink()
    report = run_doctor(result.project_root)
    assert report["status"] == "FAIL"
    failed = {item["name"] for item in report["checks"] if item["status"] == "FAIL"}
    assert "input_files" in failed


def test_full_doctor_accepts_explicit_executable_paths(tmp_path: Path) -> None:
    result = _register(tmp_path)
    tools = tmp_path / "tools"
    tools.mkdir()
    overrides: dict[str, str] = {}
    for name in ("diamond", "wgdi", "bash", "muscle", "yn00", "perl"):
        path = tools / name
        path.write_text("test placeholder\n", encoding="utf-8")
        overrides[name] = str(path)
    report = run_doctor(result.project_root, workflow="full", executables=overrides)
    assert report["status"] in {"PASS", "WARN"}
    assert report["counts"]["FAIL"] == 0
    assert report["diamond_settings"] == {"evalue": "1e-5", "sensitivity": "very-sensitive"}


def test_init_and_doctor_cli_blackbox(tmp_path: Path) -> None:
    root = tmp_path / "cli-project"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SRC_DIR)
    initialized = subprocess.run(
        [sys.executable, "-m", "pgci_preprocess", "init-project", "--project-root", str(root), "--json"],
        cwd=PROJECT_ROOT,
        env=env,
        text=True,
        capture_output=True,
        timeout=120,
    )
    assert initialized.returncode == 0, initialized.stderr
    assert json.loads(initialized.stdout)["status"] == "OK"
    doctor = subprocess.run(
        [
            sys.executable,
            "-m",
            "pgci_preprocess",
            "doctor",
            "--project-root",
            str(root),
            "--allow-empty",
            "--json",
            "--output",
            str(root / "doctor.json"),
        ],
        cwd=PROJECT_ROOT,
        env=env,
        text=True,
        capture_output=True,
        timeout=120,
    )
    assert doctor.returncode == 0, doctor.stderr
    assert json.loads(doctor.stdout)["counts"]["FAIL"] == 0
    assert json.loads((root / "doctor.json").read_text(encoding="utf-8"))["status"] in {"PASS", "WARN"}
