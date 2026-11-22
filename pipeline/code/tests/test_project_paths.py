from __future__ import annotations

import shutil
from pathlib import Path

from pgci_preprocess.project_paths import (
    downloads_root,
    is_species_directory,
    relative_to_data,
    species_by_show_id,
    species_directories,
)


def _species(root: Path, name: str) -> Path:
    path = root / name
    path.mkdir(parents=True)
    (path / "annotation.gff3").write_text("x\n", encoding="utf-8")
    return path


def test_raw_data_sits_directly_under_data(tmp_path: Path) -> None:
    data = tmp_path / "data"
    one = _species(data, "1_Testales_Testaceae(Testus primus)")
    two = _species(data, "2_Testales_Testaceae(Testus secundus)")

    assert downloads_root(tmp_path) == data
    assert species_directories(tmp_path) == sorted([one, two], key=lambda p: p.name)
    assert set(species_by_show_id(tmp_path)) == {"1", "2"}


def test_legacy_data_downloads_still_resolves(tmp_path: Path) -> None:
    legacy = tmp_path / "data" / "downloads"
    species = _species(legacy, "7_Testales_Testaceae(Testus septimus)")

    assert downloads_root(tmp_path) == legacy
    assert species_by_show_id(tmp_path) == {"7": species}


def test_empty_legacy_directory_does_not_shadow_data(tmp_path: Path) -> None:
    # An empty data/downloads (created by an older scaffold step) must not
    # capture resolution once raw species folders live in data/.
    (tmp_path / "data" / "downloads").mkdir(parents=True)
    species = _species(tmp_path / "data", "3_Testales_Testaceae(Testus tertius)")

    assert downloads_root(tmp_path) == tmp_path / "data"
    assert species_by_show_id(tmp_path) == {"3": species}


def test_new_layout_wins_when_both_are_populated(tmp_path: Path) -> None:
    # A half-migrated project resolves to data/, the canonical layout.
    new = _species(tmp_path / "data", "8_Testales_Testaceae(Testus octavus)")
    _species(tmp_path / "data" / "downloads", "9_Testales_Testaceae(Testus nonus)")

    assert downloads_root(tmp_path) == tmp_path / "data"
    assert set(species_by_show_id(tmp_path)) == {"8"}
    assert relative_to_data(new / "annotation.gff3", tmp_path) == (
        "8_Testales_Testaceae(Testus octavus)/annotation.gff3"
    )


def test_infrastructure_and_non_species_dirs_are_ignored(tmp_path: Path) -> None:
    data = tmp_path / "data"
    _species(data, "4_Testales_Testaceae(Testus quartus)")
    (data / "derived").mkdir()
    (data / "downloads").mkdir()
    (data / "notes").mkdir()  # no underscore
    (data / "gap_Assembly").mkdir()  # underscore but non-numeric prefix

    assert [p.name for p in species_directories(tmp_path)] == [
        "4_Testales_Testaceae(Testus quartus)"
    ]
    assert not is_species_directory(data / "gap_Assembly")


def test_end_to_end_run_from_flat_data(tmp_path: Path) -> None:
    """Moving data/downloads up one level is the whole migration."""

    from pgci_preprocess.pipeline import PipelineOptions, run_project
    from tests.support import build_multiformat_project

    fixture = build_multiformat_project(tmp_path)["project"]
    flat = tmp_path / "flat"
    flat.mkdir()
    shutil.move(str(fixture / "data"), str(flat / "data"))
    shutil.move(str(fixture / "resources"), str(flat / "resources"))
    # flatten: data/downloads/<species> -> data/<species>
    for entry in sorted((flat / "data" / "downloads").iterdir()):
        shutil.move(str(entry), str(flat / "data" / entry.name))
    (flat / "data" / "downloads").rmdir()

    results = run_project(
        flat,
        options=PipelineOptions(
            project_root=flat,
            min_genes_per_sequence=1,
            min_genes_per_chromosome=0,
            sequence_mode="adaptive",
        ),
    )
    assert {str(r["code"]) for r in results} == {"s1", "s2", "s3", "s4"}
    assert all(r["status"] == "OK" for r in results), results
    assert downloads_root(flat) == flat / "data"
    for code in ("s1", "s2", "s3", "s4"):
        assert (flat / "wgdi_input" / code / f"{code}.pep").is_file()


def test_relative_to_data_handles_both_layouts(tmp_path: Path) -> None:
    new = _species(tmp_path / "data", "5_Testales_Testaceae(Testus quintus)")
    assert relative_to_data(new / "annotation.gff3", tmp_path) == (
        "5_Testales_Testaceae(Testus quintus)/annotation.gff3"
    )
    old = _species(tmp_path / "data" / "downloads", "6_Testales_Testaceae(Testus sextus)")
    assert relative_to_data(old / "annotation.gff3", tmp_path) == (
        "6_Testales_Testaceae(Testus sextus)/annotation.gff3"
    )
