from __future__ import annotations

import json
from pathlib import Path

from pgci_preprocess.pipeline import PipelineOptions, discover_species, run_project
from pgci_preprocess.reports import read_tsv
from pgci_preprocess.species_names import load_species_names, resolve_show_ids
from pgci_preprocess.species_names_migration import migrate_species_names
from pgci_preprocess.validate import sync_compat_links
from tests.support import build_multiformat_project


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _write_test_mapping(project: Path) -> None:
    path = project / "resources" / "species_abbreviations.tsv"
    path.write_text(
        "show_id\tlatin_name\torder\tfamily\tabbreviation\n"
        "1\tAcorus gramineus\tAcorales\tAcoraceae\tAdi\n"
        "2\tAcorus tatarinowii\tAcorales\tAcoraceae\tAes\n"
        "3\tAlisma plantago-aquatica\tAlismatales\tAraceae\tAlb\n"
        "4\tAmorphophallus albus\tAlismatales\tAraceae\tAct\n",
        encoding="utf-8",
    )


def test_authoritative_mapping_resolves_show_id_and_legacy_aliases() -> None:
    names = load_species_names(PROJECT_ROOT)
    assert names is not None
    assert names.code_for("2720") == "Lcu"
    assert names.code_for("2725") == "Lsb"
    assert names.resolve_show_id("Lcu") == "2720"
    assert names.resolve_show_id("s2720") == "2720"
    assert resolve_show_ids(PROJECT_ROOT, {"Lcu", "2725"}) == {"2720", "2725"}


def test_discovery_outputs_and_index_use_abbreviations(tmp_path: Path) -> None:
    project = build_multiformat_project(tmp_path)["project"]
    _write_test_mapping(project)
    species = {item.show_id: item for item in discover_species(project)}
    assert {key: item.code for key, item in species.items()} == {
        "1": "Adi",
        "2": "Aes",
        "3": "Alb",
        "4": "Act",
    }

    results = run_project(
        project,
        options=PipelineOptions(project_root=project, min_genes_per_sequence=1),
    )
    assert {str(item["code"]) for item in results} == {"Adi", "Aes", "Alb", "Act"}
    assert (project / "analysis_v2" / "Adi" / "Adi.qc.json").is_file()
    assert (project / "analysis_v2" / "Adi" / "Adi.gene_map.tsv").is_file()
    assert not (project / "clean_v2").exists()
    assert not (project / "analysis_v2" / "s1").exists()

    index = read_tsv(project / "resources" / "wgdi_v2_index.tsv")
    assert {row["code"] for row in index} == {"Adi", "Aes", "Alb", "Act"}
    assert {row["show_id"] for row in index} == {"1", "2", "3", "4"}
    assert all("analysis_v2/Adi/" in row["pep" ] for row in index if row["code"] == "Adi")


def test_species_selectors_accept_abbreviation_and_legacy_code(tmp_path: Path) -> None:
    project = build_multiformat_project(tmp_path)["project"]
    _write_test_mapping(project)
    results = run_project(
        project,
        show_ids={"s2", "Alb"},
        options=PipelineOptions(project_root=project, min_genes_per_sequence=1),
    )
    assert {str(item["show_id"]) for item in results} == {"2", "3"}
    assert {str(item["code"]) for item in results} == {"Aes", "Alb"}


def test_legacy_s_outputs_are_reused_but_index_and_compat_links_use_abbreviation(
    tmp_path: Path,
) -> None:
    paths = build_multiformat_project(tmp_path)
    project = paths["project"]
    run_project(
        project,
        options=PipelineOptions(project_root=project, min_genes_per_sequence=1),
    )
    assert (project / "analysis_v2" / "s1" / "s1.qc.json").is_file()

    _write_test_mapping(project)
    results = run_project(
        project,
        options=PipelineOptions(project_root=project, min_genes_per_sequence=1),
    )
    by_show = {str(item["show_id"]): item for item in results}
    assert by_show["1"]["status"] == "SKIPPED"
    assert by_show["1"]["code"] == "Adi"
    # 第二次运行会自动把旧 s<show_id> 产物改名为缩写，artifact_code 同步更新。
    assert by_show["1"]["artifact_code"] == "Adi"

    index = {row["code"]: row for row in read_tsv(project / "resources" / "wgdi_v2_index.tsv")}
    assert set(index) == {"Adi", "Aes", "Alb", "Act"}
    assert index["Adi"]["pep"] == "analysis_v2/Adi/Adi.pep"
    assert index["Adi"]["show_id"] == "1"
    assert not (project / "analysis_v2" / "s1").exists()

    workspace = project.parent
    sync_compat_links(project, workspace)
    assert (workspace / "Adi" / "Adi.pep").is_file()
    assert not (workspace / "Adi" / "s1.pep").exists()


def test_metadata_migration_updates_codes_without_rewriting_artifacts(
    tmp_path: Path,
) -> None:
    project = build_multiformat_project(tmp_path)["project"]
    run_project(
        project,
        options=PipelineOptions(project_root=project, min_genes_per_sequence=1),
    )
    _write_test_mapping(project)

    result = migrate_species_names(project, refresh_links=True)
    assert result.unresolved == []
    assert result.backup_dir and result.backup_dir.is_dir()
    inventory = read_tsv(project / "resources" / "download_inventory.tsv")
    assert {row["species_code"] for row in inventory} == {
        "Adi",
        "Aes",
        "Alb",
        "Act",
    }
    index = read_tsv(project / "resources" / "wgdi_v2_index.tsv")
    assert {row["code"] for row in index} == {"Adi", "Aes", "Alb", "Act"}
    assert {row["artifact_code"] for row in index} == {"s1", "s2", "s3", "s4"}
    assert {row["show_id"] for row in index} == {"1", "2", "3", "4"}
    assert all(row["pep"].startswith("analysis_v2/s") for row in index)
    assert (project.parent / "Adi" / "Adi.pep").is_file()
    qc = json.loads(
        (project / "analysis_v2" / "s1" / "s1.qc.json").read_text(encoding="utf-8")
    )
    assert qc["code"] == "Adi"
    assert qc["artifact_code"] == "s1"
    assert qc["legacy_code"] == "s1"
    second = migrate_species_names(project, refresh_links=False)
    assert second.changed_files["qc_metadata"] == 0
