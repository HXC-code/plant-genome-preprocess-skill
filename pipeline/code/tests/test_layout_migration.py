from __future__ import annotations

import json
import os
from pathlib import Path

from pgci_preprocess.layout_migration import migrate_layout
from pgci_preprocess.project_paths import output_dirname, output_root
from pgci_preprocess.validate import sync_compat_links

FIELDS = ("show_id", "code", "pep", "cds", "gff", "lens", "gene_map", "qc")


def _project(tmp_path: Path, *, compat: bool = True) -> Path:
    project = tmp_path / "project"
    bundle = project / "analysis_v2" / "A"
    bundle.mkdir(parents=True)
    resources = project / "resources"
    resources.mkdir()

    for suffix in (".pep", ".cds", ".new.gff", ".lens"):
        (bundle / f"A{suffix}").write_text(f"A{suffix}\n", encoding="utf-8")
    for suffix in (".gene_map.tsv", ".qc.json"):
        (bundle / f"A{suffix}").write_text("{}\n", encoding="utf-8")

    header = "\t".join(FIELDS)
    row = "\t".join(
        [
            "1",
            "A",
            "analysis_v2/A/A.pep",
            "analysis_v2/A/A.cds",
            "analysis_v2/A/A.new.gff",
            "analysis_v2/A/A.lens",
            "analysis_v2/A/A.gene_map.tsv",
            "analysis_v2/A/A.qc.json",
        ]
    )
    (resources / "wgdi_v2_index.tsv").write_text(
        f"{header}\n{row}\n", encoding="utf-8"
    )
    (resources / "machine_paths.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "paths": {"project_root": str(project), "output_dirname": "analysis_v2"},
            }
        ),
        encoding="utf-8",
    )

    if compat:
        link_dir = project / "A"
        link_dir.mkdir()
        for suffix in (".pep", ".cds", ".new.gff", ".lens"):
            try:
                os.link(bundle / f"A{suffix}", link_dir / f"A{suffix}")
            except OSError:
                (link_dir / f"A{suffix}").write_text(
                    (bundle / f"A{suffix}").read_text(encoding="utf-8"),
                    encoding="utf-8",
                )
    return project


def test_migrate_layout_renames_container_and_index(tmp_path: Path) -> None:
    project = _project(tmp_path)

    assert output_dirname(project) == "analysis_v2"
    result = migrate_layout(project)

    assert result.ok and result.status == "RENAMED"
    assert not (project / "analysis_v2").exists()
    assert (project / "wgdi_input" / "A" / "A.pep").is_file()
    assert not (project / "A").exists()

    index = (project / "resources" / "wgdi_v2_index.tsv").read_text(encoding="utf-8")
    assert "wgdi_input/A/A.pep" in index
    assert "analysis_v2" not in index

    config = json.loads((project / "resources" / "machine_paths.json").read_text("utf-8"))
    assert config["paths"]["output_dirname"] == "wgdi_input"
    assert output_dirname(project) == "wgdi_input"


def test_migrate_layout_keeps_directories_with_extra_files(tmp_path: Path) -> None:
    project = _project(tmp_path)
    (project / "A" / "notes.txt").write_text("keep me\n", encoding="utf-8")

    result = migrate_layout(project)

    assert result.ok
    assert result.removed_links == []
    assert result.kept_directories == ["A"]
    assert (project / "A" / "notes.txt").is_file()


def test_migrate_layout_refuses_to_merge_two_containers(tmp_path: Path) -> None:
    project = _project(tmp_path)
    (project / "wgdi_input").mkdir()

    result = migrate_layout(project)

    assert result.status == "BLOCKED"
    assert not result.ok
    assert (project / "analysis_v2" / "A" / "A.pep").is_file()


def test_migrate_layout_dry_run_changes_nothing(tmp_path: Path) -> None:
    project = _project(tmp_path)

    result = migrate_layout(project, dry_run=True)

    assert result.status == "DRY_RUN"
    assert result.removed_links == ["A"]
    assert (project / "analysis_v2" / "A" / "A.pep").is_file()
    assert (project / "A" / "A.pep").is_file()


def test_new_container_needs_no_compatibility_links(tmp_path: Path) -> None:
    project = _project(tmp_path)
    migrate_layout(project)

    assert output_root(project) == project / "wgdi_input"
    assert sync_compat_links(project, project) == []
    assert not (project / "A").exists()
