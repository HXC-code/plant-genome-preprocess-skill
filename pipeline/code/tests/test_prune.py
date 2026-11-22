from __future__ import annotations

from pathlib import Path

import pytest

from pgci_preprocess.pipeline import PipelineOptions, run_project
from pgci_preprocess.prune import prune_project
from pgci_preprocess.reports import read_tsv
from pgci_preprocess.validate import ValidationOptions, validate_project

from tests.test_pipeline import _make_project


def _project_with_snapshot(tmp_path: Path) -> tuple[Path, str]:
    project, _ = _make_project(tmp_path)
    results = run_project(
        project,
        options=PipelineOptions(
            project_root=project,
            min_genes_per_sequence=0,
            write_all_sequences=True,
        ),
    )
    code = str(results[0]["code"])
    analysis = project / "wgdi_input" / code
    assert (analysis / f"{code}.all.pep.fa").is_file()
    return project, code


def test_prune_drops_snapshots_and_keeps_the_project_valid(tmp_path: Path) -> None:
    project, code = _project_with_snapshot(tmp_path)
    analysis = project / "wgdi_input" / code

    result = prune_project(project, all_snapshots=True)

    assert result["removed_files"] == 4
    assert result["removed_bytes"] > 0
    assert result["index_rows_updated"] == 1
    for suffix in (".all.pep.fa", ".all.cds.fa", ".all.gff", ".all.lens"):
        assert not (analysis / f"{code}{suffix}").exists(), suffix
    # the four-piece set is untouched
    for suffix in (".pep", ".cds", ".new.gff", ".lens"):
        assert (analysis / f"{code}{suffix}").is_file(), suffix

    # index no longer points at the deleted files
    row = next(
        item
        for item in read_tsv(project / "resources" / "wgdi_v2_index.tsv")
        if item["code"] == code
    )
    assert row["all_pep"] == "" and row["all_lens"] == ""

    # no hash manifest is produced at all
    assert not (project / "resources" / "wgdi_v2_manifest_sha256.tsv").exists()

    # and the species still validates
    validation = validate_project(
        ValidationOptions(project_root=project, require_compat_links=False)
    )
    assert all(row["status"] == "PASS" for row in validation.rows), validation.rows


def test_prune_dry_run_deletes_nothing(tmp_path: Path) -> None:
    project, code = _project_with_snapshot(tmp_path)
    analysis = project / "wgdi_input" / code

    result = prune_project(project, all_snapshots=True, dry_run=True)

    assert result["dry_run"] is True
    assert result["removed_files"] == 4
    for suffix in (".all.pep.fa", ".all.cds.fa", ".all.gff", ".all.lens"):
        assert (analysis / f"{code}{suffix}").is_file(), suffix
    row = next(
        item
        for item in read_tsv(project / "resources" / "wgdi_v2_index.tsv")
        if item["code"] == code
    )
    assert row["all_pep"].endswith(".all.pep.fa")


def test_prune_compat_links_removes_root_directories(tmp_path: Path) -> None:
    project, code = _project_with_snapshot(tmp_path)
    link_dir = project / code
    link_dir.mkdir()
    for suffix in (".pep", ".cds", ".new.gff", ".lens"):
        (link_dir / f"{code}{suffix}").write_text("x\n", encoding="utf-8")

    result = prune_project(project, compat_links=True)

    assert result["removed_dirs"] == 1
    assert not link_dir.exists()
    assert "hard link" in str(result["note"])


def test_prune_requires_an_index(tmp_path: Path) -> None:
    empty = tmp_path / "empty"
    (empty / "resources").mkdir(parents=True)
    with pytest.raises(FileNotFoundError):
        prune_project(empty, all_snapshots=True)
