from __future__ import annotations

from pathlib import Path

from pgci_preprocess.audit import audit_project
from pgci_preprocess.pipeline import PipelineOptions, load_format_exceptions, run_project
from pgci_preprocess.reports import merge_index, merge_registry, merge_summary, read_tsv
from tests.support import build_multiformat_project


def test_load_format_exceptions_reads_evidence_table(tmp_path: Path) -> None:
    project = tmp_path / "project"
    resources = project / "resources"
    resources.mkdir(parents=True)
    (resources / "format_exceptions.tsv").write_text(
        "show_id\tmerge_split_genes\tgene_span_only\tgene_cds_direct\treason\n"
        "4200\tN\tY\tN\tgene-only\n"
        "5630\tY\tN\tY\tsegmented genes\n",
        encoding="utf-8",
    )
    mapping = load_format_exceptions(project)
    assert mapping["4200"]["gene_span_only"] is True
    assert mapping["5630"]["merge_split_genes"] is True
    assert mapping["5630"]["gene_cds_direct"] is True
    assert mapping["4200"]["merge_split_genes"] is False


def test_merge_helpers_preserve_unselected_species(tmp_path: Path) -> None:
    project = tmp_path / "project"
    resources = project / "resources"
    resources.mkdir(parents=True)
    (resources / "stage_preprocess_v2_summary.tsv").write_text(
        "show_id\tcode\tspecies\taccession\tstatus\twgdi_ready\n"
        "999\ts999\tOldus testus\tGCF_1\tOK\tyes\n",
        encoding="utf-8",
    )
    (resources / "species_registry_v2.tsv").write_text(
        "show_id\tcode\tspecies\torder\tfamily\taccession\tgff\tcds\tpep\tgenome\tdirectory\n"
        "999\ts999\tOldus testus\tTestales\tTestaceae\tGCF_1\ta.gff\ta.cds\ta.pep\t\td999\n",
        encoding="utf-8",
    )
    (resources / "wgdi_v2_index.tsv").write_text(
        "code\tspecies\taccession\twgdi_ready\tpep\tcds\tgff\tlens\tall_pep\tall_cds\tall_gff\tall_lens\tgene_map\tchr_map\trejected\tqc\n"
        "s999\tOldus testus\tGCF_1\tyes\ta.pep\ta.cds\ta.gff\ta.lens\t\t\t\t\ta.gmap\tachr\ta.rej\ta.qc\n",
        encoding="utf-8",
    )
    new_result = {
        "show_id": "1",
        "code": "s1",
        "species": "Testus primus",
        "accession": "GCF_2",
        "status": "OK",
        "wgdi_ready": "yes",
    }
    merge_summary(project, [new_result], replace_codes={"s1"})
    merge_registry(
        project,
        [
            {
                "show_id": "1",
                "code": "s1",
                "species": "Testus primus",
                "order": "Testales",
                "family": "Testaceae",
                "accession": "GCF_2",
                "gff": "b.gff",
                "cds": "b.cds",
                "pep": "b.pep",
                "genome": "",
                "directory": "d1",
            }
        ],
        replace_codes={"s1"},
    )
    merge_index(project, [{**new_result, "qc": "analysis_v2/s1/s1.qc.json"}], replace_codes={"s1"})
    assert {row["code"] for row in read_tsv(resources / "stage_preprocess_v2_summary.tsv")} == {
        "s1",
        "s999",
    }
    assert {row["code"] for row in read_tsv(resources / "species_registry_v2.tsv")} == {
        "s1",
        "s999",
    }
    assert {row["code"] for row in read_tsv(resources / "wgdi_v2_index.tsv")} == {
        "s1",
        "s999",
    }


def test_audit_is_read_only(tmp_path: Path) -> None:
    project = build_multiformat_project(tmp_path)["project"]
    results = run_project(
        project,
        options=PipelineOptions(project_root=project, min_genes_per_sequence=1),
    )
    assert all(result["status"] == "OK" for result in results)

    before = {
        str(path.relative_to(project)): path.stat().st_mtime_ns
        for path in project.rglob("*")
        if path.is_file()
    }
    report = audit_project(project)
    after = {
        str(path.relative_to(project)): path.stat().st_mtime_ns
        for path in project.rglob("*")
        if path.is_file()
    }
    assert before == after
    assert report.inventory_species == 4
    assert report.index_species == 4
    assert report.summary_species == 4
    assert report.qc_ok == 4
    assert report.qc_bad == 0
    assert report.status_distribution["OK"] == 4
    assert report.ready_distribution["yes"] == 4
    assert report.annotation_levels["transcript_only"] == 1
    assert report.annotation_levels["simplified5"] == 1
    assert report.format_exceptions["transcript_only_gff"] == 1
    assert report.format_exceptions["simplified5_gff"] == 1
