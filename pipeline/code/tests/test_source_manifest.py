"""Source manifest: verified provenance, evidence-only repairs, fail-closed."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from pgci_preprocess.source_manifest import (
    SourceManifestError,
    build_source_manifest,
    load_repairs,
    verify_source_manifest,
)

STANDARD_GFF = (
    "chr1\tsrc\tregion\t1\t1000\t.\t+\t.\tID=chr1:1..1000\n"
    "chr1\tsrc\tgene\t10\t90\t.\t+\t.\tID=gene-A\n"
    "chr1\tsrc\tmRNA\t10\t90\t.\t+\t.\tID=rna-A1;Parent=gene-A\n"
    "chr1\tsrc\tCDS\t10\t90\t.\t+\t0\tID=cds-P1;Parent=rna-A1;protein_id=P1\n"
)
TRANSCRIPT_ONLY_GFF = (
    "chr1\tEVM\tmRNA\t10\t90\t.\t+\t.\tID=t1\n"
    "chr1\tEVM\tCDS\t10\t90\t.\t+\t0\tParent=t1\n"
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(path: Path) -> dict[str, object]:
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": _sha(path)}


def _make_project(tmp_path: Path, *, species: int = 1) -> tuple[Path, Path, dict[str, Path]]:
    project = tmp_path / "project"
    raw = tmp_path / "raw"
    folder = raw / f"{species}_Testales_Testaceae(Testus plantus)"
    folder.mkdir(parents=True)
    (project / "resources").mkdir(parents=True)
    (project / "analysis_v2" / f"s{species}").mkdir(parents=True)
    files = {
        "gff": folder / "test.gff",
        "cds": folder / "test.cds.fa",
        "pep": folder / "test.pep.fa",
        "genome": folder / "test.genome.fa",
    }
    files["gff"].write_text(STANDARD_GFF, encoding="utf-8")
    files["cds"].write_text(">P1\nATGAAATAA\n", encoding="utf-8")
    files["pep"].write_text(">P1\nMK\n", encoding="utf-8")
    files["genome"].write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    qc = {
        "order": "Testales",
        "family": "Testaceae",
        "inputs": {role: _record(path) for role, path in files.items()},
    }
    qc_path = project / "analysis_v2" / f"s{species}" / f"s{species}.qc.json"
    qc_path.write_text(json.dumps(qc), encoding="utf-8")
    index = project / "resources" / "wgdi_v2_index.tsv"
    index.write_text(
        "code\tspecies\taccession\twgdi_ready\tqc\n"
        f"s{species}\tTestus plantus\tGWHTEST00000000.1\tyes\tanalysis_v2/s{species}/s{species}.qc.json\n",
        encoding="utf-8",
    )
    return project, raw, files


def _write_repairs(path: Path, rows: list[dict[str, str]]) -> Path:
    fields = ["code", "role", "action", "new_path", "expected_sha256", "evidence", "note"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_build_verifies_every_role_and_writes_sidecar(tmp_path: Path) -> None:
    project, raw, files = _make_project(tmp_path)
    output = tmp_path / "manifest.tsv"
    rows, summary = build_source_manifest(
        project, raw_roots=(raw,), output_path=output
    )
    assert summary["blockers"] == []
    assert len(rows) == 4
    assert all(row["status"] == "VERIFIED" for row in rows)
    assert all(row["integrity"] == "sha256_verified" for row in rows)
    gff_row = next(row for row in rows if row["role"] == "gff")
    assert gff_row["content_format"] == "standard"
    assert gff_row["conversion_required"] == "none"
    sidecar = output.with_suffix(".tsv.sha256")
    assert sidecar.is_file()
    assert sidecar.read_text(encoding="utf-8").split()[0] == summary["manifest_sha256"]
    report = verify_source_manifest(output)
    assert report["ok"], report["problems"]
    assert report["files_rehashed"] == 4


def test_build_is_deterministic(tmp_path: Path) -> None:
    project, raw, _files = _make_project(tmp_path)
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    build_source_manifest(project, raw_roots=(raw,), output_path=first)
    build_source_manifest(project, raw_roots=(raw,), output_path=second)
    assert first.read_bytes() == second.read_bytes()


def test_dead_path_remap_requires_hash_match(tmp_path: Path) -> None:
    project, raw, files = _make_project(tmp_path)
    qc_path = project / "analysis_v2" / "s1" / "s1.qc.json"
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    qc["inputs"]["gff"]["path"] = str(tmp_path / "gone" / "test.gff")
    qc_path.write_text(json.dumps(qc), encoding="utf-8")
    rows, summary = build_source_manifest(project, raw_roots=(raw,))
    gff_row = next(row for row in rows if row["role"] == "gff")
    assert gff_row["status"] == "REMAPPED_HASH_VERIFIED"
    assert Path(gff_row["source_path"]) == files["gff"]
    assert summary["blockers"] == []


def test_dead_path_without_hash_match_is_blocker(tmp_path: Path) -> None:
    project, raw, files = _make_project(tmp_path)
    qc_path = project / "analysis_v2" / "s1" / "s1.qc.json"
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    qc["inputs"]["gff"]["path"] = str(tmp_path / "gone" / "test.gff")
    qc_path.write_text(json.dumps(qc), encoding="utf-8")
    files["gff"].write_text(STANDARD_GFF + "chr1\tsrc\tgene\t95\t99\t.\t+\t.\tID=gene-B\n", encoding="utf-8")
    rows, summary = build_source_manifest(project, raw_roots=(raw,))
    gff_row = next(row for row in rows if row["role"] == "gff")
    assert gff_row["status"] == "BLOCKER_MISSING"
    assert summary["blockers"] == ["s1/gff:BLOCKER_MISSING"]
    assert all(row["status"] == "VERIFIED" for row in rows if row["role"] != "gff")


def test_missing_role_blocks_only_that_species(tmp_path: Path) -> None:
    project, raw, files = _make_project(tmp_path)
    files["pep"].unlink()
    rows, summary = build_source_manifest(project, raw_roots=(raw,))
    pep_row = next(row for row in rows if row["role"] == "pep")
    assert pep_row["status"] == "BLOCKER_MISSING"
    assert summary["status_counts"]["VERIFIED"] == 3
    assert summary["blockers"] == ["s1/pep:BLOCKER_MISSING"]


def test_genome_absent_is_recorded_as_omitted_optional(tmp_path: Path) -> None:
    project, raw, _files = _make_project(tmp_path)
    qc_path = project / "analysis_v2" / "s1" / "s1.qc.json"
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    del qc["inputs"]["genome"]
    qc_path.write_text(json.dumps(qc), encoding="utf-8")
    rows, summary = build_source_manifest(project, raw_roots=(raw,))
    genome_row = next(row for row in rows if row["role"] == "genome")
    assert genome_row["status"] == "OMITTED_OPTIONAL"
    assert "optional" in genome_row["provenance_note"]
    assert summary["blockers"] == []


def test_repair_replaces_artifact_with_hash_verified_original(tmp_path: Path) -> None:
    project, raw, files = _make_project(tmp_path)
    original = raw / "1_Testales_Testaceae(Testus plantus)" / "original.gff"
    original.write_text(TRANSCRIPT_ONLY_GFF, encoding="utf-8")
    repairs = _write_repairs(tmp_path / "repairs.tsv", [{
        "code": "s1", "role": "gff", "action": "replace_with_original",
        "new_path": str(original), "expected_sha256": _sha(original),
        "evidence": "task1_artifact_tracing", "note": "normalized artifact rebuilt from original",
    }])
    rows, summary = build_source_manifest(
        project, repairs_path=repairs, raw_roots=(raw,)
    )
    gff_row = next(row for row in rows if row["role"] == "gff")
    assert gff_row["status"] == "REPAIRED_ORIGINAL"
    assert gff_row["content_format"] == "transcript_only"
    assert gff_row["conversion_required"] == "transcript_only"
    assert gff_row["converter"] == "pgci_preprocess.gff_formats"
    assert "qc_recorded=" in gff_row["conversion_chain"]
    assert summary["blockers"] == []


def test_repair_with_wrong_hash_is_blocker(tmp_path: Path) -> None:
    project, raw, _files = _make_project(tmp_path)
    original = raw / "1_Testales_Testaceae(Testus plantus)" / "original.gff"
    original.write_text(TRANSCRIPT_ONLY_GFF, encoding="utf-8")
    repairs = _write_repairs(tmp_path / "repairs.tsv", [{
        "code": "s1", "role": "gff", "action": "replace_with_original",
        "new_path": str(original), "expected_sha256": "0" * 64,
        "evidence": "test", "note": "",
    }])
    rows, summary = build_source_manifest(project, repairs_path=repairs, raw_roots=(raw,))
    gff_row = next(row for row in rows if row["role"] == "gff")
    assert gff_row["status"] == "BLOCKER_HASH_MISMATCH"
    assert summary["blockers"]


def test_verify_detects_source_tampering(tmp_path: Path) -> None:
    project, raw, files = _make_project(tmp_path)
    output = tmp_path / "manifest.tsv"
    build_source_manifest(project, raw_roots=(raw,), output_path=output)
    with files["cds"].open("a", encoding="utf-8") as handle:
        handle.write(">extra\nATG\n")
    report = verify_source_manifest(output)
    assert not report["ok"]
    assert any("s1/cds" in problem for problem in report["problems"])


def test_verify_detects_manifest_edit(tmp_path: Path) -> None:
    project, raw, _files = _make_project(tmp_path)
    output = tmp_path / "manifest.tsv"
    build_source_manifest(project, raw_roots=(raw,), output_path=output)
    text = output.read_text(encoding="utf-8").replace("VERIFIED", "FORGED")
    output.write_text(text, encoding="utf-8")
    report = verify_source_manifest(output, rehash=False)
    assert not report["ok"]
    assert any("manifest_content_hash_mismatch" in problem for problem in report["problems"])


def test_repair_accepts_gbff_original_for_all_extracted_roles(tmp_path: Path) -> None:
    """s486 pattern: one original GBFF supplies gff/cds/pep via the adapter."""

    project, raw, files = _make_project(tmp_path)
    gbff = raw / "1_Testales_Testaceae(Testus plantus)" / "original.gbff.gz"
    import gzip

    record = (
        "LOCUS       TEST000001              90 bp    DNA     linear   PLN 19-NOV-2026\n"
        "FEATURES             Location/Qualifiers\n"
        "     source          1..90\n"
        '                     /organism="Testus plantus"\n'
        "     CDS             1..90\n"
        '                     /protein_id="P1"\n'
        '                     /translation="MKMKMKMKMKMKMKMKMKMKMKMKMKMKM"\n'
        "ORIGIN\n"
        "        1 atgaaaatga aaatgaaaat gaaaatgaaa atgaaaatga aaatgaaaat gaaaatgaaa\n"
        "       61 atgaaaatga aaatgaaaat gaaaatgaaa\n"
        "//\n"
    )
    with gzip.open(gbff, "wt", encoding="utf-8") as handle:
        handle.write(record)
    repairs = _write_repairs(tmp_path / "repairs.tsv", [
        {
            "code": "s1", "role": role, "action": "replace_with_original",
            "new_path": str(gbff), "expected_sha256": _sha(gbff),
            "evidence": "task3 gbff chain", "note": "derived artifact replaced by original GBFF",
        }
        for role in ("gff", "cds", "pep")
    ])
    rows, summary = build_source_manifest(project, repairs_path=repairs, raw_roots=(raw,))
    assert summary["blockers"] == []
    for role in ("gff", "cds", "pep"):
        entry = next(row for row in rows if row["role"] == role)
        assert entry["status"] == "REPAIRED_ORIGINAL"
        assert entry["content_format"] == "GBFF"
        assert entry["conversion_required"] == "gbff"
        assert entry["converter"] == "pgci_preprocess.annotation_convert"


def test_invalid_repairs_are_rejected(tmp_path: Path) -> None:
    bad_action = _write_repairs(tmp_path / "bad1.tsv", [{
        "code": "s1", "role": "gff", "action": "guess",
        "new_path": "x", "expected_sha256": "0" * 64, "evidence": "e", "note": "",
    }])
    with pytest.raises(SourceManifestError):
        load_repairs(bad_action)
    no_evidence = _write_repairs(tmp_path / "bad2.tsv", [{
        "code": "s1", "role": "gff", "action": "remap",
        "new_path": "x", "expected_sha256": "0" * 64, "evidence": "", "note": "",
    }])
    with pytest.raises(SourceManifestError):
        load_repairs(no_evidence)
