"""Shared TSV I/O and project registry/report consolidation.

This module is the single implementation of the "finalize" step that used to be
copy-pasted between ``pipeline.run_project`` and the legacy
``convert_simplified_species.py`` script.  It only merges report rows; it never
rewrites per-species scientific outputs.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

from pgci_preprocess.selection import natural_key

REGISTRY_FIELDS = [
    "show_id",
    "code",
    "artifact_code",
    "species",
    "order",
    "family",
    "accession",
    "gff",
    "cds",
    "pep",
    "genome",
    "directory",
]
SUMMARY_FIELDS = [
    "show_id",
    "code",
    "species",
    "accession",
    "status",
    "assembly_build",
    "n_regions",
    "n_representative_genes",
    "n_all_genes",
    "n_wgdi_genes",
    "n_wgdi_sequences",
    "sequence_mode",
    "retain_ratio",
    "wgdi_ready",
    "nuclear_taxa",
    "n_anomalies",
    "rejection_counts",
    "reason",
    "detail",
]
INDEX_FIELDS = [
    "show_id",
    "code",
    "artifact_code",
    "species",
    "accession",
    "wgdi_ready",
    "pep",
    "cds",
    "gff",
    "lens",
    "all_pep",
    "all_cds",
    "all_gff",
    "all_lens",
    "gene_map",
    "chr_map",
    "rejected",
    "qc",
]
MANIFEST_FIELDS = ["code", "layer", "relative_path", "bytes", "sha256"]


def read_tsv(path: Path | str) -> list[dict[str, str]]:
    path = Path(path)
    if not path.is_file():
        return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_tsv_with_fields(path: Path | str) -> tuple[list[str], list[dict[str, str]]]:
    path = Path(path)
    opened = path.open(encoding="utf-8-sig", newline="")
    with opened as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    return fields, rows


def write_tsv(
    path: Path | str,
    fields: Sequence[str],
    rows: Iterable[Mapping[str, Any]],
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def _code(row: Mapping[str, Any]) -> str:
    return str(row.get("code", "")).strip()


def _merge_rows(
    path: Path,
    fields: Sequence[str],
    new_rows: Sequence[Mapping[str, Any]],
    replace_codes: set[str],
) -> list[dict[str, Any]]:
    existing = read_tsv(path)
    kept = [row for row in existing if _code(row) not in replace_codes]
    merged: list[dict[str, Any]] = [
        dict(row) for row in [*kept, *new_rows]
    ]
    merged.sort(key=lambda row: natural_key(_code(row)))
    write_tsv(path, fields, merged)
    return merged


def merge_registry(
    project_root: Path | str,
    species_rows: Sequence[Mapping[str, Any]],
    *,
    replace_codes: set[str] | None = None,
) -> list[dict[str, Any]]:
    """Incrementally refresh ``species_registry_v2.tsv``.

    ``replace_codes`` defaults to every species code in ``species_rows`` so a
    full discovery pass refreshes all known inventory species while preserving
    legacy entries that are not present in the current inventory.
    """

    path = Path(project_root) / "resources" / "species_registry_v2.tsv"
    if replace_codes is None:
        replace_codes = {_code(row) for row in species_rows}
    rows_to_write = [row for row in species_rows if _code(row) in replace_codes]
    return _merge_rows(path, REGISTRY_FIELDS, rows_to_write, replace_codes)


def merge_summary(
    project_root: Path | str,
    results: Sequence[Mapping[str, Any]],
    *,
    replace_codes: set[str] | None = None,
) -> list[dict[str, Any]]:
    path = Path(project_root) / "resources" / "stage_preprocess_v2_summary.tsv"
    if replace_codes is None:
        replace_codes = {_code(row) for row in results}
    return _merge_rows(path, SUMMARY_FIELDS, results, replace_codes)


def build_index_rows(
    results: Sequence[Mapping[str, Any]],
    *,
    analysis_dirname: str = "analysis_v2",
    clean_dirname: str = "clean_v2",
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for result in results:
        code = str(result.get("code", ""))
        artifact_code = str(result.get("artifact_code") or code)
        if result.get("status") not in {"OK", "SKIPPED"}:
            continue
        base = f"{analysis_dirname}/{artifact_code}"
        rows.append(
            {
                **result,
                "artifact_code": artifact_code,
                "pep": f"{base}/{artifact_code}.pep",
                "cds": f"{base}/{artifact_code}.cds",
                "gff": f"{base}/{artifact_code}.new.gff",
                "lens": f"{base}/{artifact_code}.lens",
                "all_pep": f"{base}/{artifact_code}.all.pep.fa",
                "all_cds": f"{base}/{artifact_code}.all.cds.fa",
                "all_gff": f"{base}/{artifact_code}.all.gff",
                "all_lens": f"{base}/{artifact_code}.all.lens",
                "gene_map": f"{clean_dirname}/{artifact_code}/{artifact_code}.gene_map.tsv",
                "chr_map": f"{clean_dirname}/{artifact_code}/{artifact_code}.chr_map.tsv",
                "rejected": f"{clean_dirname}/{artifact_code}/{artifact_code}.rejected.tsv",
                "qc": f"{base}/{artifact_code}.qc.json",
            }
        )
    return rows


def merge_index(
    project_root: Path | str,
    index_rows: Sequence[Mapping[str, Any]],
    *,
    replace_codes: set[str],
) -> list[dict[str, Any]]:
    path = Path(project_root) / "resources" / "wgdi_v2_index.tsv"
    return _merge_rows(path, INDEX_FIELDS, index_rows, replace_codes)


def manifest_sha256(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rebuild_manifest(
    project_root: Path | str,
    *,
    analysis_dirname: str = "analysis_v2",
    clean_dirname: str = "clean_v2",
    codes: Sequence[str] | None = None,
) -> Path:
    """Rebuild the SHA-256 manifest for every indexed species directory."""

    project_root = Path(project_root)
    index_path = project_root / "resources" / "wgdi_v2_index.tsv"
    if codes is None:
        codes = [_code(row) for row in read_tsv(index_path) if _code(row)]
    wanted = set(codes)
    rows: list[dict[str, object]] = []
    for layer in (clean_dirname, analysis_dirname):
        layer_root = project_root / layer
        for code in sorted(wanted, key=natural_key):
            code_dir = layer_root / code
            if not code_dir.is_dir():
                continue
            for path in sorted(code_dir.rglob("*")):
                if not path.is_file():
                    continue
                rows.append(
                    {
                        "code": code,
                        "layer": layer,
                        "relative_path": path.relative_to(project_root).as_posix(),
                        "bytes": path.stat().st_size,
                        "sha256": manifest_sha256(path),
                    }
                )
    manifest_path = project_root / "resources" / "wgdi_v2_manifest_sha256.tsv"
    write_tsv(manifest_path, MANIFEST_FIELDS, rows)
    return manifest_path


def finalize_project(
    project_root: Path | str,
    *,
    species_rows: Sequence[Mapping[str, Any]],
    results: Sequence[Mapping[str, Any]],
    processed_codes: set[str],
    analysis_dirname: str = "analysis_v2",
    clean_dirname: str = "clean_v2",
) -> dict[str, Path]:
    """Consolidate registry, summary, index and manifest after a run."""

    project_root = Path(project_root)
    merge_registry(project_root, species_rows, replace_codes=processed_codes)
    merge_summary(project_root, results, replace_codes=processed_codes)
    index_rows = build_index_rows(
        results,
        analysis_dirname=analysis_dirname,
        clean_dirname=clean_dirname,
    )
    merge_index(project_root, index_rows, replace_codes=processed_codes)
    manifest = rebuild_manifest(
        project_root,
        analysis_dirname=analysis_dirname,
        clean_dirname=clean_dirname,
    )
    return {
        "registry": project_root / "resources" / "species_registry_v2.tsv",
        "summary": project_root / "resources" / "stage_preprocess_v2_summary.tsv",
        "index": project_root / "resources" / "wgdi_v2_index.tsv",
        "manifest": manifest,
    }
