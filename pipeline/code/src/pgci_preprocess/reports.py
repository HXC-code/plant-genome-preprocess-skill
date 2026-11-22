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
    clean_dirname: str = "analysis_v2",
    project_root: Path | str | None = None,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    root = Path(project_root) if project_root is not None else None
    for result in results:
        code = str(result.get("code", ""))
        artifact_code = str(result.get("artifact_code") or code)
        if result.get("status") not in {"OK", "SKIPPED"}:
            continue
        base = f"{analysis_dirname}/{artifact_code}"
        all_paths = {
            "all_pep": f"{base}/{artifact_code}.all.pep.fa",
            "all_cds": f"{base}/{artifact_code}.all.cds.fa",
            "all_gff": f"{base}/{artifact_code}.all.gff",
            "all_lens": f"{base}/{artifact_code}.all.lens",
        }
        # ``.all.*`` is opt-in (run --keep-all). Index only files that exist so
        # validation never chases paths that were deliberately not written.
        if root is not None:
            for key, relative in all_paths.items():
                if not (root / relative).is_file():
                    all_paths[key] = ""
        rows.append(
            {
                **result,
                "artifact_code": artifact_code,
                "pep": f"{base}/{artifact_code}.pep",
                "cds": f"{base}/{artifact_code}.cds",
                "gff": f"{base}/{artifact_code}.new.gff",
                "lens": f"{base}/{artifact_code}.lens",
                **all_paths,
                "gene_map": f"{clean_dirname}/{artifact_code}/{artifact_code}.gene_map.tsv",
                "chr_map": f"{clean_dirname}/{artifact_code}/{artifact_code}.chr_map.tsv",
                "rejected": f"{clean_dirname}/{artifact_code}/{artifact_code}.removed_list.tsv",
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


def finalize_project(
    project_root: Path | str,
    *,
    species_rows: Sequence[Mapping[str, Any]],
    results: Sequence[Mapping[str, Any]],
    processed_codes: set[str],
    analysis_dirname: str = "analysis_v2",
    clean_dirname: str = "analysis_v2",
) -> dict[str, Path]:
    """Consolidate registry, summary and index after a run.

    No SHA-256 manifest is written: the project keeps one index and one QC
    record per species, and per-file provenance lives inside ``*.qc.json``.
    """

    project_root = Path(project_root)
    merge_registry(project_root, species_rows, replace_codes=processed_codes)
    merge_summary(project_root, results, replace_codes=processed_codes)
    index_rows = build_index_rows(
        results,
        analysis_dirname=analysis_dirname,
        clean_dirname=clean_dirname,
        project_root=project_root,
    )
    merge_index(project_root, index_rows, replace_codes=processed_codes)
    return {
        "registry": project_root / "resources" / "species_registry_v2.tsv",
        "summary": project_root / "resources" / "stage_preprocess_v2_summary.tsv",
        "index": project_root / "resources" / "wgdi_v2_index.tsv",
    }
