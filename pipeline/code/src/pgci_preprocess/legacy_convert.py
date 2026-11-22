"""Legacy folder-scan conversion path.

Historically, species that were absent from ``download_inventory.tsv`` were
converted by scanning the raw download tree for the largest plausible
GFF/CDS/PEP triplet.  That behaviour is retained here as an explicit,
parameterized module so it can be tested and so the old scripts can become thin
deprecation wrappers.  It is not the default path for inventory-driven runs.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from pgci_preprocess.pipeline import (
    PipelineOptions,
    SpeciesInput,
    process_species,
)
from pgci_preprocess.project_paths import downloads_root as resolve_downloads_root
from pgci_preprocess.reports import finalize_project, read_tsv
from pgci_preprocess.species_names import (
    load_species_names,
    resolve_species_code,
)

GFF_RE = re.compile(r"\.gff3?(\.(?:gz|bz2|xz|zip))?$", re.IGNORECASE)
COMPRESSION_RE = r"(?:\.(?:gz|bz2|xz|zip))?"
CDS_RE = re.compile(rf"\.(?:fa|fasta|fna|cds){COMPRESSION_RE}$", re.IGNORECASE)
PEP_RE = re.compile(rf"\.(?:fa|fasta|faa|pep){COMPRESSION_RE}$", re.IGNORECASE)
ALTERNATE_PARTS = {"alt_annotation", "manual_annotation"}


class LegacyConversionError(RuntimeError):
    pass


@dataclass(frozen=True, slots=True)
class BundleFiles:
    gff: Path
    cds: Path
    pep: Path


def _truthy(value: object) -> bool:
    return str(value or "").strip().casefold() in {"1", "true", "yes", "y"}


def default_downloads_root(project_root: Path | str) -> Path:
    """Resolve the raw download root for the legacy folder-scan converter.

    Prefers an explicit ``01_DATA_DOWNLOAD`` sibling of old lab trees, then the
    shared layout resolver (``data/`` with ``data/downloads`` compatibility).
    """

    project = Path(project_root).resolve()
    candidate = project.parents[2] / "01_DATA_DOWNLOAD" / "data" / "downloads"
    if candidate.is_dir():
        return candidate
    return resolve_downloads_root(project)


def find_bundle_files(folder: Path, *, include_alt: bool = False) -> BundleFiles | None:
    """Return the largest plausible GFF/CDS/PEP triplet in a species folder."""

    gff_candidates: list[Path] = []
    cds_candidates: list[Path] = []
    pep_candidates: list[Path] = []
    for directory_text, _subdirs, filenames in os.walk(folder):
        directory = Path(directory_text)
        if not include_alt and any(part in ALTERNATE_PARTS for part in directory.parts):
            continue
        for filename in filenames:
            lowered = filename.lower()
            if lowered.endswith(".part") or "chunks" in {part.lower() for part in directory.parts}:
                continue
            path = directory / filename
            if GFF_RE.search(lowered):
                gff_candidates.append(path)
            elif CDS_RE.search(lowered) and ("cds" in lowered or "cdna" in lowered):
                cds_candidates.append(path)
            elif PEP_RE.search(lowered) and re.search(r"pep|protein|\.pro\.", lowered):
                pep_candidates.append(path)
    if not (gff_candidates and cds_candidates and pep_candidates):
        return None
    return BundleFiles(
        gff=max(gff_candidates, key=lambda path: path.stat().st_size),
        cds=max(cds_candidates, key=lambda path: path.stat().st_size),
        pep=max(pep_candidates, key=lambda path: path.stat().st_size),
    )


def discover_species_directory(downloads_root: Path, show_id: str) -> Path | None:
    matches = [
        path
        for path in downloads_root.iterdir()
        if path.is_dir() and path.name.startswith(f"{show_id}_")
    ]
    if not matches:
        return None
    return max(matches, key=lambda path: path.stat().st_size)


def load_targets(status_table: Path | str) -> list[dict[str, str]]:
    rows = read_tsv(status_table)
    if not rows:
        raise LegacyConversionError(f"status table is empty or missing: {status_table}")
    return rows


def _row_value(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = str(row.get(name, "")).strip()
        if value:
            return value
    return ""


def build_species_input(
    project_root: Path | str,
    downloads_root: Path | str,
    row: dict[str, str],
    *,
    include_alt: bool | None = None,
) -> SpeciesInput:
    show_id = _row_value(row, "show_id", "showid")
    if not show_id:
        raise LegacyConversionError("target row has no show_id")
    folder_text = _row_value(row, "folder")
    if folder_text:
        directory = Path(folder_text)
    else:
        found = discover_species_directory(Path(downloads_root), show_id)
        if found is None:
            raise LegacyConversionError(f"no species directory for s{show_id}")
        directory = found
    if not directory.is_dir():
        raise LegacyConversionError(f"s{show_id}: species directory missing: {directory}")
    alt = _truthy(row.get("include_alt_annotation")) if include_alt is None else include_alt
    bundle = find_bundle_files(directory, include_alt=alt)
    if bundle is None:
        raise LegacyConversionError(f"s{show_id}: no GFF/CDS/PEP triplet under {directory}")
    accession = _row_value(row, "candidate_accession", "accession")
    species = _row_value(row, "Species", "species") or show_id
    order = _row_value(row, "Order", "order") or "NA"
    family = _row_value(row, "Family", "family") or "NA"
    code = resolve_species_code(
        project_root,
        show_id,
        latin_name=species,
    )
    return SpeciesInput(
        show_id=show_id,
        code=code,
        species=species,
        order=order,
        family=family,
        directory=directory,
        gff_name=bundle.gff.relative_to(directory).as_posix(),
        cds_name=bundle.cds.relative_to(directory).as_posix(),
        pep_name=bundle.pep.relative_to(directory).as_posix(),
        genome_name="",
        accession=accession,
    )


def run_legacy_conversion(
    project_root: Path | str,
    status_table: Path | str,
    *,
    downloads_root: Path | str | None = None,
    show_ids: Iterable[str] | None = None,
    exclude_ids: Iterable[str] | None = None,
    only_missing: bool = False,
    overwrite: bool = False,
    gff_format: str = "auto",
    protein_normalization: str = "auto",
    accession_policy: str = "legacy",
    min_genes_per_sequence: int = 50,
    target_gene_coverage: float = 0.90,
    max_sequences: int = 200,
) -> list[dict[str, object]]:
    project = Path(project_root).resolve()
    rows = load_targets(status_table)
    downloads = Path(downloads_root).resolve() if downloads_root else default_downloads_root(project)
    selected = {str(item).strip() for item in show_ids or [] if str(item).strip()}
    excluded = {str(item).strip() for item in exclude_ids or [] if str(item).strip()}
    names = load_species_names(project)
    indexed: set[str] = set()
    for index_row in read_tsv(project / "resources" / "wgdi_v2_index.tsv"):
        show_id = str(index_row.get("show_id", "")).strip()
        code = str(index_row.get("code", "")).strip()
        if show_id:
            indexed.add(show_id)
        elif names is not None and code:
            indexed.add(names.resolve_show_id(code))
        elif code.casefold().startswith("s") and code[1:].isdigit():
            indexed.add(code[1:])

    options = PipelineOptions(
        project_root=project,
        overwrite=overwrite,
        gff_format=gff_format,
        protein_normalization=protein_normalization,
        accession_policy=accession_policy,
        min_genes_per_sequence=min_genes_per_sequence,
        target_gene_coverage=target_gene_coverage,
        max_sequences=max_sequences,
    )
    species_rows: list[dict[str, object]] = []
    results: list[dict[str, object]] = []
    for row in rows:
        show_id = _row_value(row, "show_id", "showid")
        if not show_id:
            continue
        if selected and show_id not in selected:
            continue
        if show_id in excluded:
            continue
        if only_missing and show_id in indexed:
            continue
        try:
            species = build_species_input(project, downloads, row)
        except LegacyConversionError as exc:
            results.append(
                {
                    "show_id": show_id,
                    "code": resolve_species_code(
                        project,
                        show_id,
                        latin_name=_row_value(row, "Species", "species"),
                    ),
                    "species": _row_value(row, "Species", "species"),
                    "status": "ERROR",
                    "reason": "input_discovery_failed",
                    "detail": str(exc),
                }
            )
            continue
        species_rows.append(
            {
                "show_id": species.show_id,
                "code": species.code,
                "species": species.species,
                "order": species.order,
                "family": species.family,
                "accession": species.accession,
                "gff": species.gff_name,
                "cds": species.cds_name,
                "pep": species.pep_name,
                "genome": species.genome_name,
                "directory": species.directory.name,
            }
        )
        results.append(process_species(species, options))

    processed_codes = {str(result.get("code", "")) for result in results if result.get("code")}
    if species_rows:
        finalize_project(
            project,
            species_rows=species_rows,
            results=results,
            processed_codes=processed_codes,
        )
    return results
