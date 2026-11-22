"""Batch orchestration: prepare -> normalize-gff -> run -> validate."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from pgci_preprocess.batch_prepare import (
    BatchPreparationError,
    BatchPreparationResult,
    prepare_batch,
)
from pgci_preprocess.gff_standardize import (
    GffStandardizationError,
    normalize_gff_species,
)
from pgci_preprocess.pipeline import PipelineOptions, run_project
from pgci_preprocess.validate import (
    ValidationOptions,
    ValidationReport,
    validate_project,
)


@dataclass(slots=True)
class BatchResult:
    candidates: tuple[str, ...] = ()
    prepared_files: int = 0
    normalized: int = 0
    converted: int = 0
    extracted: int = 0
    normalization_errors: list[str] = field(default_factory=list)
    run_results: list[dict[str, object]] = field(default_factory=list)
    validation: ValidationReport | None = None
    skipped_prepare_reason: str = ""
    prepare_failed: bool = False

    @property
    def ok(self) -> bool:
        run_ok = all(
            result.get("status") in {"OK", "SKIPPED"} for result in self.run_results
        )
        validation_ok = self.validation is None or self.validation.ok
        return (
            run_ok
            and validation_ok
            and not self.normalization_errors
            and not self.prepare_failed
        )


def run_batch(
    project_root: Path | str,
    source_status: Path | str,
    *,
    expected_count: int | None = None,
    batch_tag: str = "current",
    dry_run: bool = False,
    do_prepare: bool = True,
    do_normalize: bool = True,
    do_run: bool = True,
    do_validate: bool = True,
    pipeline_options: PipelineOptions | None = None,
    sequence_mode: str = "chromosomes",
    validate_expected_count: int | None = None,
    require_compat_links: bool = False,
    species_names_path: Path | str | None = None,
) -> BatchResult:
    project_root = Path(project_root).resolve()
    options = pipeline_options or PipelineOptions(
        project_root=project_root,
        species_names_path=(
            Path(species_names_path).expanduser().resolve()
            if species_names_path
            else None
        ),
    )
    if species_names_path:
        options.species_names_path = Path(species_names_path).expanduser().resolve()
    options.sequence_mode = sequence_mode
    result = BatchResult()

    prepared: BatchPreparationResult | None = None
    if do_prepare and not dry_run:
        try:
            prepared = prepare_batch(
                project_root,
                source_status,
                expected_count=expected_count,
                dry_run=False,
                batch_tag=batch_tag,
                species_names_path=species_names_path,
            )
        except BatchPreparationError as exc:
            result.skipped_prepare_reason = str(exc)
            if "no unprocessed complete species found" not in str(exc):
                result.prepare_failed = True
    else:
        try:
            probe = prepare_batch(
                project_root,
                source_status,
                expected_count=expected_count,
                dry_run=True,
                batch_tag=batch_tag,
                species_names_path=species_names_path,
            )
        except BatchPreparationError as exc:
            result.skipped_prepare_reason = str(exc)
            if "no unprocessed complete species found" not in str(exc):
                result.prepare_failed = True
        else:
            prepared = probe

    if prepared is not None:
        result.candidates = prepared.candidate_ids
        result.prepared_files = len(prepared.prepared_files)
        result.converted = prepared.converted
        result.extracted = prepared.extracted
    if not result.candidates:
        return result

    if do_normalize and not dry_run:
        try:
            normalizations = normalize_gff_species(
                project_root,
                list(result.candidates),
                update_inventory=True,
                species_names_path=species_names_path,
            )
        except GffStandardizationError as exc:
            result.normalization_errors.append(str(exc))
        else:
            result.normalized = sum(
                1 for item in normalizations if item.status == "CONVERTED"
            )
            result.normalization_errors.extend(
                f"s{item.show_id}: {item.detail}"
                for item in normalizations
                if item.status == "ERROR"
            )

    if do_run and not dry_run:
        result.run_results = run_project(
            project_root,
            show_ids=set(result.candidates),
            options=options,
        )

    if do_validate and not dry_run:
        result.validation = validate_project(
            ValidationOptions(
                project_root=project_root,
                expected_count=validate_expected_count,
                require_compat_links=require_compat_links,
            )
        )
    return result
