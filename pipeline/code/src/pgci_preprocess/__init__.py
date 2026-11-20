"""Reproducible NCBI annotation to WGDI preprocessing."""

from pgci_preprocess.annotation_convert import (
    ConversionError,
    ConversionResult,
    convert_annotation,
    convert_gtf_to_gff3,
    detect_conversion_format,
    normalize_functional_go_table,
)
from pgci_preprocess.gbff_convert import convert_gbff_to_bundle
from pgci_preprocess.audit import AuditReport, audit_project
from pgci_preprocess.batch import BatchResult, run_batch
from pgci_preprocess.gff_standardize import (
    GffStandardizationError,
    SpeciesNormalizationResult,
    normalize_gff_species,
    standardize_transcript_only_gff,
)
from pgci_preprocess.input_router import (
    AnnotationRoute,
    InputRoutingError,
    route_annotation,
    select_role_source,
)
from pgci_preprocess.doctor import run_doctor
from pgci_preprocess.project_init import ProjectInitError, ProjectInitResult, init_project
from pgci_preprocess.pipeline import (
    PipelineOptions,
    discover_species,
    load_format_exceptions,
    process_species,
    run_project,
)
from pgci_preprocess.sequence import validate_cds_pep
from pgci_preprocess.sequence_extract import (
    SequenceExtractionError,
    SequenceExtractionResult,
    extract_cds_pep_from_gff,
)
from pgci_preprocess.species_names import (
    SpeciesName,
    SpeciesNameMap,
    SpeciesNamingError,
    legacy_species_code,
    load_species_names,
    resolve_show_ids,
    resolve_species_code,
)
from pgci_preprocess.species_names_migration import (
    SpeciesNameMigrationResult,
    migrate_species_names,
)
from pgci_preprocess.validate import (
    ValidationOptions,
    ValidationReport,
    validate_project,
)
from pgci_preprocess.version import __version__  # single source of truth

__all__ = [
    "__version__",
    "AuditReport",
    "AnnotationRoute",
    "BatchResult",
    "ConversionError",
    "ConversionResult",
    "GffStandardizationError",
    "InputRoutingError",
    "PipelineOptions",
    "ProjectInitError",
    "ProjectInitResult",
    "SpeciesNormalizationResult",
    "SequenceExtractionError",
    "SequenceExtractionResult",
    "SpeciesName",
    "SpeciesNameMigrationResult",
    "SpeciesNameMap",
    "SpeciesNamingError",
    "ValidationOptions",
    "ValidationReport",
    "audit_project",
    "convert_annotation",
    "convert_gbff_to_bundle",
    "convert_gtf_to_gff3",
    "detect_conversion_format",
    "discover_species",
    "extract_cds_pep_from_gff",
    "init_project",
    "load_format_exceptions",
    "legacy_species_code",
    "load_species_names",
    "migrate_species_names",
    "normalize_functional_go_table",
    "normalize_gff_species",
    "process_species",
    "route_annotation",
    "run_batch",
    "run_doctor",
    "run_project",
    "resolve_show_ids",
    "resolve_species_code",
    "standardize_transcript_only_gff",
    "select_role_source",
    "validate_cds_pep",
    "validate_project",
]
