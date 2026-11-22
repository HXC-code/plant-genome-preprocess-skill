"""Unified command line interface for ``python -m pgci_preprocess``."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from pgci_preprocess.annotation_convert import ConversionError, convert_annotation
from pgci_preprocess.artifact_migration import migrate_artifacts
from pgci_preprocess.audit import audit_project, render_text
from pgci_preprocess.batch import run_batch
from pgci_preprocess.batch_prepare import BatchPreparationError, prepare_batch
from pgci_preprocess.gff_standardize import GffStandardizationError, normalize_gff_species
from pgci_preprocess.legacy_convert import (
    LegacyConversionError,
    default_downloads_root,
    run_legacy_conversion,
)
from pgci_preprocess.doctor import run_doctor
from pgci_preprocess.naming_audit import audit_naming
from pgci_preprocess.naming_audit import render_text as render_naming
from pgci_preprocess.project_init import ProjectInitError, init_project
from pgci_preprocess.pipeline import PipelineOptions, run_project
from pgci_preprocess.organelle import MT_MIN_CORE_GENES
from pgci_preprocess.readiness import read_readiness_map
from pgci_preprocess.adjudication import (
    DECISIONS,
    gate_status,
    write_adjudication,
)
from pgci_preprocess.wgdi_integrity import check_comparison_outputs
from pgci_preprocess.source_manifest import (
    SourceManifestError,
    build_source_manifest,
    verify_source_manifest,
)
from pgci_preprocess.species_names_migration import migrate_species_names
from pgci_preprocess.validate import (
    ValidationOptions,
    resolve_layout,
    sync_compat_links,
    validate_project,
)
from pgci_preprocess.version import __version__

COMMANDS = {
    "prepare",
    "init-project",
    "doctor",
    "normalize-gff",
    "convert-annotation",
    "run",
    "validate",
    "readiness",
    "check-wgdi-outputs",
    "adjudicate",
    "adjudication-status",
    "chromosomes",
    "batch",
    "audit",
    "convert",
    "source-manifest",
    "migrate-species-names",
    "migrate-artifacts",
    "check-names",
}
def _default_project_root() -> Path:
    """Locate the project root without assuming code lives inside it.

    Resolution order:
      1. ``PGCI_PROJECT_ROOT`` environment variable (explicit override)
      2. sibling ``plant-genome-resources`` of the code tree (this repository's layout)
      3. the current working directory
    """
    environment = os.environ.get("PGCI_PROJECT_ROOT", "").strip()
    if environment:
        return Path(environment).expanduser().resolve()
    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / "plant-genome-resources"
        if (candidate / "resources").is_dir() or (candidate / "analysis_v2").is_dir():
            return candidate
    return Path.cwd()


PROJECT_ROOT = _default_project_root()
VERSION = __version__


def _parse_ids(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def _add_project_root(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_ROOT,
        help="Directory containing resources/ and data/downloads/ (default: this project).",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pgci_preprocess",
        description="Normalize plant genome annotations into validated WGDI input bundles.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")

    prepare = sub.add_parser(
        "prepare",
        help="Discover complete species, verify bundles, mirror files, refresh inventory.",
    )
    _add_project_root(prepare)
    prepare.add_argument(
        "--source-status",
        type=Path,
        required=True,
        help="Status table with show_id/folder/current_* columns.",
    )
    prepare.add_argument("--expected-count", type=int)
    prepare.add_argument("--batch-tag", default="current")
    prepare.add_argument("--dry-run", action="store_true")
    prepare.add_argument(
        "--species-names",
        type=Path,
        help="Canonical species abbreviation TSV (default: project resources).",
    )
    prepare.set_defaults(func=_cmd_prepare)

    initialize = sub.add_parser(
        "init-project",
        help="Create standard directories and optionally register one explicit species bundle.",
    )
    _add_project_root(initialize)
    initialize.add_argument("--show-id", default="", help="Numeric species/show ID.")
    initialize.add_argument("--scientific-name", default="", help="Species name used in the directory.")
    initialize.add_argument("--order", default="")
    initialize.add_argument("--family", default="")
    initialize.add_argument("--assembly-id", default="", help="Evidence-backed assembly identifier.")
    for role in ("genome", "gff", "cds", "pep"):
        initialize.add_argument(f"--{role}", type=Path, help=f"Explicit {role} source file.")
    initialize.add_argument("--json", action="store_true")
    initialize.add_argument(
        "--species-names",
        type=Path,
        help="Canonical species abbreviation TSV (default: project resources).",
    )
    initialize.set_defaults(func=_cmd_init_project)

    doctor = sub.add_parser(
        "doctor",
        help="Read-only checks for Python, project inputs and external workflow tools.",
    )
    _add_project_root(doctor)
    doctor.add_argument(
        "--workflow",
        default="preprocess",
        choices=["preprocess", "diamond", "blast", "wgdi", "full"],
    )
    doctor.add_argument("--require-gbff", action="store_true")
    doctor.add_argument("--allow-empty", action="store_true")
    doctor.add_argument("--json", action="store_true")
    doctor.add_argument("--output", type=Path, help="Also write the rendered report to this file.")
    for name in ("diamond", "blastp", "makeblastdb", "wgdi", "bash", "muscle", "yn00", "perl"):
        doctor.add_argument(f"--{name}", default="", help=f"Explicit {name} executable or wrapper path.")
    doctor.set_defaults(func=_cmd_doctor)

    normalize = sub.add_parser(
        "normalize-gff",
        help="Expand transcript-only GFF rows into gene/mRNA/CDS hierarchy.",
    )
    _add_project_root(normalize)
    normalize.add_argument("--species", required=True, help="Comma-separated show IDs.")
    normalize.add_argument("--overwrite", action="store_true")
    normalize.add_argument(
        "--no-update-inventory",
        action="store_true",
        help="Write normalized files without rewriting download_inventory.tsv.",
    )
    normalize.add_argument("--species-names", type=Path)
    normalize.set_defaults(func=_cmd_normalize_gff)

    convert_annotation_parser = sub.add_parser(
        "convert-annotation",
        help="Convert GTF, GBFF or functional GO tables without modifying raw inputs.",
    )
    convert_annotation_parser.add_argument("--input", type=Path, required=True, help="Source annotation file.")
    convert_annotation_parser.add_argument("--output-dir", type=Path, required=True, help="Derived output directory.")
    convert_annotation_parser.add_argument(
        "--format",
        dest="fmt",
        default="auto",
        choices=["auto", "gtf", "gbff", "functional-go"],
        help="Input content format (default: detect from content).",
    )
    convert_annotation_parser.add_argument("--prefix", default="", help="Output filename prefix.")
    convert_annotation_parser.add_argument("--overwrite", action="store_true")
    convert_annotation_parser.add_argument(
        "--no-genome",
        action="store_true",
        help="For GBFF, do not extract genome FASTA.",
    )
    convert_annotation_parser.add_argument(
        "--max-records",
        type=int,
        help="GBFF smoke test only: limit records and mark the manifest complete=false.",
    )
    convert_annotation_parser.add_argument("--json", action="store_true", help="Print a JSON result.")
    convert_annotation_parser.set_defaults(func=_cmd_convert_annotation)

    run = sub.add_parser(
        "run",
        help="Run one-gene-one-representative-transcript preprocessing.",
    )
    _add_project_root(run)
    run.add_argument("--species", default="", help="Comma-separated show IDs; empty means all.")
    run.add_argument("--min-genes-per-sequence", type=int, default=50)
    run.add_argument("--coverage", type=float, default=0.90)
    run.add_argument("--max-sequences", type=int, default=200)
    run.add_argument(
        "--min-genes-per-chromosome",
        type=int,
        default=100,
        help=(
            "Drop any WGDI chromosome with fewer than this many representative "
            "genes (default: %(default)s; 0 disables). Applied after sequence "
            "selection; exclusions are recorded in chr_map."
        ),
    )
    run.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing valid output directories (default: preserve them).",
    )
    run.add_argument(
        "--no-overwrite",
        dest="overwrite",
        action="store_false",
        help=argparse.SUPPRESS,
    )
    run.add_argument(
        "--gff-format",
        default="auto",
        choices=["auto", "standard", "transcript_only", "simplified5", "gene_span", "gene_cds_direct"],
    )
    run.add_argument(
        "--protein-normalization",
        default="auto",
        choices=["auto", "strict", "evm-dot"],
    )
    run.add_argument(
        "--accession-policy",
        default="strict",
        choices=["strict", "legacy"],
    )
    run.add_argument(
        "--no-transcript-region-backfill",
        action="store_true",
        help="Fail when a GFF has no region rows and no genome can resolve them.",
    )
    run.add_argument(
        "--organelle-pattern",
        default="",
        help="Regex for synthetic sequence IDs that must be excluded as organelles.",
    )
    run.add_argument(
        "--no-organelle-inference",
        action="store_true",
        help="Disable identifier and compact-region organelle inference.",
    )
    run.add_argument(
        "--organelle-mt-min-genes",
        type=int,
        default=MT_MIN_CORE_GENES,
        help=(
            "Distinct mitochondrial core gene families required to auto-exclude "
            "a synthetic region as mitochondrion (default: %(default)s)."
        ),
    )
    run.add_argument("--species-names", type=Path)
    run.set_defaults(func=_cmd_run)

    validate = sub.add_parser(
        "validate",
        help="Validate index, QC, ID sets, lens, SHA-256 and optional compat links.",
    )
    _add_project_root(validate)
    validate.add_argument("--expected-count", type=int)
    validate.add_argument("--no-compat-links", action="store_true")
    validate.add_argument("--skip-hashes", action="store_true")
    validate.add_argument("--refresh-links", action="store_true")
    validate.add_argument(
        "--allow-not-ready",
        action="store_true",
        help="Accept wgdi_ready=no bundles for archival/data-integrity validation.",
    )
    validate.add_argument(
        "--require-adjudication",
        action="store_true",
        help=(
            "Hard gate: every species needs a current, accepting human "
            "adjudication of its chromosome count (the lab SOP). Missing/stale/"
            "mismatched sign-offs FAIL. The skill enables this by default."
        ),
    )
    validate.add_argument("--report", type=Path)
    validate.set_defaults(func=_cmd_validate)

    readiness = sub.add_parser(
        "readiness",
        help="Read WGDI readiness for species (authoritative gate source).",
    )
    _add_project_root(readiness)
    readiness.add_argument(
        "--species",
        help="Comma-separated species codes; default: read folder_list.txt.",
    )
    readiness.add_argument(
        "--output",
        type=Path,
        help="Write a species_id<TAB>ready TSV here instead of stdout.",
    )
    readiness.set_defaults(func=_cmd_readiness)

    integrity = sub.add_parser(
        "check-wgdi-outputs",
        help="Content-integrity check of one comparison's 8 WGDI products.",
    )
    integrity.add_argument("--folder", type=Path, required=True)
    integrity.add_argument("--prefix", required=True)
    integrity.set_defaults(func=_cmd_check_wgdi_outputs)

    adjudicate = sub.add_parser(
        "adjudicate",
        help="Record a human sign-off on the final WGDI chromosome count.",
    )
    _add_project_root(adjudicate)
    adjudicate.add_argument("--species", required=True, help="Species code.")
    adjudicate.add_argument("--confirmed-n", type=int, required=True)
    adjudicate.add_argument(
        "--decision", required=True, choices=sorted(DECISIONS)
    )
    adjudicate.add_argument("--reviewer", required=True)
    adjudicate.add_argument("--rationale", default="")
    adjudicate.set_defaults(func=_cmd_adjudicate)

    adj_status = sub.add_parser(
        "adjudication-status",
        help="Report the adjudication gate status per species (authoritative).",
    )
    _add_project_root(adj_status)
    adj_status.add_argument(
        "--species", help="Comma-separated codes; default: read folder_list.txt."
    )
    adj_status.add_argument("--output", type=Path)
    adj_status.set_defaults(func=_cmd_adjudication_status)

    chromosomes = sub.add_parser(
        "chromosomes",
        help="Print the per-species chromosome inventory (count + genes each).",
    )
    _add_project_root(chromosomes)
    chromosomes.add_argument(
        "--species", help="Comma-separated codes; default: read folder_list.txt."
    )
    chromosomes.set_defaults(func=_cmd_chromosomes)

    batch = sub.add_parser(
        "batch",
        help="Chain prepare -> optional normalize-gff -> run -> validate.",
    )
    _add_project_root(batch)
    batch.add_argument("--source-status", type=Path, required=True)
    batch.add_argument("--expected-count", type=int)
    batch.add_argument("--batch-tag", default="current")
    batch.add_argument("--dry-run", action="store_true")
    batch.add_argument("--skip-prepare", action="store_true")
    batch.add_argument("--skip-normalize", action="store_true")
    batch.add_argument("--skip-run", action="store_true")
    batch.add_argument("--skip-validate", action="store_true")
    batch.add_argument("--validate-expected-count", type=int)
    batch.add_argument("--require-compat-links", action="store_true")
    batch.add_argument("--species-names", type=Path)
    batch.set_defaults(func=_cmd_batch)

    audit = sub.add_parser(
        "audit",
        help="Read-only summary of species, suffixes, ready distribution and failures.",
    )
    _add_project_root(audit)
    audit.add_argument(
        "--species",
        default="",
        help="Optional comma-separated show IDs or s-codes for detailed inspection.",
    )
    audit.add_argument("--json", action="store_true", help="Print JSON instead of text.")
    audit.add_argument("--output", type=Path, help="Optional explicit report path.")
    audit.set_defaults(func=_cmd_audit)

    convert = sub.add_parser(
        "convert",
        help="Legacy folder-scan conversion for species absent from inventory (deprecated path).",
    )
    _add_project_root(convert)
    convert.add_argument(
        "--status-table",
        type=Path,
        default=PROJECT_ROOT / "resources" / "legacy_conversion_targets.tsv",
    )
    convert.add_argument("--downloads-root", type=Path)
    convert.add_argument("--species", default="", help="Comma-separated show IDs; empty means all.")
    convert.add_argument("--exclude", default="", help="Comma-separated show IDs to skip.")
    convert.add_argument("--only-missing", action="store_true", help="Skip codes already in the index.")
    convert.add_argument("--include-alt", action="store_true")
    convert.add_argument("--overwrite", action="store_true")
    convert.add_argument("--gff-format", default="auto")
    convert.add_argument("--protein-normalization", default="auto", choices=["auto", "strict", "evm-dot"])
    convert.add_argument("--accession-policy", default="legacy", choices=["strict", "legacy"])
    convert.set_defaults(func=_cmd_convert)

    source_manifest = sub.add_parser(
        "source-manifest",
        help="Build or verify a role-level source manifest for a frozen species index.",
    )
    manifest_sub = source_manifest.add_subparsers(dest="manifest_command", required=True)
    manifest_build = manifest_sub.add_parser(
        "build", help="Build the manifest (and .sha256 sidecar) from index + QC + repairs."
    )
    _add_project_root(manifest_build)
    manifest_build.add_argument("--index", type=Path, help="Frozen index TSV (default: resources/wgdi_v2_index.tsv).")
    manifest_build.add_argument("--repairs", type=Path, help="Evidence-backed repair TSV applied before verification.")
    manifest_build.add_argument(
        "--raw-root", type=Path, action="append", default=[],
        help="Read-only raw download root searched for hash-verified remaps (repeatable).",
    )
    manifest_build.add_argument("--output", type=Path, required=True)
    manifest_build.add_argument("--json", action="store_true")
    manifest_verify = manifest_sub.add_parser(
        "verify", help="Verify the frozen manifest hash and every referenced source file."
    )
    manifest_verify.add_argument("--manifest", type=Path, required=True)
    manifest_verify.add_argument("--no-rehash", action="store_true", help="Check content hash and sizes only.")
    manifest_verify.add_argument("--json", action="store_true")
    source_manifest.set_defaults(func=_cmd_source_manifest)

    migrate_names = sub.add_parser(
        "migrate-species-names",
        help="Update legacy s<show_id> metadata codes to abbreviations (index/QC only).",
    )
    _add_project_root(migrate_names)
    migrate_names.add_argument("--dry-run", action="store_true")
    migrate_names.add_argument("--skip-links", action="store_true")
    migrate_names.set_defaults(func=_cmd_migrate_species_names)

    migrate_artifacts_parser = sub.add_parser(
        "migrate-artifacts",
        help="Rename legacy s<show_id> bundles, file names and gene IDs to abbreviations.",
    )
    _add_project_root(migrate_artifacts_parser)
    migrate_artifacts_parser.add_argument("--dry-run", action="store_true")
    migrate_artifacts_parser.add_argument(
        "--species",
        default="",
        help="Comma-separated show IDs; empty means every legacy bundle.",
    )
    migrate_artifacts_parser.add_argument(
        "--skip-links",
        action="store_true",
        help="Do not refresh the workspace compatibility directories.",
    )
    migrate_artifacts_parser.add_argument(
        "--prune-legacy-dirs",
        action="store_true",
        help="Move byte-identical workspace/s<show_id> copies to _quarantine.",
    )
    migrate_artifacts_parser.add_argument("--json", action="store_true")
    migrate_artifacts_parser.set_defaults(func=_cmd_migrate_artifacts)

    check_names = sub.add_parser(
        "check-names",
        help="Audit the abbreviation table against downloads, index and bundles.",
    )
    _add_project_root(check_names)
    check_names.add_argument("--table", type=Path, help="Species abbreviation TSV override.")
    check_names.add_argument("--json", action="store_true")
    check_names.set_defaults(func=_cmd_check_names)
    return parser


def _cmd_prepare(args: argparse.Namespace) -> int:
    try:
        result = prepare_batch(
            args.project_root,
            args.source_status,
            expected_count=args.expected_count,
            dry_run=args.dry_run,
            batch_tag=args.batch_tag,
            species_names_path=args.species_names,
        )
    except BatchPreparationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"candidates={len(result.candidate_ids)}: {','.join(result.candidate_ids)}")
    print(f"files={len(result.prepared_files)}")
    print(f"converted={result.converted} extracted={result.extracted}")
    if result.manifest_path:
        print(f"manifest={result.manifest_path}")
    if result.suffix_summary_path:
        print(f"suffix_summary={result.suffix_summary_path}")
    if result.routing_path:
        print(f"routing={result.routing_path}")
    return 0


def _cmd_init_project(args: argparse.Namespace) -> int:
    sources = {
        role: getattr(args, role)
        for role in ("genome", "gff", "cds", "pep")
        if getattr(args, role) is not None
    }
    try:
        result = init_project(
            args.project_root,
            show_id=args.show_id,
            scientific_name=args.scientific_name,
            order=args.order,
            family=args.family,
            assembly_id=args.assembly_id,
            sources=sources,
            species_names_path=args.species_names,
        )
    except (ProjectInitError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    payload = {"status": "OK", **result.as_dict()}
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        print(f"[OK] project={result.project_root}")
        print(f"inventory={result.inventory_path}")
        if result.show_id:
            print(f"registered={result.species_code} show_id={result.show_id}")
            print(f"species_directory={result.species_directory}")
        else:
            print("registered=none (empty scaffold)")
    return 0


def _cmd_doctor(args: argparse.Namespace) -> int:
    overrides = {
        name: getattr(args, name)
        for name in ("diamond", "blastp", "makeblastdb", "wgdi", "bash", "muscle", "yn00", "perl")
        if getattr(args, name)
    }
    report = run_doctor(
        args.project_root,
        workflow=args.workflow,
        require_gbff=args.require_gbff,
        allow_empty=args.allow_empty,
        executables=overrides,
    )
    if args.json:
        rendered = json.dumps(report, ensure_ascii=False, indent=2)
    else:
        checks = report.get("checks")
        checks = checks if isinstance(checks, list) else []
        lines = [
            f"[{check['status']}] {check['name']}: {check['detail']}"
            for check in checks
        ]
        counts = report.get("counts")
        counts = counts if isinstance(counts, dict) else {}
        lines.append(f"summary: {counts['PASS']} PASS / {counts['WARN']} WARN / {counts['FAIL']} FAIL")
        rendered = "\n".join(lines)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
        if not args.json:
            print(f"report={args.output}")
    return 1 if report["status"] == "FAIL" else 0


def _cmd_convert_annotation(args: argparse.Namespace) -> int:
    try:
        result = convert_annotation(
            args.input,
            args.output_dir,
            fmt=args.fmt,
            prefix=args.prefix,
            overwrite=args.overwrite,
            include_genome=not args.no_genome,
            max_records=args.max_records,
        )
    except ConversionError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    payload = {
        "status": "OK",
        "format": result.format,
        "complete": result.complete,
        "source": str(result.source_path),
        "output_dir": str(result.output_dir),
        "outputs": {role: str(path) for role, path in sorted(result.outputs.items())},
        "manifest": str(result.manifest_path),
        "counts": result.counts,
        "warnings": list(result.warnings),
    }
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        print(f"[OK] format={result.format} complete={str(result.complete).lower()}")
        for role, path in sorted(result.outputs.items()):
            print(f"{role}={path}")
        print(f"manifest={result.manifest_path}")
        for warning in result.warnings:
            print(f"WARNING: {warning}")
    return 0 if result.complete else 1


def _cmd_normalize_gff(args: argparse.Namespace) -> int:
    ids = _parse_ids(args.species)
    if not ids:
        print("ERROR: --species is empty", file=sys.stderr)
        return 2
    try:
        results = normalize_gff_species(
            args.project_root,
            ids,
            update_inventory=not args.no_update_inventory,
            overwrite=args.overwrite,
            species_names_path=args.species_names,
        )
    except GffStandardizationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    failures = 0
    for item in results:
        print(
            f"[{item.status}] {item.species_code or item.show_id} "
            f"show_id={item.show_id} {item.detail}"
        )
        if item.status == "ERROR":
            failures += 1
    return 1 if failures else 0


def _cmd_run(args: argparse.Namespace) -> int:
    if args.organelle_mt_min_genes < 1:
        print("ERROR: --organelle-mt-min-genes must be >= 1", file=sys.stderr)
        return 2
    ids = set(_parse_ids(args.species)) or None
    options = PipelineOptions(
        project_root=args.project_root,
        min_genes_per_sequence=args.min_genes_per_sequence,
        target_gene_coverage=args.coverage,
        max_sequences=args.max_sequences,
        min_genes_per_chromosome=args.min_genes_per_chromosome,
        overwrite=args.overwrite,
        gff_format=args.gff_format,
        protein_normalization=args.protein_normalization,
        accession_policy=args.accession_policy,
        backfill_regions_from_transcripts=not args.no_transcript_region_backfill,
        infer_organelles=not args.no_organelle_inference,
        organelle_pattern=args.organelle_pattern,
        organelle_mt_min_genes=args.organelle_mt_min_genes,
        species_names_path=args.species_names,
    )
    try:
        results = run_project(args.project_root, show_ids=ids, options=options)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if not results:
        print("ERROR: no species selected or discovered", file=sys.stderr)
        return 2
    failures = sum(
        1 for result in results if result.get("status") in {"FAIL", "ERROR"}
    )
    skipped = sum(1 for result in results if result.get("status") == "SKIPPED")
    print(f"processed={len(results)} skipped={skipped} failures={failures}")
    return 1 if failures else 0


def _cmd_validate(args: argparse.Namespace) -> int:
    try:
        workspace, canonical_root = resolve_layout(args.project_root)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if args.refresh_links:
        for message in sync_compat_links(canonical_root, workspace):
            print(message)
    try:
        result = validate_project(
            ValidationOptions(
                project_root=workspace,
                expected_count=args.expected_count,
                require_compat_links=not args.no_compat_links,
                verify_hashes=not args.skip_hashes,
                allow_not_ready=args.allow_not_ready,
                require_adjudication=args.require_adjudication,
                report_path=args.report,
            )
        )
    except OSError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    for row in result.rows:
        print(f"[{row['status']}] {row['code']} {row['problems']}")
    passed = sum(1 for row in result.rows if row["status"] == "PASS")
    failed = sum(1 for row in result.rows if row["status"] == "FAIL")
    print(f"合计: {passed} PASS / {failed} FAIL")
    print(f"报告: {result.report_path}")
    return 0 if result.ok else 1


def _cmd_readiness(args: argparse.Namespace) -> int:
    if args.species:
        ids = [item.strip() for item in args.species.split(",") if item.strip()]
    else:
        folder_list = args.project_root / "folder_list.txt"
        ids = (
            [
                line.strip()
                for line in folder_list.read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            if folder_list.is_file()
            else []
        )
    mapping = read_readiness_map(args.project_root, ids)
    text = "".join(f"{species_id}\t{mapping[species_id]}\n" for species_id in ids)
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
    return 0


def _cmd_check_wgdi_outputs(args: argparse.Namespace) -> int:
    problems = check_comparison_outputs(args.folder, args.prefix)
    for problem in problems:
        print(problem)
    return 1 if problems else 0


def _species_ids(args: argparse.Namespace) -> list[str]:
    if args.species:
        return [item.strip() for item in args.species.split(",") if item.strip()]
    folder_list = args.project_root / "folder_list.txt"
    if folder_list.is_file():
        return [
            line.strip()
            for line in folder_list.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    return []


def _cmd_adjudicate(args: argparse.Namespace) -> int:
    try:
        path = write_adjudication(
            args.project_root,
            args.species,
            confirmed_n=args.confirmed_n,
            decision=args.decision,
            reviewer=args.reviewer,
            rationale=args.rationale,
        )
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    status = gate_status(args.project_root, args.species)
    print(
        f"[{status}] {args.species} confirmed_n={args.confirmed_n} "
        f"decision={args.decision} reviewer={args.reviewer}"
    )
    print(f"记录: {path}")
    return 0


def _cmd_adjudication_status(args: argparse.Namespace) -> int:
    ids = _species_ids(args)
    text = "".join(
        f"{species_id}\t{gate_status(args.project_root, species_id)}\n"
        for species_id in ids
    )
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
    return 0


def _cmd_chromosomes(args: argparse.Namespace) -> int:
    import csv

    ids = _species_ids(args)
    missing = 0
    for code in ids:
        path = args.project_root / "analysis_v2" / code / f"{code}.chromosomes.tsv"
        if not path.is_file():
            path = args.project_root / "clean_v2" / code / f"{code}.chromosomes.tsv"
        if not path.is_file():
            print(f"{code}: 缺少染色体清单（先 run）")
            missing += 1
            continue
        with path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        total = sum(int(row["gene_count"]) for row in rows)
        print(f"{code}: {len(rows)} 条染色体，共 {total} 个基因")
        print("chromosome\tgene_count\tlength\tsource_seqid\tis_primary")
        for row in rows:
            print(
                f"{row['chromosome']}\t{row['gene_count']}\t{row['length']}"
                f"\t{row['source_seqid']}\t{row['is_primary']}"
            )
    return 1 if missing else 0


def _cmd_batch(args: argparse.Namespace) -> int:
    try:
        result = run_batch(
            args.project_root,
            args.source_status,
            expected_count=args.expected_count,
            batch_tag=args.batch_tag,
            dry_run=args.dry_run,
            do_prepare=not args.skip_prepare,
            do_normalize=not args.skip_normalize,
            do_run=not args.skip_run,
            do_validate=not args.skip_validate,
            validate_expected_count=args.validate_expected_count,
            require_compat_links=args.require_compat_links,
            species_names_path=args.species_names,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"candidates={len(result.candidates)}: {','.join(result.candidates)}")
    print(
        f"prepared_files={result.prepared_files} converted={result.converted} "
        f"extracted={result.extracted} normalized={result.normalized}"
    )
    for message in result.normalization_errors:
        print(f"[NORMALIZE-ERROR] {message}", file=sys.stderr)
    for item in result.run_results:
        print(
            f"[{item.get('status')}] {item.get('code', item.get('show_id', ''))} "
            f"show_id={item.get('show_id', '')} "
            f"genes={item.get('n_wgdi_genes', '')} "
            f"seqs={item.get('n_wgdi_sequences', '')} "
            f"reason={item.get('reason', '')}"
        )
    if result.validation is not None:
        passed = sum(1 for row in result.validation.rows if row["status"] == "PASS")
        failed = sum(1 for row in result.validation.rows if row["status"] == "FAIL")
        print(f"validation: {passed} PASS / {failed} FAIL")
    if result.skipped_prepare_reason:
        print(f"prepare skipped: {result.skipped_prepare_reason}")
    return 0 if result.ok else 1


def _cmd_audit(args: argparse.Namespace) -> int:
    ids = set(_parse_ids(args.species)) or None
    report = audit_project(args.project_root, species=ids)
    if args.json:
        rendered = json.dumps(report.as_dict(), ensure_ascii=False, indent=2)
    else:
        rendered = render_text(report)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
        print(f"report={args.output}")
    return 1 if report.problems else 0


def _cmd_convert(args: argparse.Namespace) -> int:
    ids = _parse_ids(args.species) or None
    excludes = set(_parse_ids(args.exclude))
    downloads = args.downloads_root or default_downloads_root(args.project_root)
    try:
        results = run_legacy_conversion(
            args.project_root,
            args.status_table,
            downloads_root=downloads,
            show_ids=ids,
            exclude_ids=excludes,
            only_missing=args.only_missing,
            overwrite=args.overwrite,
            gff_format=args.gff_format,
            protein_normalization=args.protein_normalization,
            accession_policy=args.accession_policy,
        )
    except LegacyConversionError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    failures = 0
    for item in results:
        print(
            f"[{item.get('status')}] {item.get('code', item.get('show_id', ''))} "
            f"show_id={item.get('show_id', '')} {item.get('species', '')} "
            f"genes={item.get('n_wgdi_genes', '')} seqs={item.get('n_wgdi_sequences', '')} "
            f"reason={item.get('reason', '')}"
        )
        if item.get("status") in {"FAIL", "ERROR"}:
            failures += 1
    return 1 if failures else 0


def _cmd_source_manifest(args: argparse.Namespace) -> int:
    try:
        if args.manifest_command == "build":
            _rows, summary = build_source_manifest(
                args.project_root,
                index_path=args.index,
                repairs_path=args.repairs,
                raw_roots=tuple(args.raw_root),
                output_path=args.output,
            )
            if args.json:
                print(json.dumps(summary, ensure_ascii=False, indent=2))
            else:
                print(f"manifest={summary.get('manifest_path', '')}")
                print(f"sha256={summary.get('manifest_sha256', '')}")
                status_counts = summary.get("status_counts")
                status_counts = status_counts if isinstance(status_counts, dict) else {}
                blockers = summary.get("blockers")
                blockers = blockers if isinstance(blockers, list) else []
                print(f"species={summary.get('species_count', 0)} rows={summary.get('row_count', 0)}")
                for status, count in sorted(status_counts.items()):
                    print(f"  {status}: {count}")
                for blocker in blockers:
                    print(f"  BLOCKER {blocker}")
            blockers = summary.get("blockers")
            return 2 if isinstance(blockers, list) and blockers else 0
        report = verify_source_manifest(args.manifest, rehash=not args.no_rehash)
        if args.json:
            print(json.dumps(report, ensure_ascii=False, indent=2))
        else:
            print(f"manifest={report['manifest']}")
            print(f"sha256={report['manifest_sha256']}")
            print(f"rows={report['rows']} files_rehashed={report['files_rehashed']}")
            problems = report.get("problems")
            problems = problems if isinstance(problems, list) else []
            for problem in problems:
                print(f"  PROBLEM {problem}")
            print("[OK]" if report["ok"] else "[FAIL]")
        return 0 if bool(report.get("ok")) else 2
    except SourceManifestError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


def _cmd_migrate_species_names(args: argparse.Namespace) -> int:
    try:
        result = migrate_species_names(
            args.project_root,
            dry_run=args.dry_run,
            refresh_links=not args.skip_links,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    for name, changed in sorted(result.changed_files.items()):
        print(f"{name}: rows_changed={changed}")
    if result.backup_dir:
        print(f"backup={result.backup_dir}")
    print(f"links={len(result.links)}")
    print(f"unresolved={len(result.unresolved)}")
    for problem in result.unresolved[:20]:
        print(f"UNRESOLVED: {problem}", file=sys.stderr)
    if not result.dry_run:
        print(f"report={result.report_path}")
    return 1 if result.unresolved else 0


def _cmd_migrate_artifacts(args: argparse.Namespace) -> int:
    show_ids = set(_parse_ids(args.species)) or None
    try:
        result = migrate_artifacts(
            args.project_root,
            dry_run=args.dry_run,
            refresh_links=not args.skip_links,
            prune_legacy_dirs=args.prune_legacy_dirs,
            show_ids=show_ids,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    payload = result.as_dict()
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        counts: dict[str, int] = {}
        for item in result.species:
            counts[item.status] = counts.get(item.status, 0) + 1
        print("artifact naming migration")
        for status, count in sorted(counts.items()):
            print(f"  {status}: {count}")
        print(f"  files_renamed={sum(item.files_renamed for item in result.species)}")
        print(f"  gene_ids_replaced={sum(item.ids_replaced for item in result.species)}")
        print(f"  index_rows_changed={result.index_rows_changed}")
        print(f"  links={len(result.links)}")
        if result.quarantine_dir:
            print(f"  quarantine={result.quarantine_dir}")
        for item in result.species:
            if item.status == "FAIL":
                print(
                    f"[FAIL] {item.show_id} {item.code}: {';'.join(item.problems)}",
                    file=sys.stderr,
                )
        if not result.dry_run:
            print(f"report={result.report_path}")
    return 0 if not any(item.status == "FAIL" for item in result.species) else 1


def _cmd_check_names(args: argparse.Namespace) -> int:
    try:
        report = audit_naming(args.project_root, table_path=args.table)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if args.json:
        print(json.dumps(report.as_dict(), ensure_ascii=False, indent=2))
    else:
        print(render_naming(report))
    return 0 if report.ok else 1


def main(argv: list[str] | None = None) -> int:
    # A non-UTF-8 Windows console (cp1252) must not turn a localized message
    # into a crash; degrade unencodable characters instead of failing the run.
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            try:
                reconfigure(errors="replace")
            except (ValueError, OSError):
                pass
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] not in COMMANDS and not argv[0].startswith("-"):
        print(f"ERROR: unknown command {argv[0]!r}; choose from {', '.join(sorted(COMMANDS))}", file=sys.stderr)
        return 2
    # Legacy wrapper compatibility: old entry points passed run flags directly.
    if argv and argv[0] not in {"-h", "--help", "--version"} and argv[0].startswith("-"):
        argv = ["run", *argv]
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        parser.print_help()
        return 2
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
