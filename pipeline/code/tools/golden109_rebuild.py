"""Unified Golden-109 isolated rebuild orchestrator.

One entry point chains the canonical flow without deprecated converters:

    setup    create a fresh isolated project from the frozen manifest inputs
    mirror   hardlink (or copy) exactly the verified source files, flat per species
    run      canonical ``run_project`` in show_id-sorted batches, resume-safe
    validate canonical contract/ID/lens/QC/SHA-256 validation
    all      setup + mirror + run + validate

Old ``analysis_v2`` / ``clean_v2`` / legacy outputs are never read as inputs.
Every stage writes receipts; a single species failure never stops the others.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from pgci_preprocess.annotation_convert import convert_annotation  # noqa: E402
from pgci_preprocess.pipeline import PipelineOptions, run_project  # noqa: E402
from pgci_preprocess.species_names import (  # noqa: E402
    resolve_show_ids,
    resolve_species_code,
)
from pgci_preprocess.validate import ValidationOptions, validate_project  # noqa: E402
from pgci_preprocess.version import __version__  # noqa: E402
from build_golden109_inputs import species_folder_name  # noqa: E402


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _manifest_by_species(manifest_path: Path) -> dict[str, dict[str, dict[str, str]]]:
    by: dict[str, dict[str, dict[str, str]]] = {}
    for row in _read_tsv(manifest_path):
        by.setdefault(row["show_id"], {})[row["role"]] = row
    return by


def setup(args: argparse.Namespace) -> int:
    root = args.project_root
    if (root / "analysis_v2").exists() and any((root / "analysis_v2").iterdir()):
        print(f"ERROR: refusing setup into non-empty project: {root}", file=sys.stderr)
        return 1
    (root / "resources").mkdir(parents=True, exist_ok=True)
    (root / "data" / "downloads").mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.inventory, root / "resources" / "download_inventory.tsv")
    shutil.copyfile(args.assembly_map, root / "resources" / "assembly_accessions.tsv")
    species_names_source = PACKAGE_ROOT / "resources" / "species_abbreviations.tsv"
    shutil.copyfile(species_names_source, root / "resources" / "species_abbreviations.tsv")
    if args.format_exceptions:
        shutil.copyfile(args.format_exceptions, root / "resources" / "format_exceptions.tsv")
    receipt = {
        "stage": "setup",
        "created_utc": _now(),
        "pipeline_version": __version__,
        "project_root": str(root),
        "inputs": {
            "manifest": str(args.manifest),
            "manifest_sha256": _sha256(args.manifest),
            "inventory_sha256": _sha256(args.inventory),
            "assembly_map_sha256": _sha256(args.assembly_map),
            "format_exceptions_sha256": _sha256(args.format_exceptions) if args.format_exceptions else "",
            "species_names_sha256": _sha256(species_names_source),
        },
    }
    (root / "resources" / "rebuild_setup_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"[OK] setup {root} version={__version__}")
    return 0


def mirror(args: argparse.Namespace) -> int:
    root = args.project_root
    by_species = _manifest_by_species(args.manifest)
    inventory = _read_tsv(root / "resources" / "download_inventory.tsv")
    selected = (
        resolve_show_ids(
            root,
            [item.strip() for item in args.species.split(",") if item.strip()],
        )
        if args.species
        else None
    )
    receipts_dir = root / "resources" / "rebuild_receipts"
    receipts_dir.mkdir(parents=True, exist_ok=True)
    failures = 0
    for row in inventory:
        show_id = row["showid"].strip()
        if selected and show_id not in selected:
            continue
        code = resolve_species_code(root, show_id, latin_name=row["Species"])
        roles = by_species.get(show_id, {})
        gff = roles.get("gff", {})
        if not gff.get("source_path"):
            print(f"[BLOCKER] {code}: no verified gff source in manifest")
            failures += 1
            continue
        folder = species_folder_name(gff["source_path"])
        destination_dir = root / "data" / "downloads" / folder
        destination_dir.mkdir(parents=True, exist_ok=True)
        files: dict[str, dict[str, str]] = {}
        species_ok = True
        if gff.get("content_format") == "GBFF":
            # One original GBFF supplies gff/cds/pep (and optionally genome)
            # through the canonical adapter; derived names match the inventory.
            include_genome = bool(row.get("genome", "").strip())
            expected_names = {
                "gff": f"{code}.gff3.gz",
                "cds": f"{code}.cds.fa.gz",
                "pep": f"{code}.pep.fa.gz",
            }
            if include_genome:
                expected_names["genome"] = f"{code}.genome.fa.gz"
            receipt_path = receipts_dir / f"mirror_{code}.json"
            prior = json.loads(receipt_path.read_text(encoding="utf-8")) if receipt_path.is_file() else None
            reusable = (
                prior is not None
                and prior.get("status") == "OK"
                and prior.get("conversion", {}).get("source_sha256") == gff["sha256"]
                and all(
                    (destination_dir / name).is_file()
                    and _sha256(destination_dir / name) == prior["files"][role]["sha256"]
                    for role, name in expected_names.items()
                )
            )
            if reusable:
                print(f"[OK] mirror {code} roles={len(expected_names)} (verified derived reuse)")
                continue
            result = convert_annotation(
                Path(gff["source_path"]), destination_dir,
                prefix=code, include_genome=include_genome, overwrite=True,
            )
            conversion_manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
            if not conversion_manifest.get("complete", True):
                print(f"[FAIL] {code}: incomplete GBFF conversion")
                failures += 1
                continue
            for role, name in expected_names.items():
                info = conversion_manifest["outputs"][role]
                derived = destination_dir / info["file"]
                if derived.name != name or _sha256(derived) != info["sha256"]:
                    print(f"[FAIL] {code}/{role}: derived name/hash mismatch")
                    species_ok = False
                    continue
                files[role] = {
                    "source": gff["source_path"], "mirrored": str(derived),
                    "sha256": info["sha256"], "transfer": "gbff_adapter",
                }
            receipt = {
                "code": code, "show_id": show_id, "species": row["Species"],
                "folder": folder, "mirrored_utc": _now(),
                "status": "OK" if species_ok else "FAIL", "files": files,
                "conversion": {
                    "adapter": "pgci_preprocess.annotation_convert",
                    "source_sha256": gff["sha256"],
                    "manifest_path": str(result.manifest_path),
                    "manifest_sha256": _sha256(result.manifest_path),
                    "counts": conversion_manifest.get("counts", {}),
                },
            }
            receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2), encoding="utf-8")
            if species_ok:
                print(f"[OK] mirror {code} roles={len(files)} (gbff adapter)")
            else:
                failures += 1
            continue
        for role in ("gff", "cds", "pep", "genome"):
            entry = roles.get(role, {})
            name = row.get(role, "").strip()
            if not name:
                continue
            source = Path(entry["source_path"])
            destination = destination_dir / name
            if destination.exists():
                if _sha256(destination) == entry["sha256"]:
                    files[role] = {"source": str(source), "mirrored": str(destination),
                                   "sha256": entry["sha256"], "transfer": "already_verified"}
                    continue
                print(f"[FAIL] {code}/{role}: existing mirror hash mismatch, refusing overwrite")
                species_ok = False
                continue
            if not source.is_file():
                print(f"[FAIL] {code}/{role}: source missing {source}")
                species_ok = False
                continue
            try:
                os.link(source, destination)
                transfer = "hardlink"
            except OSError:
                shutil.copy2(source, destination)
                transfer = "copy"
            actual = _sha256(destination)
            if actual != entry["sha256"]:
                print(f"[FAIL] {code}/{role}: mirror hash mismatch")
                destination.unlink(missing_ok=True)
                species_ok = False
                continue
            files[role] = {"source": str(source), "mirrored": str(destination),
                           "sha256": actual, "transfer": transfer}
        receipt = {
            "code": code, "show_id": show_id, "species": row["Species"],
            "folder": folder, "mirrored_utc": _now(),
            "status": "OK" if species_ok else "FAIL", "files": files,
        }
        (receipts_dir / f"mirror_{code}.json").write_text(
            json.dumps(receipt, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        if species_ok:
            print(f"[OK] mirror {code} roles={len(files)}")
        else:
            failures += 1
    return 1 if failures else 0


def run_batches(args: argparse.Namespace) -> int:
    root = args.project_root
    log_dir = args.log_dir or root / "resources" / "rebuild_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    inventory = _read_tsv(root / "resources" / "download_inventory.tsv")
    show_ids = sorted(
        (row["showid"].strip() for row in inventory if row["showid"].strip()),
        key=int,
    )
    if args.species:
        selected = resolve_show_ids(
            root,
            [item.strip() for item in args.species.split(",") if item.strip()],
        )
        show_ids = [sid for sid in show_ids if sid in selected]
    options = PipelineOptions(project_root=root, overwrite=args.overwrite)
    batches = [show_ids[i:i + args.batch_size] for i in range(0, len(show_ids), args.batch_size)]
    total_fail = 0
    for number, batch in enumerate(batches, 1):
        started = _now()
        clock = time.monotonic()
        print(f"=== batch {number}/{len(batches)} species={','.join(batch)} ===", flush=True)
        try:
            results = run_project(root, show_ids=set(batch), options=options)
        except Exception as exc:  # noqa: BLE001
            results = [{"show_id": sid, "status": "ERROR", "reason": f"{type(exc).__name__}: {exc}"} for sid in batch]
        elapsed = round(time.monotonic() - clock, 2)
        failures = [r for r in results if r.get("status") in {"FAIL", "ERROR"}]
        total_fail += len(failures)
        log = {
            "batch": number, "started_utc": started, "finished_utc": _now(),
            "elapsed_seconds": elapsed, "pipeline_version": __version__,
            "command": f"run_project(show_ids={batch})",
            "options": {"gff_format": options.gff_format, "accession_policy": options.accession_policy,
                        "overwrite": options.overwrite},
            "results": results,
        }
        (log_dir / f"batch_{number:03d}.json").write_text(
            json.dumps(log, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        for item in results:
            print(f"[{item.get('status')}] {item.get('code', item.get('show_id', ''))} "
                  f"show_id={item.get('show_id', '')} genes={item.get('n_wgdi_genes','')} "
                  f"seqs={item.get('n_wgdi_sequences','')} reason={item.get('reason','')}", flush=True)
        print(f"batch {number}: ok={sum(1 for r in results if r.get('status')=='OK')} "
              f"skipped={sum(1 for r in results if r.get('status')=='SKIPPED')} fail={len(failures)} "
              f"elapsed={elapsed}s", flush=True)
    return 1 if total_fail else 0


def validate_stage(args: argparse.Namespace) -> int:
    report = validate_project(ValidationOptions(
        project_root=args.project_root,
        expected_count=args.expected_count,
        require_compat_links=False,
        report_path=args.report,
    ))
    print(f"validate ok={report.ok} report={args.report}")
    return 0 if report.ok else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="stage", required=True)

    def common(p: argparse.ArgumentParser) -> None:
        p.add_argument("--project-root", type=Path, required=True)
        p.add_argument("--manifest", type=Path, required=True)

    p_setup = sub.add_parser("setup")
    common(p_setup)
    p_setup.add_argument("--inventory", type=Path, required=True)
    p_setup.add_argument("--assembly-map", type=Path, required=True)
    p_setup.add_argument("--format-exceptions", type=Path)
    p_setup.set_defaults(func=setup)

    p_mirror = sub.add_parser("mirror")
    common(p_mirror)
    p_mirror.add_argument("--species", help="Comma-separated show_ids (default: all).")
    p_mirror.set_defaults(func=mirror)

    p_run = sub.add_parser("run")
    p_run.add_argument("--project-root", type=Path, required=True)
    p_run.add_argument("--batch-size", type=int, default=10)
    p_run.add_argument("--species")
    p_run.add_argument("--log-dir", type=Path)
    p_run.add_argument("--overwrite", action="store_true",
                       help="Explicitly replace existing validated bundles for the selected species.")
    p_run.set_defaults(func=run_batches)

    p_validate = sub.add_parser("validate")
    p_validate.add_argument("--project-root", type=Path, required=True)
    p_validate.add_argument("--expected-count", type=int, default=109)
    p_validate.add_argument("--report", type=Path, required=True)
    p_validate.set_defaults(func=validate_stage)

    p_all = sub.add_parser("all")
    common(p_all)
    p_all.add_argument("--inventory", type=Path, required=True)
    p_all.add_argument("--assembly-map", type=Path, required=True)
    p_all.add_argument("--format-exceptions", type=Path)
    p_all.add_argument("--batch-size", type=int, default=10)
    p_all.add_argument("--species")
    p_all.add_argument("--log-dir", type=Path)
    p_all.add_argument("--expected-count", type=int, default=109)
    p_all.add_argument("--report", type=Path, required=True)
    p_all.set_defaults(func=lambda a: (
        setup(a) or mirror(a) or run_batches(a) or validate_stage(a)
    ))

    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
