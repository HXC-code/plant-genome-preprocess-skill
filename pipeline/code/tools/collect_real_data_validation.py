#!/usr/bin/env python3
"""Collect a reproducible real-species preprocessing and subset-smoke receipt."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import tomllib
from datetime import datetime, timezone
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--smoke-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    smoke = json.loads(args.smoke_report.resolve().read_text(encoding="utf-8"))
    index = read_tsv(root / "resources" / "wgdi_v2_index.tsv")
    validation = read_tsv(root / "data-validation-report.tsv")
    naming = read_tsv(root / "resources" / "naming_audit.tsv")
    validation_by_code = {row["code"]: row for row in validation}
    naming_by_code = {row["code"]: row for row in naming}

    species: list[dict[str, object]] = []
    for row in index:
        code = row["code"]
        qc = json.loads((root / row["qc"]).read_text(encoding="utf-8"))
        report = validation_by_code.get(code, {})
        species.append(
            {
                "show_id": row["show_id"],
                "code": code,
                "species": row["species"],
                "accession": row["accession"],
                "wgdi_ready": row["wgdi_ready"],
                "validation": report.get("status", "MISSING"),
                "naming": naming_by_code.get(code, {}).get("status", "MISSING"),
                "genes_in_annotation": qc["n_genes"],
                "transcripts": qc["n_transcripts"],
                "representative_genes": qc["n_representative_genes"],
                "wgdi_genes": int(report.get("genes") or 0),
                "wgdi_sequences": int(report.get("sequences") or 0),
                "rejection_counts": qc["rejection_counts"],
                "input_sha256": {
                    role: data["sha256"] for role, data in qc["inputs"].items()
                },
            }
        )

    preprocessing_ok = all(
        item["validation"] == "PASS"
        and item["naming"] == "PASS"
        and item["wgdi_ready"] == "yes"
        for item in species
    )
    smoke_ok = smoke.get("status") == "PASS"
    report = {
        "schema_version": "real-species-validation-v1",
        "status": "PASS" if preprocessing_ok and smoke_ok else "FAIL",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "workflow_version": tomllib.loads(
            (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text(
                encoding="utf-8"
            )
        )["project"]["version"],
        "data_source": "NCBI Datasets v2alpha public RefSeq packages",
        "scope": (
            "Public RefSeq real-species preprocessing, validation, naming, QC, "
            "and a real-derived cross-species WGDI subset."
        ),
        "not_claimed": (
            "This is not a BUSCO lineage assessment, a complete full-proteome "
            "WGDI run, or a biological WGD conclusion."
        ),
        "species": species,
        "preprocessing_ok": preprocessing_ok,
        "smoke_report_sha256": sha256(args.smoke_report.resolve()),
        "smoke": smoke,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": report["status"],
                "species": len(species),
                "smoke_comparisons": len(smoke.get("comparisons", [])),
                "output": str(args.output),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
