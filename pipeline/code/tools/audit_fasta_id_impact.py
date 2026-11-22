"""Read-only audit: quantify how a FASTA identifier rule affects real inputs.

Format adapters in :mod:`pgci_preprocess.fasta` only matter for headers that
carry no explicit ``[protein_id=...]`` or ``Protein=`` field.  Before changing
any identifier rule, run this audit over the inventory actually referenced by the
project so the blast radius is a number, not a guess.  Nothing is written.

Usage::

    python tools/audit_fasta_id_impact.py [--project-root PATH] [--limit N]

``--limit`` caps CDS files audited per role so a spot check stays fast; the
default audits every inventory CDS file.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import TextIO

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_ROOT / "src"))

from pgci_preprocess.fasta import (  # noqa: E402
    CDS_FALLBACK_RE,
    PROTEIN_FIELD_RE,
    PROTEIN_ID_RE,
    open_text,
)

# The pattern shipped in 2.3.0: ``\\s``/``\\d`` inside a raw string matched
# literal backslashes, so the fallback never fired as documented.
LEGACY_CDS_FALLBACK_RE = re.compile(
    r"(?:^|_)cds_([^_\\s]+?)(?:_\\d+)?(?:\\s|$)", re.IGNORECASE
)


def _headers(handle: TextIO):
    for line in handle:
        if line.startswith(">"):
            yield line[1:].rstrip("\r\n")


def audit_cds_file(path: Path, legacy: bool = False) -> dict[str, object]:
    """Resolve every CDS header to an identifier and collect the resulting set."""

    pattern = LEGACY_CDS_FALLBACK_RE if legacy else CDS_FALLBACK_RE
    explicit = 0
    fallback_hits = 0
    plain = 0
    identifiers: set[str] = set()
    with open_text(path) as handle:
        for header in _headers(handle):
            if PROTEIN_ID_RE.search(header) or PROTEIN_FIELD_RE.search(header):
                explicit += 1
                identifiers.add(
                    (PROTEIN_ID_RE.search(header) or PROTEIN_FIELD_RE.search(header)).group(1).strip()
                )
                continue
            token = header.split(None, 1)[0].strip()
            match = pattern.search(header if legacy else token)
            if match:
                identifiers.add(match.group(1).strip())
                fallback_hits += 1
            else:
                identifiers.add(token)
                plain += 1
    return {
        "file": path.name,
        "headers": explicit + fallback_hits + plain,
        "explicit_protein_id": explicit,
        "cds_fallback_applied": fallback_hits,
        "first_token_only": plain,
        "distinct_ids": len(identifiers),
        "_identifiers": identifiers,
    }


def compare_rules(path: Path) -> dict[str, object]:
    """Report whether the active rule resolves any header differently.

    ``id_set_identical`` compares the resolved identifier sets, not counters, so
    a rule change that renames even one gene is reported.
    """

    current = audit_cds_file(path)
    legacy = audit_cds_file(path, legacy=True)
    changed = sorted(legacy.pop("_identifiers") ^ current.pop("_identifiers"))
    return {
        **current,
        "legacy_distinct_ids": legacy["distinct_ids"],
        "id_set_identical": not changed,
        "changed_examples": changed[:5],
    }


def inventory_cds_files(project_root: Path) -> list[Path]:
    inventory = project_root / "resources" / "download_inventory.tsv"
    if not inventory.is_file():
        raise FileNotFoundError(inventory)
    downloads = project_root / "data"
    legacy = downloads / "downloads"
    if legacy.is_dir() and any(p.is_dir() and "_" in p.name for p in legacy.iterdir()):
        downloads = legacy
    folders = {
        path.name.split("_", 1)[0]: path
        for path in downloads.iterdir()
        if path.is_dir() and "_" in path.name
    }
    files: list[Path] = []
    with inventory.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            show_id = (row.get("showid") or "").strip()
            cds_name = (row.get("cds") or "").strip()
            folder = folders.get(show_id)
            if not cds_name or folder is None:
                continue
            path = folder / cds_name
            if path.is_file():
                files.append(path)
    return files


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PACKAGE_ROOT)
    parser.add_argument("--limit", type=int, help="Audit at most this many CDS files.")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)

    project_root = args.project_root.resolve()
    files = inventory_cds_files(project_root)
    if args.limit is not None:
        files = files[: args.limit]
    results = [compare_rules(path) for path in files]
    summary = {
        "project_root": str(project_root),
        "cds_files": len(results),
        "headers": sum(int(item["headers"]) for item in results),
        "explicit_protein_id": sum(int(item["explicit_protein_id"]) for item in results),
        "cds_fallback_applied": sum(int(item["cds_fallback_applied"]) for item in results),
        "id_set_changed_by_rule": sum(0 if item["id_set_identical"] else 1 for item in results),
    }
    if args.json:
        print(json.dumps({"summary": summary, "files": results}, ensure_ascii=False, indent=2))
    else:
        for key, value in summary.items():
            print(f"{key}={value}")
        for item in results:
            if not item["id_set_identical"]:
                print(f"  DIFFERS {item['file']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
