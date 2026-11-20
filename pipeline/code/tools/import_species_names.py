#!/usr/bin/env python3
"""Import the authoritative Excel abbreviation table into canonical TSV."""

from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
from pathlib import Path


def import_xlsx(source: Path, destination: Path) -> int:
    try:
        import openpyxl
    except ImportError as exc:
        raise RuntimeError(
            "Excel import requires openpyxl; run this tool with an environment "
            "that provides the spreadsheet runtime"
        ) from exc
    workbook = openpyxl.load_workbook(source, read_only=True, data_only=True)
    worksheet = workbook.active
    rows: list[dict[str, str]] = []
    for values in worksheet.iter_rows(min_row=3, values_only=True):
        show_id = str(values[0] or "").strip()
        if not show_id:
            continue
        automatic = str(values[4] or "").strip()
        manual = str(values[5] or "").strip()
        abbreviation = manual or automatic
        if not abbreviation:
            raise RuntimeError(f"show_id={show_id}: abbreviation is empty")
        rows.append(
            {
                "show_id": show_id,
                "latin_name": str(values[1] or "").strip(),
                "order": str(values[2] or "").strip(),
                "family": str(values[3] or "").strip(),
                "abbreviation": abbreviation,
                "auto_abbreviation": automatic,
                "manual_override": manual,
            }
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    show_ids = {row["show_id"] for row in rows}
    if len(show_ids) != len(rows):
        raise RuntimeError("show_id is not unique in the abbreviation table")
    abbreviations: dict[str, str] = {}
    for row in rows:
        previous = abbreviations.setdefault(row["abbreviation"], row["show_id"])
        if previous != row["show_id"]:
            raise RuntimeError(
                f"abbreviation {row['abbreviation']} is shared by "
                f"show_id={previous} and show_id={row['show_id']}"
            )
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        dir=destination.parent,
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "show_id",
                    "latin_name",
                    "order",
                    "family",
                    "abbreviation",
                    "auto_abbreviation",
                    "manual_override",
                ],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_name, destination)
    except Exception:
        Path(temporary_name).unlink(missing_ok=True)
        raise
    return len(rows)


def _find_default_source() -> Path | None:
    here = Path(__file__).resolve()
    for parent in list(here.parents)[:4]:
        for name in ("物种缩写表.xlsx", "species_abbreviations.xlsx"):
            candidate = parent / name
            if candidate.is_file():
                return candidate
    return None


def main(argv: list[str] | None = None) -> int:
    project_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="authoritative 物种缩写表.xlsx (default: searched upwards from the package)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=project_root / "resources" / "species_abbreviations.tsv",
    )
    args = parser.parse_args(argv)
    source = args.input or _find_default_source()
    if source is None:
        print(
            "ERROR: 物种缩写表.xlsx not found; pass --input <path to xlsx>",
            file=sys.stderr,
        )
        return 1
    try:
        count = import_xlsx(source.resolve(), args.output.resolve())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"[OK] source={source} output={args.output} rows={count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
