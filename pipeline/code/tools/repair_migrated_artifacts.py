#!/usr/bin/env python3
"""Repair artifacts damaged by the interrupted 2026-11-20 bulk migration.

The faulty migration replaced ``<seqid>\t<legacy_gene_id>\t`` with
``<seqid><code><gene_id>`` (one tab lost).  This tool rebuilds every GFF row and
the ``gene_map`` first column deterministically:

* ``.new.gff`` is restored from the untouched legacy copy kept at
  ``<workspace>/s<show_id>/s<show_id>.new.gff`` and the gene IDs are rewritten;
* ``.all.gff`` has no legacy copy, so the merged cell
  ``<seq>\t<code><seq>g<order><start>`` is split back on the fixed-width order.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import uuid
from pathlib import Path

INDEX = Path("resources/wgdi_v2_index.tsv")


def atomic_write(path: Path, payload: bytes) -> None:
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        temporary.write_bytes(payload)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def rewrite_legacy_gff(payload: bytes, show_id: str, code: str) -> bytes:
    pattern = re.compile(rb"(?m)^([^\t]*\t)s" + str(show_id).encode() + rb"(\d+g\d{4,})\t")
    replaced, count = pattern.subn(rb"\1" + code.encode() + rb"\2\t", payload)
    if pattern.search(replaced):
        raise RuntimeError("legacy gene IDs survived the rewrite")
    if not count:
        raise RuntimeError("no legacy gene IDs found in the source GFF")
    return replaced


def split_merged_gff(payload: bytes, code: str) -> tuple[bytes, int]:
    merged = re.compile(
        rb"(?m)^([^\t\n]*)\t" + re.escape(code.encode()) + rb"(\d+)g(\d{4})(\d+)\t(\d+)\t"
    )

    def rebuild(match: re.Match[bytes]) -> bytes:
        seq, gene_seq, order, start, end = match.groups()
        gene = code.encode() + gene_seq + b"g" + order
        return b"\t".join([seq, gene, start, end]) + b"\t"

    replaced, count = merged.subn(rebuild, payload)
    return replaced, count


def split_merged_gene_id(merged: str, code: str, order: str) -> tuple[str, str] | None:
    """Split ``code<seq>g<order><seq>`` (optionally prefixed by ``>``)."""

    body = merged[1:] if merged.startswith(">") else merged
    if not body.startswith(code):
        return None
    for cut in range(len(code) + 3, len(body)):
        head, tail = body[:cut], body[cut:]
        match = re.fullmatch(rf"{re.escape(code)}(\d+)g(\d{{4,}})", head)
        if not match or not tail.isdigit():
            continue
        sequence = match.group(1)
        if int(match.group(2)) == int(order) and tail == sequence:
            return head, sequence
    return None


def repair_gene_map(path: Path, code: str) -> int:
    text = path.read_text(encoding="utf-8").splitlines()
    header, rows = text[0], text[1:]
    fixed = 0
    output = [header]
    for row in rows:
        if not row:
            continue
        fields = row.split("\t")
        if len(fields) == 16:
            split = split_merged_gene_id(fields[0], code, fields[1])
            if split is None:
                raise RuntimeError(f"{path}: cannot split {fields[0]!r} order={fields[1]!r}")
            gene_id, sequence = split
            fields = [gene_id, sequence, fields[1]] + fields[2:]
            fixed += 1
        elif len(fields) != 17:
            raise RuntimeError(f"{path}: unexpected column count {len(fields)}")
        output.append("\t".join(fields))
    atomic_write(path, ("\n".join(output) + "\n").encode("utf-8"))
    return fixed


def count_gff_rows(path: Path) -> int:
    with path.open(encoding="utf-8-sig") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    root = args.project_root.resolve()
    if (root / "plant-genome-resources" / INDEX).is_file():
        workspace, package = root, root / "plant-genome-resources"
    else:
        workspace, package = root, root
    rows = list(csv.DictReader((package / INDEX).open(encoding="utf-8"), delimiter="\t"))
    analysis = package / "analysis_v2"
    report = {"species": 0, "new_gff": 0, "all_gff": 0, "gene_map_rows": 0, "problems": []}
    for row in rows:
        code, show_id = row["code"], row["show_id"]
        bundle = analysis / code
        if not code or not bundle.is_dir():
            continue
        report["species"] += 1
        legacy = workspace / f"s{show_id}" / f"s{show_id}.new.gff"
        target = bundle / f"{code}.new.gff"
        if not legacy.is_file():
            report["problems"].append(f"{code}: missing legacy copy {legacy}")
            continue
        payload = rewrite_legacy_gff(legacy.read_bytes(), show_id, code)
        if not args.dry_run:
            atomic_write(target, payload)
        report["new_gff"] += 1
        all_gff = bundle / f"{code}.all.gff"
        if all_gff.is_file():
            payload = all_gff.read_bytes()
            fixed, count = split_merged_gff(payload, code)
            damaged = any(
                len(line.rstrip(b"\r\n").split(b"\t")) != 7
                for line in payload.splitlines()
                if line.strip() and not line.startswith(b"#")
            )
            if damaged and not count:
                report["problems"].append(f"{code}: all.gff rows still malformed")
            elif count and not args.dry_run:
                atomic_write(all_gff, fixed)
            if count:
                report["all_gff"] += 1
        for layer in (package / "analysis_v2", package / "clean_v2"):
            gene_map = layer / code / f"{code}.gene_map.tsv"
            if not gene_map.is_file():
                continue
            merged_rows = sum(
                1
                for line in gene_map.read_text(encoding="utf-8").splitlines()[1:]
                if line and len(line.split("\t")) == 16
            )
            if merged_rows:
                if not args.dry_run:
                    repair_gene_map(gene_map, code)
                report["gene_map_rows"] += merged_rows
        for name, expected in (
            (f"{code}.new.gff", bundle / f"{code}.lens"),
            (f"{code}.all.gff", bundle / f"{code}.all.lens"),
        ):
            gff = bundle / name
            if not gff.is_file() or not expected.is_file():
                continue
            observed = count_gff_rows(gff)
            wanted = sum(
                int(line.split("\t")[2])
                for line in expected.read_text(encoding="utf-8").splitlines()
                if line.strip()
            )
            if observed != wanted:
                report["problems"].append(f"{code}: {name} rows {observed} != lens {wanted}")
    print(
        f"species={report['species']} new_gff={report['new_gff']} "
        f"all_gff={report['all_gff']} gene_map_rows={report['gene_map_rows']} "
        f"problems={len(report['problems'])}"
    )
    for problem in report["problems"][:20]:
        print(f"  {problem}")
    return 1 if report["problems"] else 0


if __name__ == "__main__":
    raise SystemExit(main())

