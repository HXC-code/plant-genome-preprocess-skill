#!/usr/bin/env python3
"""Build an isolated WGDI fixture from a real preprocessing output subset."""

from __future__ import annotations

import argparse
import json
import tomllib
from pathlib import Path

from build_wgdi_smoke_fixture import install_smoke_runtime


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current = ""
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if current:
                records[current] = "".join(chunks)
            current = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if current:
        records[current] = "".join(chunks)
    return records


def write_fasta(path: Path, ids: list[str], records: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for gene_id in ids:
            handle.write(f">{gene_id}\n{records[gene_id]}\n")


def load_bundle(bundle: Path, chromosome: str, max_genes: int) -> dict[str, object]:
    code = bundle.name
    gff_path = bundle / f"{code}.new.gff"
    lens_path = bundle / f"{code}.lens"
    rows: list[list[str]] = []
    for line in gff_path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) == 7 and fields[0] == chromosome:
            rows.append(fields)
    rows.sort(key=lambda fields: int(fields[5]))
    rows = rows[:max_genes]
    if len(rows) < 5:
        raise ValueError(f"{code}: chromosome {chromosome} has fewer than five genes")

    chromosome_lengths: dict[str, int] = {}
    for line in lens_path.read_text(encoding="utf-8").splitlines():
        sequence, length, _count = line.split("\t")
        chromosome_lengths[sequence] = int(length)
    if chromosome not in chromosome_lengths:
        raise ValueError(f"{code}: chromosome {chromosome} is absent from lens")

    gene_ids = [row[1] for row in rows]
    pep = read_fasta(bundle / f"{code}.pep")
    cds = read_fasta(bundle / f"{code}.cds")
    missing = [gene for gene in gene_ids if gene not in pep or gene not in cds]
    if missing:
        raise ValueError(f"{code}: missing sequences for {missing[:5]}")
    return {
        "code": code,
        "chromosome": chromosome,
        "chromosome_length": chromosome_lengths[chromosome],
        "rows": rows,
        "gene_ids": gene_ids,
        "pep": pep,
        "cds": cds,
    }


def write_bundle(output: Path, bundle: dict[str, object]) -> None:
    code = str(bundle["code"])
    target = output / code
    target.mkdir(parents=True)
    rows = bundle["rows"]
    gene_ids = bundle["gene_ids"]
    assert isinstance(rows, list)
    assert isinstance(gene_ids, list)
    (target / f"{code}.new.gff").write_text(
        "".join("\t".join(row) + "\n" for row in rows),
        encoding="utf-8",
        newline="\n",
    )
    (target / f"{code}.lens").write_text(
        f"{bundle['chromosome']}\t{bundle['chromosome_length']}\t{len(rows)}\n",
        encoding="utf-8",
        newline="\n",
    )
    write_fasta(target / f"{code}.pep", gene_ids, bundle["pep"])  # type: ignore[arg-type]
    write_fasta(target / f"{code}.cds", gene_ids, bundle["cds"])  # type: ignore[arg-type]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-root", type=Path, required=True)
    parser.add_argument("--left-bundle", type=Path, required=True)
    parser.add_argument("--left-chromosome", required=True)
    parser.add_argument("--right-bundle", type=Path, required=True)
    parser.add_argument("--right-chromosome", required=True)
    parser.add_argument("--max-genes", type=int, default=1000)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    if args.max_genes < 5:
        parser.error("--max-genes must be >= 5")

    workflow = args.workflow_root.resolve()
    output = args.output_root.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing fixture: {output}")
    output.mkdir(parents=True)
    left = load_bundle(args.left_bundle.resolve(), args.left_chromosome, args.max_genes)
    right = load_bundle(args.right_bundle.resolve(), args.right_chromosome, args.max_genes)
    write_bundle(output, left)
    write_bundle(output, right)
    install_smoke_runtime(workflow, output)
    manifest = {
        "status": "REAL_FIXTURE_READY_NOT_WGDI_RUN",
        "workflow_version": tomllib.loads(
            (workflow / "code" / "pyproject.toml").read_text(encoding="utf-8")
        )["project"]["version"],
        "selection": {
            "method": "real preprocessing output subset",
            "left_code": left["code"],
            "left_chromosome": left["chromosome"],
            "right_code": right["code"],
            "right_chromosome": right["chromosome"],
            "max_genes": args.max_genes,
        },
    }
    (output / "fixture_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": "READY", "root": str(output), "manifest": manifest}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
