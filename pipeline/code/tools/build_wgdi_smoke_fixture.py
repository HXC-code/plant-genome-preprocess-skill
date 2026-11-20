#!/usr/bin/env python3
"""Build an isolated, real-data WGDI smoke workspace from existing collinear blocks.

The fixture is derived from a previously validated bundle. It never changes the
source bundle or previous BLAST/WGDI results. DIAMOND must be rerun in the new
workspace; historical alignments are used only to choose a compact gene subset.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from collections import defaultdict
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_gff(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 7:
            raise ValueError(f"{path}:{line_number}: expected seven GFF columns")
        rows.append(fields)
    if not rows:
        raise ValueError(f"empty GFF: {path}")
    return rows


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current = ""
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    records[current] = "".join(chunks)
                current = line[1:].split()[0]
                if not current or current in records:
                    raise ValueError(f"empty or duplicate FASTA ID in {path}: {current}")
                chunks = []
            else:
                if not current:
                    raise ValueError(f"sequence before FASTA header in {path}")
                chunks.append(line)
    if current:
        records[current] = "".join(chunks)
    if not records or any(not sequence for sequence in records.values()):
        raise ValueError(f"empty FASTA or empty sequence: {path}")
    return records


def read_alignments(path: Path) -> list[dict[str, object]]:
    groups: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("# Alignment"):
            current = {"header": line, "pairs": []}
            groups.append(current)
        elif line.strip() and current is not None:
            fields = line.split()
            if len(fields) >= 4:
                current["pairs"].append((fields[0], fields[2]))  # type: ignore[index]
    return [group for group in groups if len(group["pairs"]) >= 4]  # type: ignore[arg-type]


def write_fasta(path: Path, ids: set[str], records: dict[str, str], order: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for gene_id in order:
            if gene_id in ids:
                handle.write(f">{gene_id}\n{records[gene_id]}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-root", type=Path, required=True)
    parser.add_argument("--species-dir", type=Path, required=True)
    parser.add_argument("--self-icl", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--max-genes", type=int, default=1600)
    parser.add_argument("--max-alignments", type=int, default=100)
    parser.add_argument("--neighbors", type=int, default=2)
    args = parser.parse_args()
    if args.max_genes < 100 or args.max_alignments < 1 or args.neighbors < 0:
        parser.error("max-genes must be >=100; max-alignments >=1; neighbors >=0")

    workflow = args.workflow_root.resolve()
    species_dir = args.species_dir.resolve()
    output = args.output_root.resolve()
    code = species_dir.name
    source_paths = {
        suffix: species_dir / f"{code}.{suffix}" for suffix in ("new.gff", "lens", "cds", "pep")
    }
    for path in [*source_paths.values(), args.self_icl]:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing smoke workspace: {output}")

    rows = read_gff(source_paths["new.gff"])
    by_id = {row[1]: row for row in rows}
    if len(by_id) != len(rows):
        raise ValueError("duplicate gene IDs in source GFF")
    by_sequence: dict[str, list[list[str]]] = defaultdict(list)
    for row in rows:
        by_sequence[row[0]].append(row)
    for sequence_rows in by_sequence.values():
        sequence_rows.sort(key=lambda row: int(row[5]))

    selected: set[str] = set()
    selected_headers: list[str] = []
    for alignment in read_alignments(args.self_icl):
        if len(selected_headers) >= args.max_alignments:
            break
        proposed: set[str] = set()
        for left, right in alignment["pairs"]:  # type: ignore[index]
            for gene_id in (left, right):
                row = by_id.get(gene_id)
                if row is None:
                    continue
                sequence_rows = by_sequence[row[0]]
                center = int(row[5])
                for neighbor in sequence_rows:
                    if abs(int(neighbor[5]) - center) <= args.neighbors:
                        proposed.add(neighbor[1])
        if selected and len(selected | proposed) > args.max_genes:
            continue
        selected.update(proposed)
        selected_headers.append(str(alignment["header"]))
        if len(selected) >= args.max_genes * 0.90:
            break
    if len(selected) < 100 or not selected_headers:
        raise ValueError(f"subset too small: {len(selected)} genes, {len(selected_headers)} alignments")

    selected_rows = [row[:] for row in rows if row[1] in selected]
    selected_rows.sort(key=lambda row: (int(row[0]), int(row[2]), int(row[5])))
    counts: dict[str, int] = defaultdict(int)
    for row in selected_rows:
        counts[row[0]] += 1
        row[5] = str(counts[row[0]])
    selected_order = [row[1] for row in selected_rows]
    cds = read_fasta(source_paths["cds"])
    pep = read_fasta(source_paths["pep"])
    if set(cds) != set(by_id) or set(pep) != set(by_id):
        raise ValueError("source GFF/CDS/PEP ID sets differ")

    lens_lengths: dict[str, str] = {}
    for line in source_paths["lens"].read_text(encoding="utf-8").splitlines():
        if line.strip():
            fields = line.split("\t")
            if len(fields) != 3:
                raise ValueError("source lens must have three columns")
            lens_lengths[fields[0]] = fields[1]
    if set(lens_lengths) != set(by_sequence):
        raise ValueError("source lens/GFF sequence sets differ")

    target_species = output / code
    target_species.mkdir(parents=True)
    with (target_species / f"{code}.new.gff").open("w", encoding="utf-8", newline="\n") as handle:
        for row in selected_rows:
            handle.write("\t".join(row) + "\n")
    with (target_species / f"{code}.lens").open("w", encoding="utf-8", newline="\n") as handle:
        for sequence_id in sorted(counts, key=int):
            handle.write(f"{sequence_id}\t{lens_lengths[sequence_id]}\t{counts[sequence_id]}\n")
    write_fasta(target_species / f"{code}.cds", selected, cds, selected_order)
    write_fasta(target_species / f"{code}.pep", selected, pep, selected_order)

    for name in (
        "0check_data.sh", "1run_blast.sh", "1run_diamond.sh", "2run_pre_wgdi.sh",
        "3blast_to_wgdi.sh", "4run_wgdi_different.sh", "5run_wgdi_same.sh",
        "6run_wgdi.sh", "run_all_pipeline.ps1",
    ):
        shutil.copy2(workflow / name, output / name)
    shutil.copytree(workflow / "bin", output / "bin")

    manifest = {
        "status": "FIXTURE_READY_NOT_WGDI_RUN",
        "species": code,
        "source_species_dir": str(species_dir),
        "source_self_icl": str(args.self_icl.resolve()),
        "source_sha256": {path.name: sha256(path) for path in [*source_paths.values(), args.self_icl]},
        "selection": {
            "method": "historical self-collinear anchors plus neighbors; DIAMOND must be rerun",
            "max_genes": args.max_genes,
            "neighbors": args.neighbors,
            "selected_genes": len(selected),
            "selected_sequences": len(counts),
            "selected_alignments": len(selected_headers),
            "alignment_headers": selected_headers,
        },
        "outputs": {path.name: sha256(path) for path in sorted(target_species.iterdir())},
        "warning": "Derived smoke fixture only; not a biological reanalysis of the full genome.",
    }
    (output / "fixture_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"output_root": str(output), **manifest["selection"]}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
