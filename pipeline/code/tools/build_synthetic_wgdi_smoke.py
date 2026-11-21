#!/usr/bin/env python3
"""Build a self-contained synthetic two-genome workspace for release validation."""

from __future__ import annotations

import argparse
import hashlib
import json
import tomllib
from pathlib import Path

from build_wgdi_smoke_fixture import install_smoke_runtime


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
CODONS = {
    "A": ("GCT", "GCC", "GCA", "GCG"),
    "C": ("TGT", "TGC"),
    "D": ("GAT", "GAC"),
    "E": ("GAA", "GAG"),
    "F": ("TTT", "TTC"),
    "G": ("GGT", "GGC", "GGA", "GGG"),
    "H": ("CAT", "CAC"),
    "I": ("ATT", "ATC", "ATA"),
    "K": ("AAA", "AAG"),
    "L": ("CTT", "CTC", "CTA", "CTG", "TTA", "TTG"),
    "M": ("ATG",),
    "N": ("AAT", "AAC"),
    "P": ("CCT", "CCC", "CCA", "CCG"),
    "Q": ("CAA", "CAG"),
    "R": ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
    "S": ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
    "T": ("ACT", "ACC", "ACA", "ACG"),
    "V": ("GTT", "GTC", "GTA", "GTG"),
    "W": ("TGG",),
    "Y": ("TAT", "TAC"),
}


def protein(seed: str, length: int) -> str:
    values = bytearray()
    counter = 0
    while len(values) < length:
        values.extend(
            hashlib.sha256(f"{seed}:{counter}".encode("utf-8")).digest()
        )
        counter += 1
    return "".join(AMINO_ACIDS[value % len(AMINO_ACIDS)] for value in values[:length])


def mutate(sequence: str, seed: str) -> str:
    position = int.from_bytes(hashlib.sha256(seed.encode("utf-8")).digest()[:2], "big") % len(sequence)
    current = sequence[position]
    alternatives = [amino for amino in AMINO_ACIDS if amino != current]
    replacement = alternatives[position % len(alternatives)]
    return sequence[:position] + replacement + sequence[position + 1 :]


def coding_sequence(
    sequence: str,
    seed: str,
    variations: tuple[tuple[str, int], ...] = (),
) -> str:
    codons: list[str] = []
    for index, amino in enumerate(sequence):
        choices = CODONS[amino]
        digest = hashlib.sha256(f"{seed}:{index}".encode("utf-8")).digest()
        codon = choices[digest[0] % len(choices)]
        for variation_seed, step in variations:
            if step > 0 and index % step == 0 and len(choices) > 1:
                alternatives = [choice for choice in choices if choice != codon]
                variation_digest = hashlib.sha256(
                    f"{variation_seed}:{index}".encode("utf-8")
                ).digest()
                codon = alternatives[variation_digest[0] % len(alternatives)]
        codons.append(codon)
    return "".join(codons)


def write_species(
    root: Path,
    code: str,
    proteins: list[str],
    chromosome_count: int,
    codon_specs: list[tuple[str, tuple[tuple[str, int], ...]]],
) -> None:
    species_dir = root / code
    species_dir.mkdir()
    genes_per_chromosome = (len(proteins) + chromosome_count - 1) // chromosome_count
    gene_rows: list[list[str]] = []
    counts = [0] * chromosome_count
    pep_lines: list[str] = []
    cds_lines: list[str] = []
    chromosome_length = 50_000

    for index, (sequence, (codon_seed, variations)) in enumerate(
        zip(proteins, codon_specs, strict=True)
    ):
        chromosome = index // genes_per_chromosome
        order = counts[chromosome] + 1
        counts[chromosome] = order
        gene_id = f"{code}{chromosome + 1}g{order:04d}"
        start = 1_000 + (order - 1) * 400
        end = start + len(sequence) * 3 - 1
        gene_rows.append([str(chromosome + 1), gene_id, str(start), str(end), "+", str(order), gene_id])
        pep_lines.append(f">{gene_id}\n{sequence}\n")
        cds = coding_sequence(sequence, codon_seed, variations)
        cds_lines.append(f">{gene_id}\n{cds}\n")

    (species_dir / f"{code}.new.gff").write_text(
        "".join("\t".join(row) + "\n" for row in gene_rows),
        encoding="utf-8",
    )
    (species_dir / f"{code}.pep").write_text("".join(pep_lines), encoding="utf-8")
    (species_dir / f"{code}.cds").write_text("".join(cds_lines), encoding="utf-8")
    (species_dir / f"{code}.lens").write_text(
        "".join(
            f"{chromosome + 1}\t{chromosome_length}\t{counts[chromosome]}\n"
            for chromosome in range(chromosome_count)
        ),
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workflow-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--genes-per-species", type=int, default=120)
    parser.add_argument("--chromosomes", type=int, default=2)
    parser.add_argument("--protein-length", type=int, default=60)
    args = parser.parse_args()
    if args.genes_per_species < 20:
        parser.error("--genes-per-species must be >= 20")
    if args.genes_per_species % 2:
        parser.error("--genes-per-species must be even so self comparisons have paralogues")
    if args.chromosomes < 1 or args.protein_length < 20:
        parser.error("--chromosomes must be >= 1 and --protein-length >= 20")

    workflow = args.workflow_root.resolve()
    output = args.output_root.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing smoke workspace: {output}")
    output.mkdir(parents=True)

    family_count = args.genes_per_species // 2
    family_proteins = [
        protein(f"synthetic-{index}", args.protein_length)
        for index in range(family_count)
    ]
    base_proteins = family_proteins + [
        mutate(sequence, f"synthetic-paralogue-{index}")
        for index, sequence in enumerate(family_proteins)
    ]
    variant_proteins = [
        mutate(sequence, f"synthetic-variant-{index}")
        for index, sequence in enumerate(base_proteins)
    ]
    family_seeds = [f"synthetic-codons-{index}" for index in range(family_count)]
    syna_specs = [
        (seed, ()) for seed in family_seeds
    ] + [
        (seed, (("synthetic-paralogue-codons", 8),)) for seed in family_seeds
    ]
    synb_specs = [
        (seed, (("synthetic-species-codons", 12),)) for seed in family_seeds
    ] + [
        (
            seed,
            (
                ("synthetic-paralogue-codons", 8),
                ("synthetic-species-codons", 12),
            ),
        )
        for seed in family_seeds
    ]
    write_species(output, "SynA", base_proteins, args.chromosomes, syna_specs)
    write_species(output, "SynB", variant_proteins, args.chromosomes, synb_specs)
    install_smoke_runtime(workflow, output)

    manifest = {
        "status": "FIXTURE_READY_NOT_WGDI_RUN",
        "workflow_version": tomllib.loads(
            (workflow / "code" / "pyproject.toml").read_text(encoding="utf-8")
        )["project"]["version"],
        "selection": {
            "method": "deterministic synthetic orthologues",
            "genes_per_species": args.genes_per_species,
            "chromosomes": args.chromosomes,
            "protein_length": args.protein_length,
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
