from __future__ import annotations

import re
from collections import Counter

from pgci_preprocess.models import Candidate, Region, SequenceSelection

DIGIT_RE = re.compile(r"(\d+)")


def natural_key(value: str) -> tuple[object, ...]:
    return tuple(int(part) if part.isdigit() else part.casefold() for part in DIGIT_RE.split(value))


def region_sort_key(region: Region) -> tuple[object, ...]:
    chromosome = region.chromosome
    if chromosome.casefold() not in {"", ".", "unknown", "na", "n/a"}:
        return (0, *natural_key(chromosome), *natural_key(region.seqid))
    return (1, *natural_key(region.seqid))


def choose_representative(candidates: list[Candidate]) -> Candidate:
    """Pick one representative transcript per gene (lab SOP rule).

    Priority 1: longest CDS.  Priority 2: stable transcript ID, ascending
    lexicographic order.  Protein length and file/order of appearance are
    never used, so the result is reproducible across runs and platforms.
    """
    if not candidates:
        raise ValueError("candidates must not be empty")
    return min(
        candidates,
        key=lambda candidate: (
            -candidate.cds_length,
            candidate.source_transcript_id,
            candidate.source_protein_id,
        ),
    )


def select_analysis_sequences(
    regions: dict[str, Region],
    gene_counts: Counter[str],
    min_genes: int,
    target_coverage: float,
    max_sequences: int,
    min_genes_per_chromosome: int = 0,
    sequence_mode: str = "adaptive",
) -> SequenceSelection:
    if sequence_mode not in {"adaptive", "chromosomes"}:
        raise ValueError(f"unknown sequence_mode: {sequence_mode}")
    total = sum(
        count
        for seqid, count in gene_counts.items()
        if seqid in regions and not regions[seqid].is_organelle
    )
    if total == 0:
        return SequenceSelection(set(), "empty", 0, 0, 0.0, "no")

    primary = {
        seqid
        for seqid, count in gene_counts.items()
        if count > 0 and seqid in regions and regions[seqid].is_primary and not regions[seqid].is_organelle
    }
    if sequence_mode == "chromosomes":
        # SOP strict mode: only chromosome-mounted sequences are eligible.
        # Unmounted scaffolds/contigs are never substituted; a species without
        # a primary chromosome yields an empty set and fails closed upstream.
        keep = set(primary)
        mode = "chromosomes_only"
    elif primary:
        keep = set(primary)
        mode = "primary_chromosome"
    else:
        eligible = {
            seqid: count
            for seqid, count in gene_counts.items()
            if count > 0 and seqid in regions and not regions[seqid].is_organelle
        }
        threshold_keep = {
            seqid for seqid, count in eligible.items() if count >= min_genes
        }
        threshold_genes = sum(eligible[seqid] for seqid in threshold_keep)
        if (
            threshold_keep
            and len(threshold_keep) <= max_sequences
            and threshold_genes / total >= 0.5
        ):
            keep = set(threshold_keep)
            mode = "threshold"
        else:
            ordered: list[str] = []
            selected = 0
            for seqid, count in sorted(
                eligible.items(),
                key=lambda item: (-item[1], *natural_key(item[0])),
            ):
                if len(ordered) >= max_sequences:
                    break
                ordered.append(seqid)
                selected += count
                if selected / total >= target_coverage:
                    break
            keep = set(ordered)
            mode = "coverage"

    # Final chromosome gene-count floor: drop any kept sequence carrying fewer
    # than min_genes_per_chromosome representative genes (e.g. tiny scaffolds or
    # fragments). Applied across all modes; recorded for auditable exclusion.
    below_floor: set[str] = set()
    if min_genes_per_chromosome > 0:
        below_floor = {
            seqid
            for seqid in keep
            if gene_counts.get(seqid, 0) < min_genes_per_chromosome
        }
        keep = keep - below_floor

    selected_genes = sum(gene_counts.get(seqid, 0) for seqid in keep)
    ratio = selected_genes / total if total else 0.0
    if ratio >= 0.7 and len(keep) <= max_sequences:
        ready = "yes"
    elif ratio >= 0.5:
        ready = "marginal"
    else:
        ready = "no"
    return SequenceSelection(
        keep,
        mode,
        selected_genes,
        total,
        round(ratio, 4),
        ready,
        below_floor,
    )

