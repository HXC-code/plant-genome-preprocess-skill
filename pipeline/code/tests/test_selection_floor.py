from __future__ import annotations

from collections import Counter

from pgci_preprocess.models import Region
from pgci_preprocess.selection import select_analysis_sequences


def _region(seqid: str, *, primary: bool = False) -> Region:
    return Region(
        seqid=seqid,
        length=1_000_000,
        genome_type="genomic",
        chromosome="1" if primary else "",
        taxon_id="",
        source_order=0,
        is_primary=primary,
        synthetic=False,
    )


def test_floor_drops_small_chromosomes_and_recomputes() -> None:
    regions = {
        "chr1": _region("chr1", primary=True),
        "chr2": _region("chr2", primary=True),
        "scaf": _region("scaf", primary=True),
    }
    counts = Counter({"chr1": 500, "chr2": 300, "scaf": 50})
    sel = select_analysis_sequences(
        regions, counts, 50, 0.9, 200, min_genes_per_chromosome=100
    )
    assert sel.keep == {"chr1", "chr2"}
    assert sel.below_floor == {"scaf"}
    assert sel.selected_genes == 800
    assert sel.mode == "primary_chromosome"
    assert sel.ready == "yes"


def test_floor_zero_keeps_everything() -> None:
    regions = {
        "chr1": _region("chr1", primary=True),
        "scaf": _region("scaf", primary=True),
    }
    counts = Counter({"chr1": 500, "scaf": 50})
    sel = select_analysis_sequences(
        regions, counts, 50, 0.9, 200, min_genes_per_chromosome=0
    )
    assert sel.keep == {"chr1", "scaf"}
    assert sel.below_floor == set()


def test_floor_can_downgrade_readiness() -> None:
    regions = {
        "chr1": _region("chr1", primary=True),
        "chr2": _region("chr2", primary=True),
    }
    counts = Counter({"chr1": 100, "chr2": 90})
    sel = select_analysis_sequences(
        regions, counts, 50, 0.9, 200, min_genes_per_chromosome=100
    )
    assert sel.keep == {"chr1"}
    assert sel.below_floor == {"chr2"}
    # 100 / 190 ~= 0.526 -> marginal, not yes
    assert sel.ready == "marginal"


def test_floor_applies_in_coverage_mode() -> None:
    regions = {f"s{i}": _region(f"s{i}") for i in range(1, 5)}
    counts = Counter({"s1": 400, "s2": 300, "s3": 200, "s4": 60})
    sel = select_analysis_sequences(
        regions, counts, 50, 0.9, 200, min_genes_per_chromosome=100
    )
    assert "s4" in sel.below_floor
    assert "s4" not in sel.keep