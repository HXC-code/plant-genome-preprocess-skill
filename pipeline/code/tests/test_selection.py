from __future__ import annotations

from collections import Counter

from pgci_preprocess.models import Candidate, Region
from pgci_preprocess.selection import choose_representative, select_analysis_sequences


def _candidate(transcript: str, protein_length: int, cds_length: int) -> Candidate:
    return Candidate(
        source_gene_id="G",
        source_gene_name="G",
        gene_biotype="protein_coding",
        source_transcript_id=transcript,
        source_protein_id=transcript + "-P",
        seqid="chr1",
        start=1,
        end=100,
        strand="+",
        cds_length=cds_length,
        protein_length=protein_length,
        translation_exact=True,
        translation_mismatches=0,
    )


def test_representative_prefers_longest_cds_not_protein() -> None:
    # T3 has the longest protein, T2 the longest CDS: the SOP rule must pick
    # T2. Protein length is never a selection criterion.
    chosen = choose_representative(
        [
            _candidate("T3", 30, 33),
            _candidate("T1", 20, 60),
            _candidate("T2", 19, 63),
        ]
    )
    assert chosen.source_transcript_id == "T2"


def test_representative_tie_breaks_by_stable_id_not_input_order() -> None:
    forward = choose_representative(
        [_candidate("T2", 20, 60), _candidate("T1", 20, 60)]
    )
    reverse = choose_representative(
        [_candidate("T1", 20, 60), _candidate("T2", 20, 60)]
    )
    assert forward.source_transcript_id == "T1"
    assert reverse.source_transcript_id == "T1"


def _region(seqid: str, *, primary: bool, organelle: bool = False) -> Region:
    return Region(
        seqid=seqid,
        length=1000,
        genome_type="mitochondrion" if organelle else ("chromosome" if primary else "genomic"),
        chromosome="1" if primary else "",
        taxon_id="1",
        source_order=0,
        is_primary=primary,
    )


def test_sequences_use_primary_chromosomes_and_ignore_scaffolds() -> None:
    regions = {
        "chr1": _region("chr1", primary=True),
        "scaffold": _region("scaffold", primary=False),
        "mt": _region("mt", primary=False, organelle=True),
    }
    selection = select_analysis_sequences(
        regions,
        Counter({"chr1": 100, "scaffold": 1000, "mt": 20}),
        min_genes=50,
        target_coverage=0.9,
        max_sequences=200,
    )
    assert selection.keep == {"chr1"}
    assert selection.mode == "primary_chromosome"
    assert selection.ratio == 0.0909
    assert selection.ready == "no"


def test_chromosomes_mode_never_substitutes_scaffolds() -> None:
    regions = {
        "chr1": _region("chr1", primary=True),
        "scaffold": _region("scaffold", primary=False),
        "mt": _region("mt", primary=False, organelle=True),
    }
    selection = select_analysis_sequences(
        regions,
        Counter({"chr1": 10, "scaffold": 1000, "mt": 20}),
        min_genes=50,
        target_coverage=0.9,
        max_sequences=200,
        sequence_mode="chromosomes",
    )
    assert selection.keep == {"chr1"}
    assert selection.mode == "chromosomes_only"


def test_chromosomes_mode_fails_closed_on_scaffold_only_assembly() -> None:
    regions = {"scaffold": _region("scaffold", primary=False)}
    selection = select_analysis_sequences(
        regions,
        Counter({"scaffold": 1000}),
        min_genes=50,
        target_coverage=0.9,
        max_sequences=200,
        sequence_mode="chromosomes",
    )
    assert selection.keep == set()
    assert selection.ready == "no"


def test_chromosomes_mode_applies_gene_floor() -> None:
    regions = {
        "chr1": _region("chr1", primary=True),
        "chr2": _region("chr2", primary=True),
    }
    selection = select_analysis_sequences(
        regions,
        Counter({"chr1": 200, "chr2": 5}),
        min_genes=50,
        target_coverage=0.9,
        max_sequences=200,
        min_genes_per_chromosome=100,
        sequence_mode="chromosomes",
    )
    assert selection.keep == {"chr1"}
    assert selection.below_floor == {"chr2"}


def test_unknown_sequence_mode_is_rejected() -> None:
    import pytest

    with pytest.raises(ValueError):
        select_analysis_sequences(
            {},
            Counter(),
            min_genes=50,
            target_coverage=0.9,
            max_sequences=200,
            sequence_mode="bogus",
        )


def test_fragmented_assembly_falls_back_to_gene_coverage() -> None:
    regions = {
        "s1": _region("s1", primary=False),
        "s2": _region("s2", primary=False),
        "s3": _region("s3", primary=False),
        "mt": _region("mt", primary=False, organelle=True),
    }
    selection = select_analysis_sequences(
        regions,
        Counter({"s1": 45, "s2": 35, "s3": 20, "mt": 10}),
        min_genes=50,
        target_coverage=0.9,
        max_sequences=200,
    )
    assert selection.keep == {"s1", "s2", "s3"}
    assert selection.mode == "coverage"
    assert selection.ratio == 1.0
    assert selection.ready == "yes"

