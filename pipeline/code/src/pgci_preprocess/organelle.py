"""Infer organelle sequences when an annotation has only synthetic regions."""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass

from pgci_preprocess.models import Anomaly, GffData, Transcript

ORGANELLE_NAME_RE = re.compile(
    r"(?i)(?:chloroplast|plastid|mitochondr(?:ion|ial)?|cpDNA|mtDNA)"
)
CHLOROPLAST_ID_RE = re.compile(
    r"(?i)cp(?=\d)|(?:^|[_\-.|])chrc(?=\d|$|[_\-.|])"
)
MITOCHONDRION_ID_RE = re.compile(
    r"(?i)mt(?=\d)|(?:^|[_\-.|])mito(?=\d|$|[_\-.|])"
)
PLASTID_GENE_RE = re.compile(r"(?i)(?:^|[_\-.|])pgp\d+")


@dataclass(frozen=True, slots=True)
class OrganelleEvidence:
    seqid: str
    genome_type: str
    reason: str
    length: int
    gene_count: int
    genes_per_100kb: float
    intronless_fraction: float

    def as_dict(self) -> dict[str, object]:
        return {
            "seqid": self.seqid,
            "genome_type": self.genome_type,
            "reason": self.reason,
            "length": self.length,
            "gene_count": self.gene_count,
            "genes_per_100kb": self.genes_per_100kb,
            "intronless_fraction": self.intronless_fraction,
        }


def _identifier_type(values: list[str]) -> str:
    joined = "\n".join(values)
    if MITOCHONDRION_ID_RE.search(joined) or re.search(
        r"(?i)mitochondr(?:ion|ial)?", joined
    ):
        return "mitochondrion"
    if (
        CHLOROPLAST_ID_RE.search(joined)
        or PLASTID_GENE_RE.search(joined)
        or re.search(r"(?i)(?:chloroplast|plastid)", joined)
    ):
        return "chloroplast"
    return ""


def _intronless_fraction(
    transcripts: list[Transcript],
    cds_lengths: dict[str, int],
) -> float:
    measurable = 0
    intronless = 0
    for transcript in transcripts:
        if len(transcript.protein_ids) != 1:
            continue
        cds_length = cds_lengths.get(next(iter(transcript.protein_ids)))
        if not cds_length:
            continue
        measurable += 1
        span = transcript.end - transcript.start + 1
        if abs(span - cds_length) <= 6:
            intronless += 1
    return intronless / measurable if measurable else 0.0


def infer_synthetic_organelles(
    data: GffData,
    cds_lengths: dict[str, int],
    *,
    custom_pattern: str = "",
    enabled: bool = True,
) -> dict[str, object]:
    """Mark synthetic organelle regions and list unresolved dense candidates.

    Authoritative ``region`` rows and genome FASTA records are never changed.
    For synthetic regions the inference requires either a direct organelle
    name/identifier signal or a dense, compact gene set whose gene spans mostly
    equal their CDS lengths. The latter captures organelle annotations with
    species-specific identifiers while avoiding gene-density-only guessing.
    """

    result: dict[str, object] = {
        "enabled": enabled,
        "custom_pattern": custom_pattern,
        "heuristic": "dense_compact_synthetic_region",
        "inferred": [],
        "review_candidates": [],
    }
    if not enabled:
        return result

    try:
        pattern = re.compile(custom_pattern) if custom_pattern else None
    except re.error as exc:
        raise ValueError(f"invalid organelle pattern: {exc}") from exc

    transcripts_by_seqid: dict[str, list[Transcript]] = defaultdict(list)
    for transcript in data.transcripts.values():
        transcripts_by_seqid[transcript.seqid].append(transcript)

    inferred: list[OrganelleEvidence] = []
    review: list[OrganelleEvidence] = []
    for seqid, region in sorted(data.regions.items()):
        if not region.synthetic or region.is_organelle:
            continue
        transcripts = transcripts_by_seqid.get(seqid, [])
        if not transcripts:
            continue
        identifiers = [seqid, *(item.source_id for item in transcripts)]
        hint_type = _identifier_type(identifiers)
        density = len(transcripts) * 100_000 / max(1, region.length)
        compact_fraction = _intronless_fraction(transcripts, cds_lengths)
        pattern_matched = bool(pattern and pattern.search(seqid))
        compact = (
            50 <= len(transcripts) <= 300
            and region.length <= 500_000
            and density >= 10
            and compact_fraction >= 0.70
        )

        if pattern_matched:
            genome_type = hint_type or (
                "mitochondrion"
                if re.search(r"(?i)mito|mt", custom_pattern)
                else "plastid"
            )
            reason = "custom_organelle_pattern"
        elif hint_type:
            genome_type = hint_type
            reason = "sequence_or_gene_identifier"
        elif compact:
            genome_type = "plastid"
            reason = "dense_compact_synthetic_region"
        else:
            if (
                30 <= len(transcripts) <= 300
                and region.length <= 500_000
                and density >= 15
            ):
                review.append(
                    OrganelleEvidence(
                        seqid=seqid,
                        genome_type="unknown",
                        reason="dense_synthetic_region_requires_review",
                        length=region.length,
                        gene_count=len(transcripts),
                        genes_per_100kb=round(density, 2),
                        intronless_fraction=round(compact_fraction, 4),
                    )
                )
            continue

        region.genome_type = genome_type
        region.chromosome = ""
        region.is_primary = False
        evidence = OrganelleEvidence(
            seqid=seqid,
            genome_type=genome_type,
            reason=reason,
            length=region.length,
            gene_count=len(transcripts),
            genes_per_100kb=round(density, 2),
            intronless_fraction=round(compact_fraction, 4),
        )
        inferred.append(evidence)
        data.anomalies.append(
            Anomaly(
                "warning",
                "organelle_sequence_inferred",
                f"{seqid}:{genome_type}:{reason}",
            )
        )

    result["inferred"] = [item.as_dict() for item in inferred]
    result["review_candidates"] = [item.as_dict() for item in review]
    return result
