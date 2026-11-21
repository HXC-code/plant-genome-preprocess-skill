"""Infer organelle sequences when an annotation has only synthetic regions.

Two complementary signals are used, both restricted to *synthetic* regions
(transcript-only / simplified adapters). Authoritative ``region`` rows and
genome FASTA records are never reclassified.

* Identifier signal: the sequence id or a gene name carries an explicit
  organelle token (chloroplast/plastid/mitochondrion, cp/chrc/mt/chrM, plastid
  ``pgp`` gene ids).
* Gene-content / shape signal:
  - Chloroplast-like: a small, dense, mostly intronless gene set (cpDNA is
    compact, ~120-160 kb, ~80 genes, few introns).
  - Mitochondrion-like: a cluster of distinct plant mitochondrial core genes
    (nad/cob/cox/atp/matR/mttB/ccm). Plant mtDNA is large, sparse and
    recombinant, so the chloroplast density/length shape does NOT apply; the
    presence of several distinct core respiratory genes is the discriminator.

Hard limit (documented, not hidden): with a 5-column simplified GFF whose gene
ids are opaque (e.g. ``OkoG000543.1``), no genome FASTA and no informative
names, an anonymous mitochondrial scaffold is *indistinguishable* from a
nuclear scaffold by these heuristics. Such cases require ``--organelle-pattern``
or sequence similarity to a reference organellar gene set; the inference will
not guess. ``review_candidates`` lists ambiguous regions for manual judgement
instead of silently excluding them.
"""

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
    r"|(?:^|[_\-.|])chrmt?(?:$|[_\-.|])|(?:^|[_\-.|])mt_(?=\w)"
)
PLASTID_GENE_RE = re.compile(r"(?i)(?:^|[_\-.|])pgp\d+")

# Plant mitochondrial core genes (respiratory complexes + maturase + Tat +
# cytochrome c maturation). A nuclear scaffold rarely carries several distinct
# *intact* copies together; NUMTs are usually fragmented pseudogenes.
MT_CORE_GENE_RE = re.compile(
    r"(?i)(?:^|[_\-.|])"
    r"(nad[1-9]|nad4l|cob|cytb|cox[1-3]|atp[14689]|matr|mttb|ccm[bc]|ccmf[cn])"
    r"(?:$|[_\-.|])"
)

# --- Chloroplast shape heuristic (cpDNA biology: small, dense, intronless) ---
CP_MIN_GENES = 50
CP_MAX_GENES = 300
CP_MAX_LENGTH = 500_000
CP_MIN_DENSITY = 10.0          # genes per 100 kb
CP_MIN_INTRONLESS = 0.70       # fraction of genes whose span ~= CDS length

# --- Dense synthetic region routed to manual review (never auto-excluded) ---
REVIEW_MIN_GENES = 30
REVIEW_MAX_GENES = 300
REVIEW_MAX_LENGTH = 500_000
REVIEW_MIN_DENSITY = 15.0

# --- Mitochondrial gene-content heuristic ---
MT_MIN_CORE_GENES = 3          # distinct core families for auto-inference
MT_MAX_GENES = 300             # mt contigs are gene-modest; guards whole chromosomes
MT_MAX_LENGTH = 2_000_000      # plant mtDNA can exceed this; larger -> review only

_LIMITS_NOTE = (
    "Anonymous mitochondrial scaffolds in opaque-id 5-column GFFs without a "
    "genome cannot be auto-detected; use --organelle-pattern or sequence "
    "similarity. review_candidates are flagged for manual judgement, not "
    "excluded."
)


@dataclass(frozen=True, slots=True)
class OrganelleEvidence:
    seqid: str
    genome_type: str
    reason: str
    length: int
    gene_count: int
    genes_per_100kb: float
    intronless_fraction: float
    mt_core_genes: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, object]:
        return {
            "seqid": self.seqid,
            "genome_type": self.genome_type,
            "reason": self.reason,
            "length": self.length,
            "gene_count": self.gene_count,
            "genes_per_100kb": self.genes_per_100kb,
            "intronless_fraction": self.intronless_fraction,
            "mt_core_genes": list(self.mt_core_genes),
        }


def _identifier_type(values: list[str]) -> str:
    for value in values:
        if MITOCHONDRION_ID_RE.search(value) or re.search(
            r"(?i)mitochondr(?:ion|ial)?", value
        ):
            return "mitochondrion"
    for value in values:
        if (
            CHLOROPLAST_ID_RE.search(value)
            or PLASTID_GENE_RE.search(value)
            or re.search(r"(?i)(?:chloroplast|plastid)", value)
        ):
            return "chloroplast"
    return ""


def _mt_core_families(transcripts: list[Transcript]) -> tuple[str, ...]:
    families: set[str] = set()
    for transcript in transcripts:
        for value in (transcript.source_id, transcript.name, transcript.gene_id):
            if not value:
                continue
            for match in MT_CORE_GENE_RE.finditer(value):
                family = match.group(1).lower()
                if family == "cytb":
                    family = "cob"
                elif family in {"ccmfc", "ccmfn"}:
                    family = "ccmf"
                families.add(family)
    return tuple(sorted(families))


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
    mt_min_core_genes: int = MT_MIN_CORE_GENES,
) -> dict[str, object]:
    """Mark synthetic organelle regions and list unresolved candidates.

    Authoritative ``region`` rows and genome FASTA records are never changed.
    High-confidence signals (identifier, chloroplast shape, a cluster of
    distinct mitochondrial core genes) auto-exclude the region; weaker signals
    are appended to ``review_candidates`` for manual judgement.
    """

    result: dict[str, object] = {
        "enabled": enabled,
        "custom_pattern": custom_pattern,
        "heuristic": "identifier+gene_content+dense_compact",
        "mt_min_core_genes": mt_min_core_genes,
        "limits": _LIMITS_NOTE,
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
        mt_families = _mt_core_families(transcripts)
        density = len(transcripts) * 100_000 / max(1, region.length)
        compact_fraction = _intronless_fraction(transcripts, cds_lengths)
        pattern_matched = bool(pattern and pattern.search(seqid))
        cp_compact = (
            CP_MIN_GENES <= len(transcripts) <= CP_MAX_GENES
            and region.length <= CP_MAX_LENGTH
            and density >= CP_MIN_DENSITY
            and compact_fraction >= CP_MIN_INTRONLESS
        )
        mt_content = (
            len(mt_families) >= mt_min_core_genes
            and len(transcripts) <= MT_MAX_GENES
            and region.length <= MT_MAX_LENGTH
        )

        genome_type = ""
        reason = ""
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
        elif cp_compact:
            genome_type = "plastid"
            reason = "dense_compact_synthetic_region"
        elif mt_content:
            genome_type = "mitochondrion"
            reason = "mitochondrial_gene_content"

        if not genome_type:
            # No high-confidence signal: route ambiguous regions to review.
            if (
                REVIEW_MIN_GENES <= len(transcripts) <= REVIEW_MAX_GENES
                and region.length <= REVIEW_MAX_LENGTH
                and density >= REVIEW_MIN_DENSITY
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
            elif mt_families:
                review.append(
                    OrganelleEvidence(
                        seqid=seqid,
                        genome_type="mitochondrion",
                        reason="mitochondrial_gene_content_requires_review",
                        length=region.length,
                        gene_count=len(transcripts),
                        genes_per_100kb=round(density, 2),
                        intronless_fraction=round(compact_fraction, 4),
                        mt_core_genes=mt_families,
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
            mt_core_genes=mt_families if genome_type == "mitochondrion" else (),
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
