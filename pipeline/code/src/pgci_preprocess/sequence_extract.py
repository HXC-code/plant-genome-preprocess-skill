"""Derive a CDS/PEP FASTA pair from genome FASTA plus GFF CDS coordinates."""

from __future__ import annotations

import gzip
import os
import uuid
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from pgci_preprocess.fasta import open_text
from pgci_preprocess.gff import parse_attributes, strip_prefix
from pgci_preprocess.gff_formats import parse_with_adapter
from pgci_preprocess.sequence import translate_cds


class SequenceExtractionError(RuntimeError):
    pass


@dataclass(frozen=True, slots=True)
class SequenceExtractionResult:
    genome_path: Path
    gff_path: Path
    cds_path: Path
    pep_path: Path
    counts: dict[str, int]
    warnings: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "genome": str(self.genome_path),
            "gff": str(self.gff_path),
            "cds": str(self.cds_path),
            "pep": str(self.pep_path),
            "counts": self.counts,
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True, slots=True)
class _Segment:
    seqid: str
    start: int
    end: int
    strand: str
    transcript_id: str
    protein_id: str
    translation_table: int


def _coordinate(text: str, line_number: int) -> int:
    try:
        return int(text)
    except ValueError as exc:
        raise SequenceExtractionError(
            f"GFF line {line_number}: invalid coordinate {text!r}"
        ) from exc


def _read_segments(
    gff_path: Path,
    parsed_transcripts,
) -> tuple[dict[str, list[_Segment]], Counter[str]]:
    counts: Counter[str] = Counter()
    segments: dict[str, list[_Segment]] = defaultdict(list)
    with open_text(gff_path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if raw_line.startswith("#") or not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise SequenceExtractionError(
                    f"{gff_path.name}:{line_number}: expected 9 columns, "
                    f"found {len(fields)}"
                )
            seqid, _source, feature, start_text, end_text, _score, strand, phase_text, attr_text = fields
            if feature != "CDS":
                continue
            counts["cds_rows"] += 1
            start = _coordinate(start_text, line_number)
            end = _coordinate(end_text, line_number)
            if start < 1 or end < start:
                raise SequenceExtractionError(
                    f"{gff_path.name}:{line_number}: invalid CDS interval {start}-{end}"
                )
            if strand not in {"+", "-"}:
                raise SequenceExtractionError(
                    f"{gff_path.name}:{line_number}: invalid CDS strand {strand!r}"
                )
            if phase_text not in {"", ".", "0", "1", "2"}:
                raise SequenceExtractionError(
                    f"{gff_path.name}:{line_number}: invalid CDS phase {phase_text!r}"
                )
            attributes = parse_attributes(attr_text)
            raw_table = attributes.get("transl_table", "").strip() or "1"
            try:
                translation_table = int(raw_table)
            except ValueError as exc:
                raise SequenceExtractionError(
                    f"{gff_path.name}:{line_number}: invalid transl_table {raw_table!r}"
                ) from exc
            parents = [
                strip_prefix(strip_prefix(item.strip(), "rna-"), "gene-")
                for item in attributes.get("Parent", "").split(",")
                if item.strip()
            ]
            if not parents:
                counts["cds_without_parent"] += 1
                continue
            explicit_protein = (
                attributes.get("protein_id", "")
                or attributes.get("Protein_Accession", "")
            ).strip()
            raw_id = attributes.get("ID", "").strip()
            if not explicit_protein and raw_id.startswith("cds-"):
                explicit_protein = raw_id[len("cds-") :]
            for parent in parents:
                transcript = parsed_transcripts.get(parent)
                if transcript is None:
                    counts["cds_parent_missing"] += 1
                    continue
                if explicit_protein:
                    protein_id = explicit_protein
                elif len(transcript.protein_ids) == 1:
                    protein_id = next(iter(transcript.protein_ids))
                elif not transcript.protein_ids:
                    protein_id = parent
                else:
                    raise SequenceExtractionError(
                        f"transcript {parent!r} has multiple protein IDs: "
                        f"{sorted(transcript.protein_ids)}"
                    )
                segments[seqid].append(
                    _Segment(
                        seqid=seqid,
                        start=start,
                        end=end,
                        strand=strand,
                        transcript_id=parent,
                        protein_id=protein_id,
                        translation_table=translation_table,
                    )
                )
                counts["cds_segments"] += 1
    return segments, counts


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN"))[::-1]


def _write_record(handle, identifier: str, sequence: str) -> None:
    handle.write(f">{identifier}\n")
    for offset in range(0, len(sequence), 80):
        handle.write(sequence[offset : offset + 80] + "\n")


def _open_staging_output(path: Path, *, gzip_output: bool):
    if gzip_output:
        return gzip.open(path, "wt", encoding="utf-8", newline="\n")
    return path.open("w", encoding="utf-8", newline="\n")


def extract_cds_pep_from_gff(
    genome_path: Path | str,
    gff_path: Path | str,
    cds_output: Path | str,
    pep_output: Path | str,
) -> SequenceExtractionResult:
    """Extract spliced CDS and translated PEP records from genome + GFF."""

    genome = Path(genome_path)
    gff = Path(gff_path)
    cds_out = Path(cds_output)
    pep_out = Path(pep_output)
    if not genome.is_file() or genome.stat().st_size == 0:
        raise SequenceExtractionError(f"genome FASTA not found or empty: {genome}")
    if not gff.is_file() or gff.stat().st_size == 0:
        raise SequenceExtractionError(f"GFF not found or empty: {gff}")
    if cds_out.exists() or pep_out.exists():
        raise SequenceExtractionError(
            "derived CDS/PEP output already exists; prepare reuses a complete "
            "derived bundle rather than overwriting it"
        )

    parsed = parse_with_adapter(gff).data
    fatal = [item for item in parsed.anomalies if item.level == "fatal"]
    if fatal:
        raise SequenceExtractionError(
            f"GFF cannot be used for extraction: {fatal[0].code}: {fatal[0].detail}"
        )
    segments_by_seqid, counts = _read_segments(gff, parsed.transcripts)
    if not segments_by_seqid:
        raise SequenceExtractionError(
            f"{gff.name}: no CDS rows with resolvable transcript parents; "
            "genome+GFF extraction requires coding annotation"
        )

    by_transcript: dict[str, list[_Segment]] = defaultdict(list)
    for segments in segments_by_seqid.values():
        for segment in segments:
            by_transcript[segment.transcript_id].append(segment)
    protein_to_transcripts: dict[str, set[str]] = defaultdict(set)
    for transcript_id, transcript_segments in by_transcript.items():
        proteins = {item.protein_id for item in transcript_segments}
        if len(proteins) != 1:
            raise SequenceExtractionError(
                f"transcript {transcript_id!r} has multiple protein IDs: "
                f"{sorted(proteins)}"
            )
        protein_to_transcripts[next(iter(proteins))].add(transcript_id)
    reused = {
        protein: transcripts
        for protein, transcripts in protein_to_transcripts.items()
        if len(transcripts) > 1
    }
    if reused:
        protein, transcripts = next(iter(reused.items()))
        raise SequenceExtractionError(
            f"protein ID {protein!r} is reused by multiple transcripts: "
            f"{sorted(transcripts)}"
        )

    extracted: dict[str, tuple[str, str, int]] = {}
    transcript_seqid: dict[str, str] = {}
    genome_lengths: dict[str, int] = {}
    warnings: list[str] = []
    current_id = ""
    current_chunks: list[str] = []

    def flush_record() -> None:
        nonlocal current_id, current_chunks
        if not current_id:
            return
        sequence = "".join(current_chunks).replace(" ", "").upper()
        genome_lengths[current_id] = len(sequence)
        counts["genome_records"] += 1
        seqid_segments = segments_by_seqid.get(current_id, [])
        for transcript_id in sorted({item.transcript_id for item in seqid_segments}):
            transcript_segments = sorted(
                by_transcript[transcript_id],
                key=lambda item: (item.start, item.end),
            )
            if any(item.start > len(sequence) or item.end > len(sequence) for item in transcript_segments):
                raise SequenceExtractionError(
                    f"transcript {transcript_id!r} CDS exceeds genome record "
                    f"{current_id!r} length {len(sequence)}"
                )
            strands = {item.strand for item in transcript_segments}
            if len(strands) != 1:
                raise SequenceExtractionError(
                    f"transcript {transcript_id!r} CDS segments mix strands"
                )
            table_values = {item.translation_table for item in transcript_segments}
            if len(table_values) != 1:
                raise SequenceExtractionError(
                    f"transcript {transcript_id!r} CDS segments mix transl_table values"
                )
            cds = "".join(
                sequence[item.start - 1 : item.end]
                for item in transcript_segments
            )
            if next(iter(strands)) == "-":
                cds = _reverse_complement(cds)
            if not cds:
                counts["empty_cds"] += 1
                continue
            if len(cds) % 3:
                counts["cds_length_not_multiple_of_3"] += 1
                continue
            translation_table = next(iter(table_values))
            protein = translate_cds(cds, translation_table)
            if "*" in protein[:-1]:
                counts["internal_stop"] += 1
                continue
            protein = protein.rstrip("*")
            if not protein:
                counts["empty_translation"] += 1
                continue
            protein_id = transcript_segments[0].protein_id
            extracted[protein_id] = (cds, protein, translation_table)
            transcript_seqid[protein_id] = current_id
        current_id = ""
        current_chunks = []

    with open_text(genome) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush_record()
                header = line[1:].strip()
                current_id = header.split(None, 1)[0] if header else ""
                if not current_id:
                    raise SequenceExtractionError(
                        f"{genome.name}:{line_number}: FASTA header without ID"
                    )
                current_chunks = []
            elif not current_id:
                raise SequenceExtractionError(
                    f"{genome.name}:{line_number}: sequence before FASTA header"
                )
            else:
                current_chunks.extend(line.split())
    flush_record()

    missing_seqids = sorted(set(segments_by_seqid) - set(genome_lengths))
    if missing_seqids:
        raise SequenceExtractionError(
            "GFF CDS seqids missing from genome FASTA: "
            + ",".join(missing_seqids[:10])
        )
    for seqid, gff_region in parsed.regions.items():
        genome_length = genome_lengths.get(seqid)
        if genome_length is None:
            continue
        if gff_region.length > genome_length:
            raise SequenceExtractionError(
                f"{seqid}: GFF region length {gff_region.length} exceeds "
                f"genome FASTA length {genome_length}"
            )
        if not gff_region.synthetic and gff_region.length != genome_length:
            warnings.append(
                f"{seqid}: GFF region length {gff_region.length} differs from "
                f"genome FASTA length {genome_length}"
            )
            counts["region_length_mismatch"] += 1

    if not extracted:
        raise SequenceExtractionError(
            "no valid CDS/PEP records could be extracted from genome + GFF"
        )

    cds_out.parent.mkdir(parents=True, exist_ok=True)
    pep_out.parent.mkdir(parents=True, exist_ok=True)
    token = uuid.uuid4().hex
    cds_staging = cds_out.parent / f".{cds_out.name}.{token}.tmp"
    pep_staging = pep_out.parent / f".{pep_out.name}.{token}.tmp"
    try:
        with _open_staging_output(
            cds_staging,
            gzip_output=cds_out.name.casefold().endswith(".gz"),
        ) as cds_handle, _open_staging_output(
            pep_staging,
            gzip_output=pep_out.name.casefold().endswith(".gz"),
        ) as pep_handle:
            for protein_id in sorted(extracted):
                cds, protein, _table = extracted[protein_id]
                _write_record(cds_handle, protein_id, cds)
                _write_record(pep_handle, protein_id, protein)
        os.replace(cds_staging, cds_out)
        os.replace(pep_staging, pep_out)
    finally:
        for path in (cds_staging, pep_staging):
            if path.exists():
                path.unlink()

    counts["transcripts"] = len(extracted)
    counts["proteins"] = len(extracted)
    counts["sequences"] = len(set(transcript_seqid.values()))
    if counts["cds_parent_missing"]:
        warnings.append(
            "CDS rows with unresolved transcript parents were skipped: "
            f"{counts['cds_parent_missing']}"
        )
    for reason, label in (
        ("cds_length_not_multiple_of_3", "CDS records skipped because length was not divisible by 3"),
        ("internal_stop", "CDS records skipped because of an internal stop codon"),
        ("empty_cds", "empty CDS records skipped"),
        ("empty_translation", "empty translations skipped"),
    ):
        if counts[reason]:
            warnings.append(f"{label}: {counts[reason]}")
    return SequenceExtractionResult(
        genome_path=genome,
        gff_path=gff,
        cds_path=cds_out,
        pep_path=pep_out,
        counts=dict(sorted(counts.items())),
        warnings=tuple(warnings),
    )
