"""Role-aware input routing for annotation and sequence files.

Suffixes are useful hints, not authority.  Annotation decisions are confirmed
from content; sequence-role selection uses both filename semantics and a bounded
FASTA content preview.  Every decision can be written to the preparation report
so a later run can reproduce why a file was mirrored, converted or extracted.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from pgci_preprocess.annotation_convert import ConversionError, detect_conversion_format
from pgci_preprocess.fasta import open_text
from pgci_preprocess.gff_formats import GffFormatError, detect_gff_format


class InputRoutingError(RuntimeError):
    pass


@dataclass(frozen=True, slots=True)
class AnnotationRoute:
    source: Path
    suffix_hint: str
    content_format: str
    route: str
    note: str = ""


def _compound_suffix(path: Path) -> str:
    suffixes = [item.casefold() for item in path.suffixes]
    if suffixes and suffixes[-1] in {".gz", ".bz2", ".xz", ".zip"}:
        suffixes = suffixes[:-1]
    return suffixes[-1] if suffixes else ""


def annotation_suffix_hint(path: Path | str) -> str:
    suffix = _compound_suffix(Path(path))
    if suffix in {".gtf"}:
        return "gtf"
    if suffix in {".gbff", ".gbk", ".gb", ".genbank"}:
        return "gbff"
    if suffix in {".gff", ".gff3"}:
        return "gff"
    return ""


def route_annotation(path: Path | str) -> AnnotationRoute:
    """Classify annotation content and choose its preprocessing route."""

    source = Path(path)
    hint = annotation_suffix_hint(source)
    try:
        detected = detect_conversion_format(source)
    except ConversionError as conversion_error:
        try:
            detection = detect_gff_format(source)
        except GffFormatError as gff_error:
            raise InputRoutingError(
                f"{source}: cannot classify annotation by content; "
                f"converter={conversion_error}; gff={gff_error}"
            ) from gff_error
        content = detection.format
        route = "direct_gff"
    else:
        content = detected
        route = {
            "gtf": "convert_gtf",
            "gbff": "convert_gbff",
            "functional-go": "metadata_only",
        }[detected]
    note = ""
    if hint and not (
        (hint == "gff" and content in {"standard", "transcript_only", "simplified5", "gene_span", "gene_cds_direct"})
        or hint == content
    ):
        note = f"suffix_hint={hint} conflicts with content={content}; content wins"
    return AnnotationRoute(
        source=source,
        suffix_hint=hint,
        content_format=content,
        route=route,
        note=note,
    )


def _candidate_paths(folder: Path, value: str, role: str, show_id: str) -> list[Path]:
    paths: list[Path] = []
    for raw in value.split("|"):
        text = raw.strip()
        if not text:
            continue
        path = Path(text)
        if not path.is_absolute():
            path = folder / path
        path = path.resolve()
        try:
            path.relative_to(folder.resolve())
        except ValueError as exc:
            raise InputRoutingError(
                f"s{show_id}: {role} candidate escapes species folder: {path}"
            ) from exc
        if not path.is_file() or path.stat().st_size == 0:
            continue
        paths.append(path)
    return paths


def _fasta_preview(path: Path) -> tuple[str, str]:
    """Return (content_kind, first_id) using at most the first FASTA record."""

    header = ""
    sequence: list[str] = []
    try:
        with open_text(path) as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header:
                        break
                    header = line[1:].strip()
                    continue
                if not header:
                    return "not_fasta", ""
                sequence.append("".join(line.split()))
                if sum(map(len, sequence)) >= 16384:
                    break
    except (OSError, EOFError, ValueError) as exc:
        raise InputRoutingError(f"{path}: compressed input cannot be read: {exc}") from exc
    if not header:
        return "not_fasta", ""
    first_id = header.split(None, 1)[0]
    letters = set("".join(sequence).upper())
    nucleotide = {"A", "C", "G", "T", "U", "N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}
    if letters and letters <= nucleotide:
        return "nucleotide_fasta", first_id
    protein = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X", "B", "Z", "J", "U", "O", "*"}
    if letters and letters <= protein:
        return "protein_fasta", first_id
    return "fasta", first_id


def _role_score(path: Path, role: str) -> tuple[int, str]:
    name = path.name.casefold()
    if role == "gff":
        try:
            route = route_annotation(path)
        except InputRoutingError:
            return (-100, "not_annotation")
        return (10, route.content_format)
    kind, _first_id = _fasta_preview(path)
    if kind not in {"nucleotide_fasta", "protein_fasta", "fasta"}:
        return (-100, kind)
    if role == "pep":
        hinted = bool(re.search(r"(?:protein|pep|\.pro(?:\.|$))", name))
        if kind not in {"protein_fasta", "fasta"} and not hinted:
            return (-100, kind)
        return ((8 if hinted else 0) + (5 if kind == "protein_fasta" else 0), kind)
    hinted = bool(re.search(r"(?:genome|chromosome|scaffold|assembly)", name))
    cds_named = bool(re.search(r"(?:^|[._-])cds(?:[._-]|$)", name))
    if role == "genome":
        if kind != "nucleotide_fasta":
            return (-100, kind)
        protein_named = bool(re.search(r"(?:protein|pep|\.pro(?:\.|$))", name))
        score = 0
        score += 8 if hinted else 0
        score += 5 if kind == "nucleotide_fasta" else 0
        score -= 20 if cds_named or protein_named else 0
        return (score, kind)
    if role == "cds":
        if kind != "nucleotide_fasta":
            return (-100, kind)
        score = 5
        score += 8 if cds_named else 0
        score -= 8 if hinted and not cds_named else 0
        score -= 20 if re.search(r"(?:protein|pep)", name) else 0
        return (score, kind)
    return (0, kind)


def select_role_source(
    folder: Path,
    value: str,
    role: str,
    show_id: str,
) -> tuple[Path, int, str]:
    """Resolve a role value, including pipe-separated download candidates."""

    if not value.strip():
        raise InputRoutingError(f"s{show_id}: missing {role} path")
    paths = _candidate_paths(folder, value, role, show_id)
    if not paths:
        raise InputRoutingError(
            f"s{show_id}: no non-empty {role} candidate exists for {value!r}"
        )
    scored = [(_role_score(path, role), path) for path in paths]
    best_score = max(score for (score, _detail), _path in scored)
    best = [
        (detail, path)
        for ((score, detail), path) in scored
        if score == best_score
    ]
    if best_score < 0:
        raise InputRoutingError(
            f"s{show_id}: no {role}-compatible file among candidates: "
            + ", ".join(path.name for path in paths)
        )
    if len(best) > 1:
        raise InputRoutingError(
            f"s{show_id}: ambiguous {role} candidates at score {best_score}: "
            + ", ".join(path.name for _detail, path in best)
        )
    detail, selected = best[0]
    return selected, len(paths), detail
