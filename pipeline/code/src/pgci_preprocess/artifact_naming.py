"""Canonical WGDI artifact naming.

One species is identified by the globally unique abbreviation of the project
species table (``resources/species_abbreviations.tsv``).  Every artifact must
therefore satisfy two rules:

* directory and file names use ``<code>`` (``analysis_v2/<code>/<code>.cds`` …);
* every canonical gene ID is ``<code><sequence>g<order>`` (``Lcu1g0001``).

The helpers in this module are deliberately cheap: preprocessing calls
:func:`bundle_name_problems` and :func:`first_gene_id` after writing a bundle,
while the full-dataset audit lives in ``pgci_preprocess.naming_audit`` and in
``pgci_preprocess.validate``.
"""

from __future__ import annotations

import re
from pathlib import Path

ANALYSIS_REQUIRED_SUFFIXES = (
    ".pep",
    ".cds",
    ".new.gff",
    ".lens",
    ".qc.json",
    ".all.pep.fa",
    ".all.cds.fa",
    ".all.gff",
    ".all.lens",
    ".gene_map.tsv",
    ".chr_map.tsv",
)
CLEAN_REQUIRED_SUFFIXES = (".gene_map.tsv", ".chr_map.tsv", ".rejected.tsv")

# Gene IDs inside a FASTA header or a GFF column.
_GENE_ID_IN_LINE = re.compile(r"^[^\s>]+")


class NamingViolation(RuntimeError):
    """Raised when an artifact cannot be produced under its canonical code."""


def canonical_gene_id_regex(code: str) -> str:
    return rf"{re.escape(code)}\d+g\d{{4,}}"


def is_canonical_gene_id(gene_id: str, code: str) -> bool:
    return bool(re.fullmatch(canonical_gene_id_regex(code), str(gene_id)))


def legacy_gene_id_regex(show_id: str) -> str:
    # Legacy IDs are ``s<show_id><clean_sequence>g<order>`` (e.g. s199767g0001
    # for show_id 199, sequence 767), so the sequence digits must be allowed
    # between the show_id and the ``g`` separator.
    return rf"s{re.escape(str(show_id))}\d+g\d{{4,}}"


def legacy_prefix(show_id: str) -> str:
    return f"s{str(show_id).strip()}"


def _open_text(path: Path):
    # utf-8-sig also tolerates a stray BOM at the start of a FASTA/GFF file.
    return path.open(encoding="utf-8-sig", errors="replace")


def first_gene_id(path: Path, kind: str) -> str:
    """Return the first gene ID of a ``gff``/``cds``/``pep`` artifact.

    ``kind`` is ``gff`` (second column) or ``fasta`` (first header token).
    The read stops after the first record, so it is cheap even for 40 MB files.
    """

    with _open_text(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if kind == "gff":
                fields = stripped.split("\t")
                return fields[1] if len(fields) > 1 else ""
            if stripped.startswith(">"):
                header = stripped[1:].split(None, 1)
                return header[0] if header else ""
            return ""
    return ""


def bundle_name_problems(
    directory: Path,
    code: str,
    *,
    required_suffixes: tuple[str, ...],
    kind_map: dict[str, str] | None = None,
    check_dir_name: bool = True,
) -> list[str]:
    """Check directory/file naming and the first gene ID of a bundle.

    ``kind_map`` maps a suffix to ``gff``/``fasta``; suffixes present in the
    mapping are also checked for a canonical first gene ID.  Staging
    directories written by preprocessing use a temporary name, so they pass
    ``check_dir_name=False``.
    """

    problems: list[str] = []
    if not directory.is_dir():
        return [f"bundle_missing:{directory.name}"]
    if check_dir_name and directory.name != code:
        problems.append(f"bundle_dir_not_code:{directory.name}!={code}")
    present = {path.name for path in directory.iterdir() if path.is_file()}
    for suffix in required_suffixes:
        expected = f"{code}{suffix}"
        if expected not in present:
            problems.append(f"missing_file:{expected}")
    foreign = sorted(
        name
        for name in present
        if not name.startswith(f"{code}.") and not name.startswith(f"{code}_")
    )
    if foreign:
        problems.append("foreign_files:" + ",".join(foreign[:5]))
    for suffix, kind in (kind_map or {}).items():
        path = directory / f"{code}{suffix}"
        if not path.is_file():
            continue
        gene_id = first_gene_id(path, kind)
        if not gene_id:
            problems.append(f"empty_first_record:{path.name}")
        elif not is_canonical_gene_id(gene_id, code):
            problems.append(f"non_canonical_gene_id:{path.name}:{gene_id}")
    return problems
