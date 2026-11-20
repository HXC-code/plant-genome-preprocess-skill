"""Shared helpers for non-destructive annotation format conversion."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import TextIO
from urllib.parse import quote


class ConversionError(RuntimeError):
    """Raised when an input cannot be converted without guessing semantics."""


def sha256_file(path: Path | str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gff_escape(value: object) -> str:
    """Escape one GFF3 attribute value without changing its decoded form."""

    return quote(str(value), safe="._:-|")


def fasta_record(handle: TextIO, record_id: str, sequence: str, width: int = 60) -> None:
    handle.write(f">{record_id}\n")
    for offset in range(0, len(sequence), width):
        handle.write(sequence[offset : offset + width] + "\n")


def default_output_prefix(source: Path | str) -> str:
    """Derive a stable output prefix from a recognized compound suffix."""

    name = Path(source).name
    lowered = name.casefold()
    for suffix in (
        ".annotations.csv.gz",
        ".gbff.gz",
        ".gbff",
        ".gtf.gz",
        ".gtf",
        ".gff3.gz",
        ".gff3",
        ".csv.gz",
        ".csv",
    ):
        if lowered.endswith(suffix):
            name = name[: -len(suffix)]
            break
    else:
        name = Path(name).stem
    name = re.sub(r"(?:_genomic|_annotations?)$", "", name, flags=re.IGNORECASE)
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._-")
    return name or "annotation"
