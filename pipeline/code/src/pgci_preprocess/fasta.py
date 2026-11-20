from __future__ import annotations

import bz2
import gzip
import io
import lzma
import re
import zipfile
from collections import Counter
from pathlib import Path
from typing import Any, TextIO

from pgci_preprocess.models import Anomaly, FastaData

PROTEIN_ID_RE = re.compile(r"\[protein_id=([^\]]+)\]")
PROTEIN_FIELD_RE = re.compile(r"(?:^|\s)Protein=([^\s]+)", re.IGNORECASE)
# Fallback for CDS headers that carry no protein_id qualifier.  NCBI
# ``cds_from_genomic`` style identifiers look like
# ``lcl|<seqid>_cds_<protein_id>_<cds_index>`` and some archives use
# ``cds_<protein_id>`` or ``cds.<protein_id>``.  The pattern is applied to the
# first FASTA token only, so ``$`` anchors at the end of the identifier and the
# optional trailing ``_<1-3 digits>`` is the CDS index, not part of the protein
# id.  A longer numeric tail is therefore preserved (for example an Ensembl- like
# ``ENSOP_00123``).  The previous pattern used ``\s``/``\d`` inside a raw string,
# which matched literal backslashes and therefore never fired as documented.
CDS_FALLBACK_RE = re.compile(
    r"(?:^|[|_])cds[._-](?P<id>.+?)(?:_\d{1,3})?$", re.IGNORECASE
)


def detect_compression(path: Path | str) -> str:
    """Detect compression from magic bytes, then use the suffix as a fallback."""

    path = Path(path)
    with path.open("rb") as handle:
        magic = handle.read(6)
    if magic.startswith(b"\x1f\x8b"):
        return "gzip"
    if magic.startswith(b"BZh"):
        return "bzip2"
    if magic.startswith(b"\xfd7zXZ\x00"):
        return "xz"
    if magic.startswith((b"PK\x03\x04", b"PK\x05\x06", b"PK\x07\x08")):
        return "zip"
    lowered = path.name.casefold()
    for suffix, compression in (
        (".gz", "gzip"),
        (".bz2", "bzip2"),
        (".xz", "xz"),
        (".zip", "zip"),
    ):
        if lowered.endswith(suffix):
            return compression
    return "plain"


class _ZipTextIO:
    """Context-managed text view over the single file in a ZIP archive."""

    def __init__(self, archive: zipfile.ZipFile, member: str) -> None:
        self._archive = archive
        self._raw = archive.open(member, "r")
        self._text = io.TextIOWrapper(
            self._raw,
            encoding="utf-8",
            errors="replace",
            newline="",
        )

    def __iter__(self):
        return iter(self._text)

    def __next__(self):
        return next(self._text)

    def read(self, *args):
        return self._text.read(*args)

    def readline(self, *args):
        return self._text.readline(*args)

    def close(self) -> None:
        try:
            self._text.close()
        finally:
            self._archive.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()


def _open_zip_member(path: Path, *, text: bool):
    archive = zipfile.ZipFile(path)
    members = [
        item
        for item in archive.infolist()
        if not item.is_dir() and item.file_size > 0
    ]
    if len(members) != 1:
        archive.close()
        raise ValueError(
            f"{path}: ZIP input must contain exactly one non-empty data file; "
            f"found {len(members)}"
        )
    if text:
        return _ZipTextIO(archive, members[0].filename)
    raw = archive.open(members[0].filename, "r")
    # Keep the archive reachable for the lifetime of the member handle.
    setattr(raw, "_pgci_archive", archive)
    return raw


def open_binary(path: Path | str) -> Any:
    """Open plain or compressed input for binary reading using content sniffing."""

    path = Path(path)
    compression = detect_compression(path)
    if compression == "gzip":
        return gzip.open(path, "rb")
    if compression == "bzip2":
        return bz2.open(path, "rb")
    if compression == "xz":
        return lzma.open(path, "rb")
    if compression == "zip":
        return _open_zip_member(path, text=False)
    return path.open("rb")


def open_text(path: Path | str) -> TextIO:
    """Open plain, gzip, bzip2, xz or single-member ZIP text by content magic."""

    path = Path(path)
    compression = detect_compression(path)
    if compression == "gzip":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    if compression == "bzip2":
        return bz2.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    if compression == "xz":
        return lzma.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    if compression == "zip":
        return _open_zip_member(path, text=True)
    return path.open("rt", encoding="utf-8", errors="replace", newline="")


def _read_fasta(path: Path | str, header_id) -> FastaData:
    path = Path(path)
    sequences: dict[str, str] = {}
    counts: Counter[str] = Counter()
    anomalies: list[Anomaly] = []
    current_id: str | None = None
    buffer: list[str] = []

    def flush() -> None:
        nonlocal current_id, buffer
        if current_id is None:
            return
        sequence = "".join(buffer).replace(" ", "").upper()
        counts["records"] += 1
        if not sequence:
            counts["empty_sequence"] += 1
            anomalies.append(Anomaly("warning", "empty_fasta_sequence", current_id))
        if current_id in sequences:
            if sequences[current_id] == sequence:
                counts["duplicate_identical"] += 1
                anomalies.append(Anomaly("warning", "duplicate_fasta_id_identical", current_id))
            else:
                counts["duplicate_conflict"] += 1
                anomalies.append(Anomaly("fatal", "duplicate_fasta_id_conflict", current_id))
        else:
            sequences[current_id] = sequence
        current_id = None
        buffer = []

    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                current_id = header_id(header)
                if not current_id:
                    counts["invalid_header"] += 1
                    anomalies.append(
                        Anomaly("warning", "fasta_header_without_id", f"line={line_number}")
                    )
                continue
            if current_id is None:
                counts["sequence_before_header"] += 1
                anomalies.append(
                    Anomaly("warning", "sequence_before_header", f"line={line_number}")
                )
                continue
            buffer.append("".join(line.split()))
    flush()
    return FastaData(sequences=sequences, counts=counts, anomalies=anomalies)


def extract_cds_id(header: str) -> str:
    """Resolve the protein identifier for one CDS FASTA header.

    Priority: an explicit ``[protein_id=...]`` qualifier, then a ``Protein=``
    field, then the ``cds`` marker in the first token, then the first token
    itself.  Nothing is guessed across records: an unmatched header simply
    yields its own identifier and the pipeline rejects the pairing later.
    """

    match = PROTEIN_ID_RE.search(header)
    if match:
        return match.group(1).strip()
    match = PROTEIN_FIELD_RE.search(header)
    if match:
        return match.group(1).strip()
    identifier = header.split(None, 1)[0].strip()
    match = CDS_FALLBACK_RE.search(identifier)
    if match:
        return match.group("id").strip()
    return identifier


def read_cds_fasta(path: Path | str) -> FastaData:
    return _read_fasta(path, extract_cds_id)


def read_pep_fasta(path: Path | str) -> FastaData:
    return _read_fasta(path, lambda header: header.split(None, 1)[0].strip())


def normalize_evm_protein_sequences(data: FastaData) -> dict[str, int]:
    """Normalize EVM protein FASTA stop-codon notation.

    EVM writes the terminal stop as ``.`` and may leave internal ``.``
    characters.  The historical converter rewrote a trailing dot to ``*`` and
    every remaining dot to ``X``.  The transformation is explicit, counted and
    recorded as an anomaly; it is never applied silently.
    """

    trailing = 0
    internal = 0
    affected: list[str] = []
    for sequence_id, sequence in list(data.sequences.items()):
        if "." not in sequence:
            continue
        affected.append(sequence_id)
        if sequence.endswith("."):
            sequence = sequence[:-1] + "*"
            trailing += 1
        if "." in sequence:
            sequence = sequence.replace(".", "X")
            internal += 1
        data.sequences[sequence_id] = sequence
    if affected:
        data.counts["dot_normalized_records"] = len(affected)
        data.counts["dot_trailing_stop"] = trailing
        data.counts["dot_internal_wildcard"] = internal
        data.anomalies.append(
            Anomaly(
                "warning",
                "protein_dot_normalized",
                f"records={len(affected)} trailing={trailing} internal={internal}",
            )
        )
    return {
        "records": len(affected),
        "trailing_stop": trailing,
        "internal_wildcard": internal,
    }

