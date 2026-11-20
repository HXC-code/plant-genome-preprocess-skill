"""Role-level source manifest generation and verification.

A frozen species index (e.g. Golden-109) gets one manifest row per input role
with a verified path, size, SHA-256, content-level format and conversion chain.
Provenance repairs (dead-path remaps, normalized-artifact -> original
replacements) are applied only from an explicit evidence-backed repair TSV;
nothing is guessed from suffixes and no original file is ever modified.

The builder is deterministic: identical inputs produce byte-identical manifests,
so the manifest content hash can be frozen alongside the scope.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path

from pgci_preprocess.fasta import detect_compression, open_binary
from pgci_preprocess.gff_formats import GffFormatError, detect_gff_format
from pgci_preprocess.species_names import load_species_names
from pgci_preprocess.version import __version__

ROLES = ("gff", "cds", "pep", "genome")
REPAIR_ACTIONS = {"remap", "replace_with_original"}

MANIFEST_FIELDS = [
    "code",
    "show_id",
    "species_code",
    "species",
    "order",
    "family",
    "assembly_accession",
    "role",
    "source_path",
    "source_origin",
    "bytes",
    "sha256",
    "compression",
    "content_format",
    "integrity",
    "conversion_required",
    "converter",
    "converter_version",
    "conversion_chain",
    "provenance_note",
    "status",
]

REPAIR_FIELDS = [
    "code",
    "role",
    "action",
    "new_path",
    "expected_sha256",
    "evidence",
    "note",
]

GENOME_OMITTED_REASON = (
    "genome not recorded in baseline QC; optional for the WGDI four-piece bundle "
    "because regions resolve from GFF region rows or transcript coordinates"
)


class SourceManifestError(RuntimeError):
    """Raised when a manifest cannot be built or verified without guessing."""


@dataclass(frozen=True)
class Repair:
    code: str
    role: str
    action: str
    new_path: Path
    expected_sha256: str
    evidence: str
    note: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_repairs(path: Path | str) -> dict[tuple[str, str], Repair]:
    path = Path(path)
    repairs: dict[tuple[str, str], Repair] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(REPAIR_FIELDS) - set(reader.fieldnames or [])
        if missing:
            raise SourceManifestError(f"repair TSV missing columns: {sorted(missing)}")
        for row in reader:
            code = (row["code"] or "").strip()
            role = (row["role"] or "").strip()
            action = (row["action"] or "").strip()
            if role not in ROLES:
                raise SourceManifestError(f"repair role must be one of {ROLES}: {role!r}")
            if action not in REPAIR_ACTIONS:
                raise SourceManifestError(f"repair action must be one of {sorted(REPAIR_ACTIONS)}: {action!r}")
            expected = (row["expected_sha256"] or "").strip().lower()
            if len(expected) != 64 or any(c not in "0123456789abcdef" for c in expected):
                raise SourceManifestError(f"repair expected_sha256 must be a 64-hex digest: {code}/{role}")
            if not (row["evidence"] or "").strip():
                raise SourceManifestError(f"repair requires evidence text: {code}/{role}")
            key = (code, role)
            if key in repairs:
                raise SourceManifestError(f"duplicate repair for {code}/{role}")
            repairs[key] = Repair(
                code=code,
                role=role,
                action=action,
                new_path=Path(row["new_path"].strip()),
                expected_sha256=expected,
                evidence=row["evidence"].strip(),
                note=(row["note"] or "").strip(),
            )
    return repairs


def _sniff_fasta(path: Path) -> tuple[str, str]:
    compression = detect_compression(path)
    with open_binary(path) as handle:
        head = handle.read(4096)
    if not head.strip():
        return compression, "EMPTY"
    if head.lstrip().startswith(b">"):
        return compression, "FASTA"
    return compression, "NOT_FASTA"


def _sniff_annotation(path: Path) -> str:
    """Content-level annotation type: GBFF, GTF or GFF-family."""
    with open_binary(path) as handle:
        head = handle.read(256 * 1024).decode("utf-8", errors="replace")
    lines = head.splitlines()
    if lines and lines[0].startswith("LOCUS"):
        return "GBFF"
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) == 9:
            attrs = fields[8]
            if (
                "ID=" not in attrs
                and any(
                    key in attrs
                    for key in ("gene_id", "transcript_id", "locus_tag")
                )
                and re.search(r'(?:gene_id|transcript_id|locus_tag|gene)\s+"', attrs)
            ):
                return "GTF"
            return "GFF"
        break
    return "OTHER"


def _find_remap(show_id: str, basename: str, raw_roots: tuple[Path, ...]) -> list[Path]:
    candidates: list[Path] = []
    for root in raw_roots:
        if not root.is_dir():
            continue
        for folder in sorted(root.iterdir()):
            if not folder.is_dir() or folder.name.split("_", 1)[0] != show_id:
                continue
            for path in sorted(folder.rglob(basename)):
                if path.is_file():
                    candidates.append(path)
    return candidates


def _origin(path: Path, project_root: Path, raw_roots: tuple[Path, ...]) -> str:
    resolved = path.resolve()
    for root in raw_roots:
        try:
            resolved.relative_to(root.resolve())
            return f"raw:{root}"
        except ValueError:
            continue
    try:
        resolved.relative_to((project_root / "data" / "downloads").resolve())
        return "canonical_mirror:data/downloads"
    except ValueError:
        return "other"


def build_source_manifest(
    project_root: Path | str,
    *,
    index_path: Path | str | None = None,
    repairs_path: Path | str | None = None,
    raw_roots: tuple[Path | str, ...] = (),
    output_path: Path | str | None = None,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    """Build the role-level manifest for every species in the frozen index.

    Blockers are kept as rows and reported; one bad species never stops the
    others.  Returns ``(rows, summary)`` and optionally writes the TSV plus a
    ``.sha256`` sidecar next to it.
    """

    project_root = Path(project_root).resolve()
    index_path = Path(index_path) if index_path else project_root / "resources" / "wgdi_v2_index.tsv"
    if not index_path.is_file():
        raise SourceManifestError(f"index not found: {index_path}")
    repairs = load_repairs(repairs_path) if repairs_path else {}
    roots = tuple(Path(root).resolve() for root in raw_roots) or (
        project_root / "data" / "downloads",
    )

    with index_path.open(encoding="utf-8", newline="") as handle:
        index_rows = list(csv.DictReader(handle, delimiter="\t"))
    if not index_rows:
        raise SourceManifestError(f"index has no rows: {index_path}")

    rows: list[dict[str, str]] = []
    blockers: list[str] = []
    status_counts: dict[str, int] = {}

    def emit(entry: dict[str, str]) -> None:
        rows.append(entry)
        status_counts[entry["status"]] = status_counts.get(entry["status"], 0) + 1
        if entry["status"].startswith("BLOCKER"):
            blockers.append(f"{entry['code']}/{entry['role']}:{entry['status']}")

    ordered = sorted(
        index_rows,
        key=lambda r: int(r["code"][1:]) if r["code"][1:].isdigit() else 10**12,
    )
    for index_row in ordered:
        code = index_row["code"]
        show_id = str(index_row.get("show_id", "")).strip()
        if not show_id:
            names = load_species_names(project_root)
            if names is not None:
                show_id = names.resolve_show_id(code)
            elif code.casefold().startswith("s") and code[1:].isdigit():
                show_id = code[1:]
            else:
                raise SourceManifestError(
                    f"index row {code!r} has no show_id and cannot be resolved"
                )
        qc_path = project_root / index_row["qc"]
        qc = json.loads(qc_path.read_text(encoding="utf-8"))
        inputs = qc.get("inputs", {})
        base = {
            "code": code,
            "species_code": code,
            "show_id": show_id,
            "species": index_row["species"],
            "order": qc.get("order", ""),
            "family": qc.get("family", ""),
            "assembly_accession": index_row["accession"],
        }
        for role in ROLES:
            entry = {field: "" for field in MANIFEST_FIELDS}
            entry.update(base)
            entry["role"] = role
            info = inputs.get(role)
            repair = repairs.get((code, role))

            if info is None and repair is None:
                if role == "genome":
                    entry["status"] = "OMITTED_OPTIONAL"
                    entry["integrity"] = "not_used"
                    entry["provenance_note"] = GENOME_OMITTED_REASON
                else:
                    entry["status"] = "BLOCKER_MISSING"
                    entry["integrity"] = "no_qc_record"
                    entry["provenance_note"] = "required role has no QC input record"
                emit(entry)
                continue

            recorded = str((info or {}).get("path", ""))
            recorded_path = Path(recorded) if recorded else Path()
            if recorded and not recorded_path.is_absolute():
                recorded_path = project_root / recorded_path
            expected_sha = str((info or {}).get("sha256", "")).lower()
            expected_bytes = (info or {}).get("bytes")

            if repair is not None:
                source = repair.new_path
                if not source.is_absolute():
                    source = project_root / source
                expected_sha = repair.expected_sha256
                expected_bytes = None
                entry["conversion_chain"] = (
                    f"qc_recorded={recorded}; repair={repair.action}; evidence={repair.evidence}"
                )
                entry["provenance_note"] = repair.note
                status_ok = "REPAIRED_ORIGINAL" if repair.action == "replace_with_original" else "REPAIRED_REMAP"
            else:
                source = recorded_path
                status_ok = "VERIFIED"
                if recorded and not source.is_file():
                    verified_candidate: Path | None = None
                    for candidate in _find_remap(show_id, Path(recorded).name, roots):
                        if expected_bytes is not None and candidate.stat().st_size != expected_bytes:
                            continue
                        if expected_sha and sha256_file(candidate) == expected_sha:
                            verified_candidate = candidate
                            break
                    if verified_candidate is None:
                        entry["source_path"] = recorded
                        entry["status"] = "BLOCKER_MISSING"
                        entry["integrity"] = "recorded_path_dead_no_hash_verified_remap"
                        emit(entry)
                        continue
                    source = verified_candidate
                    status_ok = "REMAPPED_HASH_VERIFIED"
                    entry["conversion_chain"] = f"qc_recorded={recorded}; remap_hash_verified=yes"

            if not source.is_file():
                entry["source_path"] = str(source)
                entry["status"] = "BLOCKER_MISSING"
                entry["integrity"] = "file_not_found"
                emit(entry)
                continue

            actual_bytes = source.stat().st_size
            actual_sha = sha256_file(source)
            entry["source_path"] = str(source)
            entry["source_origin"] = _origin(source, project_root, roots)
            entry["bytes"] = str(actual_bytes)
            entry["sha256"] = actual_sha
            if expected_sha and actual_sha != expected_sha:
                entry["status"] = "BLOCKER_HASH_MISMATCH"
                entry["integrity"] = f"expected={expected_sha}"
                emit(entry)
                continue
            if expected_bytes is not None and actual_bytes != expected_bytes:
                entry["status"] = "BLOCKER_SIZE_MISMATCH"
                entry["integrity"] = f"expected_bytes={expected_bytes}"
                emit(entry)
                continue
            entry["integrity"] = "sha256_verified"

            if role == "gff":
                annotation_kind = _sniff_annotation(source)
                entry["compression"] = "gzip" if source.suffix == ".gz" else "plain"
                if annotation_kind in {"GBFF", "GTF"}:
                    entry["content_format"] = annotation_kind
                    entry["conversion_required"] = annotation_kind.lower()
                    entry["converter"] = "pgci_preprocess.annotation_convert"
                    entry["converter_version"] = __version__
                    entry["status"] = status_ok
                    emit(entry)
                    continue
                try:
                    detection = detect_gff_format(source)
                except GffFormatError as exc:
                    entry["status"] = "BLOCKER_FORMAT"
                    entry["content_format"] = f"ERROR:{exc}"
                    emit(entry)
                    continue
                entry["content_format"] = detection.format
                if detection.format != "standard":
                    entry["conversion_required"] = detection.format
                    entry["converter"] = "pgci_preprocess.gff_formats"
                    entry["converter_version"] = __version__
                else:
                    entry["conversion_required"] = "none"
            else:
                compression, fmt = _sniff_fasta(source)
                entry["compression"] = compression
                entry["content_format"] = fmt
                entry["conversion_required"] = "none"
                if fmt != "FASTA":
                    # A repair may legitimately point CDS/PEP/genome at the
                    # original GBFF that the gbff adapter extracts them from.
                    if repair is not None and _sniff_annotation(source) == "GBFF":
                        entry["content_format"] = "GBFF"
                        entry["conversion_required"] = "gbff"
                        entry["converter"] = "pgci_preprocess.annotation_convert"
                        entry["converter_version"] = __version__
                    else:
                        entry["status"] = "BLOCKER_FORMAT"
                        entry["integrity"] = f"role={role} content is {fmt}, expected FASTA"
                        emit(entry)
                        continue
            entry["status"] = status_ok
            emit(entry)

    summary: dict[str, object] = {
        "index": str(index_path),
        "index_sha256": _sha256_bytes(index_path.read_bytes()),
        "species_count": len(index_rows),
        "row_count": len(rows),
        "status_counts": status_counts,
        "blockers": blockers,
        "repairs_applied": len(repairs),
        "builder_version": __version__,
    }
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        import io

        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(buffer, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        data = buffer.getvalue().encode("utf-8")
        output_path.write_bytes(data)
        digest = _sha256_bytes(data)
        output_path.with_suffix(output_path.suffix + ".sha256").write_text(
            f"{digest}  {output_path.name}\n", encoding="utf-8"
        )
        summary["manifest_path"] = str(output_path)
        summary["manifest_sha256"] = digest
    return rows, summary


def verify_source_manifest(
    manifest_path: Path | str,
    *,
    rehash: bool = True,
) -> dict[str, object]:
    """Verify the frozen content hash and (optionally) every referenced file."""

    manifest_path = Path(manifest_path)
    if not manifest_path.is_file():
        raise SourceManifestError(f"manifest not found: {manifest_path}")
    data = manifest_path.read_bytes()
    digest = _sha256_bytes(data)
    sidecar = manifest_path.with_suffix(manifest_path.suffix + ".sha256")
    problems: list[str] = []
    if sidecar.is_file():
        frozen = sidecar.read_text(encoding="utf-8").split()[0].strip().lower()
        if frozen != digest:
            problems.append(f"manifest_content_hash_mismatch frozen={frozen} actual={digest}")
    else:
        problems.append("missing_sha256_sidecar")

    import io

    reader = csv.DictReader(io.StringIO(data.decode("utf-8")), delimiter="\t")
    rows = list(reader)
    checked = 0
    for row in rows:
        if not row.get("source_path") or not row.get("sha256"):
            continue
        path = Path(row["source_path"])
        if not path.is_file():
            problems.append(f"{row['code']}/{row['role']}:file_missing:{path}")
            continue
        if str(path.stat().st_size) != row["bytes"]:
            problems.append(f"{row['code']}/{row['role']}:size_mismatch:{path}")
            continue
        if rehash:
            checked += 1
            if sha256_file(path) != row["sha256"]:
                problems.append(f"{row['code']}/{row['role']}:sha256_mismatch:{path}")
    return {
        "manifest": str(manifest_path),
        "manifest_sha256": digest,
        "rows": len(rows),
        "files_rehashed": checked,
        "ok": not problems,
        "problems": problems,
    }
