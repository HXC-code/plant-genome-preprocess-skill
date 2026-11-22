"""Machine and data discovery.

A preprocessing project must not depend on where a previous run happened to
put things.  This module turns paths into a **scanned result**:

* locate the skill/code tree through ``__file__`` (symlinks resolved), never
  through an installed link or the current working directory;
* scan arbitrary user-supplied roots for species data folders and classify
  their files by *content* (annotation vs genome/CDS/PEP FASTA);
* probe the interpreter and external tools already present on the machine;
* write a path map (``resources/machine_paths.json``) plus a pointer next to
  the code, so every later CLI command resolves the same project root without
  an explicit ``--project-root``.

The generated inventory is a *suggestion*: each role picks the best-scoring
candidate, and ``discovery_report.tsv`` records every candidate and score so a
reviewer can override mixed-assembly choices.  Nothing here mutates raw data.
"""

from __future__ import annotations

import json
import os
import platform
import re
import shutil
import sys
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from pgci_preprocess import envcheck
from pgci_preprocess.fasta import open_text
from pgci_preprocess.gff_formats import chromosome_from_seqid
from pgci_preprocess.project_paths import OUTPUT_DIRNAME, output_dirname
from pgci_preprocess.input_router import (
    InputRoutingError,
    _fasta_preview,
    _role_score,
    route_annotation,
)

CONFIG_FILENAME = "machine_paths.json"
POINTER_FILENAME = ".pgci-machine.json"
SCHEMA_VERSION = 1
DEFAULT_MAX_DEPTH = 4

UNUSABLE_SCORE = -500_000

SKIP_DIR_NAMES = frozenset(
    {
        ".git",
        ".hg",
        ".svn",
        "__pycache__",
        ".pytest_cache",
        ".ruff_cache",
        "node_modules",
        ".venv",
        "venv",
        ".wgdi-env",
        "site-packages",
        ".range_download",
        # Derived/output trees: they contain WGDI four-piece sets, which would
        # otherwise be mistaken for raw species bundles when a project root is
        # scanned.
        "wgdi_input",
        "analysis_v2",
        "clean_v2",
        "resources",
        "logs",
        "_pipeline_logs",
        "_pipeline_state",
        "validation",
        "dist",
        "build",
    }
)
SPECIES_DIR_RE = re.compile(r"^(\d+)_([^_()]+)_([^()]+)\((.+)\)$")
ROLE_ORDER = ("gff", "genome", "cds", "pep", "rna", "gbff", "gtf")
FASTA_ROLES = ("genome", "cds", "pep")


@dataclass(frozen=True, slots=True)
class FileRecord:
    path: str
    relpath: str
    role: str
    score: int
    note: str
    size: int

    def as_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "relpath": self.relpath,
            "role": self.role,
            "score": self.score,
            "note": self.note,
            "size": self.size,
        }


@dataclass(slots=True)
class SpeciesDiscovery:
    show_id: str
    directory: str
    name: str
    order: str
    family: str
    species: str
    strict_name: bool
    files: list[FileRecord] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    chromosome_evidence: bool = False
    evidence_detail: str = ""
    sequence_mode: str = "chromosomes"
    chosen: dict[str, FileRecord] = field(default_factory=dict)
    bundle_note: str = ""

    def best(self, role: str) -> FileRecord | None:
        """The role's pick from the consistent bundle, else the best score."""

        if self.chosen:
            return self.chosen.get(role)
        candidates = [item for item in self.files if item.role == role and item.score >= 0]
        if not candidates:
            return None
        return sorted(candidates, key=lambda item: (-item.score, item.relpath))[0]

    def candidates(self, role: str) -> list[FileRecord]:
        return sorted(
            (item for item in self.files if item.role == role and item.score >= 0),
            key=lambda item: (-item.score, item.relpath),
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "show_id": self.show_id,
            "directory": self.directory,
            "name": self.name,
            "order": self.order,
            "family": self.family,
            "species": self.species,
            "strict_name": self.strict_name,
            "files": [item.as_dict() for item in self.files],
            "warnings": list(self.warnings),
            "chromosome_evidence": self.chromosome_evidence,
            "evidence_detail": self.evidence_detail,
            "sequence_mode": self.sequence_mode,
            "chosen": {role: item.relpath for role, item in self.chosen.items()},
            "bundle": self.bundle_note,
        }


# ---------------------------------------------------------------------------
# code tree / pointer
# ---------------------------------------------------------------------------


def find_code_root(start: Path | str = "") -> Path:
    """Locate the ``pipeline/code`` tree via its own file, symlinks resolved.

    An installed skill may be a copy, a junction or a symlink to somewhere
    else.  ``Path(__file__).resolve()`` gives the real source location, so the
    result is stable on every machine regardless of how the skill was wired
    into the agent's skills directory.
    """

    if start:
        base = Path(start).expanduser().resolve()
        candidates: Iterable[Path] = (base, *base.parents)
    else:
        here = Path(__file__).resolve()
        candidates = here.parents
    for candidate in candidates:
        if (candidate / "pyproject.toml").is_file() and (
            candidate / "src" / "pgci_preprocess" / "__init__.py"
        ).is_file():
            return candidate
    for candidate in candidates:
        if (candidate / "src" / "pgci_preprocess" / "__init__.py").is_file():
            return candidate
    raise FileNotFoundError("cannot locate the pgci_preprocess code tree (pyproject.toml + src/)")


def find_skill_root(code_root: Path | str) -> Path | None:
    """Return the skill package root if the code tree lives inside one."""

    base = Path(code_root).resolve()
    for candidate in (base, *base.parents):
        if (candidate / "SKILL.md").is_file() and (candidate / "pipeline" / "code").is_dir():
            return candidate
    return None


def pointer_path(code_root: Path | str | None = None) -> Path:
    root = Path(code_root).resolve() if code_root else find_code_root()
    return root / POINTER_FILENAME


def home_pointer_path() -> Path:
    """Fallback pointer for read-only skill installs."""

    return Path.home() / POINTER_FILENAME


def pointer_candidates(code_root: Path | str | None = None) -> list[Path]:
    """Pointer locations in priority order (``PGCI_POINTER_DIR`` overrides)."""

    override = os.environ.get("PGCI_POINTER_DIR", "").strip()
    if override:
        return [Path(override).expanduser() / POINTER_FILENAME]
    return [pointer_path(code_root), home_pointer_path()]


def write_pointer(project_root: Path | str, config_path: Path | str, code_root: Path | str | None = None) -> Path:
    """Write the tiny pointer the CLI reads to resolve the default project root.

    The code tree is the natural place (it travels with the skill), but an
    installed skill may be read-only, so fall back to the home directory.
    """

    payload = {
        "schema_version": SCHEMA_VERSION,
        "project_root": str(Path(project_root).expanduser().resolve()),
        "config": str(Path(config_path).expanduser().resolve()),
        "scanned_at": _now(),
    }
    text = json.dumps(payload, ensure_ascii=False, indent=2) + "\n"
    errors: list[str] = []
    for target in pointer_candidates(code_root):
        try:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(text, encoding="utf-8")
            return target
        except OSError as exc:
            errors.append(f"{target}: {exc}")
    raise OSError("cannot write the machine pointer: " + "; ".join(errors))


def read_pointer(code_root: Path | str | None = None) -> dict[str, Any] | None:
    for target in pointer_candidates(code_root):
        if not target.is_file():
            continue
        try:
            payload = json.loads(target.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if isinstance(payload, dict) and payload.get("project_root"):
            return payload
    return None


# ---------------------------------------------------------------------------
# data discovery
# ---------------------------------------------------------------------------


def _now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def parse_species_dirname(name: str) -> dict[str, str] | None:
    match = SPECIES_DIR_RE.match(name)
    if not match:
        return None
    return {
        "show_id": match.group(1),
        "order": match.group(2).strip(),
        "family": match.group(3).strip(),
        "species": match.group(4).strip(),
    }


def _first_content_line(path: Path, limit: int = 12) -> str:
    """Bounded peek used to avoid parsing huge FASTA files as annotations."""

    try:
        with open_text(path) as handle:
            for _ in range(limit):
                raw = handle.readline()
                if not raw:
                    break
                line = raw.strip()
                if line and not line.startswith("#"):
                    return line
    except (OSError, EOFError, ValueError):
        return ""
    return ""


def _gff_quality_bonus(path: Path, content_format: str = "standard") -> tuple[int, str]:
    """Rank coding-gene annotations above repeats, ncRNA and organelle stubs.

    Several bundles ship several GFFs with identical parsing scores: a small
    GeSeq chloroplast stub, a RepeatMasker track, and the real nuclear gene
    annotation.  Picking by filename order alone silently selects the wrong
    one, so we inspect the feature composition instead.

    ``content_format`` matters: a five-column simplified table has coordinates
    in the third column, so counting it as a feature type would wrongly report
    "no gene/CDS features".
    """

    name = path.name.casefold()
    if re.search(r"(?:repeat|ncrna|masked|rrna|trna)", name):
        return -20, "repeat/ncRNA track"
    counts: dict[str, int] = {}
    organelle = False
    rows = 0
    try:
        with open_text(path) as handle:
            for index, line in enumerate(handle):
                if index >= 4000:
                    break
                if line.startswith("#"):
                    if re.search(r"chloroplast|plastid|mitochondrion", line, re.IGNORECASE):
                        organelle = True
                    continue
                if not line.strip():
                    continue
                parts = line.rstrip("\r\n").split("\t")
                if content_format == "simplified5":
                    if len(parts) >= 5:
                        rows += 1
                    continue
                if len(parts) >= 3:
                    if re.search(r"GeSeq", parts[1], re.IGNORECASE):
                        organelle = True
                    counts[parts[2]] = counts.get(parts[2], 0) + 1
    except (OSError, EOFError, ValueError):
        return 0, ""
    if organelle:
        return -15, "organelle annotation"
    if content_format == "simplified5":
        if rows == 0:
            return -10, "no gene models"
        return min(6, rows // 2000), f"gene_models={rows}"
    coding = sum(counts.get(feature, 0) for feature in ("gene", "mRNA", "transcript", "CDS", "exon"))
    if coding == 0:
        return -10, "no gene/CDS features"
    return min(8, coding // 50), f"coding_features={coding}"


#: Resources that annotate genes but never describe an assembly.
REFERENCE_NAME_RE = re.compile(
    r"(?:swissprot|uniprot|trembl|pfam|interpro|eggnog|kegg|orthofinder|busco"
    r"|single[_-]?copy|multialignment|multiple[_-]?alignment|blast[_-]?(?:db|p)|go[_-]?annotation)"
)

#: A sequence whose own header names an organelle genome.
ORGANELLE_HEADER_RE = re.compile(r"chloroplast|plastid|mitochond", re.IGNORECASE)
ORGANELLE_NAME_RE = re.compile(
    r"chloroplast|plastid|mitochond|(?:^|[._-])(?:cp|mt)dna(?:[._-]|$)", re.IGNORECASE
)


def _is_organelle_fasta(name: str, headers: list[str]) -> bool:
    """Distinguish an organelle assembly from a nuclear gene set.

    Nuclear CDS/PEP headers routinely mention "chloroplastic" or
    "mitochondrial" in *product descriptions* (e.g. "biotin carboxylase 1,
    chloroplastic").  Treating any such hit as organelle evidence throws away
    real annotation.  A genuine organelle FASTA labels essentially every
    sequence, so we require a majority of the sampled headers.
    """

    if ORGANELLE_NAME_RE.search(name):
        return True
    if not headers:
        return False
    hits = sum(1 for header in headers if ORGANELLE_HEADER_RE.search(header))
    return hits / len(headers) >= 0.6


def classify_file(path: Path, base: Path) -> FileRecord:
    """Classify one file by content; name is only a hint for FASTA roles."""

    resolved = path.resolve()
    try:
        rel = resolved.relative_to(base.resolve()).as_posix()
    except ValueError:
        rel = resolved.name
    first = _first_content_line(resolved)
    if first.startswith(">"):
        headers = _fasta_headers(resolved, limit=40)
        try:
            kind, _first_id = _fasta_preview(resolved)
        except InputRoutingError as exc:
            return FileRecord(str(resolved), rel, "other", -100, str(exc), resolved.stat().st_size)
        lowered = resolved.name.casefold()
        # Reference databases and multiple-sequence alignments live next to
        # the annotation but are never assembly inputs.
        if REFERENCE_NAME_RE.search(lowered):
            return FileRecord(
                str(resolved), rel, "reference", 0,
                "reference database / alignment (not an assembly input)",
                resolved.stat().st_size,
            )
        if _is_organelle_fasta(lowered, headers):
            return FileRecord(
                str(resolved), rel, "organelle", -100, "organelle genome", resolved.stat().st_size
            )
        # RNA/transcript FASTA is a nucleotide FASTA and would otherwise win
        # the generic genome score; its headers may even mention
        # "chromosome X-like", which would fake chromosome evidence.
        if (
            re.search(r"(?:rna|transcript)", lowered)
            and not re.search(r"(?:^|[._-])cds(?:[._-]|$)", lowered)
            and kind in {"nucleotide_fasta", "fasta"}
        ):
            return FileRecord(str(resolved), rel, "rna", 6, kind, resolved.stat().st_size)
        ranked = {role: _role_score(resolved, role) for role in FASTA_ROLES}
        role = max(ranked, key=lambda item: ranked[item][0])
        score, detail = ranked[role]
        if score < 0:
            lowered = resolved.name.casefold()
            if re.search(r"(?:rna|transcript|trna|rrna)", lowered) and kind in {
                "nucleotide_fasta",
                "fasta",
            }:
                return FileRecord(str(resolved), rel, "rna", 2, kind, resolved.stat().st_size)
            return FileRecord(str(resolved), rel, "other", -50, kind, resolved.stat().st_size)
        return FileRecord(str(resolved), rel, role, score, f"{kind}; {detail}", resolved.stat().st_size)

    try:
        route = route_annotation(resolved)
    except InputRoutingError as exc:
        return FileRecord(str(resolved), rel, "other", -100, str(exc)[:200], resolved.stat().st_size)
    if route.content_format == "functional-go":
        return FileRecord(
            str(resolved), rel, "metadata", 0, "functional table (not a coordinate input)",
            resolved.stat().st_size,
        )
    if route.content_format in {"gtf", "gbff"}:
        return FileRecord(
            str(resolved), rel, route.content_format, 8,
            f"convert first ({route.note or route.content_format})", resolved.stat().st_size,
        )
    bonus, quality = _gff_quality_bonus(resolved, route.content_format)
    detail = route.content_format + (f"; {route.note}" if route.note else "")
    if quality:
        detail += f"; {quality}"
    # A five-column table is a legitimate but lossy annotation: rank it just
    # below a full GFF that carries gene/mRNA/CDS structure.
    score_base = 8 if route.content_format == "simplified5" else 10
    return FileRecord(
        str(resolved), rel, "gff", score_base + bonus, detail, resolved.stat().st_size
    )


# Require a real chromosome number/label.  NCBI RNA headers such as
# "chromosome X-like" describe a *similarity*, not a mounted chromosome.
CHROMOSOME_LABEL_RE = re.compile(
    r"(?:^|\s)Chromosome\s+(\d{1,3}|X|Y|Z|MT|M)\b(?![\w-])", re.IGNORECASE
)


#: Tokens that describe a *file role*, not an assembly version.
ROLE_TOKEN_RE = re.compile(
    r"[._-](?:cds|cdna|pep|peptides?|proteins?|transcripts?|transcriptome|genes?|gene"
    r"|genome|assembly|annotation|annotations|normalized|sorted|softmasked|hardmasked"
    r"|masked|repeat|gff3?|gtf|gbff|faa|fna|fas|fasta|fa)$"
)
ASSEMBLY_FORMAT_RE = re.compile(r"\.(?:gz|bz2|xz|zip)$")


def _assembly_key(relpath: str) -> str:
    """Coarse assembly identity used to keep a bundle on one genome version.

    ``ppa_v2.asm.genes.gff`` and ``ppa_v2.asm.CDS.fasta`` must share a key,
    while ``ppa_v1.asm.transcriptome.gtf`` must not join that bundle.  Role and
    compression tokens are stripped, then trailing role words repeatedly.
    """

    name = Path(relpath).name.casefold()
    match = re.match(r"^(gc[af]_\d+\.\d+)", name)
    if match:
        return match.group(1)
    match = re.match(r"^(gwh[a-z0-9]+?)(?:[._]|$)", name)
    if match:
        return match.group(1)
    stem = ASSEMBLY_FORMAT_RE.sub("", name)
    while True:
        reduced = ROLE_TOKEN_RE.sub("", stem)
        if reduced == stem:
            break
        stem = reduced
    return stem or name


_ID_CACHE: dict[tuple[str, str], frozenset[str]] = {}


def _normalise_id(token: str) -> str:
    """Fold common gene/transcript ID spellings onto one comparable form."""

    value = token.strip().strip('"').strip("'")
    if not value:
        return ""
    value = value.split("|")[-1].split()[0] if value.split("|")[-1].split() else ""
    if not value:
        return ""
    # NCBI ``*_cds_from_genomic`` headers embed the protein: lcl|NW_x_cds_XP_1.1_2
    if "_cds_" in value:
        value = value.split("_cds_", 1)[-1]
        value = re.sub(r"_\d+$", "", value)
    # Drop feature-type prefixes: gene-, cds-, rna-, mRNA:, ID:gene:
    value = re.sub(
        r"^(?:gene|mrna|transcript|cds|rna|protein|exon|pseudogene)[-:]", "", value,
        flags=re.IGNORECASE,
    )
    # Drop isoform suffixes: AT1G01010.1, XP_001.2, NPAL_000001_R1, gene-T1.
    value = re.sub(r"(?:\.\d+|-\d+|[._-][RTP]\d+)$", "", value)
    return value.casefold()


def _gff_seqids_of(path: Path, limit: int = 5000) -> set[str]:
    """Sequence IDs from the first column of an annotation (never sequence)."""

    key = (str(path), "gff_seqid")
    cached = _ID_CACHE.get(key)
    if cached is not None:
        return set(cached)
    found: set[str] = set()
    try:
        with open_text(path) as handle:
            for index, line in enumerate(handle):
                if index >= limit:
                    break
                if line.startswith("#") or not line.strip():
                    continue
                normalised = _normalise_id(line.split("\t", 1)[0])
                if normalised:
                    found.add(normalised)
    except (OSError, EOFError, ValueError):
        return set()
    _ID_CACHE[key] = frozenset(found)
    return found


def _fasta_seqids_of(path: Path, limit: int = 5000) -> set[str]:
    """Sequence names from FASTA headers only; sequence lines are ignored."""

    key = (str(path), "fasta_seqid")
    cached = _ID_CACHE.get(key)
    if cached is not None:
        return set(cached)
    found: set[str] = set()
    try:
        with open_text(path) as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                normalised = _normalise_id(line[1:].strip().split(None, 1)[0])
                if normalised:
                    found.add(normalised)
                if len(found) >= limit:
                    break
    except (OSError, EOFError, ValueError):
        return set()
    _ID_CACHE[key] = frozenset(found)
    return found


_ATTR_ID_RE = re.compile(r"(?:^|;)\s*(?:ID|Parent|gene_id|transcript_id|protein_id)=([^;]+)")


def _feature_ids_of(path: Path, limit: int = 20000) -> set[str]:
    """Collect ID/Parent values from a GFF-like annotation."""

    key = (str(path), "feature_id")
    cached = _ID_CACHE.get(key)
    if cached is not None:
        return set(cached)
    found: set[str] = set()
    try:
        with open_text(path) as handle:
            for index, line in enumerate(handle):
                if index >= limit:
                    break
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < 8:
                    continue
                for match in _ATTR_ID_RE.finditer(parts[8]):
                    for token in match.group(1).split(","):
                        normalised = _normalise_id(token)
                        if normalised:
                            found.add(normalised)
    except (OSError, EOFError, ValueError):
        return set()
    _ID_CACHE[key] = frozenset(found)
    return found


def _overlap(left: set[str], right: set[str]) -> float:
    if not left or not right:
        return 0.0
    # Compare on a bounded sample: annotation files can hold far more IDs than
    # a FASTA preview, and dividing by the larger side would dilute a genuine
    # match into noise.
    limit = min(len(left), len(right), 5000)
    return min(1.0, len(left & right) / limit)


def _collection_ids(path: Path, limit: int = 5000) -> set[str]:
    """IDs carried by a CDS/PEP/RNA FASTA (headers only, never sequence)."""

    key = (str(path), "collection_id")
    cached = _ID_CACHE.get(key)
    if cached is not None:
        return set(cached)
    found: set[str] = set()
    try:
        with open_text(path) as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                normalised = _normalise_id(line[1:].strip().split(None, 1)[0])
                if normalised:
                    found.add(normalised)
                if len(found) >= limit:
                    break
    except (OSError, EOFError, ValueError):
        return set()
    _ID_CACHE[key] = frozenset(found)
    return found


def select_bundle(files: list[FileRecord]) -> tuple[dict[str, FileRecord], str]:
    """Pair annotation with the genome/CDS/PEP of the *same* assembly.

    Choosing the highest score per role independently can splice unrelated
    assemblies together (a nuclear GFF with a chloroplast genome, or a GCF
    annotation with a manually annotated genome).  The annotation is the
    anchor: every partner must share sequence IDs (genome) or gene/transcript
    IDs (CDS/PEP) with it.  The best-covered, highest-scoring consistent tuple
    wins; anything left unmatched is simply omitted and reported.
    """

    anchors = sorted(
        (item for item in files if item.role in {"gff", "gtf", "gbff"} and item.score >= 0),
        key=lambda item: (-item.score, item.relpath),
    )[:4]
    if not anchors:
        return {}, "no annotation candidate"

    best_total = float("-inf")
    best_choice: dict[str, FileRecord] = {}
    best_note = ""
    for anchor in anchors:
        seqids = _gff_seqids_of(Path(anchor.path))
        feature_ids = (
            _feature_ids_of(Path(anchor.path)) if anchor.role == "gff" else set()
        )
        key = _assembly_key(anchor.relpath)
        # Keep the anchor's real role: a GBFF-only bundle must land in the
        # `gbff` inventory column so `prepare` converts it, not parse it as GFF.
        choice: dict[str, FileRecord] = {anchor.role: anchor}
        notes = [f"anchor={anchor.relpath}"]
        for role in ("genome", "cds", "pep", "rna"):
            candidates = [item for item in files if item.role == role and item.score >= 0]
            if not candidates:
                continue
            scored: list[tuple[float, FileRecord, str]] = []
            for item in candidates:
                score = float(item.score)
                extra = ""
                name_match = _assembly_key(item.relpath) == key
                if name_match:
                    score += 10
                    extra += "; name-match"
                if role == "genome":
                    ratio = _overlap(seqids, _fasta_seqids_of(Path(item.path)))
                    if ratio:
                        score += 30 + 40 * ratio
                        extra += f"; seqid-overlap={ratio:.0%}"
                    else:
                        # A genome whose sequence names do not appear in the
                        # annotation cannot host its coordinates, whatever the
                        # filename says.  Mark it unusable rather than pairing
                        # two unrelated assemblies.
                        score = UNUSABLE_SCORE
                        extra += "; seqid-mismatch"
                else:
                    ratio = _overlap(feature_ids, _collection_ids(Path(item.path)))
                    if ratio:
                        score += 30 + 40 * ratio
                        extra += f"; id-overlap={ratio:.0%}"
                    elif name_match:
                        # Different ID namespace (e.g. VIGUN gene models vs
                        # GWH transcript IDs); the filename is the only link
                        # available, so accept it without punishing.
                        extra += "; ids-unverifiable"
                    elif role == "rna":
                        # RNA is optional.  A transcriptome from another
                        # assembly would inject wrong isoforms, so leave it
                        # out instead of guessing.
                        score = UNUSABLE_SCORE
                        extra += "; unrelated transcriptome"
                    else:
                        score -= 25
                        extra += "; id-mismatch"
                scored.append((score, item, extra))
            scored.sort(key=lambda entry: (-entry[0], entry[1].relpath))
            score, item, extra = scored[0]
            if score <= UNUSABLE_SCORE:
                notes.append(f"{role}=UNPAIRED{extra}")
                continue
            choice[role] = replace(item, score=int(round(score)), note=item.note + extra)
            notes.append(f"{role}={item.relpath}{extra}")
        total = sum(item.score for item in choice.values()) + 12 * (len(choice) - 1)
        if total > best_total:
            best_total = total
            best_choice = choice
            best_note = "; ".join(notes)
    return best_choice, best_note


def _fasta_headers(path: Path, limit: int = 1000) -> list[str]:
    headers: list[str] = []
    try:
        with open_text(path) as handle:
            for line in handle:
                if line.startswith(">"):
                    headers.append(line[1:].strip())
                    if len(headers) >= limit:
                        break
    except (OSError, EOFError, ValueError):
        return []
    return headers


def detect_chromosome_evidence(files: list[FileRecord]) -> tuple[bool, str]:
    """Decide from content whether a bundle has chromosome-mounted sequences.

    This is the data-dependent half of the strict chromosome mode.  The docker
    fix made the pipeline recognise ``Chr1``/``Ajul.Chr01``; here the scan
    decides, for *this* dataset on *this* machine, whether such evidence
    exists at all.  Scaffold-only assemblies get ``adaptive`` suggested, and
    the decision is written to ``resources/sequence_mode.tsv`` for review.
    """

    for record in sorted(
        (item for item in files if item.role == "genome"),
        key=lambda item: -item.score,
    )[:3]:
        for header in _fasta_headers(Path(record.path)):
            first = header.split(None, 1)[0] if header else ""
            if CHROMOSOME_LABEL_RE.search(header) or chromosome_from_seqid(first):
                return True, f"genome header: {header[:120]}"
    for record in sorted(
        (item for item in files if item.role == "gff"),
        key=lambda item: -item.score,
    )[:3]:
        try:
            with open_text(Path(record.path)) as handle:
                for _ in range(5000):
                    line = handle.readline()
                    if not line:
                        break
                    if line.startswith("##sequence-region"):
                        parts = line.split()
                        seqid = parts[1] if len(parts) > 1 else ""
                    elif line.startswith("#") or not line.strip():
                        continue
                    else:
                        seqid = line.split("\t", 1)[0]
                    if chromosome_from_seqid(seqid):
                        return True, f"GFF seqid: {seqid}"
        except (OSError, EOFError, ValueError):
            continue
    return False, "no chromosome label in genome headers or GFF seqids"


def _iter_files(directory: Path, max_depth: int) -> Iterable[Path]:
    base_depth = len(directory.resolve().parts)
    for root, dirs, files in os.walk(directory):
        dirs[:] = [item for item in dirs if item not in SKIP_DIR_NAMES]
        current = Path(root)
        if len(current.parts) - base_depth >= max_depth:
            dirs[:] = []
        for name in files:
            if name.startswith("."):
                continue
            candidate = current / name
            if candidate.is_file() and candidate.stat().st_size > 0:
                yield candidate


def _looks_like_species(files: list[FileRecord]) -> bool:
    has_annotation = any(item.role in {"gff", "gtf", "gbff"} for item in files)
    has_fasta = any(item.role in FASTA_ROLES for item in files)
    return has_annotation and has_fasta


def _has_species_descendant(directory: Path, max_depth: int = 3) -> bool:
    """True when a sub-folder already is a species bundle.

    A data root such as ``project/data/137_...(Albizia)`` also "looks like a
    species" when judged by its recursive files, which would swallow every
    child into one giant fake species.  A container must never be mistaken for
    a bundle.
    """

    base_depth = len(directory.resolve().parts)
    try:
        for root, dirs, _files in os.walk(directory):
            current = Path(root)
            depth = len(current.parts) - base_depth
            if depth >= max_depth:
                dirs[:] = []
                continue
            dirs[:] = [item for item in dirs if item not in SKIP_DIR_NAMES]
            for name in dirs:
                if parse_species_dirname(name):
                    return True
    except OSError:
        return False
    return False


def find_species_dirs(
    scan_roots: Iterable[Path | str],
    *,
    max_depth: int = DEFAULT_MAX_DEPTH,
) -> list[SpeciesDiscovery]:
    """Scan roots for species bundles, strict names first, content fallback.

    ``max_depth`` is measured *inside each candidate directory*, so a data
    root with ``<species>/figshare_download/<file>`` is covered by the default
    value of four without walking an entire drive.
    """

    seen_roots: dict[Path, None] = {}
    for raw in scan_roots:
        path = Path(raw).expanduser().resolve()
        if path.is_dir():
            seen_roots[path] = None
        elif path.is_file():
            seen_roots[path.parent] = None

    candidates: list[tuple[Path, dict[str, str] | None, list[FileRecord]]] = []
    for root in seen_roots:
        # The root itself may be a species bundle.  Judge it by its *direct*
        # files only: a plain data root must not swallow its children.
        parsed = parse_species_dirname(root.name)
        direct = [
            classify_file(item, root)
            for item in root.iterdir()
            if item.is_file() and item.stat().st_size > 0
        ]
        if parsed or (
            not _has_species_descendant(root) and _looks_like_species(direct)
        ):
            files = [classify_file(item, root) for item in _iter_files(root, max_depth)]
            candidates.append((root, parsed, files))
        for child in sorted(root.iterdir()):
            if not child.is_dir() or child.name in SKIP_DIR_NAMES:
                continue
            child_parsed = parse_species_dirname(child.name)
            child_files = [classify_file(item, child) for item in _iter_files(child, max_depth)]
            if child_parsed or (
                not _has_species_descendant(child) and _looks_like_species(child_files)
            ):
                candidates.append((child, child_parsed, child_files))
                continue
            # A container such as ``<project>/data`` is one level deeper than
            # the species folders; registering only the container would hide
            # every species when the project root is scanned instead of data/.
            for grandchild in sorted(child.iterdir()):
                if not grandchild.is_dir() or grandchild.name in SKIP_DIR_NAMES:
                    continue
                grandchild_parsed = parse_species_dirname(grandchild.name)
                if not grandchild_parsed:
                    continue
                candidates.append(
                    (
                        grandchild,
                        grandchild_parsed,
                        [
                            classify_file(item, grandchild)
                            for item in _iter_files(grandchild, max_depth)
                        ],
                    )
                )

    ordered = sorted(candidates, key=lambda item: len(item[0].parts))
    accepted: list[tuple[Path, dict[str, str] | None, list[FileRecord]]] = []
    for directory, parsed, files in ordered:
        if any(directory == kept[0] or directory.is_relative_to(kept[0]) for kept in accepted):
            continue
        accepted.append((directory, parsed, files))

    discoveries: list[SpeciesDiscovery] = []
    for directory, parsed, files in accepted:
        if parsed:
            discovery = SpeciesDiscovery(
                show_id=parsed["show_id"],
                directory=str(directory),
                name=directory.name,
                order=parsed["order"],
                family=parsed["family"],
                species=parsed["species"],
                strict_name=True,
                files=files,
            )
        else:
            # Content-only candidate: derive a synthetic show_id from a
            # numeric prefix when present, so it can still be registered.
            match = re.match(r"^(\d+)", directory.name)
            show_id = match.group(1) if match else ""
            discovery = SpeciesDiscovery(
                show_id=show_id,
                directory=str(directory),
                name=directory.name,
                order="",
                family="",
                species=directory.name,
                strict_name=False,
                files=files,
                warnings=["directory name does not follow <show_id>_<Order>_<Family>(<Species>)"],
            )
        chosen, bundle_note = select_bundle(files)
        discovery.chosen = chosen
        discovery.bundle_note = bundle_note
        has_chromosomes, detail = detect_chromosome_evidence(
            list(chosen.values()) or files
        )
        discovery.chromosome_evidence = has_chromosomes
        discovery.evidence_detail = detail
        discovery.sequence_mode = "chromosomes" if has_chromosomes else "adaptive"
        if not has_chromosomes:
            discovery.warnings.append(
                f"no chromosome evidence ({detail}); suggested sequence_mode=adaptive"
            )
        annotation = chosen.get("gff") or chosen.get("gtf") or chosen.get("gbff")
        if annotation is None:
            discovery.warnings.append("no annotation candidate (gff/gtf/gbff)")
        for role in ("genome", "cds", "pep"):
            if discovery.best(role) is None:
                discovery.warnings.append(f"no paired {role} candidate")
        discoveries.append(discovery)
    return sorted(discoveries, key=lambda item: (not item.show_id.isdigit(), item.show_id, item.name))


def infer_project_root(scan_roots: Iterable[Path | str]) -> Path:
    """Best-effort project root from the scanned locations."""

    roots = [Path(item).expanduser().resolve() for item in scan_roots]
    for root in roots:
        for candidate in (root, *root.parents):
            if (candidate / "resources" / "download_inventory.tsv").is_file():
                return candidate
    for root in roots:
        if root.name.casefold() == "data" and root.parent.is_dir():
            return root.parent
    if len(roots) == 1:
        return roots[0]
    if roots:
        common = Path(os.path.commonpath([str(item) for item in roots]))
        return common if common.is_dir() else roots[0]
    return Path.cwd()


def infer_data_root(species: Iterable[SpeciesDiscovery]) -> str:
    """Directory that directly holds the species folders (may be off-project)."""

    directories = [Path(item.directory) for item in species]
    if not directories:
        return ""
    parents = {item.parent for item in directories}
    if len(parents) == 1:
        return str(parents.pop())
    try:
        return str(Path(os.path.commonpath([str(item) for item in directories])))
    except ValueError:
        return str(directories[0].parent)


def default_scan_roots() -> list[Path]:
    """Where to look when the user gives no paths on a fresh machine.

    ``PGCI_SCAN_ROOTS`` (path-separator separated) wins; otherwise the current
    directory plus its conventional data sub-folders are scanned.  Scanning a
    whole drive is deliberately not the default: on a new machine the user
    knows which folder holds the data.
    """

    configured = os.environ.get("PGCI_SCAN_ROOTS", "").strip()
    if configured:
        return [Path(item).expanduser() for item in configured.split(os.pathsep) if item.strip()]
    cwd = Path.cwd()
    roots = [cwd]
    seen = {os.path.normcase(str(cwd))}
    for name in ("data", "DATA", "Data", "raw", "rawdata", "genomes", "downloads"):
        candidate = cwd / name
        key = os.path.normcase(str(candidate.resolve()))
        if candidate.is_dir() and key not in seen:
            seen.add(key)
            roots.append(candidate)
    return roots


# ---------------------------------------------------------------------------
# environment discovery
# ---------------------------------------------------------------------------


def detect_interpreters() -> list[dict[str, Any]]:
    candidates: list[str] = [sys.executable]
    for name in ("python", "python3", "py"):
        resolved = shutil.which(name)
        if resolved:
            candidates.append(resolved)
    extra = os.environ.get("PGCI_PYTHON_CANDIDATES", "")
    candidates.extend(item.strip() for item in extra.split(os.pathsep) if item.strip())
    seen: set[str] = set()
    result: list[dict[str, Any]] = []
    for candidate in candidates:
        try:
            resolved = str(Path(candidate).resolve())
        except OSError:
            continue
        if resolved in seen or not Path(resolved).is_file():
            continue
        seen.add(resolved)
        info = envcheck.interpreter_info(resolved)
        if not info.get("reachable"):
            continue
        info["packages"] = envcheck.installed_packages(resolved)
        result.append(info)
    return result


def detect_tools() -> dict[str, str]:
    tools = envcheck.external_tools()
    for tool in tools:
        override = os.environ.get(f"PGCI_TOOL_{tool.upper()}", "").strip()
        if override and Path(override).expanduser().is_file():
            tools[tool] = str(Path(override).expanduser().resolve())
    return tools


# ---------------------------------------------------------------------------
# config / inventory output
# ---------------------------------------------------------------------------


def suggested_inventory(species: list[SpeciesDiscovery]) -> list[dict[str, str]]:
    """One best-scoring path per role; full candidate list goes to the report."""

    fields = (
        "showid", "species_code", "Order", "Family", "Species",
        "gff", "genome", "cds", "pep", "rna", "gbff", "gtf",
    )
    rows: list[dict[str, str]] = []
    for item in species:
        if not item.show_id:
            continue
        row = {field: "" for field in fields}
        row.update(
            {
                "showid": item.show_id,
                "Order": item.order,
                "Family": item.family,
                "Species": item.species,
            }
        )
        # A real GFF wins; GTF/GBFF are valid fallbacks because prepare converts
        # them.  The discovery report still lists every candidate.  Only one
        # annotation column is filled, so the inventory never asks `prepare` to
        # merge two competing annotations of the same assembly.
        annotation = item.best("gff") or item.best("gtf") or item.best("gbff")
        for role, column in (
            ("gff", "gff"), ("gtf", "gtf"), ("gbff", "gbff"),
            ("genome", "genome"), ("cds", "cds"), ("pep", "pep"), ("rna", "rna"),
        ):
            if role in {"gff", "gtf", "gbff"}:
                if annotation is not None and annotation.role == role:
                    row[column] = annotation.relpath
                continue
            best = item.best(role)
            if best is not None:
                row[column] = best.relpath
        rows.append(row)
    return rows


def build_config(
    scan_roots: Iterable[Path | str],
    *,
    project_root: Path | str | None = None,
    max_depth: int = DEFAULT_MAX_DEPTH,
    probe_index: bool = False,
    index_url: str = "",
) -> dict[str, Any]:
    roots: list[Path] = []
    seen_roots: set[str] = set()
    for item in scan_roots:
        resolved = Path(item).expanduser().resolve()
        key = os.path.normcase(str(resolved))
        if key in seen_roots:
            continue
        seen_roots.add(key)
        roots.append(resolved)
    root = Path(project_root).expanduser().resolve() if project_root else infer_project_root(roots)
    code_root = find_code_root()
    skill_root = find_skill_root(code_root)
    species = find_species_dirs(roots, max_depth=max_depth)
    network = (
        envcheck.probe_network(index_url)
        if probe_index
        else {"reachable": False, "skipped": True, "index": "(skipped)"}
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at": _now(),
        "machine": {
            "platform": platform.platform(),
            "python": sys.executable,
            "python_version": ".".join(map(str, sys.version_info[:3])),
            "hostname": platform.node(),
        },
        "paths": {
            "code_root": str(code_root),
            "skill_root": str(skill_root) if skill_root else "",
            "project_root": str(root),
            "data_root": infer_data_root(species),
            "output_dirname": output_dirname(root),
            "scan_roots": [str(item) for item in roots],
            "writable": os.access(root, os.W_OK),
        },
        "environment": {
            "interpreters": detect_interpreters(),
            "tools": detect_tools(),
            "network": network,
        },
        "species": [item.as_dict() for item in species],
        "suggested_inventory": suggested_inventory(species),
        "warnings": [warning for item in species for warning in item.warnings],
    }


def write_config(config: dict[str, Any], path: Path | str) -> Path:
    target = Path(path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(config, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return target


def load_config(path: Path | str) -> dict[str, Any]:
    target = Path(path).expanduser().resolve()
    payload = json.loads(target.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"unsupported path-map schema in {target}")
    return payload


def write_inventory(
    config: dict[str, Any], project_root: Path | str, *, overwrite: bool = False
) -> Path:
    """Write the suggested inventory, never clobbering a reviewed one by default."""

    root = Path(project_root).expanduser().resolve()
    target = root / "resources" / "download_inventory.tsv"
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and not overwrite:
        # A header-only file is what `init-project` leaves behind; filling it
        # is not "clobbering a reviewed inventory".
        try:
            rows = [
                line for line in target.read_text(encoding="utf-8").splitlines() if line.strip()
            ]
        except OSError:
            rows = []
        if len(rows) > 1:
            return target
    fields = (
        "showid", "species_code", "Order", "Family", "Species",
        "gff", "genome", "cds", "pep", "rna", "gbff", "gtf",
    )
    lines = ["\t".join(fields)]
    for row in config.get("suggested_inventory", []):
        lines.append("\t".join(str(row.get(field, "")) for field in fields))
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def seed_resources(project_root: Path | str) -> list[str]:
    """Copy the lookup tables every project needs but nobody wants to write.

    ``init-project`` does this when it creates a project; a scanned project is
    created by ``discover`` instead, so it must seed the same files — without
    the abbreviation table the run would fall back to legacy ``s<id>`` names.
    """

    root = Path(project_root).expanduser().resolve()
    resources = root / "resources"
    resources.mkdir(parents=True, exist_ok=True)
    bundled = Path(__file__).resolve().parents[2] / "resources"
    created: list[str] = []
    for name in ("species_abbreviations.tsv", "format_exceptions.tsv"):
        source = bundled / name
        destination = resources / name
        if source.is_file() and not destination.exists():
            shutil.copy2(source, destination)
            created.append(str(destination))
    return created


def write_discovery_report(config: dict[str, Any], project_root: Path | str) -> Path:
    root = Path(project_root).expanduser().resolve()
    target = root / "resources" / "discovery_report.tsv"
    target.parent.mkdir(parents=True, exist_ok=True)
    fields = ("show_id", "role", "score", "relpath", "format_or_kind", "size")
    lines = ["\t".join(fields)]
    for item in config.get("species", []):
        for record in item.get("files", []):
            if record.get("role") in {"other", "metadata"} or record.get("score", -1) < 0:
                continue
            lines.append(
                "\t".join(
                    str(value)
                    for value in (
                        item.get("show_id", ""),
                        record.get("role", ""),
                        record.get("score", ""),
                        record.get("relpath", ""),
                        record.get("note", ""),
                        record.get("size", ""),
                    )
                )
            )
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def write_sequence_mode_map(config: dict[str, Any], project_root: Path | str) -> Path:
    """Write the per-species sequence mode so `run --sequence-mode from-map` works."""

    root = Path(project_root).expanduser().resolve()
    target = root / "resources" / "sequence_mode.tsv"
    target.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(("show_id", "sequence_mode", "evidence"))]
    for item in config.get("species", []):
        if not item.get("show_id"):
            continue
        lines.append(
            "\t".join(
                str(value)
                for value in (
                    item["show_id"],
                    item.get("sequence_mode", "chromosomes"),
                    item.get("evidence_detail", "").replace("\t", " "),
                )
            )
        )
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def read_sequence_mode_map(project_root: Path | str) -> dict[str, str]:
    """Return show_id -> mode; missing/invalid rows fall back to strict mode."""

    target = Path(project_root).expanduser().resolve() / "resources" / "sequence_mode.tsv"
    if not target.is_file():
        return {}
    result: dict[str, str] = {}
    with target.open(encoding="utf-8", newline="") as handle:
        import csv

        for row in csv.DictReader(handle, delimiter="\t"):
            show_id = str(row.get("show_id", "")).strip()
            mode = str(row.get("sequence_mode", "")).strip()
            if show_id and mode in {"chromosomes", "adaptive"}:
                result[show_id] = mode
    return result


def render(config: dict[str, Any]) -> str:
    paths = config["paths"]
    lines = [
        "== 机器与数据扫描 ==",
        f"代码根   : {paths['code_root']}",
        f"Skill 根 : {paths['skill_root'] or '(独立代码仓，无 SKILL.md)'}",
        f"项目根   : {paths['project_root']}",
        f"数据根   : {paths.get('data_root') or '(未发现物种目录)'}",
        f"产物目录 : {paths.get('output_dirname') or OUTPUT_DIRNAME}",
        f"扫描位置 : {', '.join(paths['scan_roots']) or '(无)'}",
        f"物种目录 : {len(config['species'])} 个",
    ]
    for item in config["species"]:
        chosen = item.get("chosen") or {}
        annotation = chosen.get("gff") or chosen.get("gtf") or chosen.get("gbff")
        parts = [f"annot={annotation or 'MISSING'}"]
        for role in ("genome", "cds", "pep"):
            parts.append(f"{role}={chosen.get(role, 'MISSING')}")
        mode = item.get("sequence_mode", "chromosomes")
        flag = "" if item.get("chromosome_evidence") else "  [无染色体证据→adaptive]"
        lines.append(f"  [{item['show_id'] or '?'}] {item['name']}  mode={mode}{flag}")
        lines.append("      " + "  ".join(parts))
    tools = config["environment"]["tools"]
    lines.append(
        "外部工具 : " + ", ".join(f"{name}{'✓' if path else '✗'}" for name, path in tools.items())
    )
    pythons = config["environment"]["interpreters"]
    lines.append(
        "解释器   : "
        + (", ".join(f"{item['version']} ({item['executable']})" for item in pythons) or "(无)")
    )
    if config["warnings"]:
        lines.append("提示:")
        lines.extend(f"  - {warning}" for warning in config["warnings"][:20])
    adaptive = [
        item["show_id"]
        for item in config["species"]
        if item.get("sequence_mode") == "adaptive" and item.get("show_id")
    ]
    if adaptive:
        lines.append(
            f"无染色体证据物种: {', '.join(adaptive)}（show_id）→ 已写入 "
            "resources/sequence_mode.tsv，用 run --sequence-mode from-map 逐物种应用。"
        )
    unnamed = [item["name"] for item in config["species"] if not item.get("show_id")]
    if unnamed:
        lines.append(
            "未写入清单（目录名缺少 <show_id>_<Order>_<Family>(<Species>)）: "
            + ", ".join(unnamed[:10])
            + (" ..." if len(unnamed) > 10 else "")
            + "；改名后重跑 discover --force。"
        )
    lines.append(
        "下一步   : 复核 resources/discovery_report.tsv 中的候选，"
        "必要时修改 resources/download_inventory.tsv 后运行 run。"
    )
    return "\n".join(lines)
