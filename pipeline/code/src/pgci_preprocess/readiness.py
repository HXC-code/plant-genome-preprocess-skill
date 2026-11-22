"""Single authoritative source for a species' WGDI readiness.

The readiness gate previously had three independent implementations: the
Python validator parsed QC JSON, the shell homology runner matched JSON text
with ``grep`` (fragile), and the PowerShell orchestrator parsed JSON separately.
This module is the one canonical reader; the validator, the ``readiness`` CLI
(used by the shell runner) and the PowerShell orchestrator all derive their
answer from it, so the gate cannot drift between entry points.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from pgci_preprocess.project_paths import output_root

VALID_READINESS = frozenset({"yes", "marginal", "no"})


def extract_qc_readiness(qc: dict[str, object]) -> str:
    """Return the validated ``sequence_selection.wgdi_ready`` value, or ``""``."""

    selection = qc.get("sequence_selection")
    if not isinstance(selection, dict):
        return ""
    ready = selection.get("wgdi_ready")
    if isinstance(ready, str) and ready in VALID_READINESS:
        return ready
    return ""


def _qc_candidates(project_root: Path, species_id: str) -> list[Path]:
    return [
        output_root(project_root) / species_id / f"{species_id}.qc.json",
        project_root / "analysis_v2" / species_id / f"{species_id}.qc.json",
        project_root
        / "plant-genome-resources"
        / "analysis_v2"
        / species_id
        / f"{species_id}.qc.json",
        project_root / species_id / f"{species_id}.qc.json",
    ]


def _index_candidates(project_root: Path) -> list[Path]:
    return [
        project_root / "resources" / "wgdi_v2_index.tsv",
        project_root / "plant-genome-resources" / "resources" / "wgdi_v2_index.tsv",
    ]


def read_species_readiness(project_root: Path, species_id: str) -> str:
    """Resolve readiness from QC JSON first, then the WGDI index TSV.

    Returns ``yes``/``marginal``/``no`` or ``""`` when no readable record exists.
    QC JSON is parsed structurally (never text-matched), so a stray or nested
    key cannot produce a false reading.
    """

    root = Path(project_root)
    for qc_path in _qc_candidates(root, species_id):
        if not qc_path.is_file():
            continue
        try:
            qc = json.loads(qc_path.read_text(encoding="utf-8"))
        except (OSError, ValueError, json.JSONDecodeError):
            continue
        if isinstance(qc, dict):
            ready = extract_qc_readiness(qc)
            if ready:
                return ready
    for index_path in _index_candidates(root):
        if not index_path.is_file():
            continue
        try:
            with index_path.open(encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if row.get("code") == species_id:
                        ready = (row.get("wgdi_ready") or "").strip()
                        if ready in VALID_READINESS:
                            return ready
        except OSError:
            continue
    return ""


def read_readiness_map(
    project_root: Path, species_ids: list[str]
) -> dict[str, str]:
    return {
        species_id: read_species_readiness(project_root, species_id)
        for species_id in species_ids
    }