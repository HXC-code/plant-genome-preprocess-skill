"""Human adjudication of the final WGDI chromosome count.

The program proposes a chromosome set; a human reviewer confirms the count and
signs off. Each record is bound to the exact QC it judged (pipeline version +
SHA-256), so a re-run that changes the result invalidates stale sign-offs.
``validate`` and the downstream runners gate on a current, accepting
adjudication for every species (the lab SOP: the final chromosome number is a
human judgement, not an automatic pass).
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from pgci_preprocess.readiness import _qc_candidates

DECISIONS = frozenset({"accept", "accept_with_note", "adjust", "reject"})
ACCEPTING = frozenset({"accept", "accept_with_note"})


def adjudication_path(project_root: Path, code: str) -> Path:
    return Path(project_root) / "resources" / "adjudications" / f"{code}.json"


def read_adjudication(project_root: Path, code: str) -> dict[str, object] | None:
    path = adjudication_path(project_root, code)
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError):
        return None
    return data if isinstance(data, dict) else None


def find_qc(project_root: Path, species_id: str) -> tuple[Path, dict] | None:
    for qc_path in _qc_candidates(Path(project_root), species_id):
        if not qc_path.is_file():
            continue
        try:
            qc = json.loads(qc_path.read_text(encoding="utf-8"))
        except (OSError, ValueError, json.JSONDecodeError):
            continue
        if isinstance(qc, dict):
            return qc_path, qc
    return None


def _qc_fingerprint(qc_path: Path, qc: dict) -> tuple[str, str]:
    version = str(qc.get("pipeline_version", ""))
    sha = hashlib.sha256(qc_path.read_bytes()).hexdigest()
    return version, sha


def gate_status_from_qc(
    project_root: Path, code: str, qc_path: Path, qc: dict
) -> str:
    """Return ``ok`` or a machine-readable reason the gate must block."""

    record = read_adjudication(project_root, code)
    if record is None:
        return "missing"
    version, sha = _qc_fingerprint(qc_path, qc)
    if record.get("qc_pipeline_version") != version or record.get("qc_sha256") != sha:
        return "stale"
    decision = record.get("decision")
    if decision == "reject":
        return "rejected"
    if decision == "adjust":
        return "pending_adjustment"
    if decision not in ACCEPTING:
        return "invalid_decision"
    if not str(record.get("reviewer", "")).strip():
        return "missing_reviewer"
    if record.get("confirmed_n") != qc.get("n_wgdi_sequences"):
        return "count_mismatch"
    return "ok"


def gate_status(project_root: Path, species_id: str) -> str:
    found = find_qc(project_root, species_id)
    if found is None:
        return "no_qc"
    qc_path, qc = found
    return gate_status_from_qc(project_root, species_id, qc_path, qc)


def write_adjudication(
    project_root: Path,
    code: str,
    *,
    confirmed_n: int,
    decision: str,
    reviewer: str,
    rationale: str = "",
) -> Path:
    if decision not in DECISIONS:
        raise ValueError(f"decision must be one of {sorted(DECISIONS)}")
    if not reviewer.strip():
        raise ValueError("reviewer is required for accountability")
    found = find_qc(project_root, code)
    if found is None:
        raise ValueError(f"no QC found for {code}; run preprocessing first")
    qc_path, qc = found
    version, sha = _qc_fingerprint(qc_path, qc)
    record = {
        "code": code,
        "program_n_wgdi_sequences": qc.get("n_wgdi_sequences"),
        "confirmed_n": int(confirmed_n),
        "decision": decision,
        "reviewer": reviewer.strip(),
        "rationale": rationale,
        "reviewed_at": datetime.now(timezone.utc).isoformat(),
        "qc_pipeline_version": version,
        "qc_sha256": sha,
    }
    path = adjudication_path(project_root, code)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return path