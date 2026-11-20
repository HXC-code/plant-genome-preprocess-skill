"""Read-only dataset audit.

The audit never creates, moves or rewrites species data.  It summarizes the
inventory, registry, index, QC records and rejection reasons so operators can
see what is ready, what is missing and why species failed.
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from pgci_preprocess.batch_prepare import compound_suffix
from pgci_preprocess.reports import read_tsv


@dataclass(slots=True)
class AuditReport:
    project_root: Path
    inventory_species: int = 0
    registry_species: int = 0
    index_species: int = 0
    summary_species: int = 0
    qc_ok: int = 0
    qc_bad: int = 0
    status_distribution: dict[str, int] = field(default_factory=dict)
    ready_distribution: dict[str, int] = field(default_factory=dict)
    annotation_levels: dict[str, int] = field(default_factory=dict)
    format_exceptions: dict[str, int] = field(default_factory=dict)
    rejection_reasons: dict[str, int] = field(default_factory=dict)
    failure_reasons: dict[str, int] = field(default_factory=dict)
    suffix_distribution: list[dict[str, object]] = field(default_factory=list)
    species_details: list[dict[str, object]] = field(default_factory=list)
    problems: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return {
            "project_root": str(self.project_root),
            "inventory_species": self.inventory_species,
            "registry_species": self.registry_species,
            "index_species": self.index_species,
            "summary_species": self.summary_species,
            "qc_ok": self.qc_ok,
            "qc_bad": self.qc_bad,
            "status_distribution": self.status_distribution,
            "ready_distribution": self.ready_distribution,
            "annotation_levels": self.annotation_levels,
            "format_exceptions": self.format_exceptions,
            "rejection_reasons": self.rejection_reasons,
            "failure_reasons": self.failure_reasons,
            "suffix_distribution": self.suffix_distribution,
            "species_details": self.species_details,
            "problems": self.problems,
        }


def _bool_value(value: object) -> bool:
    return bool(value) if not isinstance(value, str) else value.strip().casefold() in {
        "1",
        "true",
        "yes",
        "y",
    }


def audit_project(
    project_root: Path | str,
    *,
    species: set[str] | None = None,
) -> AuditReport:
    project = Path(project_root).resolve()
    report = AuditReport(project_root=project)
    if not project.is_dir():
        report.problems.append(f"project root not found: {project}")
        return report

    inventory = read_tsv(project / "resources" / "download_inventory.tsv")
    registry = read_tsv(project / "resources" / "species_registry_v2.tsv")
    index = read_tsv(project / "resources" / "wgdi_v2_index.tsv")
    summary = read_tsv(project / "resources" / "stage_preprocess_v2_summary.tsv")
    report.inventory_species = len(inventory)
    report.registry_species = len(registry)
    report.index_species = len(index)
    report.summary_species = len(summary)

    if not index:
        report.problems.append("wgdi_v2_index.tsv is missing or empty")

    analysis_root = project / "analysis_v2"
    clean_root = project / "clean_v2"
    if analysis_root.is_dir():
        report.status_distribution["analysis_dirs"] = sum(
            1
            for path in analysis_root.iterdir()
            if path.is_dir() and not path.name.startswith(".")
        )
    if clean_root.is_dir():
        report.status_distribution["clean_dirs"] = sum(
            1
            for path in clean_root.iterdir()
            if path.is_dir() and not path.name.startswith(".")
        )

    statuses: Counter[str] = Counter()
    ready: Counter[str] = Counter()
    levels: Counter[str] = Counter()
    exceptions: Counter[str] = Counter()
    rejections: Counter[str] = Counter()
    failures: Counter[str] = Counter()

    for row in summary:
        status = str(row.get("status", "")).strip() or "UNKNOWN"
        statuses[status] += 1
        ready[str(row.get("wgdi_ready", "")).strip() or "UNKNOWN"] += 1
        if status not in {"OK", "SKIPPED"}:
            failures[str(row.get("reason", "")).strip() or "unknown"] += 1

    suffix_counts: Counter[tuple[str, str]] = Counter()
    for row in inventory:
        for role in ("genome", "gff", "cds", "pep"):
            value = str(row.get(role, "")).strip()
            if value:
                suffix_counts[(role, compound_suffix(value))] += 1

    for row in index:
        code = str(row.get("code", "")).strip()
        show_id = str(row.get("show_id", "")).strip()
        species_tokens = {str(item).strip().casefold() for item in species or []}
        if (
            species is not None
            and show_id.casefold() not in species_tokens
            and code.casefold() not in species_tokens
            and f"s{show_id}".casefold() not in species_tokens
        ):
            continue
        qc_relative = str(row.get("qc", "")).strip()
        if not qc_relative:
            report.qc_bad += 1
            continue
        qc_path = (project / qc_relative).resolve()
        try:
            qc_path.relative_to(project)
        except ValueError:
            report.qc_bad += 1
            report.problems.append(f"{code}: QC path escapes project root")
            continue
        if not qc_path.is_file():
            report.qc_bad += 1
            continue
        try:
            qc = json.loads(qc_path.read_text(encoding="utf-8"))
        except (OSError, ValueError, json.JSONDecodeError):
            report.qc_bad += 1
            continue
        if qc.get("status") == "OK" and qc.get("code") == code:
            report.qc_ok += 1
        else:
            report.qc_bad += 1
        levels[str(qc.get("annotation_level") or "legacy_unrecorded")] += 1
        for key, value in (qc.get("format_exceptions") or {}).items():
            if _bool_value(value):
                exceptions[key] += 1
        for reason, count in (qc.get("rejection_counts") or {}).items():
            try:
                rejections[str(reason)] += int(count)
            except (TypeError, ValueError):
                rejections[str(reason)] += 0
        qc_show_id = str(qc.get("show_id", "")).strip()
        if species is not None and (
            qc_show_id.casefold() in species_tokens
            or code.casefold() in species_tokens
            or f"s{qc_show_id}".casefold() in species_tokens
        ):
            report.species_details.append(
                {
                    "code": code,
                    "show_id": qc.get("show_id", ""),
                    "species": qc.get("species", ""),
                    "status": qc.get("status", ""),
                    "annotation_level": qc.get("annotation_level", ""),
                    "format_detection": qc.get("format_detection", {}),
                    "format_exceptions": qc.get("format_exceptions", {}),
                    "rejection_counts": qc.get("rejection_counts", {}),
                    "sequence_selection": qc.get("sequence_selection", {}),
                    "wgdi_bundle": qc.get("wgdi_bundle", {}),
                }
            )

    report.status_distribution.update(dict(sorted(statuses.items())))
    report.ready_distribution = dict(sorted(ready.items()))
    report.annotation_levels = dict(sorted(levels.items()))
    report.format_exceptions = dict(sorted(exceptions.items()))
    report.rejection_reasons = dict(sorted(rejections.items(), key=lambda item: (-item[1], item[0])))
    report.failure_reasons = dict(sorted(failures.items()))
    report.suffix_distribution = [
        {"role": role, "compound_suffix": suffix, "files": count}
        for (role, suffix), count in sorted(suffix_counts.items())
    ]
    return report


def render_text(report: AuditReport) -> str:
    lines = [
        f"project_root\t{report.project_root}",
        f"inventory_species\t{report.inventory_species}",
        f"registry_species\t{report.registry_species}",
        f"index_species\t{report.index_species}",
        f"summary_species\t{report.summary_species}",
        f"qc_ok\t{report.qc_ok}",
        f"qc_bad\t{report.qc_bad}",
    ]
    for label, values in (
        ("status_distribution", report.status_distribution),
        ("ready_distribution", report.ready_distribution),
        ("annotation_levels", report.annotation_levels),
        ("format_exceptions", report.format_exceptions),
        ("rejection_reasons", report.rejection_reasons),
        ("failure_reasons", report.failure_reasons),
    ):
        lines.append(f"{label}\t{json.dumps(values, ensure_ascii=False, sort_keys=True)}")
    lines.append("suffix_distribution")
    for row in report.suffix_distribution:
        lines.append(
            f"  {row['role']}\t{row['compound_suffix']}\t{row['files']}"
        )
    if report.species_details:
        lines.append("species_details")
        for row in report.species_details:
            lines.append(json.dumps(row, ensure_ascii=False, sort_keys=True))
    for problem in report.problems:
        lines.append(f"PROBLEM\t{problem}")
    return "\n".join(lines)
