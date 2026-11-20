"""Audit species naming across table, downloads and WGDI artifacts.

Produces ``resources/naming_audit.tsv`` with one row per species and exits
non-zero when any rule is violated.  Rules:

* ``resources/species_abbreviations.tsv`` exists, IDs are numeric and
  abbreviations are unique;
* every download directory ``<show_id>_<Order>_<Family>(<Species>)`` has a
  table row (and its Latin name agrees with the table);
* every index row uses the table abbreviation as ``code`` and its artifact
  paths point at ``analysis_v2/<code>/<code>.*`` / ``clean_v2/<code>/<code>.*``;
* the first gene ID of ``<code>.new.gff/.cds/.pep`` is ``<code><seq>g<order>``;
* no legacy ``s<show_id>`` artifact directory remains under either layer.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

from pgci_preprocess.artifact_naming import (
    ANALYSIS_REQUIRED_SUFFIXES,
    bundle_name_problems,
    first_gene_id,
    is_canonical_gene_id,
)
from pgci_preprocess.reports import read_tsv, write_tsv
from pgci_preprocess.species_names import (
    SpeciesNamingError,
    load_species_names,
)
from pgci_preprocess.validate import resolve_layout

FIELDS = [
    "show_id",
    "code",
    "latin_name",
    "download_directory",
    "index_code",
    "artifact_code",
    "status",
    "problems",
]


@dataclass
class NamingAuditReport:
    rows: list[dict[str, object]] = field(default_factory=list)
    report_path: Path | None = None
    summary: dict[str, int] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return all(row.get("status") != "FAIL" for row in self.rows)

    def as_dict(self) -> dict[str, object]:
        return {
            "report_path": str(self.report_path) if self.report_path else "",
            "summary": self.summary,
            "rows": self.rows,
        }


def _download_index(project_root: Path) -> dict[str, str]:
    downloads = project_root / "data" / "downloads"
    mapping: dict[str, str] = {}
    if not downloads.is_dir():
        return mapping
    for path in downloads.iterdir():
        if path.is_dir() and "_" in path.name:
            mapping.setdefault(path.name.split("_", 1)[0], path.name)
    return mapping


def _normalize(name: str) -> str:
    return " ".join(str(name).replace("×", "x").split()).casefold()


def audit_naming(project_root: Path | str, *, table_path: Path | str | None = None) -> NamingAuditReport:
    workspace, canonical_root = resolve_layout(project_root)
    names = load_species_names(canonical_root, explicit_path=table_path)
    if names is None:
        raise SpeciesNamingError(
            "resources/species_abbreviations.tsv is required to audit naming"
        )
    report = NamingAuditReport()
    inventory_rows = read_tsv(canonical_root / "resources" / "download_inventory.tsv")
    index_rows = read_tsv(canonical_root / "resources" / "wgdi_v2_index.tsv")
    downloads = _download_index(canonical_root)
    index_by_show_id = {
        str(row.get("show_id", "")).strip(): row
        for row in index_rows
        if str(row.get("show_id", "")).strip()
    }
    legacy_dirs = sorted(
        f"{layer}/{path.name}"
        for layer in ("analysis_v2", "clean_v2")
        for path in (canonical_root / layer).glob("s*")
        if path.is_dir() and path.name[1:].isdigit()
    )
    if legacy_dirs:
        report.rows.append(
            {
                "show_id": "",
                "code": "",
                "latin_name": "",
                "download_directory": "",
                "index_code": "",
                "artifact_code": "",
                "status": "FAIL",
                "problems": "legacy_artifact_directories:" + ",".join(legacy_dirs),
            }
        )

    seen: set[str] = set()
    taken_codes = set(names.by_code)
    for row in inventory_rows:
        show_id = str(row.get("showid", "")).strip()
        if not show_id:
            continue
        seen.add(show_id)
        problems: list[str] = []
        entry = names.by_show_id.get(show_id)
        latin_name = str(row.get("Species", ""))
        if entry is None:
            # 表里没有的新物种按命名规则现算：属首字母 + 种加词两字母，
            # 占用则换组合，全部占满再扩到 4 位及更长。
            from pgci_preprocess.species_codegen import derive_code

            try:
                code = derive_code(latin_name, taken_codes)
            except ValueError as exc:
                code = ""
                problems.append(f"cannot_derive_code:{exc}")
            if code:
                taken_codes.add(code.casefold())
                inventory_code = str(row.get("species_code", "")).strip()
                if inventory_code not in {"", code}:
                    problems.append(
                        f"inventory_species_code_mismatch:{inventory_code}!={code}"
                    )
        else:
            code = entry.abbreviation
            latin_name = entry.latin_name
            if str(row.get("species_code", "")).strip() not in {"", code}:
                problems.append(
                    f"inventory_species_code_mismatch:{row.get('species_code')}!={code}"
                )
            if _normalize(latin_name) != _normalize(row.get("Species", "")):
                problems.append(f"latin_name_mismatch:{row.get('Species')}!={latin_name}")
        download = downloads.get(show_id, "")
        if not download:
            problems.append("download_directory_missing")
        else:
            dir_species = download.split("(", 1)[-1].rstrip(")")
            if _normalize(dir_species) != _normalize(latin_name):
                problems.append(f"download_dir_species_mismatch:{dir_species}!={latin_name}")
        index_row = index_by_show_id.get(show_id, {})
        index_code = str(index_row.get("code", "")).strip()
        artifact_code = str(index_row.get("artifact_code", "")).strip()
        if index_row:
            if code and index_code != code:
                problems.append(f"index_code_mismatch:{index_code}!={code}")
            if code and artifact_code != code:
                problems.append(f"artifact_code_mismatch:{artifact_code}!={code}")
            for key in ("pep", "cds", "gff", "lens"):
                relative = str(index_row.get(key, ""))
                if not relative:
                    continue
                if code and f"/{code}/" not in relative:
                    problems.append(f"index_path_not_code:{key}:{relative}")
                if code and not Path(relative).name.startswith(f"{code}."):
                    problems.append(f"index_file_not_code:{key}:{relative}")
            if code:
                analysis_dir = canonical_root / "analysis_v2" / code
                clean_dir = canonical_root / "clean_v2" / code
                if analysis_dir.is_dir():
                    problems.extend(
                        bundle_name_problems(
                            analysis_dir,
                            code,
                            required_suffixes=ANALYSIS_REQUIRED_SUFFIXES,
                            kind_map={
                                ".new.gff": "gff",
                                ".cds": "fasta",
                                ".pep": "fasta",
                            },
                        )
                    )
                    for suffix, kind in (
                        (".new.gff", "gff"),
                        (".cds", "fasta"),
                        (".pep", "fasta"),
                    ):
                        path = analysis_dir / f"{code}{suffix}"
                        if path.is_file():
                            gene_id = first_gene_id(path, kind)
                            if gene_id and not is_canonical_gene_id(gene_id, code):
                                problems.append(
                                    f"gene_id_not_code:{suffix}:{gene_id}"
                                )
                else:
                    problems.append("analysis_bundle_missing")
                if not clean_dir.is_dir():
                    problems.append("clean_bundle_missing")
        report.rows.append(
            {
                "show_id": show_id,
                "code": code,
                "latin_name": latin_name,
                "download_directory": download,
                "index_code": index_code,
                "artifact_code": artifact_code,
                "status": "PASS" if not problems else "FAIL",
                "problems": ";".join(problems),
            }
        )

    for show_id in sorted(set(names.by_show_id) - seen):
        report.summary["table_only"] = report.summary.get("table_only", 0) + 1

    report.report_path = canonical_root / "resources" / "naming_audit.tsv"
    write_tsv(report.report_path, FIELDS, report.rows)
    report.summary.update(
        {
            "rows": len(report.rows),
            "pass": sum(1 for row in report.rows if row["status"] == "PASS"),
            "fail": sum(1 for row in report.rows if row["status"] == "FAIL"),
        }
    )
    return report


def render_text(report: NamingAuditReport) -> str:
    lines = [
        f"命名核查: {report.summary.get('rows', 0)} 行 / "
        f"PASS {report.summary.get('pass', 0)} / "
        f"FAIL {report.summary.get('fail', 0)}"
    ]
    if report.summary.get("table_only"):
        lines.append(f"表内未在当前清单: {report.summary['table_only']} 条")
    for row in report.rows:
        if row["status"] == "FAIL":
            lines.append(
                f"[FAIL] {row['show_id'] or '-'} {row['code'] or '-'} {row['problems']}"
            )
    if report.report_path:
        lines.append(f"报告: {report.report_path}")
    return "\n".join(lines)
