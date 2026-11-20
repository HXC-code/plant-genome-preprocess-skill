"""Safely migrate legacy ``s<show_id>`` metadata codes to abbreviations."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from pgci_preprocess.reports import read_tsv_with_fields
from pgci_preprocess.species_names import (
    SpeciesNamingError,
    load_species_names,
)
from pgci_preprocess.validate import (
    INDEX_RELATIVE,
    resolve_layout,
    sync_compat_links,
)


@dataclass(slots=True)
class SpeciesNameMigrationResult:
    project_root: Path
    backup_dir: Path | None
    report_path: Path
    changed_files: dict[str, int]
    unresolved: list[str]
    links: list[str]
    dry_run: bool

    def as_dict(self) -> dict[str, object]:
        return {
            "project_root": str(self.project_root),
            "backup_dir": str(self.backup_dir) if self.backup_dir else "",
            "report_path": str(self.report_path),
            "changed_files": self.changed_files,
            "unresolved": self.unresolved,
            "links": self.links,
            "dry_run": self.dry_run,
        }


def _atomic_write_rows(
    path: Path,
    fields: list[str],
    rows: list[dict[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        dir=path.parent,
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fields,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in fields})
        os.replace(temporary_name, path)
    except Exception:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _migrate_qc_metadata(
    canonical_root: Path,
    *,
    names,
    dry_run: bool,
) -> tuple[int, list[str]]:
    index_path = canonical_root / INDEX_RELATIVE
    if not index_path.is_file():
        return 0, []
    _fields, rows = read_tsv_with_fields(index_path)
    changed = 0
    unresolved: list[str] = []
    manifest_updates: dict[str, tuple[int, str]] = {}
    for row in rows:
        code = str(row.get("code", "")).strip()
        show_id = str(row.get("show_id", "")).strip()
        artifact_code = str(row.get("artifact_code", "")).strip()
        qc_relative = str(row.get("qc", "")).strip()
        if not code or not qc_relative:
            continue
        qc_path = (canonical_root / qc_relative).resolve()
        try:
            qc_path.relative_to(canonical_root.resolve())
        except ValueError:
            unresolved.append(f"{code}: QC path escapes project root: {qc_relative}")
            continue
        if not qc_path.is_file():
            unresolved.append(f"{code}: QC file missing: {qc_relative}")
            continue
        try:
            qc = json.loads(qc_path.read_text(encoding="utf-8"))
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            unresolved.append(f"{code}: unreadable QC: {exc}")
            continue
        if not isinstance(qc, dict):
            unresolved.append(f"{code}: QC root is not an object")
            continue
        before = dict(qc)
        old_code = str(
            artifact_code
            or qc.get("legacy_code")
            or qc.get("code", "")
        ).strip()
        qc["code"] = code
        qc["show_id"] = show_id or qc.get("show_id", "")
        qc["artifact_code"] = artifact_code or old_code
        if old_code:
            qc["legacy_code"] = old_code
        if qc != before and not dry_run:
            descriptor, temporary_name = tempfile.mkstemp(
                prefix=f".{qc_path.name}.",
                dir=qc_path.parent,
            )
            try:
                with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
                    json.dump(qc, handle, ensure_ascii=False, indent=2)
                    handle.write("\n")
                os.replace(temporary_name, qc_path)
            except Exception:
                Path(temporary_name).unlink(missing_ok=True)
                raise
            manifest_updates[qc_relative] = (qc_path.stat().st_size, _sha256(qc_path))
        if qc != before:
            changed += 1

    manifest_path = canonical_root / "resources" / "wgdi_v2_manifest_sha256.tsv"
    if manifest_updates and manifest_path.is_file():
        fields, manifest_rows = read_tsv_with_fields(manifest_path)
        for manifest_row in manifest_rows:
            relative = str(manifest_row.get("relative_path", "")).strip()
            if relative in manifest_updates:
                size, digest = manifest_updates[relative]
                manifest_row["bytes"] = str(size)
                manifest_row["sha256"] = digest
        _atomic_write_rows(manifest_path, fields, manifest_rows)
    return changed, unresolved


def _resolve_row_show_id(
    row: dict[str, str],
    *,
    id_field: str,
    names,
) -> str:
    show_id = str(row.get(id_field, "")).strip()
    if show_id:
        return show_id
    code = str(row.get("code", "")).strip()
    return names.resolve_show_id(code)


def _update_table(
    path: Path,
    *,
    names,
    id_field: str,
    insert_columns: list[tuple[str, str]],
    update_show_id: bool,
) -> tuple[int, list[str]]:
    if not path.is_file():
        return 0, []
    fields, rows = read_tsv_with_fields(path)
    for column, _default in insert_columns:
        if column not in fields:
            if column == "species_code" and "showid" in fields:
                fields.insert(fields.index("showid") + 1, column)
            elif column == "show_id" and "code" in fields:
                fields.insert(fields.index("code"), column)
            else:
                fields.append(column)
    changed = 0
    unresolved: list[str] = []
    for row in rows:
        try:
            show_id = _resolve_row_show_id(row, id_field=id_field, names=names)
            code = names.code_for(show_id)
        except SpeciesNamingError as exc:
            unresolved.append(str(exc))
            continue
        before = dict(row)
        if update_show_id:
            row["show_id"] = show_id
        if "species_code" in fields:
            row["species_code"] = code
        if "code" in fields:
            raw_code = str(row.get("code", "")).strip()
            if "artifact_code" in fields and not str(row.get("artifact_code", "")).strip():
                row["artifact_code"] = raw_code or code
            row["code"] = code
        if row != before:
            changed += 1
    if changed:
        _atomic_write_rows(path, fields, rows)
    return changed, unresolved


def migrate_species_names(
    project_root: Path | str,
    *,
    dry_run: bool = False,
    refresh_links: bool = True,
) -> SpeciesNameMigrationResult:
    workspace, canonical_root = resolve_layout(project_root)
    names = load_species_names(canonical_root)
    if names is None:
        raise SpeciesNamingError(
            "resources/species_abbreviations.tsv is required before migration"
        )
    resources = canonical_root / "resources"
    report_path = resources / "species_names_migration.json"
    tables: list[tuple[str, str, list[tuple[str, str]], bool]] = [
        (
            "download_inventory.tsv",
            "showid",
            [("species_code", "")],
            False,
        ),
        (
            "species_registry_v2.tsv",
            "show_id",
            [],
            False,
        ),
        (
            "stage_preprocess_v2_summary.tsv",
            "show_id",
            [],
            False,
        ),
        (
            INDEX_RELATIVE.name,
            "show_id",
            [
                ("show_id", ""),
                ("artifact_code", ""),
            ],
            True,
        ),
        (
            "wgdi_v2_manifest_sha256.tsv",
            "",
            [],
            False,
        ),
    ]
    changed_files: dict[str, int] = {}
    unresolved: list[str] = []
    targets = [resources / table[0] for table in tables]
    backup_dir: Path | None = None
    if not dry_run:
        stamp = datetime.now().strftime("%Y%m%dT%H%M%S%f")
        backup_dir = resources / "backups" / f"species_names_{stamp}"
        backup_dir.mkdir(parents=True, exist_ok=False)
        for path in targets:
            if path.is_file():
                shutil.copy2(path, backup_dir / path.name)

    for table_name, id_field, insert_columns, update_show_id in tables:
        path = resources / table_name
        if not path.is_file():
            continue
        if path.name == "wgdi_v2_manifest_sha256.tsv":
            fields, rows = read_tsv_with_fields(path)
            changed = 0
            for row in rows:
                raw_code = str(row.get("code", "")).strip()
                try:
                    code = names.resolve_show_id(raw_code)
                    new_code = names.code_for(code)
                    if row["code"] != new_code:
                        row["code"] = new_code
                        changed += 1
                except SpeciesNamingError as exc:
                    unresolved.append(str(exc))
            if changed and not dry_run:
                _atomic_write_rows(path, fields, rows)
            changed_files[path.name] = changed
            continue
        if dry_run:
            # Compute the same changes in memory without replacing files.
            fields, rows = read_tsv_with_fields(path)
            changed = 0
            for row in rows:
                try:
                    show_id = _resolve_row_show_id(
                        row,
                        id_field=id_field,
                        names=names,
                    )
                    code = names.code_for(show_id)
                except SpeciesNamingError as exc:
                    unresolved.append(str(exc))
                    continue
                if update_show_id:
                    row["show_id"] = show_id
                if "code" in row:
                    row["code"] = code
                changed += 1
            changed_files[path.name] = changed
        else:
            changed, table_unresolved = _update_table(
                path,
                names=names,
                id_field=id_field,
                insert_columns=insert_columns,
                update_show_id=update_show_id,
            )
            changed_files[path.name] = changed
            unresolved.extend(table_unresolved)

    qc_changed, qc_unresolved = _migrate_qc_metadata(
        canonical_root,
        names=names,
        dry_run=dry_run,
    )
    changed_files["qc_metadata"] = qc_changed
    unresolved.extend(qc_unresolved)

    links: list[str] = []
    if refresh_links:
        if dry_run:
            links.append("dry-run: would refresh compatibility links")
        else:
            links = sync_compat_links(canonical_root, workspace)
    result = SpeciesNameMigrationResult(
        project_root=canonical_root,
        backup_dir=backup_dir,
        report_path=report_path,
        changed_files=changed_files,
        unresolved=sorted(set(unresolved)),
        links=links,
        dry_run=dry_run,
    )
    if not dry_run:
        report_path.write_text(
            json.dumps(result.as_dict(), ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
    return result
