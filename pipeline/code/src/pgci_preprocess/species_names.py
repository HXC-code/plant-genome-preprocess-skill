"""Canonical species abbreviations and legacy ``s<show_id>`` aliases."""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

SPECIES_NAME_PATH = Path("resources") / "species_abbreviations.tsv"
SPECIES_NAME_ENV = "PGCI_SPECIES_ABBREVIATIONS"
REQUIRED_COLUMNS = {
    "show_id",
    "latin_name",
    "order",
    "family",
    "abbreviation",
}
SAFE_CODE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,63}$")


class SpeciesNamingError(RuntimeError):
    pass


@dataclass(frozen=True, slots=True)
class SpeciesName:
    show_id: str
    latin_name: str
    order: str
    family: str
    abbreviation: str


@dataclass(frozen=True, slots=True)
class SpeciesNameMap:
    path: Path
    by_show_id: dict[str, SpeciesName]
    by_code: dict[str, SpeciesName]

    def code_for(
        self,
        show_id: str,
        *,
        latin_name: str = "",
    ) -> str:
        show_id = str(show_id).strip()
        entry = self.by_show_id.get(show_id)
        if entry is None:
            raise SpeciesNamingError(
                f"{self.path}: no abbreviation row for show_id={show_id}"
                + (f" ({latin_name})" if latin_name else "")
            )
        return entry.abbreviation

    def resolve_show_id(self, token: str) -> str:
        value = str(token).strip()
        if not value:
            raise SpeciesNamingError("empty species selector")
        if value in self.by_show_id:
            return value
        if value.casefold().startswith("s") and value[1:].isdigit():
            legacy_id = value[1:]
            if legacy_id in self.by_show_id:
                return legacy_id
        entry = self.by_code.get(value.casefold())
        if entry is not None:
            return entry.show_id
        raise SpeciesNamingError(
            f"unknown species selector {value!r}; use show_id, s<show_id>, "
            "or an abbreviation from the species table"
        )

    def entry_for_show_id(self, show_id: str) -> SpeciesName:
        entry = self.by_show_id.get(str(show_id).strip())
        if entry is None:
            raise SpeciesNamingError(
                f"{self.path}: no abbreviation row for show_id={show_id}"
            )
        return entry


def species_name_path(
    project_root: Path | str,
    *,
    explicit_path: Path | str | None = None,
) -> Path | None:
    if explicit_path:
        return Path(explicit_path).expanduser().resolve()
    environment = os.environ.get(SPECIES_NAME_ENV, "").strip()
    if environment:
        return Path(environment).expanduser().resolve()
    candidate = Path(project_root).resolve() / SPECIES_NAME_PATH
    return candidate if candidate.is_file() else None


@lru_cache(maxsize=32)
def _load_cached(path_text: str, mtime_ns: int, size: int) -> SpeciesNameMap:
    del mtime_ns, size
    path = Path(path_text)
    by_id: dict[str, SpeciesName] = {}
    by_code: dict[str, SpeciesName] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(REQUIRED_COLUMNS - set(reader.fieldnames or []))
        if missing:
            raise SpeciesNamingError(
                f"{path}: missing columns: {', '.join(missing)}"
            )
        for line_number, row in enumerate(reader, 2):
            show_id = (row.get("show_id") or "").strip()
            code = (row.get("abbreviation") or "").strip()
            if not show_id or not code:
                continue
            if not show_id.isdigit():
                raise SpeciesNamingError(
                    f"{path}:{line_number}: show_id must contain digits only"
                )
            if not SAFE_CODE_RE.fullmatch(code):
                raise SpeciesNamingError(
                    f"{path}:{line_number}: unsafe abbreviation {code!r}"
                )
            if show_id in by_id:
                raise SpeciesNamingError(
                    f"{path}:{line_number}: duplicate show_id {show_id}"
                )
            folded = code.casefold()
            if folded in by_code:
                raise SpeciesNamingError(
                    f"{path}:{line_number}: duplicate abbreviation {code!r}"
                )
            entry = SpeciesName(
                show_id=show_id,
                latin_name=(row.get("latin_name") or "").strip(),
                order=(row.get("order") or "").strip(),
                family=(row.get("family") or "").strip(),
                abbreviation=code,
            )
            by_id[show_id] = entry
            by_code[folded] = entry
    if not by_id:
        raise SpeciesNamingError(f"{path}: no species rows")
    return SpeciesNameMap(path=path, by_show_id=by_id, by_code=by_code)


def load_species_names(
    project_root: Path | str,
    *,
    explicit_path: Path | str | None = None,
) -> SpeciesNameMap | None:
    path = species_name_path(project_root, explicit_path=explicit_path)
    if path is None:
        return None
    if not path.is_file():
        raise SpeciesNamingError(f"species abbreviation table not found: {path}")
    stat = path.stat()
    return _load_cached(str(path), stat.st_mtime_ns, stat.st_size)


def resolve_species_code(
    project_root: Path | str,
    show_id: str,
    *,
    latin_name: str = "",
    explicit_path: Path | str | None = None,
) -> str:
    names = load_species_names(project_root, explicit_path=explicit_path)
    if names is None:
        return f"s{show_id}"
    return names.code_for(show_id, latin_name=latin_name)


def resolve_show_ids(
    project_root: Path | str,
    tokens: set[str] | list[str] | tuple[str, ...],
    *,
    explicit_path: Path | str | None = None,
) -> set[str]:
    names = load_species_names(project_root, explicit_path=explicit_path)
    result: set[str] = set()
    for raw in tokens:
        token = str(raw).strip()
        if not token:
            continue
        if names is not None:
            result.add(names.resolve_show_id(token))
            continue
        if token.isdigit():
            result.add(token)
        elif token.casefold().startswith("s") and token[1:].isdigit():
            result.add(token[1:])
        else:
            raise SpeciesNamingError(
                f"cannot resolve species selector {token!r} without "
                f"{SPECIES_NAME_PATH.as_posix()}"
            )
    return result


def legacy_species_code(show_id: str) -> str:
    return f"s{show_id}"
