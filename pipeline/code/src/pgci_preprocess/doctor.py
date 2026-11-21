"""Read-only environment and project preflight checks."""

from __future__ import annotations

import csv
import importlib.util
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from pgci_preprocess.pipeline import discover_species
from pgci_preprocess.species_names import SpeciesNamingError, load_species_names


@dataclass(frozen=True)
class DoctorCheck:
    name: str
    status: str
    detail: str

    def as_dict(self) -> dict[str, str]:
        return {"name": self.name, "status": self.status, "detail": self.detail}


def _resolve_executable(
    requested: str = "",
    *,
    env_name: str = "",
    candidates: tuple[Path, ...] = (),
) -> str:
    value = requested.strip() or (os.environ.get(env_name, "").strip() if env_name else "")
    if value:
        path = Path(value).expanduser()
        if path.is_file():
            return str(path.resolve())
        resolved = shutil.which(value)
        return resolved or ""
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate.resolve())
    return ""


def _executable_check(
    checks: list[DoctorCheck],
    name: str,
    requested: str,
    env_name: str,
    candidates: tuple[Path, ...],
) -> str:
    resolved = _resolve_executable(requested, env_name=env_name, candidates=candidates)
    if not resolved:
        resolved = shutil.which(name) or ""
    if resolved:
        checks.append(DoctorCheck(name, "PASS", resolved))
    else:
        hint = f"; set {env_name} or pass an explicit path" if env_name else ""
        checks.append(DoctorCheck(name, "FAIL", f"not found{hint}"))
    return resolved


def _wgdi_probe(check_name: str, executable: str) -> DoctorCheck:
    """Run the actual launcher so import and startup failures fail doctor."""

    path = Path(executable)
    command = [executable, "--version"]
    if path.is_file() and path.suffix.casefold() not in {".exe", ".cmd", ".bat"}:
        bash = shutil.which("bash")
        if not bash:
            return DoctorCheck(
                check_name,
                "FAIL",
                "bash is required to run the WGDI launcher",
            )
        command = [bash, executable, "--version"]
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return DoctorCheck(check_name, "FAIL", str(exc))
    output = (completed.stdout or completed.stderr).strip().splitlines()
    detail = output[-1] if output else "no version output"
    if completed.returncode == 0:
        return DoctorCheck(check_name, "PASS", detail)
    return DoctorCheck(
        check_name,
        "FAIL",
        f"exit={completed.returncode}: {detail}",
    )


def _biopython_interpreters() -> list[str]:
    """Interpreters on this machine that can import Biopython (best effort)."""

    candidates: list[Path] = [Path(sys.executable)]
    for name in ("python", "python3", "py"):
        resolved = shutil.which(name)
        if resolved:
            candidates.append(Path(resolved))
    extra = os.environ.get("PGCI_PYTHON_CANDIDATES", "")
    candidates.extend(Path(item) for item in extra.split(os.pathsep) if item.strip())
    found: list[str] = []
    seen: set[str] = set()
    for candidate in candidates:
        if not candidate.is_file():
            continue
        resolved = str(candidate.resolve())
        if resolved in seen:
            continue
        seen.add(resolved)
        try:
            completed = subprocess.run(
                [resolved, "-c", "import Bio, sys; print(Bio.__version__)"],
                capture_output=True,
                text=True,
                timeout=20,
                check=False,
            )
        except (OSError, subprocess.SubprocessError):
            continue
        version = completed.stdout.strip()
        if completed.returncode == 0 and version:
            found.append(f"{resolved} (Bio {version})")
    return found


def run_doctor(
    project_root: Path | str,
    *,
    workflow: str = "preprocess",
    require_gbff: bool = False,
    allow_empty: bool = False,
    executables: dict[str, str] | None = None,
) -> dict[str, object]:
    """Inspect a project without creating, installing or modifying anything."""

    if workflow not in {"preprocess", "diamond", "blast", "wgdi", "full"}:
        raise ValueError("workflow must be preprocess, diamond, blast, wgdi or full")
    root = Path(project_root).expanduser().resolve()
    overrides = executables or {}
    checks: list[DoctorCheck] = []
    python_ok = sys.version_info >= (3, 11)
    checks.append(
        DoctorCheck(
            "python>=3.11",
            "PASS" if python_ok else "FAIL",
            f"{sys.version.split()[0]} ({sys.executable})",
        )
    )
    biopython = importlib.util.find_spec("Bio")
    if biopython:
        checks.append(DoctorCheck("biopython", "PASS", "available"))
    else:
        interpreters = _biopython_interpreters()
        detail = "required only for GBFF conversion" if not require_gbff else "required for GBFF conversion"
        if interpreters:
            detail += "; interpreters with Biopython available: " + ", ".join(interpreters)
            detail += " (run GBFF conversion with one of these, or install biopython>=1.83 here)"
        else:
            detail += "; install biopython>=1.83 (e.g. pip install 'biopython>=1.83')"
        checks.append(
            DoctorCheck(
                "biopython",
                "FAIL" if require_gbff else "WARN",
                detail,
            )
        )

    species_count = 0
    inventory_path = root / "resources" / "download_inventory.tsv"
    downloads = root / "data" / "downloads"
    if not root.is_dir():
        checks.append(DoctorCheck("project_root", "FAIL", str(root)))
    else:
        checks.append(DoctorCheck("project_root", "PASS", str(root)))
        checks.append(
            DoctorCheck(
                "downloads_directory",
                "PASS" if downloads.is_dir() else "FAIL",
                str(downloads),
            )
        )
        if not inventory_path.is_file():
            checks.append(DoctorCheck("inventory", "FAIL", f"missing: {inventory_path}"))
        else:
            try:
                with inventory_path.open(encoding="utf-8-sig", newline="") as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    fields = reader.fieldnames or []
                    rows = list(reader)
                required = {"showid", "Order", "Family", "Species", "gff", "genome", "cds", "pep"}
                missing_columns = sorted(required - set(fields))
                if missing_columns:
                    checks.append(DoctorCheck("inventory_header", "FAIL", ", ".join(missing_columns)))
                else:
                    checks.append(DoctorCheck("inventory_header", "PASS", str(inventory_path)))
                ids = [str(row.get("showid", "")).strip() for row in rows if str(row.get("showid", "")).strip()]
                names = load_species_names(root)
                if names is None:
                    checks.append(
                        DoctorCheck(
                            "species_abbreviations",
                            "WARN",
                            "resources/species_abbreviations.tsv not found; "
                            "artifact names would fall back to legacy s<show_id> codes. "
                            "Copy the table shipped with the package or run init-project.",
                        )
                    )
                else:
                    missing_codes = [show_id for show_id in ids if show_id not in names.by_show_id]
                    checks.append(
                        DoctorCheck(
                            "species_abbreviations",
                            "FAIL" if missing_codes else "PASS",
                            (
                                "missing show_id: " + ",".join(missing_codes[:10])
                                if missing_codes
                                else f"{len(ids)} inventory IDs mapped from {names.path}"
                            ),
                        )
                    )
                if len(ids) != len(set(ids)):
                    checks.append(DoctorCheck("inventory_unique_ids", "FAIL", "duplicate showid values"))
                else:
                    checks.append(DoctorCheck("inventory_unique_ids", "PASS", str(len(ids))))
                species = discover_species(root)
                species_count = len(species)
                if len(species) != len(ids):
                    checks.append(
                        DoctorCheck(
                            "inventory_species_directories",
                            "FAIL",
                            f"inventory IDs={len(ids)}, discovered species={len(species)}",
                        )
                    )
                else:
                    checks.append(DoctorCheck("inventory_species_directories", "PASS", str(len(species))))
                if not species:
                    checks.append(
                        DoctorCheck(
                            "species_count",
                            "WARN" if allow_empty else "FAIL",
                            "0 registered species",
                        )
                    )
                else:
                    checks.append(DoctorCheck("species_count", "PASS", str(len(species))))
                missing_files: list[str] = []
                escaping_files: list[str] = []
                for item in species:
                    for role, name in (
                        ("gff", item.gff_name),
                        ("cds", item.cds_name),
                        ("pep", item.pep_name),
                        ("genome", item.genome_name),
                    ):
                        if role == "genome" and not name:
                            continue
                        if not name:
                            missing_files.append(f"{item.code}:{role}=<empty>")
                            continue
                        path = (item.directory / name).resolve()
                        try:
                            path.relative_to(item.directory.resolve())
                        except ValueError:
                            escaping_files.append(f"{item.code}:{role}={name}")
                            continue
                        if not path.is_file() or path.stat().st_size == 0:
                            missing_files.append(f"{item.code}:{role}={name}")
                checks.append(
                    DoctorCheck(
                        "inventory_paths_inside_species",
                        "FAIL" if escaping_files else "PASS",
                        ", ".join(escaping_files[:10]) if escaping_files else "all paths are local",
                    )
                )
                checks.append(
                    DoctorCheck(
                        "input_files",
                        "FAIL" if missing_files else "PASS",
                        ", ".join(missing_files[:10]) if missing_files else f"checked {species_count} species",
                    )
                )
            except (OSError, csv.Error, ValueError, SpeciesNamingError) as exc:
                checks.append(DoctorCheck("inventory_parse", "FAIL", str(exc)))

    if workflow != "preprocess":
        bash_candidates = (
            root.parent / "bin" / "bash",
            root / "bin" / "bash",
        )
        _executable_check(checks, "bash", overrides.get("bash", ""), "BASH", bash_candidates)
    if workflow in {"diamond", "full"}:
        _executable_check(
            checks,
            "diamond",
            overrides.get("diamond", ""),
            "DIAMOND",
            (),
        )
    if workflow == "blast":
        _executable_check(
            checks,
            "blastp",
            overrides.get("blastp", ""),
            "BLASTP",
            (),
        )
        _executable_check(
            checks,
            "makeblastdb",
            overrides.get("makeblastdb", ""),
            "MAKEBLASTDB",
            (),
        )
    if workflow in {"wgdi", "full"}:
        wgdi = _executable_check(
            checks,
            "wgdi",
            overrides.get("wgdi", ""),
            "WGDI",
            (root / "bin" / "wgdi", root.parent / "bin" / "wgdi"),
        )
        if wgdi:
            checks.append(_wgdi_probe("wgdi_probe", wgdi))
        for name in ("muscle", "yn00"):
            _executable_check(
                checks,
                name,
                overrides.get(name, ""),
                name.upper(),
                (),
            )
        _executable_check(
            checks,
            "perl",
            overrides.get("perl", ""),
            "PERL",
            (),
        )

    counts = {
        status: sum(1 for check in checks if check.status == status)
        for status in ("PASS", "WARN", "FAIL")
    }
    return {
        "status": "FAIL" if counts["FAIL"] else "WARN" if counts["WARN"] else "PASS",
        "project_root": str(root),
        "workflow": workflow,
        "species_count": species_count,
        "counts": counts,
        "diamond_settings": {
            "evalue": os.environ.get("DIAMOND_EVALUE", "1e-5"),
            "sensitivity": os.environ.get("DIAMOND_SENSITIVITY", "very-sensitive"),
        },
        "checks": [check.as_dict() for check in checks],
    }
