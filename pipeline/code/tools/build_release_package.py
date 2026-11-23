#!/usr/bin/env python3
"""Build the self-contained pgci_preprocess + 0-3 WGDI workflow source release.

Source layout (this repository):

    wgdi_pipeline/                  pipeline entry scripts and bin/
    wgdi_pipeline/code/             preprocessing source (src/, scripts/, tools/, tests/)

Release layout (what the other person receives):

    wgdi_pipeline/
        README.md, LICENSE, QUICKSTART.md, VALIDATION.md, *.sh, *.ps1,
        requirements*.txt, bin/
        code/
            src/pgci_preprocess/...
            scripts/ tools/ tests/ docs/
            resources/species_abbreviations.tsv
            pyproject.toml
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

CODE_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_ROOT = CODE_ROOT.parent
sys.path.insert(0, str(CODE_ROOT / "src"))
from pgci_preprocess.version import __version__ as VERSION  # noqa: E402

TOP_LEVEL_FILES = (
    ".gitattributes",
    ".gitignore",
    "README.md",
    "QUICKSTART.md",
    "VALIDATION.md",
    "CHANGELOG.md",
    "CONTRIBUTING.md",
    "LICENSE",
    "THIRD_PARTY_NOTICES.md",
    "run_all_pipeline.ps1",
    "run_homology.sh",
    "bootstrap_wgdi.ps1",
    "bootstrap_wgdi.sh",
    "requirements-wgdi-validated.txt",
    "discover.py",
    "0check_data.sh",
    "_pgci_env.sh",
    "1run_blast.sh",
    "1run_diamond.sh",
    "2gen_wgdi_conf.sh",
    "3run_wgdi.sh",
    "bin/wgdi",
    "bin/wgdi_entry.py",
)
CODE_FILES = (
    "README.md",
    "pyproject.toml",
    "resources/species_abbreviations.tsv",
    "resources/format_exceptions.tsv",
)
CODE_GLOBS = (
    "src/pgci_preprocess/*.py",
    "pgci_preprocess/*.py",
    "scripts/*.py",
    "tests/*.py",
    "tools/*.py",
    "tools/*.sh",
    "docs/**/*.md",
)
PIPELINE_GLOBS = ("validation/*.json",)

# bash on Linux/WSL rejects CRLF shell scripts. Windows checkouts (autocrlf)
# and Windows-built archives must therefore never carry CRLF into a release.
_SHELL_SUFFIXES = (".sh",)
_SHELL_NAMES = ("bin/wgdi",)


def _is_shell_entry(arc: str) -> bool:
    return arc.endswith(_SHELL_SUFFIXES) or arc in _SHELL_NAMES


def release_bytes(path: Path, arc: str) -> bytes:
    """Return the exact bytes stored in the archive for one entry.

    Shell scripts are normalized to LF so a release built on a CRLF Windows
    checkout still runs under Linux/WSL bash. All other files are stored
    byte-for-byte.
    """
    data = path.read_bytes()
    if _is_shell_entry(arc):
        data = data.replace(b"\r\n", b"\n")
    return data


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def collect_files() -> list[tuple[Path, str]]:
    """Return (source path, archive-relative path) pairs."""
    files: dict[str, Path] = {}
    for name in TOP_LEVEL_FILES:
        path = PIPELINE_ROOT / name
        if not path.is_file():
            raise FileNotFoundError(path)
        files[name] = path
    for name in CODE_FILES:
        path = CODE_ROOT / name
        if not path.is_file():
            raise FileNotFoundError(path)
        files[f"code/{name}"] = path
    for pattern in CODE_GLOBS:
        for path in CODE_ROOT.glob(pattern):
            if path.is_file():
                files[f"code/{path.relative_to(CODE_ROOT).as_posix()}"] = path
    for pattern in PIPELINE_GLOBS:
        for path in PIPELINE_ROOT.glob(pattern):
            if path.is_file():
                files[path.relative_to(PIPELINE_ROOT).as_posix()] = path
    return [(path, arc) for arc, path in sorted(files.items())]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"refusing to overwrite {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    entries = collect_files()
    manifest_lines = ["path\tsha256\tbytes"]
    payloads: list[tuple[str, bytes, int]] = []
    for path, arc in entries:
        data = release_bytes(path, arc)
        digest = hashlib.sha256(data).hexdigest()
        manifest_lines.append(f"{arc}\t{digest}\t{len(data)}")
        mode = 0o755 if _is_shell_entry(arc) else 0o644
        payloads.append((arc, data, mode))
    info = {
        "name": "pgci-preprocess-wgdi-workflow",
        "version": VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": (
            "Canonical preprocessing source, tests, 0-3 execution scripts, wrappers, "
            "release checks and installation guides; no raw data or generated results."
        ),
        "file_count": len(entries),
        "entry_layout": "wgdi_pipeline/code/<source tree>; launchers at wgdi_pipeline/<name>",
        "release_gate": "VALIDATION.md",
    }
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for (arc, data, mode), (path, _arc) in zip(payloads, entries, strict=True):
            zip_info = zipfile.ZipInfo(
                "wgdi_pipeline/" + arc,
                date_time=time.localtime(path.stat().st_mtime)[:6],
            )
            zip_info.create_system = 3  # Unix: preserve the executable bit on extraction
            zip_info.external_attr = mode << 16
            zip_info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(zip_info, data, compresslevel=9)
        archive.writestr("wgdi_pipeline/PACKAGE_SHA256.tsv", "\n".join(manifest_lines) + "\n")
        archive.writestr("wgdi_pipeline/PACKAGE_INFO.json", json.dumps(info, ensure_ascii=False, indent=2) + "\n")
    with zipfile.ZipFile(output) as archive:
        bad = archive.testzip()
        if bad:
            raise ValueError(f"ZIP CRC validation failed: {bad}")
    print(json.dumps({
        "output": str(output),
        "version": VERSION,
        "source_files": len(entries),
        "zip_entries": len(entries) + 2,
        "bytes": output.stat().st_size,
        "sha256": sha256(output),
        "crc_check": "PASS",
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
