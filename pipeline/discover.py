#!/usr/bin/env python3
"""Standalone entry point for ``discover`` on a fresh machine.

A new computer may have nothing but the skill folder and the raw data: no
installed package, no virtualenv, no environment variables.  This script adds
``code/src`` to ``sys.path`` and runs the scanner with the interpreter that
invoked it, so:

    python /path/to/skill/pipeline/discover.py /path/to/raw-data

is enough to produce the project metadata and the path map.  Every later
command then finds the project without ``--project-root``.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "code" / "src"
if not SRC.is_dir():
    raise SystemExit(f"ERROR: cannot find the code tree at {SRC}")
sys.path.insert(0, str(SRC))

from pgci_preprocess.cli import main as cli_main  # noqa: E402

from pgci_preprocess.machine_scan import pointer_candidates  # noqa: E402


def tidy_metadata_directory(project_root: Path) -> None:
    """Name generated scan metadata by content and keep a legacy alias.

    The code tree historically reads ``resources/``; it contains only scan
    metadata, so the canonical directory is renamed ``metadata/`` and a
    zero-copy link (symlink, or a Windows directory junction) keeps the old
    name working.  Failure is non-fatal: the plain directory remains valid.
    """

    resources = project_root / "resources"
    metadata = project_root / "metadata"
    if not resources.is_dir() or resources.is_symlink() or metadata.exists():
        return
    try:
        resources.rename(metadata)
    except OSError:
        return
    try:
        os.symlink("metadata", resources, target_is_directory=True)
        return
    except (OSError, NotImplementedError):
        pass
    if os.name == "nt":
        subprocess.run(
            ["cmd", "/c", f'mklink /J "{resources}" "{metadata}"'],
            check=False,
            capture_output=True,
        )
    # Roll back if no alias could be created: code must keep working.
    if not resources.exists():
        try:
            metadata.rename(resources)
        except OSError:
            pass


def main(argv: list[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    status = cli_main(["discover", *arguments])
    if status == 0:
        for candidate in pointer_candidates():
            try:
                root = Path(json.loads(candidate.read_text(encoding="utf-8"))["project_root"])
            except (OSError, ValueError, KeyError):
                continue
            tidy_metadata_directory(root)
            break
    return status


if __name__ == "__main__":
    raise SystemExit(main())
