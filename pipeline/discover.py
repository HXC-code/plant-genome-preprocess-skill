#!/usr/bin/env python3
"""Standalone entry point for ``discover`` on a fresh machine.

A new computer may have nothing but the skill folder and the raw data: no
installed package, no virtualenv, no environment variables.  This script adds
``code/src`` to ``sys.path`` and runs the scanner with the interpreter that
invoked it, so:

    python /path/to/skill/pipeline/discover.py /path/to/raw-data

is enough to produce ``resources/machine_paths.json``,
``resources/download_inventory.tsv`` and ``resources/sequence_mode.tsv``.
Every later command then finds the project without ``--project-root``.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "code" / "src"
if not SRC.is_dir():
    raise SystemExit(f"ERROR: cannot find the code tree at {SRC}")
sys.path.insert(0, str(SRC))

from pgci_preprocess.cli import main as cli_main  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    return cli_main(["discover", *arguments])


if __name__ == "__main__":
    raise SystemExit(main())
