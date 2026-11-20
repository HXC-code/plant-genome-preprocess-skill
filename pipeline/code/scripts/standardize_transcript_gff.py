#!/usr/bin/env python3
"""Deprecated compatibility wrapper for ``python -m pgci_preprocess normalize-gff``."""

from __future__ import annotations

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR.parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pgci_preprocess.cli import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main(["normalize-gff", *sys.argv[1:]]))
