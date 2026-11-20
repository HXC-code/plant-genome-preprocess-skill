#!/usr/bin/env python3
"""Deprecated compatibility wrapper for ``python -m pgci_preprocess prepare``."""

from __future__ import annotations

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pgci_preprocess.cli import main  # noqa: E402


def _legacy_argv() -> list[str]:
    argv = list(sys.argv[1:])
    if any(item == "--source-status" or item.startswith("--source-status=") for item in argv):
        return ["prepare", *argv]
    default_status = PROJECT_ROOT.parents[2] / "01_DATA_DOWNLOAD" / "resources" / "current_gap_jobs.tsv"
    return ["prepare", "--source-status", str(default_status), *argv]


if __name__ == "__main__":
    raise SystemExit(main(_legacy_argv()))
