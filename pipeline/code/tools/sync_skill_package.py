#!/usr/bin/env python3
"""Synchronize packaged skill copies with the canonical source package.

The canonical implementation is ``src/pgci_preprocess``.  Release packages
must not be hand-edited copies; run this tool after every canonical change.

Usage::

    python tools/sync_skill_package.py --check
    python tools/sync_skill_package.py
    python tools/sync_skill_package.py --target /path/to/scripts/pgci_preprocess

``--check`` exits 1 on drift; ``--sync`` (default) copies canonical ``*.py``
files and tests, and removes stale modules from the target directories.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CANONICAL = PROJECT_ROOT / "src" / "pgci_preprocess"
CANONICAL_TESTS = PROJECT_ROOT / "tests"
DEFAULT_TARGET = (
    PROJECT_ROOT.parent
    / "skill-packages"
    / "plant-genome-wgdi-preprocess"
    / "scripts"
    / "pgci_preprocess"
)


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _python_files(directory: Path) -> dict[str, Path]:
    if not directory.is_dir():
        return {}
    return {
        path.name: path
        for path in directory.glob("*.py")
        if path.is_file() and "__pycache__" not in path.parts
    }


def _canonical_files() -> dict[str, Path]:
    return _python_files(CANONICAL)


def _canonical_test_files() -> dict[str, Path]:
    return _python_files(CANONICAL_TESTS)


def _test_target(target: Path) -> Path:
    return target.parent.parent / "tests"


def _target_files(target: Path) -> dict[str, Path]:
    return _python_files(target)


def compare(target: Path) -> list[str]:
    problems: list[str] = []
    canonical = _canonical_files()
    existing = _target_files(target)
    for name, source in sorted(canonical.items()):
        destination = existing.get(name)
        if destination is None:
            problems.append(f"missing: {name}")
        elif _digest(source) != _digest(destination):
            problems.append(f"drift: {name}")
    for name in sorted(set(existing) - set(canonical)):
        problems.append(f"stale: {name}")

    canonical_tests = _canonical_test_files()
    existing_tests = _target_files(_test_target(target))
    for name, source in sorted(canonical_tests.items()):
        destination = existing_tests.get(name)
        if destination is None:
            problems.append(f"tests/missing: {name}")
        elif _digest(source) != _digest(destination):
            problems.append(f"tests/drift: {name}")
    for name in sorted(set(existing_tests) - set(canonical_tests)):
        problems.append(f"tests/stale: {name}")
    return problems


def sync(target: Path) -> list[str]:
    canonical = _canonical_files()
    existing = _target_files(target)
    target.mkdir(parents=True, exist_ok=True)
    removed: list[str] = []
    for name in sorted(set(existing) - set(canonical)):
        existing[name].unlink()
        removed.append(name)
    for name, source in sorted(canonical.items()):
        shutil.copy2(source, target / name)

    test_target = _test_target(target)
    test_target.mkdir(parents=True, exist_ok=True)
    canonical_tests = _canonical_test_files()
    existing_tests = _target_files(test_target)
    for name in sorted(set(existing_tests) - set(canonical_tests)):
        existing_tests[name].unlink()
        removed.append(f"tests/{name}")
    for name, source in sorted(canonical_tests.items()):
        shutil.copy2(source, test_target / name)
    return removed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, action="append")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    targets = args.target or [DEFAULT_TARGET]
    status = 0
    for target in targets:
        target = target.resolve()
        if args.check:
            problems = compare(target)
            if problems:
                status = 1
                print(f"[DRIFT] {target}")
                for problem in problems:
                    print(f"  {problem}")
            else:
                print(f"[OK] {target}")
        else:
            removed = sync(target)
            print(f"[SYNCED] {target} files={len(_canonical_files())} tests={len(_canonical_test_files())} removed={len(removed)}")
            for name in removed:
                print(f"  removed {name}")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
