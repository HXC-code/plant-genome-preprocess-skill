from __future__ import annotations

import getpass
import os
import sys
import tempfile
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for candidate in (PROJECT_ROOT / "src", PROJECT_ROOT / "scripts"):
    if (candidate / "pgci_preprocess").is_dir():
        sys.path.insert(0, str(candidate))
        break

_BASETEMP_LIMIT = 90


def _usable_temp_dir(directory: Path) -> bool:
    """Return True when pytest can list and write inside ``directory``."""
    try:
        if directory.exists():
            next(directory.iterdir(), None)
        directory.mkdir(parents=True, exist_ok=True)
        probe = directory / f".pgci-probe-{os.getpid()}"
        probe.write_text("ok", encoding="ascii")
        probe.unlink()
        return True
    except OSError:
        return False


def pytest_configure(config: pytest.Config) -> None:
    """Keep pytest temp projects writable and inside the Windows MAX_PATH budget.

    A deep ``--basetemp`` pushes ``tmp_path`` projects past the legacy 260-character
    limit: directories are created but files cannot be written, producing dozens
    of misleading failures.  A stale or unreadable default ``pytest-of-<user>``
    directory causes the same class of mass failures, so fall back to a fresh
    short temp root whenever the default one is not usable.
    """

    basetemp = getattr(config.option, "basetemp", None)
    if basetemp:
        if len(str(Path(basetemp).resolve())) <= _BASETEMP_LIMIT:
            return
        try:
            short = Path(tempfile.mkdtemp(prefix="pgci_pytest_"))
        except OSError:
            return
        config.option.basetemp = str(short)
        config.issue_config_time_warning(
            pytest.PytestWarning(
                f"--basetemp path is too deep for Windows MAX_PATH 260: {basetemp}; "
                f"redirected to {short}"
            ),
            stacklevel=2,
        )
        return

    default_root = Path(tempfile.gettempdir()) / f"pytest-of-{getpass.getuser()}"
    if _usable_temp_dir(default_root):
        return
    try:
        short = Path(tempfile.mkdtemp(prefix="pgci_pytest_"))
    except OSError:
        return
    config.option.basetemp = str(short)
    config.issue_config_time_warning(
        pytest.PytestWarning(
            f"default pytest temp root is not writable: {default_root}; "
            f"redirected to {short}"
        ),
        stacklevel=2,
    )
