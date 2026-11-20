from __future__ import annotations

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


def pytest_configure(config: pytest.Config) -> None:
    """Keep pytest temp projects inside the Windows MAX_PATH budget.

    A deep ``--basetemp`` (or the default, when the repository itself lives in
    a deep folder) pushes ``tmp_path`` projects past the legacy 260-character
    limit: directories are created but files cannot be written, producing
    dozens of misleading failures.  Redirect to a short temp root instead.
    """

    basetemp = getattr(config.option, "basetemp", None)
    if not basetemp:
        return
    if len(str(Path(basetemp).resolve())) <= _BASETEMP_LIMIT:
        return
    short = Path(tempfile.mkdtemp(prefix="pgci_pytest_"))
    config.option.basetemp = str(short)
    config.issue_config_time_warning(
        pytest.PytestWarning(
            f"--basetemp path is too deep for Windows MAX_PATH 260: {basetemp}; "
            f"redirected to {short}"
        ),
        stacklevel=2,
    )
