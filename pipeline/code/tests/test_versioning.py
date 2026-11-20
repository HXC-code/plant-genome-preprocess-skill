"""Version provenance must be single-valued.

The 2.3.0 release recorded three different version strings: the package reported
2.3.0 while ``*.qc.json`` said 2.1.0 and ``*.conversion.json`` said 2.2.0.  That
made it impossible to tell which implementation produced an artefact, so all
implementation versions now derive from ``pgci_preprocess.version``.

The tests run in both layouts: the canonical project (``src/pgci_preprocess``)
and the packaged skill (``scripts/pgci_preprocess``).
"""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import pytest

from pgci_preprocess import __version__
from pgci_preprocess.annotation_convert import CONVERTER_VERSION
from pgci_preprocess.cli import VERSION as CLI_VERSION
from pgci_preprocess.pipeline import PIPELINE_VERSION

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_DIR = next(
    root
    for root in (
        PROJECT_ROOT / "src" / "pgci_preprocess",
        PROJECT_ROOT / "scripts" / "pgci_preprocess",
    )
    if root.is_dir()
)


def test_package_version_is_the_single_source() -> None:
    assert PIPELINE_VERSION == __version__
    assert CONVERTER_VERSION == __version__
    assert CLI_VERSION == __version__


def test_package_public_api_is_reachable_from_the_package_object() -> None:
    """The project-root launcher stub must expose the canonical public API.

    pytest puts the project root on ``sys.path`` before ``src``, so importing the
    package there used to yield a stub that only extended ``__path__``: submodule
    imports worked but package-level names were missing.
    """

    import pgci_preprocess

    assert pgci_preprocess.__version__ == __version__
    for name in pgci_preprocess.__all__:
        assert hasattr(pgci_preprocess, name), name


def test_version_module_has_no_package_imports() -> None:
    """``version.py`` must stay importable from anywhere without a cycle."""

    text = (PACKAGE_DIR / "version.py").read_text(encoding="utf-8")
    imports = [
        line
        for line in re.findall(r"^(?:from|import)\s.*", text, re.MULTILINE)
        if not line.startswith("from __future__")
    ]
    assert imports == []


def test_stale_version_constants_are_gone() -> None:
    """No module may hard-code an implementation version again."""

    offenders: list[str] = []
    for path in sorted(PACKAGE_DIR.glob("*.py")):
        if path.name == "version.py":
            continue
        text = path.read_text(encoding="utf-8")
        for match in re.finditer(
            r"^(?:PIPELINE_VERSION|CONVERTER_VERSION|VERSION)\s*=\s*[\"']",
            text,
            re.MULTILINE,
        ):
            offenders.append(f"{path.name}:{match.group(0).strip()}")
    assert offenders == []


def test_pyproject_version_matches_package_when_declared() -> None:
    pyproject = PROJECT_ROOT / "pyproject.toml"
    if not pyproject.is_file():
        pytest.skip("no pyproject.toml in this layout")
    text = pyproject.read_text(encoding="utf-8")
    match = re.search(r"^version\s*=\s*\"([^\"]+)\"", text, re.MULTILINE)
    if match is None:
        pytest.skip("packaged skill pyproject declares no project version")
    assert match.group(1) == __version__


def test_conversion_manifest_reports_the_package_version(tmp_path: Path) -> None:
    from pgci_preprocess.annotation_convert import convert_annotation

    source = tmp_path / "x.annotations.csv.gz"
    with gzip.open(source, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write("P1\tGO:0000001\t1e-5\t100.0\tdesc\tTestus plantus\tGO\n")

    result = convert_annotation(source, tmp_path / "derived", prefix="demo")
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest["converter_version"] == __version__
