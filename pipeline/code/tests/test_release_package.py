from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
import zipfile
from pathlib import Path

from pgci_preprocess import __version__


PROJECT_ROOT = Path(__file__).resolve().parents[1]
WORKFLOW_ROOT = PROJECT_ROOT.parent
BUILDER = PROJECT_ROOT / "tools" / "build_release_package.py"
SYNC = PROJECT_ROOT / "tools" / "sync_skill_package.py"
RECEIPT_NAME = f"smoke-report-{__version__}.json"


def _load_builder():
    spec = importlib.util.spec_from_file_location("build_release_package", BUILDER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_release_archive_has_handoff_material_and_no_vendored_wgdi(
    tmp_path: Path,
) -> None:
    output = tmp_path / "release.zip"
    result = subprocess.run(
        [sys.executable, str(BUILDER), "--output", str(output)],
        text=True,
        encoding="utf-8",
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    with zipfile.ZipFile(output) as archive:
        names = set(archive.namelist())
        assert not archive.testzip()
        assert "wgdi_pipeline/README.md" in names
        assert "wgdi_pipeline/.gitignore" in names
        assert "wgdi_pipeline/LICENSE" in names
        assert "wgdi_pipeline/VALIDATION.md" in names
        assert "wgdi_pipeline/THIRD_PARTY_NOTICES.md" in names
        assert "wgdi_pipeline/run_homology.sh" in names
        assert "wgdi_pipeline/code/tools/generate_wgdi_configs.py" in names
        assert "wgdi_pipeline/code/tools/pal2nal.v14/pal2nal.pl" in names
        assert "wgdi_pipeline/PACKAGE_SHA256.tsv" in names
        assert "wgdi_pipeline/PACKAGE_INFO.json" in names
        assert (
            f"wgdi_pipeline/validation/{RECEIPT_NAME}"
            in names
        )
        assert not any("/vendor/" in name for name in names)
        for name in names:
            arc = name.removeprefix("wgdi_pipeline/")
            if arc.endswith(".sh") or arc == "bin/wgdi":
                assert b"\r\n" not in archive.read(name), (
                    f"CRLF in release shell entry: {name}"
                )
        info = json.loads(
            archive.read("wgdi_pipeline/PACKAGE_INFO.json").decode("utf-8")
        )
        assert info["file_count"] > 100
        assert info["release_gate"] == "VALIDATION.md"

    recovered = tmp_path / "recovered"
    with zipfile.ZipFile(output) as archive:
        archive.extractall(recovered)
    env = {
        **os.environ,
        "PYTHONPATH": str(recovered / "wgdi_pipeline" / "code" / "src"),
        "PYTHONDONTWRITEBYTECODE": "1",
    }
    cli = subprocess.run(
        [sys.executable, "-m", "pgci_preprocess", "--version"],
        text=True,
        encoding="utf-8",
        capture_output=True,
        check=False,
        env=env,
    )
    assert cli.returncode == 0, cli.stdout + cli.stderr
    assert cli.stdout.strip()


def test_skill_mirror_uses_the_same_shared_wrappers(tmp_path: Path) -> None:
    pipeline = tmp_path / "skill" / "pipeline"
    result = subprocess.run(
        [sys.executable, str(SYNC), "--pipeline-target", str(pipeline)],
        text=True,
        encoding="utf-8",
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert (pipeline / "run_homology.sh").is_file()
    assert "--mode different" in (pipeline / "4run_wgdi_different.sh").read_text(
        encoding="utf-8"
    )
    assert (pipeline / "validation" / RECEIPT_NAME).is_file()
    assert (pipeline / "code" / "docs" / "HANDOFF.md").is_file()
    assert not (pipeline / "code" / "docs" / "REVIEW-README.md").exists()
    assert (pipeline / "code" / "resources" / "species_abbreviations.tsv").is_file()
    assert f'version = "{__version__}"' in (
        pipeline / "code" / "pyproject.toml"
    ).read_text(encoding="utf-8")
    assert not (pipeline / "vendor").exists()


def test_release_normalizes_crlf_shell_scripts_even_from_windows_checkout(
    tmp_path: Path,
) -> None:
    builder = _load_builder()
    script = tmp_path / "0check_data.sh"
    script.write_bytes(b"#!/usr/bin/env bash\r\necho ok\r\n")
    assert builder.release_bytes(script, "0check_data.sh") == (
        b"#!/usr/bin/env bash\necho ok\n"
    )
    launcher = tmp_path / "wgdi"
    launcher.write_bytes(b"#!/usr/bin/env bash\r\nexec python\r\n")
    assert builder.release_bytes(launcher, "bin/wgdi") == (
        b"#!/usr/bin/env bash\nexec python\n"
    )
    third_party = tmp_path / "pal2nal.pl"
    third_party.write_bytes(b"#!/usr/bin/perl\r\n# upstream bytes\r\n")
    assert builder.release_bytes(third_party, "code/tools/pal2nal.v14/pal2nal.pl") == (
        b"#!/usr/bin/perl\r\n# upstream bytes\r\n"
    )
