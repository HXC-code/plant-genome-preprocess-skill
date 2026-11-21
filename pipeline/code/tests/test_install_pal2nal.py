from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import tarfile
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INSTALLER = PROJECT_ROOT / "tools" / "install_pal2nal.py"


def _installer_module():
    spec = importlib.util.spec_from_file_location("install_pal2nal", INSTALLER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _archive() -> bytes:
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w:gz") as archive:
        for name, content in (
            ("pal2nal.v14/pal2nal.pl", "#!/usr/bin/env perl\n"),
            ("pal2nal.v14/README", "PAL2NAL test fixture\n"),
            ("pal2nal.v14/unused.txt", "not installed\n"),
        ):
            data = content.encode()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
    return payload.getvalue()


def test_installer_verifies_checksum_and_extracts_only_runtime_files(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = _installer_module()
    payload = _archive()
    monkeypatch.setattr(module, "download", lambda _url: payload)
    digest = hashlib.sha256(payload).hexdigest()
    output = module.install("https://example.test/pal2nal.tar.gz", digest, tmp_path / "pal2nal", False)
    assert (output / "pal2nal.pl").read_text(encoding="utf-8").startswith("#!/usr/bin")
    assert (output / "README").is_file()
    assert not (output / "unused.txt").exists()
    metadata = json.loads((output / "SOURCE.json").read_text(encoding="utf-8"))
    assert metadata["archive_sha256"] == digest
    assert "does not bundle" in metadata["redistribution"]


def test_installer_rejects_checksum_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = _installer_module()
    monkeypatch.setattr(module, "download", lambda _url: _archive())
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        module.install("https://example.test/pal2nal.tar.gz", "0" * 64, tmp_path / "pal2nal", False)
