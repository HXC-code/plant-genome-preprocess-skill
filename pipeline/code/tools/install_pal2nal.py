#!/usr/bin/env python3
"""Download PAL2NAL v14 from its official site without bundling it in releases."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import shutil
import tarfile
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


DEFAULT_URL = "https://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz"
DEFAULT_SHA256 = "0e458d38571ab6d61078047a9e80d843bb2e733414eeff8b1144d1c6bc13848f"
MEMBERS = ("pal2nal.v14/pal2nal.pl", "pal2nal.v14/README")


def download(url: str) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "plant-genome-preprocess/install-pal2nal"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read()


def install(url: str, expected_sha256: str, output: Path, force: bool) -> Path:
    if output.exists() and not force:
        raise FileExistsError(f"refusing to overwrite {output}; use --force")
    payload = download(url)
    actual = hashlib.sha256(payload).hexdigest()
    if actual.lower() != expected_sha256.lower():
        raise ValueError(
            f"SHA-256 mismatch: expected {expected_sha256}, received {actual}"
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}-", dir=output.parent))
    try:
        extracted: dict[str, bytes] = {}
        with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
            for name in MEMBERS:
                try:
                    member = archive.getmember(name)
                except KeyError as exc:
                    raise ValueError(f"archive is missing member: {name}") from exc
                source = archive.extractfile(member)
                if source is None:
                    raise ValueError(f"archive member is not a file: {name}")
                extracted[Path(name).name] = source.read()
        for name, data in extracted.items():
            (staging / name).write_bytes(data)
        metadata = {
            "name": "PAL2NAL",
            "version": "v14",
            "source_url": url,
            "archive_sha256": actual,
            "installed_at": datetime.now(timezone.utc).isoformat(),
            "redistribution": (
                "Downloaded directly from the upstream site by the user. "
                "The workflow release does not bundle this file."
            ),
            "license_note": (
                "The v14 upstream archive contains no explicit license file. "
                "Do not assume the MIT license of the newer GitHub port applies "
                "to this historical v14 distribution."
            ),
        }
        (staging / "SOURCE.json").write_text(
            json.dumps(metadata, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        if output.exists():
            shutil.rmtree(output)
        os.replace(staging, output)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "external" / "pal2nal.v14",
    )
    parser.add_argument("--url", default=DEFAULT_URL)
    parser.add_argument("--sha256", default=DEFAULT_SHA256)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    try:
        output = install(args.url, args.sha256, args.output.resolve(), args.force)
    except (OSError, ValueError, tarfile.TarError) as exc:
        parser.exit(1, f"错误：{exc}\n")
    print(f"PAL2NAL v14 installed: {output / 'pal2nal.pl'}")
    print("Set WGDI_PAL2NAL only if you want to use a different installation.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
