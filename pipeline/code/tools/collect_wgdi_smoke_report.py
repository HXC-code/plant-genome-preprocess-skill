#!/usr/bin/env python3
"""Collect a fail-closed JSON receipt for an isolated WGDI smoke workspace."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_key_value(path: Path, *, skip_header: bool = True) -> dict[str, str]:
    result: dict[str, str] = {}
    lines = path.read_text(encoding="utf-8").splitlines()
    for line in lines[1:] if skip_header else lines:
        if line.strip():
            key, value = line.split("\t", 1)
            result[key] = value
    return result


def command_output(argv: list[str]) -> str:
    try:
        result = subprocess.run(
            argv,
            text=True,
            encoding="utf-8",
            errors="replace",
            capture_output=True,
            timeout=15,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    text = (result.stdout + "\n" + result.stderr).strip()
    return text.splitlines()[0] if text else ""


def binary_version_hint(path: str) -> str:
    try:
        result = subprocess.run(
            ["strings", path],
            text=True,
            encoding="utf-8",
            errors="replace",
            capture_output=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    match = re.search(r"YN00 in paml version[^\r\n]*", result.stdout)
    return match.group(0).strip() if match else ""


def runtime_environment(
    settings: dict[str, str], *, public: bool
) -> dict[str, object]:
    wgdi_python = os.environ.get("WGDI_PYTHON") or sys.executable
    wgdi_version = command_output(
        [
            wgdi_python,
            "-c",
            "import importlib.metadata; print(importlib.metadata.version('wgdi'))",
        ]
    )
    wgdi_python_version = command_output(
        [wgdi_python, "-c", "import platform; print(platform.python_version())"]
    )
    tools: dict[str, object] = {}
    workflow_root = Path(__file__).resolve().parents[2]
    perl_bin = os.environ.get("WGDI_PERL_BIN", "")
    perl = os.environ.get("WGDI_PERL") or (
        str(Path(perl_bin) / "perl") if perl_bin else shutil.which("perl")
    )
    pal2nal = (
        os.environ.get("WGDI_PAL2NAL")
        or os.environ.get("PAL2NAL")
        or str(workflow_root / "external" / "pal2nal.v14" / "pal2nal.pl")
    )
    probes = {
        "diamond": (
            settings.get("executable") or shutil.which("diamond") or "",
            ["version"],
        ),
        "blastp": (shutil.which("blastp") or "", ["-version"]),
        "makeblastdb": (shutil.which("makeblastdb") or "", ["-version"]),
        "muscle": (
            os.environ.get("WGDI_MUSCLE")
            or os.environ.get("MUSCLE")
            or shutil.which("muscle")
            or "",
            ["-version"],
        ),
        "yn00": (
            os.environ.get("WGDI_YN00")
            or os.environ.get("YN00")
            or shutil.which("yn00")
            or "",
            [],
        ),
        "perl": (perl or "", ["-v"]),
        "pal2nal": (pal2nal if Path(pal2nal).is_file() else "", []),
    }
    for name, (path, arguments) in probes.items():
        version = os.environ.get(f"WGDI_{name.upper()}_VERSION", "")
        if not version:
            version = command_output([path, *arguments]) if path else ""
        if name == "diamond" and not version:
            version = settings.get("version", "")
        if name == "pal2nal" and path:
            version = version or "v14 (downloaded from upstream)"
        if name == "yn00" and not version and path:
            version = binary_version_hint(path)
        tools[name] = {
            "path": Path(path).name if public and path else path,
            "version": version,
        }
    return {
        "os": platform.platform(),
        "machine": platform.machine(),
        "python": platform.python_version(),
        "wgdi_python": Path(wgdi_python).name if public else wgdi_python,
        "wgdi_python_version": wgdi_python_version,
        "wgdi": wgdi_version,
        "tools": tools,
    }


def comparison_pairs(root: Path) -> list[tuple[str, str]]:
    """Derive (left, right) pairs from the homology output layout.

    ``wgdi_input/blast/<left>_pep/<left>_<right>.blast`` is the authoritative
    BLAST location; results live in the single ``wgdi_output/results`` folder.
    """

    pairs: list[tuple[str, str]] = []
    for pep_dir in sorted((root / "wgdi_input" / "blast").glob("*_pep")):
        if not pep_dir.is_dir():
            continue
        left = pep_dir.name[: -len("_pep")]
        for blast in sorted(pep_dir.glob("*.blast")):
            stem = blast.stem
            right = stem[len(left) + 1 :] if stem.startswith(f"{left}_") else ""
            if right:
                pairs.append((left, right))
    return pairs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("smoke_root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--public",
        action="store_true",
        help="Remove host-specific absolute paths from the receipt.",
    )
    args = parser.parse_args()
    root = args.smoke_root.resolve()
    output = (args.output or root / "smoke_report.json").resolve()
    fixture = json.loads((root / "fixture_manifest.json").read_text(encoding="utf-8"))
    settings = read_key_value(root / "wgdi_input" / "blast" / "diamond_run_settings.tsv")
    status = read_key_value(root / "smoke_status.tsv", skip_header=False)

    comparisons: list[dict[str, object]] = []
    all_ok = status.get("status") == "PASS"
    ks_enabled = status.get("ks_enabled", "true") == "true"
    for left, right in comparison_pairs(root):
        stem = left if left == right else f"{left}_{right}"
        comparison = root / "wgdi_output" / "results"
        state = root / "wgdi_output" / "state"
        names = {
            "dotplot": f"{stem}.png",
            "collinearity": f"{stem}.icl",
        }
        markers = [".done.dotplot", ".done.collinearity"]
        if ks_enabled:
            names.update({
                "ks": f"{stem}.ks",
                "blockinfo": f"{stem}.bi.csv",
                "correspondence": f"{stem}.bi.c.csv",
                "blockks": f"{stem}.bk.png",
                "kspeaks": f"{stem}.kp.png",
                "peaksfit": f"{stem}.pf.png",
            })
            markers += [
                ".done.ks",
                ".done.blockinfo",
                ".done.correspondence",
                ".done.blockks",
                ".done.kspeaks",
                ".done.peaksfit",
            ]
        outputs: dict[str, object] = {}
        comparison_ok = True
        for stage, name in names.items():
            path = comparison / name
            ok = path.is_file() and path.stat().st_size > 0
            comparison_ok = comparison_ok and ok
            outputs[stage] = {
                "path": str(path.relative_to(root)),
                "exists_nonempty": ok,
                "bytes": path.stat().st_size if path.is_file() else 0,
                "sha256": sha256(path) if ok else "",
            }
        # Single authoritative BLAST table under wgdi_input/blast.
        source_blast = root / "wgdi_input" / "blast" / f"{left}_pep" / f"{left}_{right}.blast"
        blast_ok = source_blast.is_file() and source_blast.stat().st_size > 0
        markers_ok = all(
            (state / f"{stem}{marker}").is_file() for marker in markers
        )
        comparison_ok = comparison_ok and blast_ok and markers_ok
        comparisons.append({
            "left": left,
            "right": right,
            "status": "PASS" if comparison_ok else "FAIL",
            "auto_pep_blast_present": blast_ok,
            "stage_markers_present": markers_ok,
            "outputs": outputs,
        })
        all_ok = all_ok and comparison_ok
    all_ok = all_ok and bool(comparisons)
    all_ok = all_ok and settings.get("evalue") == status.get("diamond_evalue")
    all_ok = all_ok and settings.get("sensitivity") == status.get("diamond_sensitivity")
    synthetic = "synthetic" in fixture["selection"].get("method", "").lower()
    report = {
        "schema_version": "wgdi-smoke-v4",
        "status": "PASS" if all_ok else "FAIL",
        "created_at": datetime.now().astimezone().isoformat(),
        "workflow_version": fixture.get("workflow_version", "unknown"),
        "scope": (
            "Isolated synthetic fixture" if synthetic else "Derived real-data fixture"
        )
        + "; DIAMOND, 2gen_wgdi_conf.sh and 3run_wgdi.sh actually executed. Step 0 is covered separately by preprocessing tests/validation.",
        "not_claimed": "This is not a full WGDI rerun of the 109-species dataset or a new biological WGD inference.",
        "smoke_root": "<isolated-workspace>" if args.public else str(root),
        "ks_enabled": ks_enabled,
        "fixture": fixture["selection"],
        "diamond": settings,
        "environment": runtime_environment(settings, public=args.public),
        "comparisons": comparisons,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "comparisons": len(comparisons),
        "key_outputs_per_comparison": 8 if ks_enabled else 2,
        "report": str(output),
    }, ensure_ascii=False, indent=2))
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
