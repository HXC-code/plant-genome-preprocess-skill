#!/usr/bin/env python3
"""Collect a fail-closed JSON receipt for an isolated WGDI smoke workspace."""
from __future__ import annotations

import argparse
import hashlib
import json
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("smoke_root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.smoke_root.resolve()
    output = (args.output or root / "smoke_report.json").resolve()
    fixture = json.loads((root / "fixture_manifest.json").read_text(encoding="utf-8"))
    settings = read_key_value(root / "Auto_pep" / "diamond_run_settings.tsv")
    status = read_key_value(root / "smoke_status.tsv", skip_header=False)

    comparisons: list[dict[str, object]] = []
    all_ok = status.get("status") == "PASS"
    for base in sorted((root / "Auto_wgdi").glob("*_base")):
        left = base.name[:-5]
        for comparison in sorted(path for path in base.iterdir() if path.is_dir()):
            right = comparison.name[len(left) + 1:]
            stem = left if left == right else f"{left}_{right}"
            names = {
                "dotplot": f"{stem}.png",
                "collinearity": f"{stem}.icl",
                "ks": f"{stem}.ks",
                "blockinfo": f"{stem}.bi.csv",
                "correspondence": f"{stem}.bi.c.csv",
                "blockks": f"{stem}.bk.png",
                "kspeaks": f"{stem}.kp.png",
                "peaksfit": f"{stem}.pf.png",
            }
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
            source_blast = root / "Auto_pep" / f"{left}_pep" / f"{left}_{right}.blast"
            target_blast = comparison / f"{left}_{right}.blast"
            preserved = source_blast.is_file() and target_blast.is_file()
            if preserved:
                preserved = sha256(source_blast) == sha256(target_blast)
            comparison_ok = comparison_ok and preserved and (comparison / ".wgdi.done").is_file()
            comparisons.append({
                "left": left,
                "right": right,
                "status": "PASS" if comparison_ok else "FAIL",
                "auto_pep_blast_preserved_and_identical": preserved,
                "outputs": outputs,
            })
            all_ok = all_ok and comparison_ok
    all_ok = all_ok and bool(comparisons)
    all_ok = all_ok and settings.get("evalue") == status.get("diamond_evalue")
    all_ok = all_ok and settings.get("sensitivity") == status.get("diamond_sensitivity")
    report = {
        "schema_version": "wgdi-smoke-v1",
        "status": "PASS" if all_ok else "FAIL",
        "created_at": datetime.now().astimezone().isoformat(),
        "scope": "Derived real-data fixture; DIAMOND and scripts 2-6 actually executed. Step 0 is covered separately by preprocessing tests/validation.",
        "not_claimed": "This is not a full WGDI rerun of the 109-species dataset or a new biological WGD inference.",
        "smoke_root": str(root),
        "fixture": fixture["selection"],
        "diamond": settings,
        "comparisons": comparisons,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "comparisons": len(comparisons),
        "key_outputs_per_comparison": 8,
        "report": str(output),
    }, ensure_ascii=False, indent=2))
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
