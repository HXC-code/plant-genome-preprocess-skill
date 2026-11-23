#!/usr/bin/env python3
"""Generate flat WGDI configs from the homology scan map.

Layout (v5 -- two data directories, no per-comparison folders):

* ``wgdi_input/`` holds everything WGDI reads: species four-piece bundles
  (``<code>/``), DIAMOND tables (``blast/<code>_pep/``) and one flat config per
  comparison (``config/<left>_<right>.conf``).
* ``wgdi_output/`` holds everything produced: ``results/`` (PNG dotplots and
  ICL collinearity blocks), ``logs/`` and ``state/`` (resume markers, Ks
  intermediates).

All references are relative to the project root, so a config stays valid after
the project moves to another machine or drive.  No input file is copied.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

INPUT_DIR = "wgdi_input"
BLAST_DIR = "blast"
CONFIG_DIR = "config"
OUTPUT_DIR = "wgdi_output"
RESULTS_DIR = "results"
LOGS_DIR = "logs"
STATE_DIR = "state"


def read_species_map(project_root: Path) -> dict[str, dict[str, str]]:
    """Parse ``species_map.tsv`` (homology output): code -> scan-resolved paths.

    Columns: ``<code>\t<bundle_dir>\t<pep>\t<cds>\t<gff>\t<lens>``.
    """

    candidates = (
        project_root / "metadata" / "species_map.tsv",
        project_root / "resources" / "species_map.tsv",
        project_root / "species_map.tsv",
    )
    path = next((candidate for candidate in candidates if candidate.is_file()), candidates[0])
    if not path.is_file() or not path.read_text(encoding="utf-8").strip():
        raise FileNotFoundError("species_map.tsv 不存在或为空（先运行第 3 步同源比对）")
    mapping: dict[str, dict[str, str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = [field.strip() for field in line.split("\t")]
        if len(fields) < 6 or not fields[0]:
            continue
        code, directory, pep, cds, gff, lens = fields[:6]
        if code in mapping:
            raise ValueError(f"species_map.tsv 物种重复：{code}")
        if not all((directory, pep, cds, gff, lens)):
            raise ValueError(f"species_map.tsv 行不完整：{line!r}")
        mapping[code] = {
            "dir": directory,
            "pep": pep,
            "cds": cds,
            "gff": gff,
            "lens": lens,
        }
    if not mapping:
        raise ValueError("species_map.tsv 没有可用行")
    return mapping


def bundle_dir(project_root: Path, row: dict[str, str]) -> Path:
    directory = Path(row["dir"])
    resolved = directory if directory.is_absolute() else project_root / directory
    if not resolved.is_dir():
        raise FileNotFoundError(f"物种目录不存在：{resolved}")
    return resolved


def bundle_file(project_root: Path, row: dict[str, str], key: str) -> Path:
    return bundle_dir(project_root, row) / row[key]


def blast_file(project_root: Path, left: str, right: str) -> Path:
    return (
        project_root
        / INPUT_DIR
        / BLAST_DIR
        / f"{left}_pep"
        / f"{left}_{right}.blast"
    )


def rel_ref(project_root: Path, target: Path) -> str:
    """Forward-slash path relative to the project root (the WGDI working dir)."""

    return Path(os.path.relpath(target, project_root)).as_posix()


def _config_text(
    left: str,
    right: str,
    *,
    blast: str,
    gff1: str,
    gff2: str,
    lens1: str,
    lens2: str,
    cds_file: str,
    pep_file: str,
    png: str,
    icl: str,
    ks: str,
    bi: str,
    bic: str,
    bk: str,
    kp: str,
    kpp: str,
    pf: str,
    evalue: str,
    tandem: bool,
    ks_column: str,
    ks_enabled: bool,
) -> str:
    tandem_value = "true" if tandem else "false"
    core = f"""[dotplot]
blast = {blast}
gff1 = {gff1}
gff2 = {gff2}
lens1 = {lens1}
lens2 = {lens2}
genome1_name = {left}
genome2_name = {right}
multiple = 1
score = 100
evalue = {evalue}
repeat_number = 10
position = order
blast_reverse = false
ancestor_left = none
ancestor_top = none
markersize = 0.5
figsize = 10,10
savefig = {png}

[collinearity]
blast = {blast}
gff1 = {gff1}
gff2 = {gff2}
lens1 = {lens1}
lens2 = {lens2}
blast_reverse = false
multiple = 1
process = 8
evalue = {evalue}
score = 100
grading = 50,25,10
mg = 25,25
pvalue = 1
repeat_number = 10
position = order
savefile = {icl}
"""
    ks_part = f"""
[ks]
cds_file = {cds_file}
pep_file = {pep_file}
align_software = muscle
pairs_file = {icl}
ks_file = {ks}

[blockinfo]
blast = {blast}
gff1 = {gff1}
gff2 = {gff2}
lens1 = {lens1}
lens2 = {lens2}
collinearity = {icl}
score = 100
evalue = {evalue}
repeat_number = 20
position = order
ks = {ks}
ks_col = {ks_column}
savefile = {bi}

[correspondence]
blockinfo = {bi}
lens1 = {lens1}
lens2 = {lens2}
tandem = {tandem_value}
tandem_length = 200
pvalue = 0.2
block_length = 5
tandem_ratio = 0.5
multiple = 1
homo = -1,1
savefile = {bic}

[blockks]
lens1 = {lens1}
lens2 = {lens2}
genome1_name = {left}
genome2_name = {right}
blockinfo = {bi}
pvalue = 0.2
tandem = {tandem_value}
tandem_length = 200
markersize = 1
area = 0,2
block_length = 5
figsize = 8,8
savefig = {bk}

[kspeaks]
blockinfo = {bi}
pvalue = 0.2
tandem = {tandem_value}
block_length = 5
ks_area = 0,10
multiple = 1
homo = 0,1
fontsize = 9
area = 0,3
figsize = 10,6.18
savefig = {kp}
savefile = {kpp}

[peaksfit]
blockinfo = {bi}
mode = median
bins_number = 200
ks_area = 0,10
fontsize = 9
area = 0,3
figsize = 10,6.18
shadow = true
savefig = {pf}
"""
    return core + ks_part if ks_enabled else core


def generate(project_root: Path, mode: str) -> int:
    root = project_root.resolve()
    mapping = read_species_map(root)
    species = list(mapping)

    tandem_raw = os.environ.get("WGDI_TANDEM", "false").strip().lower()
    if tandem_raw not in {"true", "false"}:
        raise ValueError("WGDI_TANDEM 必须是 true 或 false")
    tandem = tandem_raw == "true"
    evalue = os.environ.get("WGDI_EVALUE", "1e-5")
    ks_column = os.environ.get("WGDI_KS_COL", "ks_NG86")
    if ks_column not in {"ks_NG86", "ks_YN00"}:
        raise ValueError("WGDI_KS_COL 必须是 ks_NG86 或 ks_YN00")
    ks_raw = os.environ.get("WGDI_KS", "false").strip().lower()
    if ks_raw not in {"true", "false"}:
        raise ValueError("WGDI_KS 必须是 true 或 false")
    ks_enabled = ks_raw == "true"

    config_root = root / INPUT_DIR / CONFIG_DIR
    config_root.mkdir(parents=True, exist_ok=True)
    for name in (RESULTS_DIR, LOGS_DIR, STATE_DIR):
        (root / OUTPUT_DIR / name).mkdir(parents=True, exist_ok=True)

    comparisons: list[tuple[str, str]] = []
    for left in species:
        for right in species:
            is_same = left == right
            if mode == "same" and not is_same:
                continue
            if mode == "different" and is_same:
                continue
            comparisons.append((left, right))

    configured = 0
    for left, right in comparisons:
        prefix = left if left == right else f"{left}_{right}"
        blast = blast_file(root, left, right)
        gff1 = bundle_file(root, mapping[left], "gff")
        gff2 = bundle_file(root, mapping[right], "gff")
        lens1 = bundle_file(root, mapping[left], "lens")
        lens2 = bundle_file(root, mapping[right], "lens")
        cds1 = bundle_file(root, mapping[left], "cds")
        pep1 = bundle_file(root, mapping[left], "pep")
        cds2 = bundle_file(root, mapping[right], "cds")
        pep2 = bundle_file(root, mapping[right], "pep")

        required = [blast, gff1, gff2, lens1, lens2]
        if ks_enabled:
            required += [cds1, pep1]
            if left != right:
                required += [cds2, pep2]
        for path in required:
            if not path.is_file():
                raise FileNotFoundError(f"缺少输入文件：{path}")

        result_root = root / OUTPUT_DIR / RESULTS_DIR
        if left == right:
            cds_file = rel_ref(root, cds1)
            pep_file = rel_ref(root, pep1)
        else:
            # Cross-species Ks needs combined CDS/PEP, built by the executor
            # in wgdi_output/state/<prefix>/ on demand.
            cds_file = rel_ref(root, root / OUTPUT_DIR / STATE_DIR / prefix / "all.cds")
            pep_file = rel_ref(root, root / OUTPUT_DIR / STATE_DIR / prefix / "all.pep")

        (config_root / f"{prefix}.conf").write_text(
            _config_text(
                left,
                right,
                blast=rel_ref(root, blast),
                gff1=rel_ref(root, gff1),
                gff2=rel_ref(root, gff2),
                lens1=rel_ref(root, lens1),
                lens2=rel_ref(root, lens2),
                cds_file=cds_file,
                pep_file=pep_file,
                png=rel_ref(root, result_root / f"{prefix}.png"),
                icl=rel_ref(root, result_root / f"{prefix}.icl"),
                ks=rel_ref(root, result_root / f"{prefix}.ks"),
                bi=rel_ref(root, root / OUTPUT_DIR / RESULTS_DIR / f"{prefix}.bi.csv"),
                bic=rel_ref(root, root / OUTPUT_DIR / RESULTS_DIR / f"{prefix}.bi.c.csv"),
                bk=rel_ref(root, root / OUTPUT_DIR / RESULTS_DIR / f"{prefix}.bk.png"),
                kp=rel_ref(root, root / OUTPUT_DIR / RESULTS_DIR / f"{prefix}.kp.png"),
                kpp=rel_ref(root, root / OUTPUT_DIR / RESULTS_DIR / f"{prefix}.kp"),
                pf=rel_ref(root, root / OUTPUT_DIR / RESULTS_DIR / f"{prefix}.pf.png"),
                evalue=evalue,
                tandem=tandem,
                ks_column=ks_column,
                ks_enabled=ks_enabled,
            ),
            encoding="utf-8",
            newline="\n",
        )
        configured += 1

    label = "同物种" if mode == "same" else "不同物种" if mode == "different" else "全部"
    print(f"{label} WGDI 配置完成（两目录布局）：{configured} 个。")
    return configured


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--mode", choices=("same", "different", "all"), required=True)
    args = parser.parse_args()
    try:
        generate(args.project_root, args.mode)
    except (FileNotFoundError, ValueError, OSError) as exc:
        parser.exit(1, f"错误：{exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
