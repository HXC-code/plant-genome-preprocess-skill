#!/usr/bin/env python3
"""Generate WGDI comparison directories and configuration from one implementation."""

from __future__ import annotations

import argparse
import os
import re
import shlex
import shutil
import stat
from pathlib import Path


def read_species(project_root: Path) -> list[str]:
    folder_list = project_root / "folder_list.txt"
    if not folder_list.is_file() or not folder_list.read_text(encoding="utf-8").strip():
        raise ValueError("folder_list.txt 不存在或为空")
    species: list[str] = []
    seen: set[str] = set()
    for line in folder_list.read_text(encoding="utf-8").splitlines():
        code = line.strip()
        if code and code not in seen:
            seen.add(code)
            species.append(code)
    return species


def resolve_wgdi_entry(project_root: Path) -> Path | str:
    requested = os.environ.get("WGDI", "")
    if requested:
        requested_path = Path(requested)
        if requested_path.is_file():
            return requested_path.resolve()
        if shutil.which(requested):
            return requested
        raise FileNotFoundError(f"找不到 WGDI：{requested}（可用 WGDI 环境变量指定）")
    elif (project_root / "bin" / "wgdi").is_file():
        return (project_root / "bin" / "wgdi").resolve()
    else:
        if shutil.which("wgdi"):
            return "wgdi"
        raise FileNotFoundError("找不到 WGDI：wgdi（可用 WGDI 环境变量指定）")


def _bash_path(value: Path | str) -> str:
    """Render an absolute Windows path for Bash without backslash escapes."""

    text = value.as_posix() if isinstance(value, Path) else str(value)
    match = re.fullmatch(r"([A-Za-z]):/(.+)", text.replace("\\", "/"))
    if match:
        return f"/{match.group(1).lower()}/{match.group(2)}"
    return text


def _run_script(
    left: str,
    right: str,
    wgdi_entry: Path | str,
    tandem: bool,
    ks_enabled: bool,
) -> str:
    common_products = left != right
    prep = ""
    if common_products:
        prep = f'''
cat "{left}.cds" "{right}.cds" > all.cds
cat "{left}.pep" "{right}.pep" > all.pep
'''

    remove_tandem = (
        not tandem
        and isinstance(wgdi_entry, Path)
        and wgdi_entry.name == "wgdi"
    )
    tandem_setup = (
        "WGDI_TANDEM_ARGS=(--remove-tandem)"
        if remove_tandem
        else "WGDI_TANDEM_ARGS=()"
    )
    entry = shlex.quote(_bash_path(wgdi_entry))
    prefix = left if left == right else f"{left}_{right}"
    # SOP step 4 core: dotplot + collinearity blocks. The six Ks-dependent
    # stages (ks/blockinfo/correspondence/blockks/kspeaks/peaksfit) run only
    # when the caller opts in with WGDI_KS=true: upstream blockinfo and
    # correspondence consume Ks values, and with ks=none correspondence
    # filters out every block, so they must not run in the default mode.
    stages = [
        f'run_stage dotplot "{prefix}.png" "${{WGDI_CMD[@]}}" -d total.conf',
        f'run_stage collinearity "{prefix}.icl" "${{WGDI_CMD[@]}}" -icl total.conf',
    ]
    if ks_enabled:
        stages += [
            f'run_stage ks "{prefix}.ks" "${{WGDI_CMD[@]}}" -ks total.conf',
            f'run_stage blockinfo "{prefix}.bi.csv" "${{WGDI_CMD[@]}}" -bi total.conf',
            f'run_stage correspondence "{prefix}.bi.c.csv" "${{WGDI_CMD[@]}}" "${{WGDI_TANDEM_ARGS[@]}}" -c total.conf',
            f'run_stage blockks "{prefix}.bk.png" "${{WGDI_CMD[@]}}" "${{WGDI_TANDEM_ARGS[@]}}" -bk total.conf',
            f'run_stage kspeaks "{prefix}.kp.png" "${{WGDI_CMD[@]}}" "${{WGDI_TANDEM_ARGS[@]}}" -kp total.conf',
            f'run_stage peaksfit "{prefix}.pf.png" "${{WGDI_CMD[@]}}" -pf total.conf',
        ]
    stage_lines = "\n".join(stages)
    return f"""#!/usr/bin/env bash
set -euo pipefail
export MPLBACKEND=Agg

mkdir -p .wgdi_steps
{prep}
WGDI_ENTRY={entry}
if [[ -f "$WGDI_ENTRY" ]]; then
    WGDI_CMD=(bash "$WGDI_ENTRY")
else
    WGDI_CMD=("$WGDI_ENTRY")
fi
{tandem_setup}

run_stage() {{
    local stage_name="$1"
    local output_file="$2"
    shift 2
    local marker=".wgdi_steps/${{stage_name}}.done"
    if [[ -f "$marker" && -f "$output_file" ]]; then
        echo "跳过已完成阶段：$stage_name"
        return 0
    fi
    "$@"
    [[ -f "$output_file" ]] || {{
        echo "错误：阶段 $stage_name 未生成 $output_file"
        return 1
    }}
    touch "$marker"
}}

{stage_lines}
"""


def _config_text(
    left: str,
    right: str,
    *,
    evalue: str,
    tandem: bool,
    ks_column: str,
    ks_enabled: bool,
) -> str:
    same = left == right
    output_prefix = left if same else f"{left}_{right}"
    blast_file = f"{left}_{right}.blast"
    cds_file = f"{left}.cds" if same else "all.cds"
    pep_file = f"{left}.pep" if same else "all.pep"
    tandem_value = "true" if tandem else "false"
    core = f"""[dotplot]
blast = {blast_file}
gff1 = {left}.new.gff
gff2 = {right}.new.gff
lens1 = {left}.lens
lens2 = {right}.lens
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
savefig = {output_prefix}.png

[collinearity]
blast = {blast_file}
gff1 = {left}.new.gff
gff2 = {right}.new.gff
lens1 = {left}.lens
lens2 = {right}.lens
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
savefile = {output_prefix}.icl
"""
    ks_part = f"""
[ks]
cds_file = {cds_file}
pep_file = {pep_file}
align_software = muscle
pairs_file = {output_prefix}.icl
ks_file = {output_prefix}.ks

[blockinfo]
blast = {blast_file}
gff1 = {left}.new.gff
gff2 = {right}.new.gff
lens1 = {left}.lens
lens2 = {right}.lens
collinearity = {output_prefix}.icl
score = 100
evalue = {evalue}
repeat_number = 20
position = order
ks = {output_prefix}.ks
ks_col = {ks_column}
savefile = {output_prefix}.bi.csv

[correspondence]
blockinfo = {output_prefix}.bi.csv
lens1 = {left}.lens
lens2 = {right}.lens
tandem = {tandem_value}
tandem_length = 200
pvalue = 0.2
block_length = 5
tandem_ratio = 0.5
multiple = 1
homo = -1,1
savefile = {output_prefix}.bi.c.csv

[blockks]
lens1 = {left}.lens
lens2 = {right}.lens
genome1_name = {left}
genome2_name = {right}
blockinfo = {output_prefix}.bi.csv
pvalue = 0.2
tandem = {tandem_value}
tandem_length = 200
markersize = 1
area = 0,2
block_length = 5
figsize = 8,8
savefig = {output_prefix}.bk.png

[kspeaks]
blockinfo = {output_prefix}.bi.csv
pvalue = 0.2
tandem = {tandem_value}
block_length = 5
ks_area = 0,10
multiple = 1
homo = 0,1
fontsize = 9
area = 0,3
figsize = 10,6.18
savefig = {output_prefix}.kp.png
savefile = {output_prefix}.kp

[peaksfit]
blockinfo = {output_prefix}.bi.csv
mode = median
bins_number = 200
ks_area = 0,10
fontsize = 9
area = 0,3
figsize = 10,6.18
shadow = true
savefig = {output_prefix}.pf.png

[ksfigure]
ksfit = ks_fit_result.csv
labelfontsize = 15
legendfontsize = 15
xlabel = none
ylabel = none
title = none
area = 0,2
figsize = 10,6.18
shadow = false
savefig = {output_prefix}.kf.png
"""
    return core + ks_part if ks_enabled else core


def generate(project_root: Path, mode: str) -> int:
    root = project_root.resolve()
    species = read_species(root)
    wgdi_entry = resolve_wgdi_entry(root)
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
        folder = root / "Auto_wgdi" / f"{left}_base" / f"{left}_{right}"
        if not folder.is_dir():
            raise FileNotFoundError(f"缺少目录：{folder}")
        (folder / "total.conf").write_text(
            _config_text(
                left,
                right,
                evalue=evalue,
                tandem=tandem,
                ks_column=ks_column,
                ks_enabled=ks_enabled,
            ),
            encoding="utf-8",
            newline="\n",
        )
        run_path = folder / "run.sh"
        run_path.write_text(
            _run_script(left, right, wgdi_entry, tandem, ks_enabled),
            encoding="utf-8",
            newline="\n",
        )
        run_path.chmod(run_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        steps = folder / ".wgdi_steps"
        if steps.is_dir():
            for marker in steps.glob("*.done"):
                marker.unlink()
        done_marker = folder / ".wgdi.done"
        if done_marker.exists():
            done_marker.unlink()
        configured += 1

    label = "同物种" if mode == "same" else "不同物种"
    print(f"{label} WGDI 配置完成：{configured} 个。")
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
