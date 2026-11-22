#!/usr/bin/env bash
# Run an isolated fixture through DIAMOND and the existing 2-6 WGDI scripts.
set -euo pipefail
export PYTHONUTF8=1
export PYTHONIOENCODING=utf-8

if [[ $# -ne 1 ]]; then
    echo "用法：bash tools/run_wgdi_smoke.sh /absolute/path/to/smoke_wgdi" >&2
    exit 2
fi
root="$(cd "$1" && pwd -P)"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
workflow_root="$(cd "$script_dir/../.." && pwd -P)"
python_bin="${PYTHON:-python}"
cd "$root"
[[ -s fixture_manifest.json ]] || { echo "错误：缺少 fixture_manifest.json"; exit 1; }
shopt -s nullglob
pep_inputs=(*/*.pep)
[[ ${#pep_inputs[@]} -ge 1 ]] || { echo "错误：缺少 smoke 蛋白文件"; exit 1; }

mkdir -p logs
export PATH="$root/bin:${PATH}"
export MPLBACKEND=Agg
export HOMOLOGY_ENGINE=diamond
if [[ -z "${WGDI_PYTHON:-}" ]]; then
    for candidate in \
        "$workflow_root/.wgdi-env/Scripts/python.exe" \
        "$workflow_root/.wgdi-env/bin/python"; do
        if [[ -f "$candidate" ]]; then
            export WGDI_PYTHON="$candidate"
            break
        fi
    done
fi
if [[ -z "${WGDI_KS_COL:-}" ]]; then
    selection_method="$("$python_bin" -c \
        'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["selection"]["method"])' \
        "$root/fixture_manifest.json")"
    if [[ "$selection_method" == "deterministic synthetic orthologues" ]]; then
        export WGDI_KS_COL=ks_YN00
    else
        export WGDI_KS_COL=ks_NG86
    fi
fi
export DIAMOND_EVALUE="${DIAMOND_EVALUE:-1e-5}"
export DIAMOND_SENSITIVITY="${DIAMOND_SENSITIVITY:-very-sensitive}"
export DIAMOND_THREADS="${DIAMOND_THREADS:-4}"
export WGDI_EVALUE="${WGDI_EVALUE:-$DIAMOND_EVALUE}"
# The release gate exercises the full chain including the optional Ks branch;
# set WGDI_KS=false to smoke only the SOP core (dotplot + collinearity).
export WGDI_KS="${WGDI_KS:-true}"

run_stage() {
    local script="$1"
    local name="$2"
    local log="logs/${name}.log"
    echo "[SMOKE] $script"
    if bash "$script" >"$log" 2>&1; then
        echo "[PASS] $script (log=$log)"
    else
        local code=$?
        echo "[FAIL] $script (exit=$code, log=$log)" >&2
        tail -n 80 "$log" >&2 || true
        exit "$code"
    fi
}

run_stage run_homology.sh 01_diamond
run_stage 2run_pre_wgdi.sh 02_pre_wgdi
run_stage 3blast_to_wgdi.sh 03_blast_to_wgdi
run_stage 4run_wgdi_different.sh 04_config_different
run_stage 5run_wgdi_same.sh 05_config_same
run_stage 6run_wgdi.sh 06_wgdi

printf 'status\tPASS\ndiamond_evalue\t%s\ndiamond_sensitivity\t%s\nks_enabled\t%s\nscope\tDIAMOND plus scripts 2-6; step 0 covered separately by preprocess validation\n' \
    "$DIAMOND_EVALUE" "$DIAMOND_SENSITIVITY" "$WGDI_KS" > smoke_status.tsv
"$python_bin" "$script_dir/collect_wgdi_smoke_report.py" "$root"
# Capture the public receipt in the same environment as the run: the collector
# probes PATH for tool versions, so collecting later from a different shell
# would record a different (wrong) environment.
"$python_bin" "$script_dir/collect_wgdi_smoke_report.py" "$root" --public --output "$root/smoke_report.public.json"
echo "[PASS] smoke completed; status=smoke_status.tsv; report=smoke_report.json; public=smoke_report.public.json"
