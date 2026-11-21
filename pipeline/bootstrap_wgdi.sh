#!/usr/bin/env bash
set -euo pipefail
export PYTHONUTF8=1
export PYTHONIOENCODING=utf-8
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
python_bin="${PYTHON:-python}"
args=("$project_root/code/tools/bootstrap_wgdi.py" --project-root "$project_root")
if [[ $# -gt 1 ]]; then
    echo "用法：bash bootstrap_wgdi.sh [ENV_DIR]" >&2
    exit 2
fi
if [[ $# -eq 1 ]]; then
    args+=(--env-dir "$1")
elif [[ -n "${WGDI_BOOTSTRAP_ENV_DIR:-}" ]]; then
    args+=(--env-dir "$WGDI_BOOTSTRAP_ENV_DIR")
fi
if [[ -n "${WGDI_BOOTSTRAP_PYTHON:-}" ]]; then
    args+=(--python-exe "$WGDI_BOOTSTRAP_PYTHON")
fi
if [[ "${WGDI_BOOTSTRAP_FORCE:-false}" == "true" ]]; then
    args+=(--force)
fi
if [[ "${WGDI_BOOTSTRAP_SKIP_EXTERNAL_CHECK:-false}" == "true" ]]; then
    args+=(--skip-external-check)
fi
if [[ "${WGDI_BOOTSTRAP_SKIP_PAL2NAL:-false}" == "true" ]]; then
    args+=(--skip-pal2nal)
fi
if [[ "${WGDI_BOOTSTRAP_ALLOW_UNPINNED_FALLBACK:-false}" == "true" ]]; then
    args+=(--allow-unpinned-fallback)
fi
exec "$python_bin" "${args[@]}"
