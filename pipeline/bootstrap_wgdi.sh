#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
python_bin="${PYTHON:-python}"
args=("$project_root/code/tools/bootstrap_wgdi.py" --project-root "$project_root")
if [[ -n "${WGDI_BOOTSTRAP_ENV_DIR:-}" ]]; then
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
if [[ "${WGDI_BOOTSTRAP_ALLOW_UNPINNED_FALLBACK:-false}" == "true" ]]; then
    args+=(--allow-unpinned-fallback)
fi
exec "$python_bin" "${args[@]}"
