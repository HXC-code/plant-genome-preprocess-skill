#!/usr/bin/env bash
set -euo pipefail
export PYTHONUTF8=1
export PYTHONIOENCODING=utf-8
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
python_bin="${PYTHON:-python}"
exec "$python_bin" "$script_dir/code/tools/generate_wgdi_configs.py" \
    --project-root "$script_dir" --mode different "$@"
