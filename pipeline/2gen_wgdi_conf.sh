#!/usr/bin/env bash
# SOP step 4a: create lean comparison folders + WGDI configs (scan-driven).
#
# Replaces the old copy scripts (2run_pre_wgdi.sh / 3blast_to_wgdi.sh) and the
# two mode wrappers (4run / 5run): one entry, one shared implementation.
# One flat <pair>.conf per comparison in wgdi_input/config; bundle and BLAST
# inputs are referenced by project-relative paths from species_map.tsv (the
# homology scan map), never copied. Results go to wgdi_output/results.
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=_pgci_env.sh
if [[ -f "$script_dir/_pgci_env.sh" ]]; then
    source "$script_dir/_pgci_env.sh"
else
    PGCI_SKILL_DIR="$script_dir"
    PGCI_PROJECT_ROOT="$script_dir"
    PGCI_CODE_SRC="$script_dir/code/src"
    export PGCI_SKILL_DIR PGCI_PROJECT_ROOT PGCI_CODE_SRC
fi

mode="${1:-all}"
case "$mode" in
    same|different|all) ;;
    *) echo "用法：$0 [same|different|all]" >&2; exit 2 ;;
esac

python_bin="${PYTHON:-$(command -v python || command -v python3 || true)}"
[[ -n "$python_bin" ]] || { echo "错误：找不到 Python。" >&2; exit 2; }

exec "$python_bin" "$PGCI_SKILL_DIR/code/tools/generate_wgdi_configs.py" \
    --project-root "$PGCI_PROJECT_ROOT" --mode "$mode"
