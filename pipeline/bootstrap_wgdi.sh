#!/usr/bin/env bash
# Create a local WGDI runtime for this workflow repository.
#
#   bash bootstrap_wgdi.sh [ENV_DIR]
#
# Creates .wgdi-env (or ENV_DIR) and installs the pinned WGDI Python
# dependencies.  bin/wgdi detects the environment automatically, so the
# repository can be copied to another machine without editing paths.
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
env_dir="${1:-$project_root/.wgdi-env}"
requirements="$project_root/requirements-wgdi-validated.txt"

python_bin="${PYTHON:-}"
if [[ -z "$python_bin" ]]; then
    for candidate in python3 python; do
        if command -v "$candidate" >/dev/null 2>&1; then
            python_bin="$(command -v "$candidate")"
            break
        fi
    done
fi
[[ -n "$python_bin" ]] || { echo "错误：找不到 Python 3.11+。" >&2; exit 1; }
echo "Python: $python_bin"

"$python_bin" -m venv "$env_dir"
venv_python="$env_dir/bin/python"
[[ -f "$venv_python" ]] || venv_python="$env_dir/Scripts/python.exe"
[[ -f "$venv_python" ]] || { echo "错误：虚拟环境创建失败：$env_dir" >&2; exit 1; }

"$venv_python" -m pip install --upgrade pip
installed=false
if [[ -f "$requirements" ]]; then
    echo "安装固定版本依赖: $requirements"
    if "$venv_python" -m pip install -r "$requirements"; then
        installed=true
    fi
fi
if [[ "$installed" != true ]]; then
    echo "固定版本安装失败，改用 WGDI 官方依赖（不锁定版本）。"
    "$venv_python" -m pip install wgdi matplotlib numpy scipy pandas biopython tabulate
fi
"$venv_python" -c "import wgdi, wgdi.run; print('wgdi import OK')"

echo
echo "Python 运行时已就绪: $venv_python"
echo "启动方式: bin/wgdi <wgdi 参数>（自动使用 $env_dir）"
echo
echo "外部程序（DIAMOND/BLAST/MUSCLE/yn00）需另行安装并在 PATH 中可用："
for tool in diamond blastp makeblastdb muscle yn00 perl; do
    if command -v "$tool" >/dev/null 2>&1; then
        echo "  [OK]   $tool -> $(command -v "$tool")"
    else
        echo "  [MISS] $tool"
    fi
done
