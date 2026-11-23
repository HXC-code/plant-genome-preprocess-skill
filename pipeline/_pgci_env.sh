#!/usr/bin/env bash
# Shared environment resolver for the pipeline shell scripts.
#
# New-machine contract:
#   - Skill-internal resources (code/, bin/, tools/) are located relative to
#     this script, so the skill can be placed anywhere.
#   - The project root (DATA/, metadata/, wgdi_input/, wgdi_output/)
#     is resolved from the pointer written by `discover`, or from the explicit
#     PGCI_PROJECT_ROOT override.  No absolute path is ever hard-coded.
#
# Source this file from a pipeline script, then cd "$PGCI_PROJECT_ROOT".

_pgci_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

PGCI_SKILL_DIR="$_pgci_here"
PGCI_CODE_SRC="${_pgci_here}/code/src"

pgci_is_wsl() {
    [[ -r /proc/version ]] && grep -qi microsoft /proc/version
}

pgci_to_windows_path() {
    local value="$1"
    if command -v cygpath >/dev/null 2>&1 && [[ "$value" != [A-Za-z]:[\\/]* ]]; then
        cygpath -w "$value"
    elif pgci_is_wsl && command -v wslpath >/dev/null 2>&1; then
        wslpath -w "$value"
    else
        printf '%s\n' "$value"
    fi
}

pgci_to_shell_path() {
    local value="$1"
    if [[ "$value" != [A-Za-z]:[\\/]* ]]; then
        printf '%s\n' "$value"
    elif command -v cygpath >/dev/null 2>&1; then
        cygpath -u "$value"
    elif pgci_is_wsl && command -v wslpath >/dev/null 2>&1; then
        wslpath -u "$value"
    else
        printf '%s\n' "$value"
    fi
}

# The pinned WGDI launcher travels with the skill, not with the project data.
# Expose it via WGDI so generate_wgdi_configs.py resolves it regardless of
# where the project root lives.
if [[ -z "${WGDI:-}" && -f "${_pgci_here}/bin/wgdi" ]]; then
    WGDI="${_pgci_here}/bin/wgdi"
    export WGDI
fi

_pgci_root=""
if [[ -n "${PGCI_PROJECT_ROOT:-}" ]]; then
    _pgci_root="$PGCI_PROJECT_ROOT"
else
    _pgci_pointer="${_pgci_here}/code/.pgci-machine.json"
    if [[ -f "$_pgci_pointer" ]]; then
        # New machines may expose only python3; probe rather than assume.
        _pgci_python="$(command -v python || command -v python3 || true)"
        if [[ -n "$_pgci_python" ]]; then
            _pgci_root="$("$_pgci_python" -c 'import json,sys;print(json.load(open(sys.argv[1],encoding="utf-8"))["project_root"])' \
                "$_pgci_pointer" 2>/dev/null || true)"
        fi
    fi
fi
# Classic layout fallback: skill directory itself is the project.
if [[ -z "$_pgci_root" ]]; then
    _pgci_root="$_pgci_here"
fi
# Normalize a Windows native path (D:\...) to POSIX for bash, when needed.
_pgci_root="$(pgci_to_shell_path "$_pgci_root")"
PGCI_PROJECT_ROOT="$(cd "$_pgci_root" 2>/dev/null && pwd -P || printf '%s\n' "$_pgci_root")"

export PGCI_SKILL_DIR PGCI_CODE_SRC PGCI_PROJECT_ROOT
