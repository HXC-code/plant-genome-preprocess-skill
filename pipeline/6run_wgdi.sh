#!/usr/bin/env bash

set -euo pipefail
cd "$(dirname "$0")"
project_root="$(pwd -P)"
export MPLBACKEND=Agg

wgdi_ks="${WGDI_KS:-false}"
case "$wgdi_ks" in
    true|false) ;;
    *) echo "错误：WGDI_KS 必须是 true 或 false，实际为：$wgdi_ks" >&2; exit 2 ;;
esac

[[ -s folder_list.txt ]] || { echo "错误：folder_list.txt 不存在或为空。"; exit 1; }
mapfile -t species_ids < folder_list.txt

resolve_python() {
    if [[ -n "${PYTHON:-}" ]]; then
        printf '%s\n' "$PYTHON"
        return 0
    fi
    command -v python 2>/dev/null || command -v python3 2>/dev/null || true
}

check_outputs() {
    local folder="$1"
    local base_id="$2"
    local target_id="$3"
    local prefix
    local required=()

    if [[ "$base_id" == "$target_id" ]]; then
        prefix="$base_id"
    else
        prefix="${base_id}_${target_id}"
    fi
    required=(
        "$prefix.png"
        "$prefix.icl"
    )
    if [[ "$wgdi_ks" == "true" ]]; then
        required+=(
            "$prefix.ks"
            "$prefix.bi.csv"
            "$prefix.bi.c.csv"
            "$prefix.bk.png"
            "$prefix.kp.png"
            "$prefix.pf.png"
        )
    fi

    local missing=()
    local name
    for name in "${required[@]}"; do
        [[ -s "$folder/$name" ]] || missing+=("$name")
    done
    if [[ ${#missing[@]} -gt 0 ]]; then
        printf '%s\n' "${missing[@]}" > "$folder/missing_outputs.txt"
        return 1
    fi
    # Content-integrity check via the installed package; a truncated PNG or a
    # garbage table is not a finished result. Degrade to presence-only when
    # Python/the package is unavailable so the runner never hard-depends on it.
    local py integrity_out
    local ks_args=()
    if [[ "$wgdi_ks" == "true" ]]; then
        ks_args=(--ks)
    fi
    py="$(resolve_python)"
    if [[ -n "$py" ]]; then
        if integrity_out="$(PYTHONPATH="$project_root/code/src${PYTHONPATH:+:$PYTHONPATH}" \
            "$py" -m pgci_preprocess check-wgdi-outputs \
            --folder "$folder" --prefix "$prefix" "${ks_args[@]}" 2>/dev/null)"; then
            :
        elif [[ -n "$integrity_out" ]]; then
            printf '%s\n' "$integrity_out" > "$folder/missing_outputs.txt"
            return 1
        fi
    fi
    rm -f "$folder/missing_outputs.txt"
    return 0
}

ran=0
skipped=0
failed=0

for base_id in "${species_ids[@]}"; do
    for target_id in "${species_ids[@]}"; do
        folder="Auto_wgdi/${base_id}_base/${base_id}_${target_id}"

        if [[ ! -f "$folder/run.sh" || ! -f "$folder/total.conf" ]]; then
            echo "错误：配置缺失：$folder"
            failed=$((failed + 1))
            continue
        fi

        if [[ -f "$folder/.wgdi.done" ]] && check_outputs "$folder" "$base_id" "$target_id"; then
            echo "WGDI 已完成，跳过：$base_id / $target_id"
            skipped=$((skipped + 1))
            continue
        fi

        rm -f "$folder/.wgdi.done"
        echo "运行 WGDI：$base_id / $target_id"

        if ! (
            cd "$folder"
            bash run.sh
        ) > "$folder/run.log" 2>&1; then
            echo "错误：WGDI 运行失败：$base_id / $target_id"
            tail -n 10 "$folder/run.log" || true
            failed=$((failed + 1))
            continue
        fi

        if ! check_outputs "$folder" "$base_id" "$target_id"; then
            echo "错误：WGDI 返回成功，但结果文件不完整：$base_id / $target_id"
            cat "$folder/missing_outputs.txt"
            failed=$((failed + 1))
            continue
        fi

        touch "$folder/.wgdi.done"
        ran=$((ran + 1))
        echo "WGDI 完成：$base_id / $target_id"
    done
done

echo "WGDI 汇总：运行 $ran 个，已跳过 $skipped 个，失败 $failed 个。"
[[ $failed -eq 0 ]]
