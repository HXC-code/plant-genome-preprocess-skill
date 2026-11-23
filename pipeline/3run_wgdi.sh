#!/usr/bin/env bash
# SOP step 4b: run all WGDI comparisons from flat configs in wgdi_input/config.
#
# Two-directory layout: wgdi_input holds every input (bundles, BLAST, configs);
# wgdi_output holds every product (results/, logs/, state/).  WGDI runs from the
# project root, so all paths inside a config are project-relative and portable.
# Resume is stage-level and fingerprint-keyed: a stage is skipped only when its
# input fingerprint still matches the marker.
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
cd "$PGCI_PROJECT_ROOT"
export MPLBACKEND=Agg

mode="${1:-all}"
case "$mode" in
    same|different|all) ;;
    *) echo "用法：$0 [same|different|all]" >&2; exit 2 ;;
esac

wgdi_ks="${WGDI_KS:-false}"
wgdi_tandem="${WGDI_TANDEM:-false}"
wgdi_ks_col="${WGDI_KS_COL:-ks_NG86}"
case "$wgdi_ks" in true|false) ;; *) echo "错误：WGDI_KS 必须是 true 或 false。" >&2; exit 2 ;; esac
case "$wgdi_tandem" in true|false) ;; *) echo "错误：WGDI_TANDEM 必须是 true 或 false。" >&2; exit 2 ;; esac
case "$wgdi_ks_col" in ks_NG86|ks_YN00) ;; *) echo "错误：WGDI_KS_COL 必须是 ks_NG86 或 ks_YN00。" >&2; exit 2 ;; esac

index_dir="metadata"
if [[ ! -d "$index_dir" ]]; then
    if [[ -d resources ]]; then
        index_dir="resources"
    else
        index_dir="."
    fi
fi
[[ -s "$index_dir/folder_list.txt" ]] || { echo "错误：folder_list.txt 不存在或为空。" >&2; exit 1; }
mapfile -t species_ids < "$index_dir/folder_list.txt"
[[ -s "$index_dir/species_map.tsv" ]] || { echo "错误：species_map.tsv 不存在或为空。" >&2; exit 1; }

config_root="wgdi_input/config"
result_root="wgdi_output/results"
log_root="wgdi_output/logs"
state_root="wgdi_output/state"
[[ -d "$config_root" ]] || { echo "错误：配置目录不存在：$config_root（先运行 2gen_wgdi_conf.sh）" >&2; exit 1; }
mkdir -p "$result_root" "$log_root" "$state_root"

# Bundle locations come from the homology scan map (code -> dir + file names).
declare -A bundle_dir_of=() cds_name_of=() pep_name_of=() gff_name_of=() lens_name_of=()
while IFS=$'\t' read -r code dir pep cds gff lens; do
    [[ -z "${code:-}" ]] && continue
    bundle_dir_of["$code"]="$dir"
    cds_name_of["$code"]="$cds"
    pep_name_of["$code"]="$pep"
    gff_name_of["$code"]="$gff"
    lens_name_of["$code"]="$lens"
done < "$index_dir/species_map.tsv"

wgdi_entry="${WGDI:-}"
if [[ -z "$wgdi_entry" ]]; then
    wgdi_entry="$(command -v wgdi || true)"
fi
if [[ -z "$wgdi_entry" ]]; then
    echo "错误：找不到 WGDI（用 WGDI 环境变量指定，或加入 PATH）。" >&2
    exit 2
fi
if [[ -f "$wgdi_entry" ]]; then
    WGDI_CMD=(bash "$wgdi_entry")
else
    WGDI_CMD=("$wgdi_entry")
fi
if [[ "$wgdi_tandem" == false && "$(basename -- "$wgdi_entry")" == "wgdi" ]]; then
    tandem_args=(--remove-tandem)
else
    tandem_args=()
fi

resolve_python() {
    if [[ -n "${PYTHON:-}" ]]; then
        printf '%s\n' "$PYTHON"
        return 0
    fi
    command -v python 2>/dev/null || command -v python3 2>/dev/null || true
}

fingerprint() {
    local path="$1"
    [[ -s "$path" ]] || { printf 'missing\n'; return 0; }
    if stat -c '%s-%y' "$path" 2>/dev/null; then
        return 0
    fi
    if stat -f '%z-%m' "$path" 2>/dev/null; then
        return 0
    fi
    cksum "$path" | awk '{print $1 "-" $2}'
}

content_digest() {
    local path="$1"
    [[ -f "$path" ]] || { printf 'missing\n'; return 0; }
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$path" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$path" | awk '{print $1}'
    else
        cksum "$path" | awk '{print $1 "-" $2}'
    fi
}

marker_path() {
    printf '%s/%s.done.%s\n' "$state_root" "$1" "$2"
}

check_outputs() {
    local prefix="$1"
    local required=("$prefix.png" "$prefix.icl")
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
        [[ -s "$result_root/$name" ]] || missing+=("$name")
    done
    if [[ ${#missing[@]} -gt 0 ]]; then
        printf '%s\n' "${missing[@]}" > "$state_root/${prefix}.missing_outputs.txt"
        return 1
    fi
    # Content-integrity check via the Python package; degrade to presence-only
    # when Python/the package is unavailable so the runner never hard-depends.
    local py integrity_out
    local ks_args=()
    if [[ "$wgdi_ks" == "true" ]]; then
        ks_args=(--ks)
    fi
    py="$(resolve_python)"
    if [[ -n "$py" ]]; then
        if integrity_out="$(PYTHONPATH="$PGCI_CODE_SRC${PYTHONPATH:+:$PYTHONPATH}" \
            "$py" -m pgci_preprocess check-wgdi-outputs \
            --folder "$result_root" --prefix "$prefix" "${ks_args[@]}" 2>/dev/null)"; then
            :
        elif [[ -n "$integrity_out" ]]; then
            printf '%s\n' "$integrity_out" > "$state_root/${prefix}.missing_outputs.txt"
            return 1
        fi
    fi
    rm -f "$state_root/${prefix}.missing_outputs.txt"
    return 0
}

run_stage() {
    local prefix="$1"
    local stage_name="$2"
    local output_file="$3"
    local key="$4"
    shift 4
    local marker
    marker="$(marker_path "$prefix" "$stage_name")"
    if [[ -f "$marker" && -s "$output_file" && "$(cat "$marker" 2>/dev/null || true)" == "$key" ]]; then
        echo "跳过已完成阶段：$stage_name"
        return 0
    fi
    echo "执行阶段：$stage_name"
    "$@"
    [[ -s "$output_file" ]] || {
        echo "错误：阶段 $stage_name 未生成 $output_file"
        return 1
    }
    printf '%s\n' "$key" > "$marker"
}

ran=0
skipped=0
failed=0

for left in "${species_ids[@]}"; do
    for right in "${species_ids[@]}"; do
        is_same=false
        [[ "$left" == "$right" ]] && is_same=true
        case "$mode" in
            same) [[ "$is_same" == false ]] && continue ;;
            different) [[ "$is_same" == true ]] && continue ;;
        esac

        prefix="$left"
        [[ "$left" != "$right" ]] && prefix="${left}_${right}"
        config="$config_root/$prefix.conf"
        if [[ ! -f "$config" ]]; then
            echo "错误：配置缺失：$config（先运行 2gen_wgdi_conf.sh）"
            failed=$((failed + 1))
            continue
        fi

        blast="wgdi_input/blast/${left}_pep/${left}_${right}.blast"
        log="$log_root/$prefix.log"

        if ! (
            blast_key="$(fingerprint "$blast")"
            config_key="$(content_digest "$config")"
            launcher_key="$(content_digest "$wgdi_entry")"
            gff_l="${bundle_dir_of[$left]}/${gff_name_of[$left]}"
            gff_r="${bundle_dir_of[$right]}/${gff_name_of[$right]}"
            lens_l="${bundle_dir_of[$left]}/${lens_name_of[$left]}"
            lens_r="${bundle_dir_of[$right]}/${lens_name_of[$right]}"
            cds_l="${bundle_dir_of[$left]}/${cds_name_of[$left]}"
            cds_r="${bundle_dir_of[$right]}/${cds_name_of[$right]}"
            pep_l="${bundle_dir_of[$left]}/${pep_name_of[$left]}"
            pep_r="${bundle_dir_of[$right]}/${pep_name_of[$right]}"
            core_key="core;config=$config_key;launcher=$launcher_key;blast=$blast_key"
            core_key="$core_key;gff1=$(fingerprint "$gff_l");gff2=$(fingerprint "$gff_r")"
            core_key="$core_key;lens1=$(fingerprint "$lens_l");lens2=$(fingerprint "$lens_r")"
            if [[ "$wgdi_ks" == true && "$left" != "$right" ]]; then
                state_dir="$state_root/$prefix"
                mkdir -p "$state_dir"
                if [[ ! -s "$state_dir/all.cds" || "$cds_l" -nt "$state_dir/all.cds" || "$cds_r" -nt "$state_dir/all.cds" ]]; then
                    cat "$cds_l" "$cds_r" > "$state_dir/all.cds"
                fi
                if [[ ! -s "$state_dir/all.pep" || "$pep_l" -nt "$state_dir/all.pep" || "$pep_r" -nt "$state_dir/all.pep" ]]; then
                    cat "$pep_l" "$pep_r" > "$state_dir/all.pep"
                fi
            fi
            run_stage "$prefix" dotplot "$result_root/$prefix.png" "$core_key" \
                "${WGDI_CMD[@]}" -d "$config"
            run_stage "$prefix" collinearity "$result_root/$prefix.icl" "$core_key" \
                "${WGDI_CMD[@]}" -icl "$config"
            if [[ "$wgdi_ks" == true ]]; then
                if [[ "$left" == "$right" ]]; then
                    cds_key="$(fingerprint "$cds_l")"
                    pep_key="$(fingerprint "$pep_l")"
                else
                    cds_key="$(fingerprint "$state_root/$prefix/all.cds")"
                    pep_key="$(fingerprint "$state_root/$prefix/all.pep")"
                fi
                ks_key="ks;$core_key;icl=$(fingerprint "$result_root/$prefix.icl")"
                ks_key="$ks_key;cds=$cds_key;pep=$pep_key;col=$wgdi_ks_col;tandem=$wgdi_tandem"
                run_stage "$prefix" ks "$result_root/$prefix.ks" "$ks_key" \
                    "${WGDI_CMD[@]}" -ks "$config"
                run_stage "$prefix" blockinfo "$result_root/$prefix.bi.csv" "$ks_key" \
                    "${WGDI_CMD[@]}" -bi "$config"
                run_stage "$prefix" correspondence "$result_root/$prefix.bi.c.csv" "$ks_key" \
                    "${WGDI_CMD[@]}" ${tandem_args[@]+"${tandem_args[@]}"} -c "$config"
                run_stage "$prefix" blockks "$result_root/$prefix.bk.png" "$ks_key" \
                    "${WGDI_CMD[@]}" ${tandem_args[@]+"${tandem_args[@]}"} -bk "$config"
                run_stage "$prefix" kspeaks "$result_root/$prefix.kp.png" "$ks_key" \
                    "${WGDI_CMD[@]}" ${tandem_args[@]+"${tandem_args[@]}"} -kp "$config"
                run_stage "$prefix" peaksfit "$result_root/$prefix.pf.png" "$ks_key" \
                    "${WGDI_CMD[@]}" -pf "$config"
            fi
        ) > "$log" 2>&1; then
            echo "错误：WGDI 运行失败：$left / $right"
            tail -n 10 "$log" || true
            failed=$((failed + 1))
            continue
        fi

        if ! check_outputs "$prefix"; then
            echo "错误：WGDI 返回成功，但结果文件不完整：$left / $right"
            cat "$state_root/${prefix}.missing_outputs.txt" 2>/dev/null || true
            failed=$((failed + 1))
            continue
        fi

        if grep -q '^执行阶段' "$log"; then
            ran=$((ran + 1))
            echo "WGDI 完成：$left / $right"
        else
            skipped=$((skipped + 1))
            echo "WGDI 已完成，跳过：$left / $right"
        fi
    done
done

echo "WGDI 汇总：运行 $ran 个，已跳过 $skipped 个，失败 $failed 个。"
[[ $failed -eq 0 ]]
