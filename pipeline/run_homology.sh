#!/usr/bin/env bash
# Shared DIAMOND/BLAST implementation behind the 1run_*.sh compatibility wrappers.
set -euo pipefail
shopt -s nullglob

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=_pgci_env.sh
if [[ -f "$script_dir/_pgci_env.sh" ]]; then
    source "$script_dir/_pgci_env.sh"
else
    # Copied standalone (tests, fixtures): the script sits at the project root.
    PGCI_SKILL_DIR="$script_dir"
    PGCI_PROJECT_ROOT="$script_dir"
    PGCI_CODE_SRC="$script_dir/code/src"
    export PGCI_SKILL_DIR PGCI_PROJECT_ROOT PGCI_CODE_SRC

    pgci_to_windows_path() {
        local value="$1"
        if command -v cygpath >/dev/null 2>&1 && [[ "$value" != [A-Za-z]:[\\/]* ]]; then
            cygpath -w "$value"
        elif command -v wslpath >/dev/null 2>&1 && grep -qi microsoft /proc/version 2>/dev/null; then
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
        elif command -v wslpath >/dev/null 2>&1 && grep -qi microsoft /proc/version 2>/dev/null; then
            wslpath -u "$value"
        else
            printf '%s\n' "$value"
        fi
    }
fi
cd "$PGCI_PROJECT_ROOT"
project_root="$(pwd -P)"
# Small scan indexes are named by content; old projects keep root-level files.
index_dir="metadata"
if [[ ! -d "$index_dir" ]]; then
    if [[ -d resources ]]; then
        index_dir="resources"
    else
        index_dir="."
    fi
fi
engine="${HOMOLOGY_ENGINE:-diamond}"
allow_not_ready="${HOMOLOGY_ALLOW_NOT_READY:-false}"

case "$engine" in
    diamond|blast) ;;
    *)
        echo "错误：HOMOLOGY_ENGINE 必须是 diamond 或 blast，实际为：$engine"
        exit 2
        ;;
esac
case "${allow_not_ready,,}" in
    true|1|yes) allow_not_ready=true ;;
    false|0|no|"") allow_not_ready=false ;;
    *)
        echo "错误：HOMOLOGY_ALLOW_NOT_READY 必须是 true 或 false，实际为：$allow_not_ready"
        exit 2
        ;;
esac
require_adjudication="${HOMOLOGY_REQUIRE_ADJUDICATION:-false}"
case "${require_adjudication,,}" in
    true|1|yes) require_adjudication=true ;;
    false|0|no|"") require_adjudication=false ;;
    *)
        echo "错误：HOMOLOGY_REQUIRE_ADJUDICATION 必须是 true 或 false，实际为：$require_adjudication"
        exit 2
        ;;
esac

if [[ -n "${HOMOLOGY_EVALUE:-}" ]]; then
    homology_evalue="$HOMOLOGY_EVALUE"
elif [[ "$engine" == "diamond" ]]; then
    homology_evalue="${DIAMOND_EVALUE:-1e-5}"
else
    homology_evalue="${BLAST_EVALUE:-1e-5}"
fi
[[ "$homology_evalue" =~ ^[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$ ]] || {
    echo "错误：HOMOLOGY_EVALUE 不是有效数值：$homology_evalue"
    exit 2
}

if [[ "$engine" == "diamond" ]]; then
    engine_label="DIAMOND"
    tool_path="${DIAMOND:-$(command -v diamond || true)}"
    threads="${DIAMOND_THREADS:-8}"
    blocks="${DIAMOND_BLOCKS:-4}"
    max_targets="${DIAMOND_MAX_TARGETS:-20}"
    sensitivity="${DIAMOND_SENSITIVITY:-very-sensitive}"
    case "$sensitivity" in
        sensitive|more-sensitive|very-sensitive|ultra-sensitive) ;;
        *)
            echo "错误：不支持的 DIAMOND_SENSITIVITY：$sensitivity"
            exit 2
            ;;
    esac
else
    engine_label="BLAST"
    blastp_path="${BLASTP:-$(command -v blastp || true)}"
    makeblastdb_path="${MAKEBLASTDB:-$(command -v makeblastdb || true)}"
    threads="${BLAST_THREADS:-16}"
    max_targets="${BLAST_MAX_TARGETS:-20}"
fi

[[ "$threads" =~ ^[1-9][0-9]*$ ]] || {
    echo "错误：线程数必须是正整数：$threads"
    exit 2
}
[[ "$max_targets" =~ ^[1-9][0-9]*$ ]] || {
    echo "错误：最大目标数必须是正整数：$max_targets"
    exit 2
}

resolve_tool() {
    local path="$1"
    local label="$2"
    if [[ -z "$path" ]] || { [[ ! -x "$path" ]] && ! command -v "$path" >/dev/null 2>&1; }; then
        echo "错误：找不到 $label：${path:-<empty>}"
        exit 1
    fi
}

resolve_python() {
    if [[ -n "${PYTHON:-}" ]]; then
        printf '%s\n' "$PYTHON"
        return 0
    fi
    command -v python 2>/dev/null || command -v python3 2>/dev/null || true
}

# Git Bash pwd -P yields POSIX paths (/d/...) that a Windows Python cannot
# resolve. WSL has the same issue when WGDI uses a Windows interpreter.
python_is_windows() {
    case "${1,,}" in
        *.exe|/mnt/*.exe) return 0 ;;
        *) return 1 ;;
    esac
}

to_python_path() {
    local py="$1"
    local value="$2"
    if python_is_windows "$py"; then
        pgci_to_windows_path "$value"
    else
        printf '%s\n' "$value"
    fi
}

python_path_sep() {
    if python_is_windows "$1"; then printf ';'; else printf ':'; fi
}

species_readiness_via_python() {
    local species_id="$1" py out proot psrc sep
    py="$(resolve_python)"
    [[ -n "$py" ]] || return 1
    proot="$(to_python_path "$py" "$project_root")"
    psrc="$(to_python_path "$py" "$PGCI_CODE_SRC")"
    sep="$(python_path_sep "$py")"
    out="$(PYTHONPATH="$psrc${PYTHONPATH:+$sep$PYTHONPATH}" \
        "$py" -m pgci_preprocess readiness \
        --project-root "$proot" --species "$species_id" 2>/dev/null)" || return 1
    printf '%s\n' "${out#*$'\t'}"
}

species_readiness_grep() {
    local species_id="$1"
    local ready=""
    local qc_path
    local index_path
    for qc_path in \
        "$project_root/wgdi_input/$species_id/$species_id.qc.json" \
        "$project_root/analysis_v2/$species_id/$species_id.qc.json" \
        "$project_root/plant-genome-resources/analysis_v2/$species_id/$species_id.qc.json" \
        "$project_root/$species_id/$species_id.qc.json"; do
        if [[ -f "$qc_path" ]]; then
            if grep -Eq '"wgdi_ready"[[:space:]]*:[[:space:]]*"yes"' "$qc_path"; then
                ready="yes"
            elif grep -Eq '"wgdi_ready"[[:space:]]*:[[:space:]]*"marginal"' "$qc_path"; then
                ready="marginal"
            elif grep -Eq '"wgdi_ready"[[:space:]]*:[[:space:]]*"no"' "$qc_path"; then
                ready="no"
            fi
            [[ -n "$ready" ]] && break
        fi
    done
    if [[ -z "$ready" ]]; then
        for index_path in \
            "$project_root/resources/wgdi_v2_index.tsv" \
            "$project_root/plant-genome-resources/resources/wgdi_v2_index.tsv"; do
            if [[ -f "$index_path" ]]; then
                ready="$(awk -F '\t' -v code="$species_id" '
                    NR == 1 {
                        for (i = 1; i <= NF; i++) {
                            if ($i == "code") code_col = i
                            if ($i == "wgdi_ready") ready_col = i
                        }
                        next
                    }
                    code_col && ready_col && $code_col == code {
                        print $ready_col
                        exit
                    }
                ' "$index_path")"
                [[ -n "$ready" ]] && break
            fi
        done
    fi
    printf '%s\n' "$ready"
}

species_readiness() {
    local species_id="$1" ready
    if ready="$(species_readiness_via_python "$species_id")"; then
        printf '%s\n' "$ready"
        return 0
    fi
    species_readiness_grep "$species_id"
}

species_adjudication_status() {
    local species_id="$1" py out proot psrc sep
    py="$(resolve_python)"
    [[ -n "$py" ]] || { printf 'unknown\n'; return 0; }
    proot="$(to_python_path "$py" "$project_root")"
    psrc="$(to_python_path "$py" "$PGCI_CODE_SRC")"
    sep="$(python_path_sep "$py")"
    out="$(PYTHONPATH="$psrc${PYTHONPATH:+$sep$PYTHONPATH}" \
        "$py" -m pgci_preprocess adjudication-status \
        --project-root "$proot" --species "$species_id" 2>/dev/null)" \
        || { printf 'unknown\n'; return 0; }
    printf '%s\n' "${out#*$'\t'}"
}

tool_path="$(pgci_to_shell_path "$tool_path")"
if [[ "$engine" == "blast" ]]; then
    blastp_path="$(pgci_to_shell_path "$blastp_path")"
    makeblastdb_path="$(pgci_to_shell_path "$makeblastdb_path")"
fi

if [[ "$engine" == "diamond" ]]; then
    resolve_tool "$tool_path" "DIAMOND"
    active_tool="$tool_path"
else
    resolve_tool "$blastp_path" "blastp"
    resolve_tool "$makeblastdb_path" "makeblastdb"
    active_tool="$blastp_path"
fi

tool_uses_windows_paths() {
    case "${active_tool,,}" in
        *.exe|/mnt/*.exe) return 0 ;;
        *) return 1 ;;
    esac
}

tool_path_arg() {
    if tool_uses_windows_paths; then
        pgci_to_windows_path "$1"
    else
        printf '%s\n' "$1"
    fi
}

mkdir -p "wgdi_input/blast/db"
if [[ "$engine" == "diamond" ]]; then
    diamond_version="$("$tool_path" version 2>&1 | head -n 1)"
    printf 'key\tvalue\nengine\tdiamond\nexecutable\t%s\nversion\t%s\nevalue\t%s\nsensitivity\t%s\nmax_target_seqs\t%s\nthreads\t%s\nblock_size\t%s\ncomp_based_stats\t0\n' \
        "$tool_path" "$diamond_version" "$homology_evalue" "$sensitivity" \
        "$max_targets" "$threads" "$blocks" > "wgdi_input/blast/diamond_run_settings.tsv"
else
    printf 'key\tvalue\nengine\tblast\nblastp\t%s\nmakeblastdb\t%s\nevalue\t%s\nmax_target_seqs\t%s\nthreads\t%s\n' \
        "$blastp_path" "$makeblastdb_path" "$homology_evalue" "$max_targets" \
        "$threads" > "wgdi_input/blast/blast_run_settings.tsv"
fi

# Remove stale maps from interrupted runs, then use this shell's PID so the
# EXIT trap always removes the file it created.
command rm -f -- "$index_dir"/species_map.tsv.tmp.* 2>/dev/null || true
map_tmp="$index_dir/species_map.tsv.tmp.${BASHPID:-$$}"
cleanup_tmp_map() {
    local status=$?
    command rm -f -- "$map_tmp" 2>/dev/null || true
    return "$status"
}
trap cleanup_tmp_map EXIT
skipped_file="$index_dir/skipped_species.txt"
: > "$skipped_file"

declare -A used_ids=()
declare -A source_dir_by_id=()
declare -A pep_path_by_id=()
declare -A pep_name_by_id=()
complete_count=0

# Preprocessing writes one directory per species under wgdi_input/; older
# projects kept those directories at the project root.
species_root="."
if [[ -d wgdi_input ]]; then
    # wgdi_input also holds blast/ and config/; only treat it as the bundle
    # root when it actually contains species four-piece sets.
    for candidate in wgdi_input/*/; do
        candidate="${candidate%/}"
        case "$(basename "$candidate")" in blast|config) continue ;; esac
        if compgen -G "$candidate/*.pep" > /dev/null; then
            species_root="wgdi_input"
            break
        fi
    done
fi
for species_dir in "$species_root"/*/; do
    species_dir="${species_dir%/}"
    species_dir="${species_dir#./}"
    case "$species_dir" in
        Auto_*|_pipeline_*|analysis_v2|bin|blast|clean_v2|code|config|data|dist|external|metadata|resources|validation|wgdi_input|wgdi_run|.*) continue ;;
    esac

    pep_files=("$species_dir"/*.pep)
    if [[ ${#pep_files[@]} -ne 1 ]]; then
        echo "$species_dir：需要且只能有一个 .pep 文件，实际 ${#pep_files[@]} 个" >> "$skipped_file"
        continue
    fi

    pep_path="${pep_files[0]}"
    pep_name="$(basename "$pep_path")"
    species_id="${pep_name%.pep}"
    if [[ ! "$species_id" =~ ^[A-Za-z0-9._-]+$ ]]; then
        echo "$species_dir：物种代码含有不安全字符：$species_id" >> "$skipped_file"
        continue
    fi
    if [[ -n "${used_ids[$species_id]:-}" ]]; then
        echo "$species_dir：物种代码 $species_id 与 ${used_ids[$species_id]} 重复" >> "$skipped_file"
        continue
    fi

    cds_name="${species_id}.cds"
    gff_name="${species_id}.new.gff"
    lens_name="${species_id}.lens"
    missing=()
    for required in "$pep_name" "$cds_name" "$gff_name" "$lens_name"; do
        [[ -s "$species_dir/$required" ]] || missing+=("$required")
    done
    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "$species_dir：缺少 ${missing[*]}" >> "$skipped_file"
        continue
    fi
    readiness="$(species_readiness "$species_id")"
    if [[ "$readiness" == "no" && "$allow_not_ready" != true ]]; then
        echo "$species_dir：wgdi_ready=no，已跳过；如需显式复核可设 HOMOLOGY_ALLOW_NOT_READY=true" >> "$skipped_file"
        continue
    fi
    if [[ -z "$readiness" && "$allow_not_ready" != true ]]; then
        echo "$species_dir：缺少可读取的 wgdi_ready 记录，已跳过；如需显式复核可设 HOMOLOGY_ALLOW_NOT_READY=true" >> "$skipped_file"
        continue
    fi
    if [[ "$require_adjudication" == true ]]; then
        adj_status="$(species_adjudication_status "$species_id")"
        if [[ "$adj_status" != "ok" ]]; then
            echo "$species_dir：染色体数人工裁决未通过（$adj_status），已跳过；先用 pgci-preprocess adjudicate 签字" >> "$skipped_file"
            continue
        fi
    fi

    used_ids[$species_id]="$species_dir"
    source_dir_by_id[$species_id]="$species_dir"
    pep_name_by_id[$species_id]="$pep_name"
    pep_path_by_id[$species_id]="$species_dir/$pep_name"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$species_id" "$species_dir" "$pep_name" "$cds_name" "$gff_name" "$lens_name" >> "$map_tmp"
    complete_count=$((complete_count + 1))
done

if [[ $complete_count -eq 0 ]]; then
    echo "错误：没有任何物种同时具备 .pep、.cds、.new.gff 和 .lens。"
    echo "详见：$skipped_file"
    exit 1
fi

LC_ALL=C sort -t $'\t' -k1,1 "$map_tmp" > "$index_dir/species_map.tsv"
cut -f1 "$index_dir/species_map.tsv" > "$index_dir/folder_list.txt"
mapfile -t species_ids < "$index_dir/folder_list.txt"

echo "准备执行 ${#species_ids[@]} 个完整物种的 ${engine_label} 比对。"
echo "不完整物种已写入：$skipped_file"

for species_id in "${species_ids[@]}"; do
    db_prefix="wgdi_input/blast/db/$species_id"
    db_prefix_abs="$project_root/$db_prefix"
    if [[ "$engine" == "diamond" ]]; then
        echo "建立 DIAMOND 数据库：$species_id"
        "$tool_path" makedb \
            --in "$(tool_path_arg "${pep_path_by_id[$species_id]}")" \
            --db "$(tool_path_arg "$db_prefix")" \
            --threads "$threads" \
            --log > "wgdi_input/blast/db/${species_id}.makedb.log" 2>&1
    else
        echo "建立 BLAST 数据库：$species_id"
        species_dir="${source_dir_by_id[$species_id]}"
        pep_name="${pep_name_by_id[$species_id]}"
        (
            cd "$species_dir"
            "$makeblastdb_path" \
                -in "$pep_name" \
                -dbtype prot \
                -out "$(tool_path_arg "$db_prefix_abs")" \
                -logfile "$(tool_path_arg "$project_root/wgdi_input/blast/db/${species_id}.makeblastdb.log")"
        )
    fi
done

pair_count=0
for base_id in "${species_ids[@]}"; do
    base_pep="${pep_path_by_id[$base_id]}"
    mkdir -p "wgdi_input/blast/${base_id}_pep"
    for target_id in "${species_ids[@]}"; do
        output="wgdi_input/blast/${base_id}_pep/${base_id}_${target_id}.blast"
        partial="${output}.part.${BASHPID:-$$}"
        log="wgdi_input/blast/${base_id}_pep/${base_id}_${target_id}.${engine}.log"
        db_prefix="$project_root/wgdi_input/blast/db/$target_id"
        echo "${engine_label}：$base_id -> $target_id"
        if [[ "$engine" == "diamond" ]]; then
            "$tool_path" blastp \
                --query "$(tool_path_arg "$base_pep")" \
                --db "$(tool_path_arg "$db_prefix")" \
                --out "$(tool_path_arg "$partial")" \
                --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
                --evalue "$homology_evalue" \
                --max-target-seqs "$max_targets" \
                --threads "$threads" \
                --block-size "$blocks" \
                --comp-based-stats 0 \
                "--$sensitivity" > "$log" 2>&1
        else
            "$blastp_path" \
                -query "$(tool_path_arg "$base_pep")" \
                -db "$(tool_path_arg "$db_prefix")" \
                -out "$(tool_path_arg "$partial")" \
                -outfmt 6 \
                -evalue "$homology_evalue" \
                -max_target_seqs "$max_targets" \
                -num_threads "$threads" > "$log" 2>&1
        fi
        [[ -s "$partial" ]] || {
            echo "错误：$engine 未生成有效结果：$base_id -> $target_id"
            tail -n 20 "$log" || true
            exit 1
        }
        mv -f "$partial" "$output"
        pair_count=$((pair_count + 1))
    done
done

echo "${engine_label} 完成：本轮重新生成 $pair_count 个结果（evalue=$homology_evalue）。"
echo "完整物种数：${#species_ids[@]}；物种映射：$index_dir/species_map.tsv"
