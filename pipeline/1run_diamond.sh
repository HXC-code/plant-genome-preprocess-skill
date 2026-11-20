#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob
cd "$(dirname "$0")"

diamond_path="${DIAMOND:-$(command -v diamond || true)}"
diamond_threads="${DIAMOND_THREADS:-8}"
diamond_blocks="${DIAMOND_BLOCKS:-4}"
diamond_evalue="${DIAMOND_EVALUE:-1e-5}"
diamond_sensitivity="${DIAMOND_SENSITIVITY:-very-sensitive}"
diamond_max_targets="${DIAMOND_MAX_TARGETS:-20}"

case "$diamond_sensitivity" in
    sensitive|more-sensitive|very-sensitive|ultra-sensitive) ;;
    *) echo "错误：不支持的 DIAMOND_SENSITIVITY：$diamond_sensitivity"; exit 2 ;;
esac
[[ "$diamond_evalue" =~ ^[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$ ]] || {
    echo "错误：DIAMOND_EVALUE 不是有效数值：$diamond_evalue"; exit 2;
}
[[ "$diamond_threads" =~ ^[1-9][0-9]*$ && "$diamond_max_targets" =~ ^[1-9][0-9]*$ ]] || {
    echo "错误：DIAMOND_THREADS 和 DIAMOND_MAX_TARGETS 必须是正整数"; exit 2;
}

if [[ ! -x "$diamond_path" ]] && ! command -v "$diamond_path" >/dev/null 2>&1; then
    echo "错误：找不到 DIAMOND：$diamond_path"
    exit 1
fi

mkdir -p "Auto_pep/db"
diamond_version="$("$diamond_path" version 2>&1 | head -n 1)"
printf 'key\tvalue\nexecutable\t%s\nversion\t%s\nevalue\t%s\nsensitivity\t%s\nmax_target_seqs\t%s\nthreads\t%s\nblock_size\t%s\ncomp_based_stats\t0\n' \
    "$diamond_path" "$diamond_version" "$diamond_evalue" "$diamond_sensitivity" \
    "$diamond_max_targets" "$diamond_threads" "$diamond_blocks" > "Auto_pep/diamond_run_settings.tsv"
map_tmp="$(mktemp)"
skipped_file="skipped_species.txt"
: > "$skipped_file"

declare -A used_ids=()
complete_count=0

for species_dir in */; do
    species_dir="${species_dir%/}"
    case "$species_dir" in
        Auto_*|_pipeline_*|.*) continue ;;
    esac

    pep_files=("$species_dir"/*.pep)
    if [[ ${#pep_files[@]} -ne 1 ]]; then
        echo "$species_dir：需要且只能有一个 .pep 文件，实际 ${#pep_files[@]} 个" >> "$skipped_file"
        continue
    fi

    pep_name="$(basename "${pep_files[0]}")"
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

    used_ids[$species_id]="$species_dir"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$species_id" "$species_dir" "$pep_name" "$cds_name" "$gff_name" "$lens_name" >> "$map_tmp"
    complete_count=$((complete_count + 1))
done

if [[ $complete_count -eq 0 ]]; then
    echo "错误：没有任何物种同时具备 .pep、.cds、.new.gff 和 .lens。"
    rm -f "$map_tmp"
    exit 1
fi

LC_ALL=C sort -t $'\t' -k1,1 "$map_tmp" > species_map.tsv
cut -f1 species_map.tsv > folder_list.txt
rm -f "$map_tmp"
mapfile -t species_ids < folder_list.txt

declare -A source_dir_by_id=()
declare -A pep_name_by_id=()
while IFS=$'\t' read -r species_id species_dir pep_name cds_name gff_name lens_name; do
    source_dir_by_id[$species_id]="$species_dir"
    pep_name_by_id[$species_id]="$pep_name"
done < species_map.tsv

echo "准备执行 ${#species_ids[@]} 个完整物种的 DIAMOND 比对。"
echo "不完整物种已写入：$skipped_file"

for species_id in "${species_ids[@]}"; do
    species_dir="${source_dir_by_id[$species_id]}"
    pep_name="${pep_name_by_id[$species_id]}"
    db_prefix="Auto_pep/db/$species_id"
    echo "建立 DIAMOND 数据库：$species_id"
    "$diamond_path" makedb \
        --in "$species_dir/$pep_name" \
        --db "$db_prefix" \
        --threads "$diamond_threads" \
        --log > "Auto_pep/db/${species_id}.makedb.log" 2>&1
done

pair_count=0
for base_id in "${species_ids[@]}"; do
    base_dir="${source_dir_by_id[$base_id]}"
    base_pep="$base_dir/${pep_name_by_id[$base_id]}"
    mkdir -p "Auto_pep/${base_id}_pep"

    for target_id in "${species_ids[@]}"; do
        output="Auto_pep/${base_id}_pep/${base_id}_${target_id}.blast"
        partial="${output}.part"
        echo "DIAMOND：$base_id -> $target_id"
        rm -f "$partial"
        "$diamond_path" blastp \
            --query "$base_pep" \
            --db "Auto_pep/db/$target_id" \
            --out "$partial" \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --evalue "$diamond_evalue" \
            --max-target-seqs "$diamond_max_targets" \
            --threads "$diamond_threads" \
            --block-size "$diamond_blocks" \
            --comp-based-stats 0 \
            "--$diamond_sensitivity" > "Auto_pep/${base_id}_pep/${base_id}_${target_id}.diamond.log" 2>&1
        [[ -s "$partial" ]] || {
            echo "错误：DIAMOND 未生成有效结果：$base_id -> $target_id"
            exit 1
        }
        mv -f "$partial" "$output"
        pair_count=$((pair_count + 1))
    done
done

echo "DIAMOND 完成：本轮重新生成 $pair_count 个结果（evalue=$diamond_evalue, $diamond_sensitivity）。"
echo "完整物种数：${#species_ids[@]}；物种映射：species_map.tsv"
