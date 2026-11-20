#!/usr/bin/env bash

set -euo pipefail
cd "$(dirname "$0")"
project_root="$(pwd -P)"
wgdi_evalue="${WGDI_EVALUE:-1e-5}"
wgdi_tandem="${WGDI_TANDEM:-false}"
wgdi_entry="${WGDI:-}"
if [[ -z "$wgdi_entry" ]]; then
    if [[ -f "$project_root/bin/wgdi" ]]; then
        wgdi_entry="$project_root/bin/wgdi"
    else
        wgdi_entry="wgdi"
    fi
fi
if [[ ! -f "$wgdi_entry" ]] && ! command -v "$wgdi_entry" >/dev/null 2>&1; then
    echo "错误：找不到 WGDI：$wgdi_entry（可用 WGDI 环境变量指定）"
    exit 1
fi
printf -v wgdi_entry_q '%q' "$wgdi_entry"
if [[ "$wgdi_tandem" == "false" && "$wgdi_entry" == "$project_root/bin/wgdi" ]]; then
    tandem_setup="WGDI_TANDEM_ARGS=(-r)"
else
    tandem_setup="WGDI_TANDEM_ARGS=()"
fi

[[ -s folder_list.txt ]] || { echo "错误：folder_list.txt 不存在或为空。"; exit 1; }
mapfile -t species_ids < folder_list.txt

configured=0
for species1 in "${species_ids[@]}"; do
    base_folder="Auto_wgdi/${species1}_base"

    for species2 in "${species_ids[@]}"; do
        [[ "$species1" == "$species2" ]] && continue

        sub_folder="$base_folder/${species1}_${species2}"
        [[ -d "$sub_folder" ]] || { echo "错误：缺少目录 $sub_folder"; exit 1; }

        cat > "$sub_folder/run.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export MPLBACKEND=Agg

mkdir -p .wgdi_steps
cat "${species1}.cds" "${species2}.cds" > all.cds
cat "${species1}.pep" "${species2}.pep" > all.pep
WGDI_ENTRY=${wgdi_entry_q}
if [[ -f "\$WGDI_ENTRY" ]]; then
    WGDI_CMD=(bash "\$WGDI_ENTRY")
else
    WGDI_CMD=("\$WGDI_ENTRY")
fi
${tandem_setup}

run_stage() {
    local stage_name="\$1"
    local output_file="\$2"
    shift 2
    local marker=".wgdi_steps/\${stage_name}.done"
    if [[ -f "\$marker" && -f "\$output_file" ]]; then
        echo "跳过已完成阶段：\$stage_name"
        return 0
    fi
    "\$@"
    [[ -f "\$output_file" ]] || {
        echo "错误：阶段 \$stage_name 未生成 \$output_file"
        return 1
    }
    touch "\$marker"
}

run_stage dotplot "${species1}_${species2}.png" "\${WGDI_CMD[@]}" -d total.conf
run_stage collinearity "${species1}_${species2}.icl" "\${WGDI_CMD[@]}" -icl total.conf
run_stage ks "${species1}_${species2}.ks" "\${WGDI_CMD[@]}" -ks total.conf
run_stage blockinfo "${species1}_${species2}.bi.csv" "\${WGDI_CMD[@]}" -bi total.conf
run_stage correspondence "${species1}_${species2}.bi.c.csv" "\${WGDI_CMD[@]}" "\${WGDI_TANDEM_ARGS[@]}" -c total.conf
run_stage blockks "${species1}_${species2}.bk.png" "\${WGDI_CMD[@]}" "\${WGDI_TANDEM_ARGS[@]}" -bk total.conf
run_stage kspeaks "${species1}_${species2}.kp.png" "\${WGDI_CMD[@]}" "\${WGDI_TANDEM_ARGS[@]}" -kp total.conf
run_stage peaksfit "${species1}_${species2}.pf.png" "\${WGDI_CMD[@]}" -pf total.conf
EOF

        cat > "$sub_folder/total.conf" <<EOF
[dotplot]
blast = ${species1}_${species2}.blast
gff1 = ${species1}.new.gff
gff2 = ${species2}.new.gff
lens1 = ${species1}.lens
lens2 = ${species2}.lens
genome1_name = ${species1}
genome2_name = ${species2}
multiple = 1
score = 100
evalue = ${wgdi_evalue}
repeat_number = 10
position = order
blast_reverse = false
ancestor_left = none
ancestor_top = none
markersize = 0.5
figsize = 10,10
savefig = ${species1}_${species2}.png

[collinearity]
blast = ${species1}_${species2}.blast
gff1 = ${species1}.new.gff
gff2 = ${species2}.new.gff
lens1 = ${species1}.lens
lens2 = ${species2}.lens
blast_reverse = false
multiple = 1
process = 8
evalue = ${wgdi_evalue}
score = 100
grading = 50,25,10
mg = 25,25
pvalue = 1
repeat_number = 10
position = order
savefile = ${species1}_${species2}.icl

[ks]
cds_file = all.cds
pep_file = all.pep
align_software = muscle
pairs_file = ${species1}_${species2}.icl
ks_file = ${species1}_${species2}.ks

[blockinfo]
blast = ${species1}_${species2}.blast
gff1 = ${species1}.new.gff
gff2 = ${species2}.new.gff
lens1 = ${species1}.lens
lens2 = ${species2}.lens
collinearity = ${species1}_${species2}.icl
score = 100
evalue = ${wgdi_evalue}
repeat_number = 20
position = order
ks = ${species1}_${species2}.ks
ks_col = ks_NG86
savefile = ${species1}_${species2}.bi.csv

[correspondence]
blockinfo = ${species1}_${species2}.bi.csv
lens1 = ${species1}.lens
lens2 = ${species2}.lens
tandem = ${wgdi_tandem}
tandem_length = 200
pvalue = 0.2
block_length = 5
tandem_ratio = 0.5
multiple = 1
homo = -1,1
savefile = ${species1}_${species2}.bi.c.csv

[blockks]
lens1 = ${species1}.lens
lens2 = ${species2}.lens
genome1_name = ${species1}
genome2_name = ${species2}
blockinfo = ${species1}_${species2}.bi.csv
pvalue = 0.2
tandem = ${wgdi_tandem}
tandem_length = 200
markersize = 1
area = 0,2
block_length = 5
figsize = 8,8
savefig = ${species1}_${species2}.bk.png

[kspeaks]
blockinfo = ${species1}_${species2}.bi.csv
pvalue = 0.2
tandem = ${wgdi_tandem}
block_length = 5
ks_area = 0,10
multiple = 1
homo = 0,1
fontsize = 9
area = 0,3
figsize = 10,6.18
savefig = ${species1}_${species2}.kp.png
savefile = ${species1}_${species2}.kp

[peaksfit]
blockinfo = ${species1}_${species2}.bi.csv
mode = median
bins_number = 200
ks_area = 0,10
fontsize = 9
area = 0,3
figsize = 10,6.18
shadow = true
savefig = ${species1}_${species2}.pf.png

[ksfigure]
ksfit = ks_fit_result.csv
labelfontsize = 15
legendfontsize = 15
xlabel = none
ylabel = none
title = none
area = 0,2
figsize = 10,6.18
shadow = false
savefig = ${species1}_${species2}.kf.png
EOF

        find "$sub_folder/.wgdi_steps" -maxdepth 1 -type f -name '*.done' -delete 2>/dev/null || true
        rm -f "$sub_folder/.wgdi.done"
        configured=$((configured + 1))
        echo "已生成不同物种配置：$species1 / $species2"
    done
done

echo "不同物种 WGDI 配置完成：$configured 个。"
