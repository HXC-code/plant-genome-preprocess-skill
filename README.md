# Plant WGDI Preprocess

植物基因组 WGDI 与 DIAMOND 分析 skill。把本地下载的基因组注释整理成可验证、
可追溯的 WGDI 输入，并继续完成 DIAMOND 同源比对、WGDI 点图和共线性块。

当前版本：**3.8.0**

## 覆盖范围

```text
1. 本地数据整理、内容识别与同套组装配对
2. WGDI 四件套预处理、染色体筛选、人工裁决与 QC
3. DIAMOND 物种内和物种间两两同源比对
4. WGDI 点图与共线性块
   └── Ks 可选分支（默认关闭）
```

## 主要能力

- 跨机器扫描数据、工具和解释器
- 按文件内容识别 GFF/GFF3、GTF、GBFF、genome、CDS 和蛋白。
- 自动配对同一套组装，避免混用注释和序列版本。
- 统一注释层级、ID、region、坐标、链和序列长度。
- 按“CDS 最长优先 + 稳定转录本 ID”选择代表转录本。
- 默认启用严格染色体模式，为碎片化组装保留显式 adaptive 逃生门。
- 输出可审计的基因映射、染色体映射、剔除记录和 QC。
- 使用项目相对路径生成 DIAMOND、WGDI 配置并集中保存结果。
- 阶段级断点续跑，输入或配置变化时自动重建对应阶段。

## 项目布局

```text
<PROJECT_ROOT>/
├── data/
│   └── <show_id>_<Order>_<Family>(<Species>)/
│       ├── genome.fa*
│       ├── annotation.gff*
│       ├── cds.fa*
│       └── protein.fa*
├── metadata/
│   ├── machine_paths.json
│   ├── download_inventory.tsv
│   ├── sequence_mode.tsv
│   ├── discovery_report.tsv
│   ├── species_map.tsv
│   └── adjudications/
├── resources/                   # metadata/ 的兼容别名
├── wgdi_input/
│   ├── <Abbr>/                  # 四件套和追溯文件
│   ├── blast/                   # DIAMOND 数据库、结果和设置
│   └── config/                  # 扁平 <pair>.conf
└── wgdi_output/
    ├── results/                 # .png、.icl、可选 Ks 产物
    ├── logs/
    └── state/
```

## 安装

要求 Python 3.11 或更高版本。

```bash
cd pipeline/code
python3 -m pip install -e ".[gbff]"
```

普通预处理只需要 Python 标准库。GTF/GBFF 转换和部分测试需要 Biopython。

下游工具：

| 功能 | 工具 |
| --- | --- |
| 同源比对 | DIAMOND，或 BLAST+ |
| WGDI | `wgdi==0.75` |
| Ks 分支 | MUSCLE、PAML `yn00`、Perl、PAL2NAL v14 |
| Shell 入口 | Linux/WSL Bash；Windows 使用 Git Bash |

## 1. 新机器先扫描

新机器上先确认 Python 和 WGDI 环境。skill 不扫描其他用户的 Linux 主目录；
请激活目标 conda/venv、使用项目内环境，或显式提供解释器：

```bash
export WGDI_PYTHON=/confirmed/path/to/python
```

扫描项目：

```bash
python3 /path/to/plant-wgdi-preprocess/pipeline/discover.py \
  /path/to/raw-data \
  --project-root /path/to/project
```

数据已经在项目 `data/` 下时可省略数据参数：

```bash
python3 /path/to/plant-wgdi-preprocess/pipeline/discover.py
```

Windows PowerShell：

```powershell
python <skill>\pipeline\discover.py E:\raw_data --project-root D:\plant_project
```

扫描结果写入 `metadata/`，并记住项目根。之后命令可以省略 `--project-root`。

可用以下命令检查环境：

```bash
python3 -m pgci_preprocess doctor
bash /path/to/plant-wgdi-preprocess/pipeline/bin/wgdi --version
```

## 2. 建立项目和登记数据

```bash
python3 -m pgci_preprocess init-project --project-root /path/to/project
```

登记单个物种：

```bash
python3 -m pgci_preprocess init-project \
  --project-root /path/to/project \
  --show-id 2848 \
  --scientific-name "Morus notabilis" \
  --order Rosales \
  --family Moraceae \
  --genome genome.fasta.gz \
  --gff annotation.gff3.gz \
  --cds cds.fa.gz \
  --pep pep.fa.gz
```

批量数据放入：

```text
data/<show_id>_<Order>_<Family>(<Species>)/
```

然后检查 `metadata/download_inventory.tsv`，使用 `prepare` / `batch` 批量路由。

## 3. 预处理和人工裁决

```bash
python3 -m pgci_preprocess run --sequence-mode from-map
python3 -m pgci_preprocess chromosomes --species <code>
```

确认最终染色体数：

```bash
python3 -m pgci_preprocess adjudicate \
  --species <code> \
  --confirmed-n <N> \
  --decision accept \
  --reviewer "<实际复核人>" \
  --rationale "<核型或文献依据>"
```

严格验收：

```bash
python3 -m pgci_preprocess validate --require-adjudication
python3 -m pgci_preprocess check-names
python3 -m pgci_preprocess audit
```

每个物种产出：

```text
wgdi_input/<Abbr>/<Abbr>.pep
wgdi_input/<Abbr>/<Abbr>.cds
wgdi_input/<Abbr>/<Abbr>.new.gff
wgdi_input/<Abbr>/<Abbr>.lens
```

同时保留 `gene_map.tsv`、`chromosome_map.tsv`、`chr_map.tsv`、
`removed_list.tsv`、`qc.json` 和 `qc_record.md`。

## 4. DIAMOND 和 WGDI

```bash
cd /path/to/plant-wgdi-preprocess/pipeline
bash run_homology.sh
bash 2gen_wgdi_conf.sh all
bash 3run_wgdi.sh
```

默认产物：

```text
wgdi_input/blast/<A>_pep/<A>_<B>.blast
wgdi_output/results/<pair>.png
wgdi_output/results/<pair>.icl
```

配置和结果全部使用项目相对路径。比较目录不再复制四件套、BLAST 或 `run.sh`，
项目可以整体搬移。

### 可选 Ks

只有明确需要 Ks 时才执行：

```bash
WGDI_KS=true bash 2gen_wgdi_conf.sh all
WGDI_KS=true bash 3run_wgdi.sh
```

Windows 总控：

```powershell
powershell -ExecutionPolicy Bypass -File <skill>\pipeline\run_all_pipeline.ps1
powershell -ExecutionPolicy Bypass -File <skill>\pipeline\run_all_pipeline.ps1 -WithKs
```

## 输入格式

| 输入 | 支持范围 |
| --- | --- |
| 注释 | GFF/GFF3、GTF、GBFF/GenBank |
| 序列 | genome、CDS、蛋白 FASTA |
| 压缩 | gzip、bzip2、xz、单成员 ZIP |
| 辅助信息 | 七列功能注释 CSV/TSV |

非标准后缀只要能按内容识别，也可进入流程。`.lens` 是输出文件，不是输入。

## 异常行为

| 情况 | 结果 |
| --- | --- |
| 无染色体证据 | 严格模式报 `no_chromosome_sequences` |
| 碎片化组装 | 人工确认后使用 `--sequence-mode adaptive` |
| 只有 genome、没有注释 | 报 `annotation_missing` |
| 没有 genome | 允许；由 CDS 和转录本坐标回填 |
| GTF/GBFF | 自动或显式转换 |
| 多套候选 | 评分选择，歧义时停止 |
| 细胞器混入 | 推断排除并记录复核候选 |
| 旧布局 | `migrate-layout` / `migrate-artifacts` 迁移 |

## 验证

```bash
cd pipeline/code
python3 -m pytest -q
python3 -m ruff check src tests tools ../discover.py
python3 -m mypy src
```

当前测试为 259 项。完整验收记录和发布说明见
[`pipeline/VALIDATION.md`](pipeline/VALIDATION.md) 与
[`CHANGELOG.md`](CHANGELOG.md)。

## 文档

- [Agent 执行规则](SKILL.md)
- [快速上手](pipeline/QUICKSTART.md)
- [实验室流程映射](pipeline/code/docs/WORKFLOW.md)
- [数据校验说明](pipeline/code/docs/data-validation.md)
- [程序仓说明](pipeline/code/README.md)
