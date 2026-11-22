# 植物基因组 WGDI + DIAMOND 比较分析流程

把本地植物基因组注释与序列整理成 WGDI 可用的四件套，并继续完成 **DIAMOND 同源比对、WGDI 点图与共线性块鉴定**。Ks 分析作为可选分支，由使用者明确开启。

当前版本：**3.7.0**

- [快速上手](QUICKSTART.md)
- [实验室整体流程](code/docs/WORKFLOW.md)
- [版本记录](CHANGELOG.md)
- [验证与发布门禁](VALIDATION.md)
- [v3.7.0 Release](https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.7.0)
- [Agent Skill](https://github.com/HXC-code/plant-genome-preprocess-skill)

## 流程范围

本版本按实验室《整体流程》交付 **第 1–4 步**，并把每一步的结果整理成可追溯文件。

```text
1 数据收集与整理
      ↓
2 预处理：pep + cds + new.gff + lens
      ↓
3 DIAMOND 物种内与物种间同源比对
      ↓
4 WGDI 点图与共线性块
      └── Ks 分支（使用者按需开启）
```

| 流程步骤 | 本版本交付 |
| --- | --- |
| 1 数据收集与整理 | 用户按推荐来源下载；工具注册、识别、路由并校验本地文件；不可用物种记录原因 |
| 1 染色体数量确认 | 输出 `chromosome_map.tsv`，由复核人确认最终数量并签字 |
| 2 数据预处理 | 一基因一代表转录本、序列质控、染色体筛选、统一命名和追溯 |
| 3 DIAMOND 比对 | 物种内与物种间两两比对，默认 `very-sensitive`、e-value `1e-5` |
| 4 WGDI 共线性 | 默认产出点图 `.png` 与共线性块 `.icl` |
| 4′ Ks 可选分支 | 明确开启后追加 `.ks`、blockinfo、correspondence、block-Ks 和 Ks 峰 |

第 4 步的点图、共线性块和统计表可以直接交给后续系统发育、文献核对、多重共线性和圈图分析使用。后续步骤由对应分析流程或人工判断接管。

## 3.7.0 重点变化

- 代表转录本改为 **CDS 最长优先**；CDS 等长时按稳定 transcript ID 字典序确定。
- CLI 默认采用严格染色体模式，只保留有染色体挂载证据的序列。
- 无染色体证据时以 `no_chromosome_sequences` 失败关闭；碎片化组装经确认后可显式使用 `--sequence-mode adaptive`。
- Ks 分析默认关闭，通过 `WGDI_KS=true` 或 Windows `-WithKs` 开启。
- SOP 追溯文件采用 `chromosome_map.tsv`、`removed_list.tsv` 和 `qc_record.md`。
- Skill 已更名为 `plant-wgdi-preprocess`，范围与实验室流程第 1–4 步一致。

## 安装环境

预处理要求：

```text
Python >= 3.11
GBFF/GenBank 转换需要 biopython >= 1.83
```

安装程序包：

```powershell
cd code
python -m pip install -e ".[gbff]"
```

下游按需准备：

| 用途 | 工具 |
| --- | --- |
| 同源比对 | DIAMOND，或 BLAST+ |
| WGDI | `wgdi==0.75` |
| Ks 分支 | MUSCLE、PAML `yn00`、PAL2NAL v14 |
| Shell 脚本 | Bash；Windows 使用 Git Bash |

WGDI 环境：

```powershell
powershell -ExecutionPolicy Bypass -File .\bootstrap_wgdi.ps1
```

环境检查：

```powershell
python -m pgci_preprocess doctor `
  --project-root D:\plant_project `
  --workflow full
```

## 输入

单物种常规输入：

```text
GFF/GFF3 注释 + CDS FASTA + 蛋白 FASTA
```

genome FASTA 是常用输入，可用来核对 region 和真实序列长度。缺少 genome 时，程序可由 CDS 与转录本坐标回填 region。

| 数据 | 支持格式 |
| --- | --- |
| 注释 | GFF/GFF3、GTF、GBFF/GenBank |
| 序列 | genome、CDS、蛋白 FASTA |
| 压缩 | gzip、bzip2、xz、单成员 ZIP |
| 辅助数据 | 功能注释 CSV/TSV |

GFF/GFF3 支持：

```text
standard
transcript_only
simplified5
gene_span
gene_cds_direct
```

文件内容负责最终格式和角色判定，后缀作为提示。同一角色有多个候选时，在清单中写 `A|B`，程序选择唯一兼容文件。

## 快速开始

### 1. 创建项目

```powershell
python -m pgci_preprocess init-project --project-root D:\plant_project
python -m pgci_preprocess doctor `
  --project-root D:\plant_project `
  --allow-empty
```

### 2. 注册物种

```powershell
python -m pgci_preprocess init-project `
  --project-root D:\plant_project `
  --show-id 2848 `
  --scientific-name "Morus notabilis" `
  --order Rosales `
  --family Moraceae `
  --genome D:\raw\Morus.genome.fasta.gz `
  --gff D:\raw\Morus.gff3.gz `
  --cds D:\raw\Morus.CDS.fasta.gz `
  --pep D:\raw\Morus.Protein.faa.gz
```

`--genome` 可省略。GTF/GBFF 可先转换，或交给 `prepare/batch` 自动路由：

```powershell
python -m pgci_preprocess convert-annotation `
  --input D:\raw\species.gtf.gz `
  --output-dir D:\plant_project\data\derived `
  --prefix species
```

### 3. 运行预处理

```powershell
python -m pgci_preprocess run --project-root D:\plant_project
```

默认严格染色体模式：

```text
只保留有染色体挂载证据的核基因组序列
剔除细胞器
无染色体证据时报告 no_chromosome_sequences
```

碎片化组装经复核后可改为：

```powershell
python -m pgci_preprocess run `
  --project-root D:\plant_project `
  --sequence-mode adaptive `
  --overwrite
```

运行完成后会直接打印：

```text
[OK] Mun show_id=2848 Morus notabilis genes=... seqs=10 ready=yes reason=
    chromosomes 10 (id:genes) 1:4210  2:3980  3:4521 ...
    inventory: D:\plant_project\analysis_v2\Mun\Mun.chromosome_map.tsv
```

### 4. 复核染色体

```powershell
python -m pgci_preprocess chromosomes `
  --project-root D:\plant_project `
  --species Mun
```

重点查看：

```text
analysis_v2\Mun\Mun.chromosome_map.tsv
analysis_v2\Mun\Mun.chr_map.tsv
analysis_v2\Mun\Mun.qc_record.md
analysis_v2\Mun\Mun.qc.json
```

确认：

- 最终保留多少条染色体；
- 每条染色体有多少代表基因；
- 原始序列 ID 与物种核型是否对应；
- 是否出现未挂载 scaffold；
- 少于 100 个基因的序列是否需要人工复核。

### 5. 人工裁决

```powershell
python -m pgci_preprocess adjudicate `
  --project-root D:\plant_project `
  --species Mun `
  --confirmed-n 10 `
  --decision accept `
  --reviewer "<实际复核人>" `
  --rationale "<核型或文献依据>"
```

裁决记录写入：

```text
resources\adjudications\Mun.json
```

记录绑定 QC 版本与 SHA-256。重跑后旧签字自动变为 `stale`，需要重新复核。

### 6. 严格验收

```powershell
python -m pgci_preprocess validate `
  --project-root D:\plant_project `
  --refresh-links `
  --require-adjudication

python -m pgci_preprocess check-names --project-root D:\plant_project
python -m pgci_preprocess audit --project-root D:\plant_project
```

验收通过后，四件套位于：

```text
analysis_v2\<Abbr>\<Abbr>.pep
analysis_v2\<Abbr>\<Abbr>.cds
analysis_v2\<Abbr>\<Abbr>.new.gff
analysis_v2\<Abbr>\<Abbr>.lens
```

## 关键 SOP 规则

### 代表转录本

每个基因保留一个代表转录本：

1. CDS 最长优先；
2. CDS 等长时按稳定 transcript ID 字典序；
3. 再按 protein ID 稳定排序。

蛋白长度和文件出现顺序不参与选择，确保跨机器结果可复现。

### 染色体筛选

CLI 默认：

```text
--sequence-mode chromosomes
```

严格模式只保留注释明确标为染色体的序列。未挂载 scaffold、contig 和细胞器进入剔除记录。

碎片化组装选项：

```text
--sequence-mode adaptive
```

启用前需要确认物种组装与分析目标允许保留 scaffold。

### 每染色体基因数下限

| 参数 | CLI 默认值 |
| --- | ---: |
| `--min-genes-per-sequence` | 50 |
| `--coverage` | 0.90 |
| `--max-sequences` | 200 |
| `--min-genes-per-chromosome` | 100 |

`--min-genes-per-chromosome 100` 表示：

```text
99 个代表基因  移除
100 个代表基因 保留
101 个代表基因 保留
```

过滤发生在序列选择之后；移除原因写入 `chr_map.tsv` 和 QC，保留率与 `wgdi_ready` 按过滤后结果重新计算。

CLI 默认 100；Python API 的 `PipelineOptions` 默认 0，需要显式传入 100。传入 0 可关闭过滤。

### `wgdi_ready`

| 状态 | 基因保留率 | 用途 |
| --- | --- | --- |
| `yes` | ≥ 70%，且序列数不超过上限 | 进入标准下游 |
| `marginal` | ≥ 50% | 结合组装与筛选记录复核 |
| `no` | < 50% | 归档或专项复核；默认不进入下游 |

`wgdi_ready` 是当前筛选集的基因保留率，论文分析应同时记录 BUSCO、组装连续性和倍性信息。

## 输出

### `analysis_v2/<Abbr>/`

```text
<Abbr>.pep
<Abbr>.cds
<Abbr>.new.gff
<Abbr>.lens
<Abbr>.all.pep.fa
<Abbr>.all.cds.fa
<Abbr>.all.gff
<Abbr>.all.lens
<Abbr>.gene_map.tsv
<Abbr>.chr_map.tsv
<Abbr>.chromosome_map.tsv
<Abbr>.qc.json
<Abbr>.qc_record.md
```

### `clean_v2/<Abbr>/`

```text
<Abbr>.gene_map.tsv
<Abbr>.chr_map.tsv
<Abbr>.chromosome_map.tsv
<Abbr>.removed_list.tsv
```

### 文件格式

`new.gff` 固定七列：

```text
序列号  唯一基因ID  起点  终点  链方向  序列内顺序号  原始基因ID
```

`lens` 固定三列：

```text
序列号  序列长度  代表基因数
```

GFF、CDS、PEP 的基因 ID 集合一致；序列编号、坐标、顺序号、长度和基因数经过闭合校验。

## 批量处理

状态表必需列：

```text
show_id
Order
Family
Species
status
folder
current_genome
current_gff_or_gtf
current_cds
current_pep
```

`status` 使用 `complete`：

```powershell
python -m pgci_preprocess batch `
  --project-root D:\plant_project `
  --source-status D:\plant_project\resources\source_status.tsv `
  --expected-count 20
```

批量顺序：

```text
prepare → 可选 GFF 规范化 → run → validate
```

批量运行完成后仍需逐物种复核染色体清单、完成人工裁决，再执行严格验收。

## 下游分析

### 完成

```powershell
bash 0check_data.sh --preprocess --refresh-links
```

### DIAMOND

```bash
bash run_homology.sh
```

默认参数：

```text
DIAMOND very-sensitive
e-value 1e-5
物种内 + 物种间两两比对
```

兼容入口：

```bash
bash 1run_diamond.sh
bash 1run_blast.sh
```

### WGDI 核心模式

```bash
bash 2run_pre_wgdi.sh
bash 3blast_to_wgdi.sh
python code/tools/generate_wgdi_configs.py --mode all
bash 6run_wgdi.sh
```

核心模式产出：

```text
<prefix>.png
<prefix>.icl
```

### Ks 可选分支

```bash
WGDI_KS=true python code/tools/generate_wgdi_configs.py --mode all
WGDI_KS=true bash 6run_wgdi.sh
```

追加产物：

```text
<prefix>.ks
<prefix>.bi.csv
<prefix>.bi.c.csv
<prefix>.bk.png
<prefix>.kp.png
<prefix>.pf.png
```

Windows 总控：

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all_pipeline.ps1
powershell -ExecutionPolicy Bypass -File .\run_all_pipeline.ps1 -WithKs
```

下游按 `wgdi_ready` 与人工裁决状态选择物种；跳过记录写入 `skipped_species.txt`。

## 数据检查

程序检查：

- CDS 长度是否为 3 的倍数；
- CDS 与蛋白翻译是否一致；
- 内部终止、非法字符和缺失序列；
- 坐标、链方向、region 和序列长度；
- gene、mRNA 和 CDS 层级；
- 基因、蛋白和 FASTA ID 重复或冲突；
- 核基因组 taxon 污染；
- 四件套的 ID 与序列集合闭合；
- 输入 SHA-256、来源清单与转换链。

所有剔除记录保存在：

```text
clean_v2/<Abbr>/<Abbr>.removed_list.tsv
```

`qc_record.md` 提供人类可读摘要；`qc.json` 保存完整机器可读结果。

## 命名

物种缩写从以下表读取：

```text
code/resources/species_abbreviations.tsv
```

新物种按拉丁名生成全局唯一缩写。基因 ID 形如：

```text
Mun1g0001
```

旧 `s<show_id>` 项目可迁移：

```powershell
python -m pgci_preprocess migrate-artifacts --project-root D:\plant_project --dry-run
python -m pgci_preprocess migrate-artifacts --project-root D:\plant_project
```

迁移后运行 `check-names`。

## 使用注意

- 严格染色体模式适用于具有明确染色体挂载信息的 GWH、NCBI 和高质量组装数据。
- 简化 GFF 使用裸数字 seqid、且缺少 region 行时，先做小样本验证；经确认后再选择 adaptive。
- 100 基因下限是机械规则，可能移除基因数少的真实小染色体，需要结合 `chr_map.tsv` 复核。
- 多倍体和高杂合数据应使用能区分等位 contig、homeolog 和亚基因组的注释来源。
- 细胞器判断同时查看 `organelle_inference.inferred` 与 `review_candidates`。
- 大规模 Ks 分析先在代表性小样本验证参数和运行时间。
- WGD 定年、Ks 峰解释、系统发育拓扑和物种演化结论由研究团队结合背景证据判断。

## 版本与 Release

| 版本 | 状态 | 说明 |
| --- | --- | --- |
| **v3.7.0** | 当前主页版本 | 对齐整体流程第 1–4 步；CDS 优先代表转录本；严格染色体模式；Ks 可选 |
| [v3.6.2](https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.6.2) | 历史 Release | 修复非 UTF-8 Windows 控制台打印染色体清单异常 |
| [v3.6.1](https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.6.1) | 历史 Release | `run` 调用时打印染色体清单 |
| [v3.6.0](https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.6.0) | 历史 Release | 引入染色体清单、100 基因下限和人工裁决门禁 |
| [v3.5.5](https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.5.5) | 历史 Release | 线粒体推断、readiness 和产物完整性检查 |

历史 Release 保持原有 tag 与内容，不改写、不移动。当前首页、开发和后续修正以 `main` 分支的 v3.7.0 系列为主。

## 测试与验证

```powershell
cd code
python -m pytest -q
python -m ruff check src tests scripts tools
python -m mypy src --ignore-missing-imports
```

当前验证：

```text
212 项 Python 测试
Ruff / Mypy 通过
Skill pipeline 镜像一致性检查通过
```

冒烟回执：

```text
validation/smoke-report-3.7.0.json        Ks 模式：4 组比较 × 8 产物
validation/smoke-report-3.7.0-core.json   核心模式：4 组比较 × 2 产物
```

详细环境与发布门禁见 [VALIDATION.md](VALIDATION.md)。

## 仓库关系

- `plant-genome-preprocess`：程序实现、测试、发布材料和主文档；
- `plant-genome-preprocess-skill`：Agent 调用入口 `plant-wgdi-preprocess` 及同步的 `pipeline/` 运行副本。

两仓由镜像检查保持一致。

## 许可

项目代码使用 MIT License。WGDI、DIAMOND、BLAST+、MUSCLE、PAML、Perl 和 PAL2NAL 按其上游许可使用，详见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
