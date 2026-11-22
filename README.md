# 植物基因组 WGDI 预处理技能

**Plant Genome WGDI Preprocessing Skill**

`plant-wgdi-preprocess` 是面向植物比较基因组研究的 Agent Skill。它扫描本地下载的植物基因组数据，自动识别注释、genome、CDS 和蛋白文件，清洗并规范化为 WGDI 可直接使用的四件套：

```text
<Abbr>.pep
<Abbr>.cds
<Abbr>.new.gff
<Abbr>.lens
```

预处理完成后，同一目录还保留染色体清单、基因映射、剔除记录和 QC，用于回答：

- 这个物种最终保留了多少条染色体/序列；
- 每条染色体有多少个代表基因；
- 哪些序列或基因为什么被剔除；
- 当前结果是否适合交给后续 DIAMOND 和 WGDI 分析。

当前代码版本：**3.8.0**

- [程序仓](https://github.com/HXC-code/plant-genome-preprocess)
- [完整快速上手](pipeline/QUICKSTART.md)
- [数据格式与验证](pipeline/code/docs/data-validation.md)
- [验证门禁](pipeline/VALIDATION.md)
- [版本记录](CHANGELOG.md)

## 这个 Skill 主要做什么

| 能力 | 结果 |
| --- | --- |
| 跨机器数据扫描 | 自动找出数据根、代码根、项目根和可用工具 |
| 文件角色识别 | 从文件内容识别 GFF/GFF3、GTF、GBFF、genome、CDS 和蛋白 |
| 同一套组装配对 | 以注释为锚点匹配 genome、CDS 和 PEP，避免混用不同装配版本 |
| 注释规范化 | 修复 gene/mRNA/CDS 层级、补充 region、统一 ID |
| 代表转录本 | 每个基因保留一个代表转录本 |
| 染色体筛选 | 默认保留有染色体证据的序列，也支持碎片化组装的 adaptive 模式 |
| 基因数下限 | 默认移除代表基因数少于 100 的 WGDI 序列 |
| 数据质控 | 检查 CDS、翻译、坐标、链、ID、taxon 和重复冲突 |
| 可追溯输出 | 输出 gene map、chr map、chromosome map、removed list 和 QC |
| 人工裁决 | 复核最终染色体数，记录复核人和判断依据 |

DIAMOND、WGDI 点图和共线性分析是四件套的下游用途。本 Skill 的首页和实现重点是数据预处理与可交付输入构建。

## 新项目目录

v3.8 起推荐使用单层物种目录和单一输出目录：

```text
<PROJECT_ROOT>/
├── data/
│   ├── <show_id>_<Order>_<Family>(<Species>)/
│   │   ├── <genome>
│   │   ├── <gff>
│   │   ├── <cds>
│   │   └── <pep>
│   └── ...
│
├── resources/
│   ├── machine_paths.json
│   ├── download_inventory.tsv
│   ├── sequence_mode.tsv
│   ├── discovery_report.tsv
│   ├── species_abbreviations.tsv
│   ├── format_exceptions.tsv
│   ├── assembly_accessions.tsv
│   ├── adjudications/
│   │   └── <Abbr>.json
│   └── wgdi_v2_index.tsv
│
└── wgdi_input/
    ├── <Abbr1>/
    │   ├── <Abbr1>.pep
    │   ├── <Abbr1>.cds
    │   ├── <Abbr1>.new.gff
    │   ├── <Abbr1>.lens
    │   ├── <Abbr1>.gene_map.tsv
    │   ├── <Abbr1>.chr_map.tsv
    │   ├── <Abbr1>.chromosome_map.tsv
    │   ├── <Abbr1>.removed_list.tsv
    │   ├── <Abbr1>.qc.json
    │   └── <Abbr1>.qc_record.md
    └── <Abbr2>/
        └── ...
```

一个物种只有一个权威目录 `wgdi_input/<Abbr>/`。WGDI 可以直接使用该目录中的四件套。

旧项目中的以下布局仍可读取：

```text
data/downloads/
analysis_v2/
clean_v2/
```

旧布局可以原地迁移：

```powershell
python -m pgci_preprocess migrate-layout --dry-run
python -m pgci_preprocess migrate-layout
```

迁移会更新索引并清理经确认的旧兼容目录，不修改四件套内容。

## 安装

要求 Python 3.11 或更高版本。

```powershell
cd pipeline\code
python -m pip install -e ".[gbff]"
```

Core 预处理只使用 Python 标准库；处理 GBFF/GenBank 时需要 Biopython。

## 快速开始

### 1. 换机器后先扫描数据

`discover.py` 不需要先安装项目，使用任意 Python 3.11+ 即可运行。

数据位于项目 `data/` 下：

```powershell
python <skill>\pipeline\discover.py
```

数据位于其他位置：

```powershell
python <skill>\pipeline\discover.py E:\raw_genomes --project-root D:\plant_project
```

扫描会生成：

```text
resources/machine_paths.json       本机路径地图
resources/download_inventory.tsv   配对后的物种输入清单
resources/sequence_mode.tsv        逐物种染色体模式
resources/discovery_report.tsv     候选文件、角色和得分
```

扫描后，后续命令不再需要每次传 `--project-root`。

### 2. 检查环境与输入

```powershell
python -m pgci_preprocess doctor
```

检查项目输入与可用工具：

```powershell
python -m pgci_preprocess doctor --workflow full
```

### 3. 运行预处理

按 `sequence_mode.tsv` 为每个物种选择模式：

```powershell
python -m pgci_preprocess run --sequence-mode from-map
```

含义：

```text
有染色体证据的物种 → chromosomes
仅 scaffold 的碎片化组装 → adaptive
```

也可以统一强制严格染色体模式：

```powershell
python -m pgci_preprocess run --sequence-mode chromosomes
```

统一使用 adaptive：

```powershell
python -m pgci_preprocess run --sequence-mode adaptive
```

无染色体证据时，严格模式会报告：

```text
no_chromosome_sequences
```

经确认属于碎片化组装后，再使用 adaptive 模式。

### 4. 查看染色体数与基因数

```powershell
python -m pgci_preprocess chromosomes --species <code>
```

也可以直接查看：

```text
wgdi_input/<Abbr>/<Abbr>.lens
wgdi_input/<Abbr>/<Abbr>.chromosome_map.tsv
```

`lens` 固定三列：

```text
序列号    序列长度    代表基因数
```

因此：

```text
行数              = 最终保留的染色体/序列数
第三列             = 每条序列的代表基因数
第三列求和         = 进入 WGDI 的总基因数
```

示例：

```text
1	45753714	4210
2	39802311	3980
3	45120982	4521
```

表示保留 3 条序列，分别有 4210、3980、4521 个代表基因。

严格染色体模式下，这些序列为染色体级序列；adaptive 模式下应称为 WGDI 序列，可能包含 scaffold。

### 5. 人工确认染色体数

```powershell
python -m pgci_preprocess adjudicate `
  --species <code> `
  --confirmed-n <N> `
  --decision accept `
  --reviewer "<实际复核人>" `
  --rationale "<核型或文献依据>"
```

查看裁决状态：

```powershell
python -m pgci_preprocess adjudication-status --species <code>
```

裁决记录绑定 QC 版本和 SHA-256。重跑产生新 QC 后，旧签字自动变为 `stale`。

### 6. 校验

```powershell
python -m pgci_preprocess validate
python -m pgci_preprocess check-names
python -m pgci_preprocess audit
```

需要人工裁决门禁时：

```powershell
python -m pgci_preprocess validate --require-adjudication
```

## 单物种显式注册

没有自动扫描结果时，可以显式注册：

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

`--genome` 可以省略。显式注册要求 GFF/GFF3、CDS 和蛋白三个不同且非空的源文件。

GTF 或 GBFF 可以先转换：

```powershell
python -m pgci_preprocess convert-annotation `
  --input D:\raw\species.gtf.gz `
  --output-dir D:\plant_project\data\derived `
  --prefix species
```

批量处理也可以自动路由 GTF/GBFF，或在已有 genome+GFF 时提取缺失的 CDS 和蛋白。

## 输入格式

| 类型 | 支持 |
| --- | --- |
| GFF/GFF3 | standard、transcript-only、simplified5、gene-span、gene-CDS-direct |
| GTF | 转换为 GFF3 后进入流程 |
| GBFF/GenBank | 提取编码 GFF3、CDS、蛋白，可导出 genome |
| genome/CDS/蛋白 | FASTA 格式 |
| 压缩 | gzip、bzip2、xz、单成员 ZIP |
| 功能注释 | CSV/TSV 作为辅助元数据 |

文件后缀作为提示，格式和角色最终由内容决定。

## 预处理规则

### 代表转录本

每个基因按以下顺序选择一个代表转录本：

1. CDS 最长；
2. CDS 等长时按 transcript ID 字典序；
3. 再按 protein ID 稳定排序。

### 染色体模式

| 模式 | 行为 |
| --- | --- |
| `chromosomes` | 只保留有染色体证据的序列；无证据时失败关闭 |
| `adaptive` | 没有明确染色体证据时，按基因数和覆盖率选择序列 |
| `from-map` | 读取 `sequence_mode.tsv`，逐物种选择以上模式 |

### 序列与基因下限

| 参数 | CLI 默认 |
| --- | ---: |
| `--min-genes-per-sequence` | 50 |
| `--coverage` | 0.90 |
| `--max-sequences` | 200 |
| `--min-genes-per-chromosome` | 100 |

少于 100 个代表基因的序列会被移除，记录写入 `chr_map.tsv` 和 QC。

### 可选全量快照

默认不写筛选前的 `.all.*`。需要时：

```powershell
python -m pgci_preprocess run --keep-all
```

已有 `.all.*` 且筛选参数已经定稿时，可以安全瘦身：

```powershell
python -m pgci_preprocess prune --all-snapshots --dry-run
python -m pgci_preprocess prune --all-snapshots
```

`prune` 会同步清理索引中的失效路径。

## 关键输出

### 四件套

```text
<Abbr>.pep
<Abbr>.cds
<Abbr>.new.gff
<Abbr>.lens
```

### 追溯文件

| 文件 | 用途 |
| --- | --- |
| `<Abbr>.gene_map.tsv` | 新基因 ID 与原始基因、转录本、蛋白 ID |
| `<Abbr>.chr_map.tsv` | 原始序列、内部编号和保留/排除原因 |
| `<Abbr>.chromosome_map.tsv` | 最终染色体数、长度和每条基因数 |
| `<Abbr>.removed_list.tsv` | 全部剔除记录与原因 |
| `<Abbr>.qc.json` | 机器可读 QC、输入指纹、阈值和异常 |
| `<Abbr>.qc_record.md` | 人类可读 QC 摘要 |

### 最小交付

只交给 WGDI 使用时，至少需要：

```text
<Abbr>.pep
<Abbr>.cds
<Abbr>.new.gff
<Abbr>.lens
```

其中 `.lens` 同时提供最终染色体数和每条染色体的代表基因数。

## 命名

物种缩写由属名首字母和种加词字母组成，优先三位，冲突后扩展。基因 ID 示例：

```text
Mun1g0001
```

旧 `s<show_id>` 项目可以迁移：

```powershell
python -m pgci_preprocess migrate-artifacts --dry-run
python -m pgci_preprocess migrate-artifacts
```

## 常见问题

| 现象 | 处理 |
| --- | --- |
| 换电脑后找不到项目 | 先运行 `pipeline/discover.py` |
| 文件配对错误 | 查看 `resources/discovery_report.tsv` 和 `download_inventory.tsv` |
| 报 `no_chromosome_sequences` | 检查 GFF/genome 的染色体 ID；确认碎片化后再用 adaptive |
| 某条染色体消失 | 查看 `chr_map.tsv` 中的 `below_min_genes_per_chromosome` |
| 想解释基因数为什么减少 | 查看 `removed_list.tsv` 和 `qc_record.md` |
| 重跑后旧签字失效 | QC 已更新，需要重新执行 `adjudicate` |
| 旧项目存在 `analysis_v2`/`clean_v2` | 使用 `migrate-layout --dry-run` 预演迁移 |
| 需要全量候选快照 | 使用 `run --keep-all` |

## 下游使用

预处理完成后，四件套可直接交给 DIAMOND/BLAST 和 WGDI。后续点图、共线性及可选 Ks 命令见：

```text
pipeline/QUICKSTART.md
pipeline/code/docs/WORKFLOW.md
```

预处理结果就绪表示四件套和追溯文件完成；下游共线性和生物学结论由后续分析继续。

## 测试

```powershell
cd pipeline\code
python -m pytest -q
python -m ruff check src tests scripts tools
python -m mypy src --ignore-missing-imports
```

## 版本与仓库

- 当前代码版本：v3.8.0
- 本仓库包含 Skill 调用说明、完整 `pipeline/` 运行副本和最新预处理目录约定
- 程序实现以 `pipeline/code/` 为准

## 许可

MIT License。第三方工具按其上游许可使用，详见 [pipeline/THIRD_PARTY_NOTICES.md](pipeline/THIRD_PARTY_NOTICES.md)。
