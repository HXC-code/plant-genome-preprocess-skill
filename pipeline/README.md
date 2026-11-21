# HXC 植物基因组规范化与质控工具

本仓库将植物基因组的原始下载文件规范化为可审计的 WGDI 输入，并提供可选的
DIAMOND/BLAST、共线性、block-Ks 和 Ks-peak 执行链。

核心预处理输出为：

```text
<Abbr>.pep  <Abbr>.cds  <Abbr>.new.gff  <Abbr>.lens
```

同时生成来源映射、染色体映射、剔除记录、QC、物种索引、哈希清单和命名审计。
当前版本为 **3.5.5**。

- [发布包](https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.5.5)
- [完整快速上手](QUICKSTART.md)
- [验证与发布门禁](VALIDATION.md)
- [Agent Skill 仓库](https://github.com/HXC-code/plant-genome-preprocess-skill)

## 项目定位

本工具把人工整理和审核后的植物基因组预处理规则，固化为可执行、可复现、可审计的
标准流程，为 DIAMOND/BLAST 和 WGDI 准备一致的输入。

> **一句话定位：** 实验室植物基因组预处理的数据底座与 SOP 执行工具；它保证输入
> 干净、一致、可追溯，不产出 WGD 生物学结论。

| 层面 | 定位 |
| --- | --- |
| 实验室 | **标准预处理流程（SOP）的执行载体**：研究者维护规则，程序负责批量执行、QC 和留痕，人工负责例外审核 |
| 科学 | **下游分析的数据基础设施**：统一注释、序列、命名和筛选口径，不替代生物学判断 |
| Agent | `plant-genome-preprocess-skill` 是 Agent 调用与分发层，`pipeline/` 是程序仓的版本化镜像 |
| 软件类别 | 科学软件中的 **data preprocessing / QC / provenance tool**，不是新算法、WGDI 分析平台或结论生成器 |

若作为论文或社区工具推广，较合适的形态是 software/tool 类工作。当前仍包含实验室
目录和命名资源，并需要更广泛的物种与平台验证，才能进一步成为通用社区工具。

## 3.5.5 更新

- 新增植物线粒体核心基因识别：对合成 region 检测多个不同的
  `nad/cob/cox/atp/matR/mttB/ccm` 基因家族，不再套用只适合叶绿体的
  “小、密、无内含子”假设。
- 新参数 `--organelle-mt-min-genes` 控制自动排除所需的线粒体核心基因家族数。
- 只有部分线粒体证据的 region 进入 `review_candidates`，不会被静默排除。
- 对于基因 ID 完全匿名、没有参考相似度的大而稀疏线粒体 scaffold，工具仍不猜测；
  QC 的 `organelle_inference.limits` 会明确记录该边界。
- 新增 `pgci_preprocess readiness` 作为 readiness 的唯一权威读取入口；
  `validate`、`run_homology.sh` 和 Windows 总控共享同一结构解析逻辑。
- `6run_wgdi.sh` 新增 8 个 WGDI 产物的内容完整性检查：PNG 魔数、Ks 表头、
  共线性块和 block-info CSV 表头；截断但非空的文件不再视为完成。

## 3.5.4 修复的问题

| # | 原问题 | 3.5.4 修复结果 |
|---|---|---|
| 1 | 无 genome、无显式 region 的简化注释无法识别细胞器，叶绿体/线粒体序列进入核基因组 WGDI 集 | 对合成 region 使用标识符信号和紧凑基因集推断；证据写入 QC `organelle_inference`，模糊候选写入 `review_candidates`；支持 `--organelle-pattern` 和 `--no-organelle-inference` |
| 2 | `wgdi_ready=no` 仍可通过 `validate`，下游仅检查四件套是否存在 | `validate` 默认报 `wgdi_not_ready_low_coverage`；归档检查可显式使用 `--allow-not-ready`；`run_homology.sh` 和 `run_all_pipeline.ps1` 跳过未就绪或无准备度记录的物种 |
| 3 | `init-project` 强制要求 genome 和 assembly ID，三件套物种无法正常注册 | GFF + CDS + PEP 三件套可直接注册，`--genome` 与 `--assembly-id` 均可省略 |
| 4 | `run_homology.sh` 将 `bin/`、`code/`、`resources/` 等基础设施目录误认为物种目录 | 下游扫描明确排除基础设施目录，`skipped_species.txt` 不再产生虚假的 `bin` 记录 |
| 5 | 生成的 `run.sh` 使用 Windows 反斜杠路径，Git Bash/Linux 无法正确执行 | WGDI 入口统一转换为 Bash 路径，例如 `D:\path` 转换为 `/d/path` |
| 6 | 临时文件清理失败会使成功的比对任务返回退出码 1 | 清理 trap 容错；中间结果使用按进程区分的 `.part.<pid>` 文件，成功结果不再被清理错误污染 |
| 7 | Windows `TMPDIR` 和 WGDI 环境错误只能在正式运行时暴露 | `bin/wgdi` 规范化 Windows 临时路径并在不可写时回退；`doctor --workflow wgdi|full` 新增实际启动入口的 `wgdi_probe` 检查 |

## 项目可以完成的工作

| 工作 | 具体能力 |
| --- | --- |
| 单物种规范化 | 将一个物种的注释与 genome/CDS/蛋白 FASTA 转换为 WGDI 四件套 |
| 批量物种清洗 | 从下载目录或清单发现物种，解析文件角色，逐物种运行、校验并汇总结果 |
| 注释转换 | 将 GTF、GBFF/GenBank 转为标准 GFF3，并从 GBFF 提取 CDS、蛋白和 genome |
| 基因模型整理 | 修复缺失 region、缺失 mRNA、gene/mRNA/CDS 层级不一致及拆分基因问题 |
| 代表转录本选择 | 每个基因保留一个代表转录本，形成严格的一基因一转录本数据集 |
| 主序列选择 | 优先保留主染色体；无主染色体标记时按基因数、累计覆盖率和序列数上限选取 |
| 序列一致性检查 | 检查 CDS 长度、翻译结果、内部终止、非法字符、ID 重复和蛋白一致性 |
| 稳定命名 | 生成全局唯一的物种缩写、序列编号、基因 ID 和旧 ID 来源映射 |
| 旧数据迁移 | 将 `s<show_id>` 目录、文件名和基因 ID 迁移为物种缩写命名 |
| 数据验收 | 校验 index、QC、ID 集合、lens、SHA-256、兼容链接和命名一致性 |
| 来源追踪 | 记录每个代表基因对应的原始文件、原始基因 ID、原始转录本和蛋白 ID |
| 冻结清单 | 生成并验证 role-level source manifest 及其 SHA-256 sidecar |
| 可选比较分析 | 执行 DIAMOND/BLAST、WGDI dotplot、共线性、Ks、blockinfo、correspondence、block-Ks 和 Ks 峰分析 |

## 两种运行层级

```text
原始注释 + genome/CDS/蛋白
          │
          ▼
  prepare / init-project
          │
          ▼
        run
          │
          ├── analysis_v2/<Abbr>/   一基因一代表转录本 + WGDI 主序列集
          └── clean_v2/<Abbr>/      映射、剔除记录和清洗审计
          │
          ▼
 validate + check-names + audit

          │ 可选下游
          ▼
 DIAMOND 或 BLAST
          │
          ▼
 WGDI dotplot / collinearity / Ks / block-Ks / peaks
```

`pgci_preprocess` 负责预处理、检查、命名、迁移和审计。仓库根目录的
`run_homology.sh`、`2run_pre_wgdi.sh`、`3blast_to_wgdi.sh`、
`6run_wgdi.sh` 与 `run_all_pipeline.ps1` 组成可选比较层。

## 可直接处理的项目形态

| 已有数据 | 处理结果 |
| --- | --- |
| GFF/GFF3 + genome + CDS + 蛋白 FASTA | 完整清洗、序列核对、代表转录本选择和 WGDI 输出 |
| GTF + genome + CDS/蛋白 | GTF 先转 GFF3，再进入标准流程 |
| GBFF/GenBank | 提取编码 GFF3、CDS、蛋白，并可导出 genome 后继续处理 |
| GFF + CDS + 蛋白，无 genome | 使用 CDS 计算长度；region 缺失时由转录本坐标回填 |
| 只有 genome、无注释 | 报告 `annotation_missing`，不生成虚构注释 |
| 下载目录中有多套候选文件 | 在清单中用 `A|B` 指定候选，程序按角色兼容性打分并处理平局歧义 |
| 已有 `s<show_id>` 旧产物 | 提供 dry-run、迁移、兼容链接刷新和旧目录隔离 |
| 已有完整项目产物 | 使用 `validate`、`check-names`、`audit` 和 `source-manifest` 进行验收与冻结 |

## 输入格式

文件后缀只作为提示。注释结构、FASTA 角色和压缩类型均结合内容识别；内容无法
唯一判定时，程序报歧义，不做静默猜测。

### 注释

| 输入 | 处理 |
| --- | --- |
| `.gff`、`.gff3` | 识别 standard、transcript-only、simplified 5 列、gene-span、gene-CDS-direct 五种结构 |
| `.gtf` | 转换 gene、transcript、exon、CDS、start_codon、stop_codon 后输出标准 GFF3 |
| `.gbff`、`.gbk`、`.gb`、`.genbank` | 提取编码区 GFF3、CDS 和蛋白；可选导出 genome，需要 `biopython>=1.83` |
| `.csv`、`.csv.gz`、TSV 功能注释表 | 规范化为功能元数据，不承担坐标注释或 WGDI 输入角色 |

### 序列

支持的 FASTA 后缀包括 `.fa`、`.fasta`、`.fna`、`.faa`、`.fas`，以及
`.genome.fasta.gz`、`_genomic.fna.gz`、`.CDS.fasta.gz`、
`.Protein.faa.gz`、`.pep.fa`、`.pro.fa`、`.RNA.fasta.gz` 等组合名称。

genome、CDS、蛋白三个角色由文件名语义和 FASTA 首条记录共同判定。具有 FASTA
内容的其它扩展名也可读取。

### 压缩

`.gz`、`.bz2`、`.xz`、`.zip` 由压缩魔数识别，不依赖后缀是否正确。

## 项目目录

```text
my_project/
├── data/
│   └── downloads/
│       └── {show_id}_{Order}_{Family}({Species})/
├── resources/
│   ├── download_inventory.tsv
│   ├── species_abbreviations.tsv
│   ├── species_registry_v2.tsv
│   ├── wgdi_v2_index.tsv
│   ├── wgdi_v2_manifest_sha256.tsv
│   ├── naming_audit.tsv
│   └── preprocess_batches/
├── analysis_v2/<Abbr>/
├── clean_v2/<Abbr>/
├── <Abbr>/                       # 可选兼容链接目录
└── data-validation-report.tsv
```

原始数据和生成结果位于项目目录，代码与文档位于本仓库。`init-project` 会创建
标准子目录、清单表头和资源文件；单物种注册必须给出显式的 `--gff`、`--cds`、
`--pep` 路径，`--genome` 与 `--assembly-id` 可选（只有三件套也能注册）。

## 安装

要求 Python 3.11 或更新版本。核心预处理以标准库为主；GBFF/GenBank 转换需要
Biopython extra。

```powershell
cd code
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[gbff]"
```

Linux/WSL：

```bash
cd code
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[gbff]'
```

可选比较层还需要 Bash、DIAMOND 或 BLAST+、WGDI 0.75、MUSCLE、PAML `yn00`
和 Perl。`doctor` 对 Python 依赖、项目输入和这些外部程序做只读检查；
`--workflow wgdi|full` 还会实际启动 WGDI 入口（`wgdi_probe`），暴露解释器
不兼容等“存在但跑不起来”的问题。

## 单物种工作流

```powershell
cd code

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

python -m pgci_preprocess run --project-root D:\plant_project
python -m pgci_preprocess validate --project-root D:\plant_project --refresh-links
python -m pgci_preprocess check-names --project-root D:\plant_project
python -m pgci_preprocess audit --project-root D:\plant_project
```

## 批量工作流

`prepare` 从 status 表和项目清单发现完整物种、校验候选文件、镜像输入并刷新
inventory。`batch` 将 `prepare → 可选 normalize-gff → run → validate` 串成一条
命令。

```powershell
cd code

python -m pgci_preprocess prepare `
  --project-root D:\plant_project `
  --source-status D:\plant_project\resources\source_status.tsv `
  --expected-count 20

python -m pgci_preprocess batch `
  --project-root D:\plant_project `
  --source-status D:\plant_project\resources\source_status.tsv `
  --expected-count 20
```

批处理的文件路由决策写入
`resources/preprocess_batches/<batch_tag>_routing.tsv`。

## CLI 命令

| 命令 | 用途 |
| --- | --- |
| `init-project` | 创建项目目录和标准表头，可注册一个显式物种 bundle |
| `prepare` | 发现完整物种、识别输入角色、镜像文件并刷新 inventory |
| `doctor` | 检查 Python、项目输入、清单以及 DIAMOND/BLAST/WGDI 等外部程序 |
| `normalize-gff` | 将 transcript-only GFF 展开为 gene/mRNA/CDS 层级 |
| `convert-annotation` | 单独转换 GTF、GBFF 或功能 GO 表，不修改原始输入 |
| `run` | 执行一基因一代表转录本预处理和主序列选择 |
| `validate` | 校验 index、QC、ID 集合、lens、SHA-256 和兼容链接；`wgdi_ready=no` 默认判 FAIL |
| `readiness` | 从 QC/index 结构化读取物种准备度，作为各入口共用的权威门禁 |
| `check-wgdi-outputs` | 检查一个比较目录中 8 个 WGDI 产物的结构和截断 |
| `batch` | 串联 prepare、可选 normalize 和 run、validate |
| `audit` | 汇总物种、输入后缀、`wgdi_ready` 分布和失败项目 |
| `source-manifest build` | 从 index、QC 和修复表生成冻结的来源清单 |
| `source-manifest verify` | 验证清单哈希及其引用的全部源文件 |
| `migrate-species-names` | 将 index/QC 中的旧 `s<show_id>` 元数据改为缩写 |
| `migrate-artifacts` | 重命名旧 bundle、文件名和基因 ID |
| `check-names` | 核对缩写表、清单、下载目录、index 和产物 bundle |
| `convert` | 保留的旧目录扫描兼容入口 |

完整参数和平台示例见 [QUICKSTART.md](QUICKSTART.md)。

## 清洗与规范化规则

一次 `run` 按以下顺序处理每个物种：

1. 结合后缀、内容预览和结构检测识别注释格式及 FASTA 角色。
2. 解析 GFF/GTF/GBFF，补 region，补 mRNA，并修复 gene/mRNA/CDS 层级。
3. 按 `protein_id` / transcript 关系匹配 CDS 与蛋白序列。
4. 检查 CDS 长度是否为 3 的倍数、翻译是否一致、是否存在内部终止或非法字符。
5. 检查坐标、链方向、重复 ID、同一蛋白对应多基因以及拆分基因同名问题。
6. 排除细胞器序列进入核基因组 WGDI 集，并保留排除原因。有 genome 或显式
   region 时按其判定；无 genome 的简化注释对合成 region 做细胞器推断
   （标识符信号 + 高密度紧凑基因集），可用 `--organelle-pattern` 显式指定、
   `--no-organelle-inference` 关闭；证据写入 QC 的 `organelle_inference`。
7. 按蛋白长度、CDS 长度和稳定 ID 为一个基因选择一个代表转录本。
8. 选择核基因组主序列，生成新的序列号、基因 ID 和原始 ID 映射。
9. 写出 `.pep`、`.cds`、`.new.gff`、`.lens` 及全部审计文件。

被剔除的转录本写入 `<Abbr>.rejected.tsv`，字段包括 `scope`、`id`、`reason`
和 `detail`。

## 代表转录本和主序列选择

每个基因的候选转录本按以下顺序选择：

1. 最长蛋白；
2. 蛋白等长时选择最长 CDS；
3. 仍相同时按 transcript ID 和 protein ID 的稳定顺序确定。

主序列选择由下列参数控制：

| 参数 | 默认值 | 作用 |
| --- | ---: | --- |
| `--min-genes-per-sequence` | 50 | 无主染色体标记时保留基因数足够多的序列 |
| `--coverage` | 0.90 | 兜底选择达到累计基因覆盖率所需的序列 |
| `--max-sequences` | 200 | 限制进入 WGDI 集的序列数 |

选择顺序为：

1. 存在 primary chromosome 标记时，使用全部非细胞器主序列；
2. 否则保留基因数不少于阈值的序列；
3. 碎片组装或缺少主序列标记时，按基因数从高到低累计选取，直到达到目标覆盖率；
4. 输出 `wgdi_ready=yes|marginal|no` 作为输入准备度标记。

`wgdi_ready` 的判定阈值：

- `yes`：保留率至少 70%，且序列数不超过 `--max-sequences`；
- `marginal`：保留率至少 50%；
- `no`：保留率低于 50%。

## 命名与 ID

物种缩写由属名首字母和种加词字母组成。默认生成三位缩写，冲突时依次使用种加词
中的其他两位字母组合，三位组合耗尽后扩展为四位及更长。随机数字后缀只作为最后
兜底。仓库携带的 `species_abbreviations.tsv` 保留已发布物种的既有缩写。

代表基因 ID 格式：

```text
<Abbr><sequence_number>g<four_or_more_digit_order>
```

例如：

```text
Mun1g0001
```

`.pep` 和 `.cds` 头行同时记录原始基因 ID：

```text
>Mun1g0001 source_gene=<original_gene_id>
```

旧项目可先 dry-run，再执行完整迁移：

```powershell
python -m pgci_preprocess migrate-artifacts --project-root D:\plant_project --dry-run
python -m pgci_preprocess migrate-artifacts --project-root D:\plant_project --prune-legacy-dirs
```

## 输出文件

### WGDI 输入集

`analysis_v2/<Abbr>/` 是可交给 WGDI 的规范输出目录：

| 文件 | 内容 |
| --- | --- |
| `<Abbr>.pep` | 主序列集、一基因一代表转录本的蛋白 FASTA |
| `<Abbr>.cds` | 与蛋白对应的编码序列 FASTA |
| `<Abbr>.new.gff` | 七列简化 GFF，包含新序列号、唯一基因 ID、坐标、链、序号和原基因 ID |
| `<Abbr>.lens` | 所选序列长度表 |
| `<Abbr>.all.pep.fa` | 全部有效代表基因的蛋白 FASTA |
| `<Abbr>.all.cds.fa` | 全部有效代表基因的 CDS FASTA |
| `<Abbr>.all.gff` | 全部有效代表基因的简化 GFF |
| `<Abbr>.all.lens` | 全部有效代表基因涉及的序列长度表 |
| `<Abbr>.gene_map.tsv` | 新基因 ID 与原始基因、转录本、蛋白 ID 的映射 |
| `<Abbr>.chr_map.tsv` | 新序列名与原始序列名的映射 |
| `<Abbr>.qc.json` | 输入哈希、格式检测、计数、阈值、异常和版本信息 |

`.new.gff` 的列定义为：

```text
简化序列号  唯一基因ID  起点  终点  链方向  序列内序号  原始基因ID
```

### 清洗与审计

`clean_v2/<Abbr>/` 保存 gene map、chr map 和 rejected 记录。

| 文件 | 内容 |
| --- | --- |
| `<Abbr>.rejected.tsv` | 被剔除记录及机器可读原因 |
| `<Abbr>.gene_map.tsv` | 完整基因来源映射 |
| `<Abbr>.chr_map.tsv` | 完整序列映射 |
| `resources/wgdi_v2_index.tsv` | 进入分析集的物种和产物索引 |
| `resources/species_registry_v2.tsv` | 物种元数据注册表 |
| `resources/wgdi_v2_manifest_sha256.tsv` | index 和产物的 SHA-256 清单 |
| `resources/naming_audit.tsv` | 命名五方一致性审计 |
| `data-validation-report.tsv` | 逐物种验证结果 |

## 校验和审计

`validate` 检查：

- 四件套是否存在且基因数一致；
- GFF 坐标、链方向和序列长度是否合法；
- `.lens` 与 GFF 序列集合是否一致；
- QC 与产物中的基因 ID 集合是否一致；
- 文件 SHA-256 是否与 index/manifest 一致；
- 可选兼容链接是否指向同一内容；
- QC 的 `wgdi_ready=no` 默认判 FAIL（`wgdi_not_ready_low_coverage`），
  防止碎片化物种被“文件齐全”误送进 WGDI；仅做数据完整性归档时显式加
  `--allow-not-ready`。

`check-names` 核对：

```text
物种缩写表 ↔ 下载 inventory ↔ 物种目录 ↔ WGDI index ↔ 产物 bundle
```

`audit` 提供只读的项目级概览，包括物种数、文件后缀、准备度分布和失败原因。
`doctor` 提供依赖与输入就绪检查，不安装外部工具，也不执行 WGDI 分析
（`wgdi_probe` 只运行入口的版本探针）。

冻结项目数据时，`source-manifest build` 生成 role-level 来源清单和 sidecar；
`source-manifest verify` 可重新验证清单及其引用的源文件。

## 可选 DIAMOND/BLAST + WGDI

比较层从预处理后的物种 bundle 开始，默认使用：

| 参数 | 默认值 |
| --- | --- |
| 同源引擎 | DIAMOND |
| 灵敏度 | `very-sensitive` |
| E-value | `1e-5` |
| tandem blocks | 从 correspondence、block-Ks 和 Ks-peak 中移除 |

`run_homology.sh` 与 `run_all_pipeline.ps1` 会读取 QC/index 的 `wgdi_ready`，
跳过 `ready=no` 和没有可读准备度记录的物种（写入 `skipped_species.txt`），
并且不把 `bin/`、`code/`、`resources/` 等基础设施目录当作物种扫描。
确需复核未就绪物种时显式设置 `HOMOLOGY_ALLOW_NOT_READY=true`。

Windows 总控执行步骤 1–6：

```powershell
bash 0check_data.sh --preprocess --refresh-links
powershell -ExecutionPolicy Bypass -File .\run_all_pipeline.ps1
```

可切换 BLAST 或指定外部工具路径：

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all_pipeline.ps1 `
  -HomologyEngine blast `
  -HomologyEvalue 1e-5 `
  -BioinfoBin D:\bioinfo\bin
```

逐步执行：

```bash
bash 0check_data.sh --preprocess --refresh-links
bash run_homology.sh
bash 2run_pre_wgdi.sh
bash 3blast_to_wgdi.sh
python code/tools/generate_wgdi_configs.py --mode all
bash 6run_wgdi.sh
```

`generate_wgdi_configs.py --mode all` 同时生成同物种和种间比较配置。WGDI 阶段
可输出 dotplot、`.icl` 共线性、Ks、blockinfo、correspondence、block-Ks、
Ks peaks 和 peaks fitting。

WGDI 环境可使用：

```powershell
powershell -ExecutionPolicy Bypass -File .\bootstrap_wgdi.ps1
```

或：

```bash
bash bootstrap_wgdi.sh
```

PAL2NAL v14 从上游安装并进行 SHA-256 校验：

```powershell
python code/tools/install_pal2nal.py
```

发布包不分发 PAL2NAL。外部工具来源与许可见
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。

## 验证记录

3.5.5 验证覆盖：

- Python 自动化测试 **188** 项通过；
- 植物线粒体核心基因识别、readiness 单一权威源和 WGDI 产物完整性检查的测试通过；
- Windows 合成数据完整 WGDI 链通过，四组比较每组检查 8 个关键产物
  （`validation/smoke-report-3.5.5.json`）；
- 运行环境：Python 3.13.9、pandas 2.3.3、WGDI 0.75；
- Ruff、Mypy 通过。

3.5.4 验证覆盖：

- Python 自动化测试 **172** 项通过；
- Ruff、Mypy、Git Bash 语法检查和严格 LF 的 Linux/WSL `bash -n` 通过；
- Windows 合成数据完整 WGDI 链通过，四组比较每组检查 8 个关键产物
  （`validation/smoke-report-3.5.4.json`）；
- 真实数据抽查（41 物种试点审查）：无 genome 三件套注册成功；橡胶树
  `Hbr` 叶绿体 `NC_015308.1`（83 基因）与银白杨 `Pab` 叶绿体
  `NC_008235.1`（83 基因）被推断排除，不再进入 WGDI 集；
- 试点项目严格验收 36 PASS / 5 FAIL（5 个 `wgdi_ready=no` 碎片化物种被
  门禁拦截），`--allow-not-ready` 归档模式 41 PASS；
- `doctor` 的 `wgdi_probe` 在专用环境（Python 3.13 + wgdi 0.75）实际启动
  入口并通过。

3.5.3 及更早验证覆盖：

- Windows 与 Linux/WSL 的合成数据完整 WGDI 链均通过，四组比较每组检查 8 个关键产物；
- NCBI RefSeq 真实物种 `Arabidopsis thaliana` TAIR10.1 和
  `Arabidopsis lyrata` v1.0 完成清洗、校验、命名和 QC；
- 真实物种子集 `Ath` chromosome 1 对 `Aly` chromosome 364 完成双方向 WGDI
  比较，每个方向 8 个关键产物通过。

回执位于：

```text
validation/smoke-report-3.5.5.json
validation/smoke-report-3.5.4.json
validation/smoke-report-3.5.3.json
validation/smoke-report-3.5.2-linux.json
validation/real-species-cross-3.5.2.json
validation/real-species-3.5.2.json
```

发布门禁和运行环境记录见 [VALIDATION.md](VALIDATION.md)。

## 当前不足与适用边界

以下限制是当前版本的已知边界。它们不否定已经通过的数据处理和执行验证，但意味着
`PASS`、`ready=yes` 或四件套齐全不能直接解释为组装完整或 WGD 结论可信。

| 类别 | 当前不足 |
| --- | --- |
| 细胞器识别 | 叶绿体、明确命名或含多个线粒体核心基因的合成 region 可推断；基因 ID 完全匿名、缺少参考相似度的大而稀疏线粒体仍可能漏判。此类序列需要 `--organelle-pattern` 或人工复核 |
| 代表转录本 | 当前按最长蛋白、最长 CDS 和稳定 ID 选择，尚未优先使用注释自带的 primary/canonical 标记，也没有表达或 NMD 证据 |
| 多倍体与杂合 | 尚未建模等位 contig、homeolog 和亚基因组；多倍体或高杂合数据可能被通用去重和拒绝规则错误压缩 |
| readiness 与组装质量 | `wgdi_ready` 只衡量当前筛选集的基因保留率，不包含 BUSCO 或独立组装完整性证据；`--max-sequences` 还可能限制碎片化组装的序列覆盖 |
| 阈值可配置性 | 部分 readiness 和细胞器阈值仍是内部默认值，缺少统一的配置来源、版本和依据说明 |
| 续跑与版本迁移 | 3.5.5 已校验 8 个 WGDI 产物的常见截断/结构错误；仍无阶段级哈希链，QC 版本变化也可能触发大范围重算而非增量迁移 |
| 规模与超时 | 批处理和 WGDI 阶段以串行为主；真实全蛋白组 Ks 可能超过当前超时，尚无集群编排和细粒度恢复 |
| 输入覆盖 | 工具从注释开始，不负责对仅有 assembly 的物种运行 Liftoff/BRAKER/Helixer；41/77 的试点覆盖只代表有完整注释的部分 |
| 验证深度 | 自动化测试以合成数据和有限真实物种为主；多倍体、高杂合、复杂线粒体等硬案例仍缺少带 ground truth 的生物学验收 |
| 灰度状态 | 当前使用 rejected、`review_candidates` 和 readiness 等级，但尚未提供经过校准的置信度或“存疑但纳入”证据层 |
| 科学解释 | 工具止于数据准备、执行和证据记录；Ks 峰解释、WGD 定年、速率异质性、饱和校正和 NG86/YN00 取舍仍由专家判断 |

合成 smoke 验证执行链和工具兼容性；真实物种子集验证预处理和跨物种执行。当前真实
WGDI 验证为 chromosome-level 子集，完整全蛋白组生物学 WGD 结论需要在目标物种全集上
单独运行并由专家解释。DIAMOND/BLAST、WGDI、MUSCLE、PAML、Perl 和 PAL2NAL 均作为
独立外部工具安装，其许可独立适用；项目源代码使用 [MIT License](LICENSE)。

## Skill 包

Agent Skill 单独发布于：

<https://github.com/HXC-code/plant-genome-preprocess-skill>

Skill 仓库包含 `SKILL.md`、完整 `pipeline/` 镜像、CI 和验证回执。安装与调用方式
见该仓库 README。

## 测试

```powershell
cd code
python -m pytest -q
python -m ruff check src tests scripts tools
python -m mypy src --ignore-missing-imports
```
