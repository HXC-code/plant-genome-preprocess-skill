---
name: plant-genome-preprocess
description: 植物基因组预处理与可选 WGDI 比较工作流。把 GFF/GFF3/GTF/GBFF/功能注释 CSV 与 genome/CDS/蛋白 FASTA 清洗成 WGDI 四件套，完成命名校验、代表转录本选择、QC、细胞器审核和结果留痕；也可继续执行 DIAMOND/BLAST、共线性、block-Ks 和 Ks-peak 分析。
---

# 植物基因组预处理

`plant-genome-preprocess` 是植物基因组规范化与质控工具的 Agent 调用入口。`pipeline/` 提供与程序仓同步的完整运行副本，包含源码、测试、资源表、下游脚本和 WGDI 启动器。

预处理核心结束于 `pgci_preprocess run`，生成 WGDI 输入四件套。DIAMOND/BLAST、共线性和 Ks 分析作为可选下游模式执行。

## 适用任务

- 将植物注释与序列整理为 WGDI 的 `pep/cds/new.gff/lens`；
- 复核已有 WGDI 输入的命名、长度、坐标、重复基因和链接；
- 批量处理几十到几百个物种，并逐物种记录 QC 和失败原因；
- 将 `s<show_id>` 旧命名迁移为物种缩写；
- 继续执行 DIAMOND/BLAST、WGDI 共线性和 Ks 分析。

## 运行目录

```text
pipeline/
├── QUICKSTART.md               教师使用说明
├── README.md                   仓库说明
├── run_homology.sh             DIAMOND/BLAST 共用实现
├── 0check_data.sh … 6run_wgdi.sh
├── run_all_pipeline.ps1        Windows 分析总控
├── bootstrap_wgdi.ps1 / .sh    WGDI 环境准备
├── bin/wgdi                    可移植启动器
├── requirements-wgdi-validated.txt
├── VALIDATION.md               发布门禁
├── validation/                 验证回执
└── code/
    ├── src/pgci_preprocess/
    ├── scripts/ tools/ tests/ docs/
    ├── resources/species_abbreviations.tsv
    └── pyproject.toml
```

## 物种缩写

- 默认编码为属名首字母加种加词前两个字母，例如 `Cercis chinensis` → `Cch`；
- 编码冲突时使用种加词内的其他字母组合，例如 `Cci`、`Ccn`；
- 三位编码全部占用后扩展为四位或更长；
- `resources/species_abbreviations.tsv` 保存已有物种编码；
- 新物种按拉丁名生成编码，并避开表中全部已占用编码；
- 已注册物种编码保持稳定，历史目录和 gene ID 继续沿用。

## 标准工作流

### 1. 确认本地环境

```powershell
cd pipeline\code
python -m pgci_preprocess --version
```

项目创建后检查依赖和输入：

```powershell
python -m pgci_preprocess doctor `
  --project-root D:\my_project `
  --workflow preprocess
```

### 2. 创建项目

```powershell
python -m pgci_preprocess init-project --project-root D:\my_project
```

程序创建 `data/downloads`、`resources`、库存表、物种缩写表和格式例外表。

### 3. 注册物种

显式 GFF/GFF3、CDS 和蛋白文件可直接注册：

```powershell
python -m pgci_preprocess init-project `
  --project-root D:\my_project `
  --show-id 2848 `
  --scientific-name "Morus notabilis" `
  --order Rosales `
  --family Moraceae `
  --genome D:\raw\genome.fa.gz `
  --gff D:\raw\annotation.gff3.gz `
  --cds D:\raw\cds.fa.gz `
  --pep D:\raw\pep.fa.gz
```

`--genome` 可省略。`init-project` 将三个或四个显式文件镜像到项目目录。

GTF 和 GBFF 数据通过 `convert-annotation` 转换，或通过批量 `prepare` / `batch` 自动路由。

### 4. 预处理

```powershell
python -m pgci_preprocess run --project-root D:\my_project
```

已有且版本一致的有效结果默认保留并标记为 `SKIPPED`。指定物种重算：

```powershell
python -m pgci_preprocess run `
  --project-root D:\my_project `
  --species 2848 `
  --overwrite
```

升级版本或修改筛选参数时，先备份结果，再对相关物种使用 `--overwrite`。

### 5. 验收

```powershell
python -m pgci_preprocess validate --project-root D:\my_project --refresh-links
python -m pgci_preprocess check-names --project-root D:\my_project
python -m pgci_preprocess audit --project-root D:\my_project
```

验收通过后，`analysis_v2/<Abbr>/` 下的四件套作为 WGDI 输入。

### 6. 染色体清单与人工裁决（3.6.0，必做）

最终染色体数由人工判断，这是硬门禁，不是可选步骤：

```powershell
# 1) run 输出已直接打印清单；如需单独重看：
python -m pgci_preprocess chromosomes --project-root D:\my_project --species <code>
# 2) 老师复核后签字（confirmed_n 须等于程序染色体数；记录绑定 QC 版本+哈希）
python -m pgci_preprocess adjudicate --project-root D:\my_project --species <code> `
  --confirmed-n <N> --decision accept --reviewer "<姓名>" --rationale "<依据>"
# 3) 带门禁验收
python -m pgci_preprocess validate --project-root D:\my_project --require-adjudication
```

`run` 默认 `--min-genes-per-chromosome 100`，丢弃基因数不足 100 的染色体（记入 chr_map/QC）；该阈值按基因数机械过滤，也可能移除基因数少的真实染色体，需在 chr_map 中核对。重跑后 QC 哈希变化，旧签字转为 `stale`，需重新裁决。

两个门禁开关默认关闭，**必须由本 skill 显式传入**：验收加 `--require-adjudication`，下游设 `HOMOLOGY_REQUIRE_ADJUDICATION=true`。不传时程序不检查裁决记录，因此不得省略。`reviewer` 填实际复核人，不得代签；仓库与发布包不包含任何裁决记录。

## 输入识别

注释格式：

| 格式 | 处理 |
| --- | --- |
| `.gff`、`.gff3` | 支持 standard、transcript_only、simplified5、gene_span 和 gene_cds_direct 五种结构 |
| `.gtf` | 转换 gene、transcript、exon、CDS、start_codon、stop_codon 后输出标准 GFF3 |
| `.gbff`、`.gbk`、`.gb`、`.genbank` | 提取编码区 GFF3、CDS、蛋白，可导出 genome |
| `.annotations.csv.gz` | 识别 GO/NR/SP/TR 等功能注释并规范化 |
| `.csv`、`.csv.gz` | 功能注释规范化为辅助元数据 |

序列格式：

```text
.fa .fasta .fna .faa .fas
.genome.fasta.gz .genome.fa.gz _genomic.fna.gz
.CDS.fasta.gz .cds.fa.gz .Protein.faa.gz .pep.fa .pro.fa
.RNA.fasta.gz
```

程序结合文件名、FASTA 首条记录和内容特征识别 genome、CDS、蛋白和 RNA。压缩格式支持 gzip、bzip2、xz 与单成员 ZIP。

同一角色存在多个候选时，清单使用 `A|B`。程序按角色兼容性选择唯一文件，并将选择结果写入 routing 表。

## 清洗内容

一次 `run` 完成以下步骤：

1. **解析注释**：识别注释结构，补充 region 和 mRNA，修复 gene/mRNA/CDS 层级；
2. **识别细胞器**：根据标识符、基因家族和区域特征识别叶绿体与线粒体，并写入 QC；
3. **检查基因模型**：剔除 CDS 长度异常、内部终止、翻译不一致、坐标越界、链方向错误和缺少必要 ID 的候选；
4. **选择代表转录本**：按最长蛋白、最长 CDS 和稳定 ID 的顺序选择；
5. **选择主序列**：根据主染色体、基因数、覆盖率和序列上限生成 WGDI 序列集；
6. **统一命名**：生成物种缩写、序列编号、基因 ID 和来源映射；
7. **生成审计记录**：输出 QC、gene map、chr map、rejected、index、manifest 和 naming audit。

## 关键结果

```text
analysis_v2/<Abbr>/
├── <Abbr>.pep
├── <Abbr>.cds
├── <Abbr>.new.gff
├── <Abbr>.lens
├── <Abbr>.all.pep.fa
├── <Abbr>.all.cds.fa
├── <Abbr>.all.gff
├── <Abbr>.all.lens
├── <Abbr>.gene_map.tsv
├── <Abbr>.chr_map.tsv
└── <Abbr>.qc.json
```

```text
clean_v2/<Abbr>/
├── <Abbr>.gene_map.tsv
├── <Abbr>.chr_map.tsv
└── <Abbr>.rejected.tsv
```

`.new.gff` 使用七列：

```text
序列编号  唯一基因ID  起点  终点  链方向  序列内序号  原始基因ID
```

## 批量模式

状态表必需列：

```text
show_id Order Family Species status folder
current_genome current_gff_or_gtf current_cds current_pep
```

执行：

```powershell
python -m pgci_preprocess batch `
  --project-root D:\my_project `
  --source-status D:\my_project\resources\source_status.tsv `
  --expected-count 20
```

`. batch` 的顺序为：

```text
prepare → 可选 normalize-gff → run → validate
```

路由结果写入 `resources/preprocess_batches/<batch_tag>_routing.tsv`。已有且版本一致的有效物种默认保留，便于分次完成整批任务。

## 检查命令

| 命令 | 检查内容 |
| --- | --- |
| `validate` | 四件套、基因数、坐标、链、长度、ID 集合、SHA-256、兼容链接和 readiness |
| `check-names` | 缩写表、库存表、索引、下载目录和产物目录 |
| `readiness` | 从 QC/index 读取物种 `wgdi_ready` |
| `chromosomes` | 物种染色体清单：染色体数、每条基因数/长度/来源 |
| `adjudicate` | 记录人工对染色体数的签字（绑定 QC 版本+哈希） |
| `adjudication-status` | 逐物种裁决门禁状态（ok/missing/stale/count_mismatch/rejected/pending_adjustment） |
| `audit` | 物种数、后缀、准备度分布和失败原因 |
| `doctor` | Python、输入文件、清单、缩写表和外部工具 |
| `check-wgdi-outputs` | 8 个 WGDI 产物的 PNG、Ks、共线性和 CSV 结构 |

## 特殊情况

| 情况 | 处理 |
| --- | --- |
| 只有 genome | `run` 报告 `annotation_missing` |
| 只有注释、CDS 和蛋白 | 使用 CDS 和转录本坐标完成长度与 region 处理 |
| 简化 GFF 混入细胞器序列 | 记录识别证据；有一定证据的候选写入 `review_candidates` |
| 匿名且稀疏的线粒体序列 | 使用 `--organelle-pattern` 指定，或根据参考信息人工复核 |
| 碎片化组装 | `wgdi_ready` 反映当前筛选集保留率；结合 BUSCO 和组装质量判断 |
| 同一目录存在 GCF 与 GWH 多套文件 | 清单使用 `A|B` 或写入唯一文件名 |
| 文件名缺少 GCA/GCF 编号 | 使用 `assembly_accessions.tsv` 记录装配来源 |
| 旧 `s<show_id>` 产物 | 使用 `migrate-artifacts` 迁移目录、文件名和基因 ID |

## 执行与汇报规则

处理任务时按以下顺序继续：

1. **先检查环境**：运行 `doctor`，确认 Python、输入和所需外部工具；
2. **报告输入路由**：列出所选文件、转换步骤和候选选择记录；
3. **报告准备度**：运行后读取每个物种的 `wgdi_ready` 和保留率；
4. **报告细胞器结果**：列出 `inferred` 和 `review_candidates`；
5. **染色体清单与人工裁决**：`run` 输出已逐物种打印「染色体 N 条（编号:基因数）」，**必须把这份清单展示给调用方**（不得只报 seqs 总数或只字不提）；**必须取得人工签字**（`adjudicate`，reviewer 填实际复核人）后才能进入下游；不得代签或静默跳过；
6. **启用硬门禁**：验收用 `validate --require-adjudication`，下游设 `HOMOLOGY_REQUIRE_ADJUDICATION=true`；未裁决/过期/数量不符的物种会被 FAIL/跳过，如实报告；
7. **报告跳过物种**：下游运行后读取 `skipped_species.txt`；
8. **复核旧版本结果**：升级后备份并对指定物种使用 `--overwrite`（重跑会使旧裁决失效，需重新签字）；
9. **检查产物结构**：运行 `check-wgdi-outputs`，确认 8 个关键产物结构与内容完整；
10. **区分两个结论**：预处理完成代表输入就绪；WGDI 分析完成代表下游执行结束；
11. **归档依据**：使用 `--allow-not-ready` 时，在报告中注明其用于完整性归档，并同时报告未就绪原因。

## 测试

```powershell
cd pipeline\code
python -m pytest -q
```

测试覆盖多格式项目、GBFF 往返、GTF 转换、命名迁移、兼容链接、失败路径、版本一致性、共享脚本契约、发布包、Skill 镜像、染色体下限与人工裁决门禁（当前 201 项）。

## 交接材料

- 使用说明：`pipeline/QUICKSTART.md`
- 发布门禁：`pipeline/VALIDATION.md`
- 当前 Windows 验证回执：`pipeline/validation/smoke-report-3.6.1.json`
- Linux 验证回执：`pipeline/validation/smoke-report-3.5.2-linux.json`
- 真实物种验证：`pipeline/validation/real-species-3.5.2.json`

`analysis_v2/<Abbr>/` 是权威预处理输出，`clean_v2` 保存清洗记录，项目根 `<Abbr>/` 是下游兼容链接目录。
