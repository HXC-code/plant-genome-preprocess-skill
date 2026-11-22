# HXC 植物基因组规范化与质控工具

HXC 植物基因组规范化与质控工具是一套面向植物比较基因组学研究的标准化数据预处理软件。它把原始下载的注释、基因组和编码序列整理为命名统一、结构一致、来源可追溯的 WGDI 输入，并可继续执行 DIAMOND/BLAST、共线性、Ks 和 Ks-peak 分析。

当前版本：**3.6.0**

- [快速上手](QUICKSTART.md)
- [发布包](https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.6.0)
- [验证说明](VALIDATION.md)
- [Agent Skill 仓库](https://github.com/HXC-code/plant-genome-preprocess-skill)

## 工具定位

本工具是实验室植物基因组预处理 SOP 的执行载体。研究者负责确定物种、组装版本和筛选规则；程序负责批量执行、数据检查、结果留痕和失败原因汇总。

它主要用于以下工作：

- 将单物种或多物种的 GFF/GFF3、GTF、GBFF、genome、CDS 和蛋白序列整理为标准数据；
- 为每个基因选择一个代表转录本，形成一基因一转录本数据集；
- 选择适合 WGDI 分析的主序列，并生成统一命名和来源映射；
- 检查坐标、链方向、序列长度、CDS 翻译、ID 重复和文件哈希；
- 生成 QC、剔除记录、染色体映射、基因映射、物种索引和命名审计；
- 按需执行 DIAMOND/BLAST、WGDI 共线性、block-Ks 和 Ks-peak 分析。

## 适用对象

| 用户 | 使用方式 |
| --- | --- |
| 老师、研究生和技术人员 | 按标准流程整理并验收多个植物基因组的 WGDI 输入 |
| 课题组数据负责人 | 维护批次清单、物种缩写和版本记录，批量复现数据准备过程 |
| 下游分析人员 | 从统一四件套开始开展共线性、Ks 和 WGD 相关分析 |
| Agent 用户 | 通过 `plant-genome-preprocess-skill` 调用同一套工程和校验规则 |

## 输入数据

单物种预处理的常规输入为：

```text
GFF/GFF3 注释 + CDS FASTA + 蛋白 FASTA
```

genome FASTA 是常用输入。缺少 genome 时，程序可根据 CDS 和转录本坐标计算长度并回填 region。

支持的数据形式：

| 数据 | 支持格式 | 处理方式 |
| --- | --- | --- |
| 注释 | GFF/GFF3 | 识别标准层级、transcript-only、5 列简化、gene-span 和 gene-CDS-direct 结构 |
| 注释 | GTF | 转换为标准 GFF3 后进入流程 |
| 注释 | GBFF/GenBank | 提取编码区 GFF3、CDS 和蛋白，可同时导出 genome |
| 序列 | FASTA | 按文件名语义、首条记录和内容特征识别 genome、CDS 与蛋白角色 |
| 压缩文件 | `.gz`、`.bz2`、`.xz`、`.zip` | 按压缩文件内容识别，支持压缩包内单个文件 |
| 功能注释 | CSV/TSV | 规范化为辅助元数据，与 GFF/CDS/PEP 一起归档 |

文件后缀作为识别线索，最终格式和角色由文件内容确认。同一角色有多个候选时，可在清单中写 `A|B`，程序按角色兼容性选择唯一文件。

## 输出结果

每个物种在 `analysis_v2/<Abbr>/` 下生成 WGDI 四件套：

```text
<Abbr>.pep
<Abbr>.cds
<Abbr>.new.gff
<Abbr>.lens
```

同时生成：

```text
<Abbr>.all.pep.fa       全部有效代表基因蛋白
<Abbr>.all.cds.fa       全部有效代表基因 CDS
<Abbr>.all.gff          全部有效代表基因注释
<Abbr>.all.lens         全部有效代表基因涉及的序列长度
<Abbr>.gene_map.tsv     新基因 ID 与原始 ID 映射
<Abbr>.chr_map.tsv      新序列名与原始序列名映射
<Abbr>.qc.json          计数、阈值、异常、来源哈希和准备度
```

`clean_v2/<Abbr>/` 保存完整清洗记录：

```text
<Abbr>.gene_map.tsv     完整基因来源映射
<Abbr>.chr_map.tsv      完整染色体来源映射
<Abbr>.rejected.tsv     被剔除转录本、原因和详情
```

项目级结果还包括物种索引、物种注册表、SHA-256 清单、命名审计和数据验收报告。

## 使用流程

```text
原始注释与序列
      │
      ▼
项目初始化 / 输入识别 / 格式转换
      │
      ▼
一基因一代表转录本 + 主序列筛选
      │
      ├── analysis_v2/<Abbr>/   WGDI 四件套
      └── clean_v2/<Abbr>/      清洗与来源记录
      │
      ▼
validate + check-names + audit
      │
      ▼
DIAMOND/BLAST → WGDI → Ks / block-Ks / peaks
```

## 快速开始

以下流程适用于已经完成软件部署的本地环境。

### 1. 检查程序

```powershell
cd code
python -m pgci_preprocess --version
```

### 2. 创建项目并检查环境

```powershell
python -m pgci_preprocess init-project --project-root D:\plant_project
python -m pgci_preprocess doctor `
  --project-root D:\plant_project `
  --allow-empty
```

`init-project` 会创建标准目录、库存表、装配信息表和资源文件，并自动带上物种缩写表与格式例外表。

### 3. 注册一个物种

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

`--genome` 可省略。`init-project` 注册显式 GFF/GFF3、CDS 和蛋白文件，并将源文件镜像到项目目录。

GTF 或 GBFF 数据先通过以下任一方式转换：

```powershell
# 单文件转换
python -m pgci_preprocess convert-annotation `
  --input D:\raw\species.gtf.gz `
  --output-dir D:\plant_project\data\derived `
  --prefix species

# 或使用批量路由
python -m pgci_preprocess prepare `
  --project-root D:\plant_project `
  --source-status D:\plant_project\resources\source_status.tsv
```

### 4. 运行预处理

```powershell
python -m pgci_preprocess run --project-root D:\plant_project
```

已存在且版本一致的有效结果默认保留，并显示为 `SKIPPED`。需要重算指定物种时使用：

```powershell
python -m pgci_preprocess run `
  --project-root D:\plant_project `
  --species 2848 `
  --overwrite
```

更新到新版本或修改筛选参数时，先备份项目结果，再对相关物种使用 `--overwrite` 重算。

### 5. 验收结果

```powershell
python -m pgci_preprocess validate `
  --project-root D:\plant_project `
  --refresh-links

python -m pgci_preprocess check-names --project-root D:\plant_project
python -m pgci_preprocess audit --project-root D:\plant_project
```

验收通过后，`analysis_v2/<Abbr>/` 中的四件套可直接作为 WGDI 输入。

## 批量处理

批量模式使用一张状态表描述待处理物种。状态表至少包含：

```text
show_id  Order  Family  Species  status  folder
current_genome  current_gff_or_gtf  current_cds  current_pep
```

将 `status` 写为 `complete`，并在 `folder` 和四列 `current_*` 中填写实际文件。然后执行：

```powershell
python -m pgci_preprocess batch `
  --project-root D:\plant_project `
  --source-status D:\plant_project\resources\source_status.tsv `
  --expected-count 20
```

`batch` 按以下顺序执行：

```text
prepare → 可选 GFF 规范化 → run → validate
```

文件识别和候选选择记录在：

```text
resources/preprocess_batches/<batch_tag>_routing.tsv
```

批处理可以重复运行。已经生成且版本一致的有效物种默认保留，便于中断后继续处理其余物种。

## 数据检查与准备度

### 基因模型检查

程序会检查：

- CDS 长度是否为 3 的倍数；
- CDS 翻译与蛋白序列是否一致；
- 是否存在内部终止密码子和非法氨基酸字符；
- 坐标、链方向、region 和序列长度是否一致；
- gene、mRNA 和 CDS 层级是否匹配；
- 基因 ID、蛋白 ID 和 FASTA ID 是否重复或冲突。

被剔除的候选记录在 `<Abbr>.rejected.tsv`，每个条目包含作用域、ID、原因和详情。

### 代表转录本

每个基因按以下顺序选择一个代表转录本：

1. 蛋白最长；
2. 蛋白等长时选择 CDS 最长；
3. 仍相同时按 transcript ID 和 protein ID 的稳定顺序确定。

原始基因、转录本和蛋白 ID 均保存在 `<Abbr>.gene_map.tsv`。

### 主序列与 `wgdi_ready`

主序列选择优先使用注释中明确标记的主染色体。缺少主染色体标记时，程序按基因数、累计覆盖率和序列数上限选择。

| 参数 | 默认值 | 作用 |
| --- | ---: | --- |
| `--min-genes-per-sequence` | 50 | 保留基因数足够的序列 |
| `--coverage` | 0.90 | 达到目标累计基因覆盖率 |
| `--max-sequences` | 200 | 限制进入 WGDI 的序列数量 |

`wgdi_ready` 分为：

| 状态 | 条件 | 使用方式 |
| --- | --- | --- |
| `yes` | 基因保留率 ≥ 70%，序列数不超过上限 | 进入标准下游分析 |
| `marginal` | 基因保留率 ≥ 50% | 结合物种组装和研究目标复核后使用 |
| `no` | 基因保留率 < 50% | 默认停止下游分析；仅做完整性归档时使用 `--allow-not-ready` |

`wgdi_ready` 描述当前筛选集的基因保留率。论文分析时建议与 BUSCO、组装连续性和物种倍性信息联合判断。

### 细胞器序列

无 genome、无显式 region 的简化注释会进行叶绿体和线粒体推断。证据写入 `<Abbr>.qc.json` 的 `organelle_inference`：

- `inferred`：被识别并排除的细胞器序列；
- `review_candidates`：证据不足、建议人工复核的候选序列。

线粒体识别使用 `nad`、`cob`、`cox`、`atp`、`matR`、`mttB`、`ccm` 等核心基因家族。可用 `--organelle-mt-min-genes` 调整阈值，用 `--organelle-pattern` 指定必须排除的合成序列。

人工复核 `review_candidates` 后，将确认的序列加入 `--organelle-pattern` 并重跑。

## 命名规则

物种缩写采用“属名首字母 + 种加词字母”的全局唯一编码。例如 `Morus notabilis` 使用 `Mun`。

```text
基因 ID：Mun1g0001
含义：Mun = 物种缩写，1 = 序列编号，g0001 = 基因顺序号
```

已有缩写从 `code/resources/species_abbreviations.tsv` 读取。新物种按规则生成缩写，并避开表中所有已占用编码。

旧 `s<show_id>` 项目可以迁移：

```powershell
python -m pgci_preprocess migrate-artifacts `
  --project-root D:\plant_project `
  --dry-run

python -m pgci_preprocess migrate-artifacts `
  --project-root D:\plant_project
```

迁移完成后运行 `check-names`，核对缩写表、下载目录、库存表、索引和产物目录。

## 可选下游分析

预处理完成后，可继续执行 DIAMOND/BLAST 和 WGDI：

```powershell
bash 0check_data.sh --preprocess --refresh-links
bash run_homology.sh
bash 2run_pre_wgdi.sh
bash 3blast_to_wgdi.sh
python code/tools/generate_wgdi_configs.py --mode all
bash 6run_wgdi.sh
```

Windows 总控：

```powershell
bash 0check_data.sh --preprocess --refresh-links
powershell -ExecutionPolicy Bypass -File .\run_all_pipeline.ps1
```

下游运行使用 QC 和索引中的 `wgdi_ready`。未就绪物种记录在 `skipped_species.txt`，结题报告应同时列出成功物种和跳过物种。

## 使用建议

- **多倍体和高杂合物种**：选用能够区分等位 contig、homeolog 或亚基因组的注释来源，并在代表转录本筛选后进行人工抽查。
- **复杂细胞器**：查看 `organelle_inference.review_candidates`，对匿名、稀疏或证据不足的序列进行人工复核。
- **大规模 Ks 分析**：先在小样本验证参数和运行时间，再扩展到全蛋白组。
- **结果解释**：预处理记录输入是否干净、一致和可追溯；WGD 定年、Ks 峰解释和生物学结论由研究团队结合物种背景判断。

## 仓库关系

- `plant-genome-preprocess`：程序实现、测试、发布包和使用文档；
- `plant-genome-preprocess-skill`：Agent 调用入口和与程序仓同步的 `pipeline/` 运行副本。

两套仓库使用同一套代码、资源表和验证规则。

## 3.6.0 新增：染色体清单与人工裁决门禁

把“老师看 QC 和 chr_map 判断染色体数”的口头 SOP 变成程序里可审计、版本绑定的正式裁决步骤。

- **染色体清单**：`run` 后每个物种生成 `<Abbr>.chromosomes.tsv`（analysis_v2 与 clean_v2 各一份），列出保留了几条染色体、每条的基因数、长度、原始 seqid、是否 primary；QC 增加 `chromosome_summary`。用 `python -m pgci_preprocess chromosomes --project-root <root> --species <code>` 直接查看，供人工复核。
- **染色体基因数下限**：`run --min-genes-per-chromosome`（CLI 默认 **100**，0 关闭）在序列选择后丢弃基因数不足 100 的染色体（多为小 scaffold/片段），排除原因 `below_min_genes_per_chromosome` 记入 chr_map，丢弃数记入 QC。该下限可能降低保留率甚至 readiness，属显式可审计的取舍。
- **人工裁决硬门禁**：`adjudicate` 记录老师对最终染色体数的确认（confirmed_n、decision、署名 reviewer、理由 rationale），并**绑定被裁决 QC 的版本与 SHA-256**——重跑改变结果后旧签字自动失效（stale）。`validate --require-adjudication` 对缺失/过期/数量不符/拒绝的物种判 FAIL；下游 `run_homology.sh` / `run_all_pipeline.ps1` 设 `HOMOLOGY_REQUIRE_ADJUDICATION=true` 后跳过未裁决物种。**Skill 调用默认启用该门禁**，即“全部物种必须人工签字才能进入下游”。

```powershell
python -m pgci_preprocess chromosomes --project-root D:\proj --species Hbr
python -m pgci_preprocess adjudicate --project-root D:\proj --species Hbr `
  --confirmed-n 21 --decision accept --reviewer "刘涛" --rationale "与核型一致，叶绿体已排除"
python -m pgci_preprocess adjudication-status --project-root D:\proj --species Hbr
python -m pgci_preprocess validate --project-root D:\proj --require-adjudication
```

裁决决定：`accept` / `accept_with_note`（confirmed_n 须等于程序染色体数）、`adjust`（待调整，门禁拦截）、`reject`（拒绝，门禁拦截）。门禁默认关闭以便归档复核，由 SOP/Skill 显式开启。

## 测试与验证

```powershell
cd code
python -m pytest -q
python -m ruff check src tests scripts tools
python -m mypy src --ignore-missing-imports
```

当前版本包含 200 项 Python 自动化测试，并通过 Ruff、Mypy、Shell 语法检查、Windows/Linux 合成 WGDI 流程和真实物种数据验收。详细回执见 [VALIDATION.md](VALIDATION.md) 和 `validation/`。

## 许可

本项目使用 MIT License。WGDI、DIAMOND、BLAST+、MUSCLE、PAML、Perl 和 PAL2NAL 按其各自上游许可使用，具体说明见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
