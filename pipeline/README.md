# HXC 植物基因组规范化与质控工具

HXC 植物基因组规范化与质控工具把植物注释、genome、CDS 和蛋白序列整理为命名统一、结构一致、来源可追溯的 WGDI 输入，并提供染色体清单、QC、人工裁决和可选的下游比较分析。

当前版本：**3.7.0**

- [快速上手](QUICKSTART.md)
- [版本记录](CHANGELOG.md)
- [验证说明](VALIDATION.md)
- [整体流程（权威范围文档）](code/docs/WORKFLOW.md)
- [3.7.0 发布页](https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.7.0)
- [Agent Skill 仓库](https://github.com/HXC-code/plant-genome-preprocess-skill)

**范围**：按实验室《整体流程》实现第 1–4 步（数据收集整理 → 预处理四件套 →
DIAMOND 同源比对 → WGDI 共线性）；Ks 为第 4 步内可选分支（`WGDI_KS=true`）。
第 5–8 步（系统发育树、文献核对、多重共线性、圈图）不在范围内，分析前必须声明。

## 3.7.0 重点

- **代表转录本 SOP 规则**：CDS 最长优先，等长按稳定转录本 ID 字典序；不用蛋白长度或出现顺序。
- **严格染色体筛选**：CLI 默认 `--sequence-mode chromosomes`，无染色体证据的物种 fail-closed（`no_chromosome_sequences`）；碎片化组装可用 `adaptive` 逃生门。
- **Ks 可选**：默认只跑 dotplot + collinearity（`.png`/`.icl`）；`WGDI_KS=true`（Windows 总控 `-WithKs`）追加六个 Ks 阶段。
- **SOP 追溯文件名**：`chromosome_map.tsv`、`removed_list.tsv`、新增人类可读 `qc_record.md`；审计链（哈希/manifest/裁决绑定）保留。

## 3.6.x 重点

- **染色体清单**：`run` 输出每个物种保留了多条序列、序号和每条序列的代表基因数，并生成 `<Abbr>.chromosome_map.tsv`。
- **染色体基因数下限**：CLI 默认丢弃代表基因数少于 100 的 WGDI 序列，排除记录写入 `chr_map` 和 QC。
- **最终数量人工裁决**：复核人确认染色体数并签字；记录绑定 QC 版本和 SHA-256，重跑后旧签字自动失效。
- **严格验收门禁**：`validate --require-adjudication` 拦截缺失、过期、拒绝或数量不符的裁决记录。
- **Windows 控制台兼容**：3.6.2 修复了 3.6.1 在非 UTF-8 控制台打印染色体清单时的异常。

## 一个合格的预处理结果

一个物种进入下游分析前，应同时具备：

```text
来源文件和版本可追溯
GFF/GTF/GBFF 路由与转换正确
每个基因一个代表转录本
CDS 与蛋白翻译一致
细胞器和异常 taxon 已处理
WGDI 四件套结构闭合
命名审计通过
wgdi_ready=yes 或经复核的 marginal
染色体清单已人工确认
裁决状态为 ok
```

## 主要能力

| 能力 | 结果 |
| --- | --- |
| 输入识别 | 识别 GFF/GFF3、GTF、GBFF、genome、CDS、蛋白和压缩格式 |
| 格式转换 | GTF/GBFF 转为标准 GFF3，并可提取 CDS、蛋白和 genome |
| 基因模型整理 | 补充 region、mRNA，处理 gene/mRNA/CDS 层级和拆分基因 |
| 代表转录本 | 每个基因保留一个代表转录本（CDS 最长优先 + 稳定 ID 字典序） |
| 序列筛选 | CLI 默认严格染色体模式；adaptive 模式按基因数/覆盖率/序列上限回退 |
| 染色体清单 | 统计最终保留数、每条序列长度、代表基因数和来源 |
| 数据 QC | 检查坐标、链、长度、CDS 翻译、ID 重复、taxon 和哈希 |
| 来源追踪 | 生成 gene map、chromosome map、剔除记录和输入哈希 |
| 人工裁决 | 记录最终染色体数、复核人、依据和 QC 指纹 |
| 可选下游 | 执行 DIAMOND/BLAST、WGDI 共线性；Ks/block-Ks/peaks 需 `WGDI_KS=true` 显式启用 |

## 输入数据

常规输入：

```text
注释 + CDS FASTA + 蛋白 FASTA
```

genome FASTA 是常用输入。缺少 genome 时，程序可使用 CDS 和转录本坐标计算长度并回填 region。

| 数据 | 支持格式 |
| --- | --- |
| GFF/GFF3 | standard、transcript-only、5 列简化、gene-span、gene-CDS-direct |
| GTF | 转换为标准 GFF3 后进入流程 |
| GBFF/GenBank | 提取编码区 GFF3、CDS、蛋白，可同时导出 genome |
| genome/CDS/蛋白 | FASTA；由文件名语义、首条记录和内容特征识别角色 |
| 压缩文件 | gzip、bzip2、xz、单成员 ZIP；按文件内容识别 |
| 功能注释 | CSV/TSV；规范化为辅助元数据 |

输入规则：

- 文件后缀作为识别线索，最终格式和角色由内容确认。
- 同一角色有多个候选时，可在清单中写 `A|B`，程序必须选出一个唯一兼容文件。
- 候选文件必须位于对应物种目录内。
- 显式注册要求 GFF/GFF3、CDS 和蛋白三个不同且非空的源文件。
- GTF/GBFF 先经过 `convert-annotation`，或通过 `prepare/batch` 自动路由。

## 输出结果

`analysis_v2/<Abbr>/` 包含：

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

`clean_v2/<Abbr>/` 包含：

```text
<Abbr>.gene_map.tsv
<Abbr>.chr_map.tsv
<Abbr>.removed_list.tsv
<Abbr>.chromosome_map.tsv
```

### WGDI 四件套

```text
<Abbr>.pep
<Abbr>.cds
<Abbr>.new.gff
<Abbr>.lens
```

`new.gff` 固定七列：

```text
数字序列号  唯一基因ID  起点  终点  链方向  序列内顺序号  原始基因ID
```

`lens` 固定三列：

```text
序列号  序列长度  代表基因数
```

### 染色体清单

`<Abbr>.chromosome_map.tsv` 记录最终保留的 WGDI 序列：

| 字段 | 含义 |
| --- | --- |
| `chromosome` | 产物内部序列编号 |
| `source_seqid` | 原始序列 ID |
| `length` | 序列长度 |
| `gene_count` | 该序列保留的代表基因数 |
| `is_primary` | 注释是否将其标为主染色体 |

程序统计的是最终 WGDI 序列数。主染色体、未挂载 scaffold 和碎片化 contig 都可能进入该清单，因此最终数量需要人工结合 `source_seqid`、`is_primary` 和物种核型确认。

## 快速开始

以下命令适用于已经部署好的本地环境。

### 1. 检查程序

```powershell
cd code
python -m pgci_preprocess --version
```

当前版本应输出：

```text
pgci_preprocess 3.6.2
```

### 2. 创建项目

```powershell
python -m pgci_preprocess init-project --project-root D:\plant_project
python -m pgci_preprocess doctor `
  --project-root D:\plant_project `
  --allow-empty
```

`init-project` 创建标准目录、库存表、装配信息表、物种缩写表和格式例外表。

### 3. 注册物种

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

GTF/GBFF 的单文件转换示例：

```powershell
python -m pgci_preprocess convert-annotation `
  --input D:\raw\species.gtf.gz `
  --output-dir D:\plant_project\data\derived `
  --prefix species
```

### 4. 运行预处理

```powershell
python -m pgci_preprocess run --project-root D:\plant_project
```

每个物种处理完成后，命令直接显示类似：

```text
[OK] Mun show_id=2848 Morus notabilis genes=... seqs=... ready=yes reason=
    chromosomes 10 (id:genes) 1:4210  2:3980  3:4521 ...
    inventory: D:\plant_project\analysis_v2\Mun\Mun.chromosome_map.tsv
```

其中 `chromosomes` 是最终 WGDI 序列数，括号中是对应序列的代表基因数。

已有且版本一致的有效结果默认保留。需要重算指定物种时：

```powershell
python -m pgci_preprocess run `
  --project-root D:\plant_project `
  --species 2848 `
  --overwrite
```

### 5. 复核染色体清单

```powershell
python -m pgci_preprocess chromosomes `
  --project-root D:\plant_project `
  --species Mun
```

同时查看：

```text
analysis_v2\Mun\Mun.chromosome_map.tsv
analysis_v2\Mun\Mun.chr_map.tsv
analysis_v2\Mun\Mun.qc.json
```

重点核对：

- 总条数是否符合该物种核型；
- `source_seqid` 是否为真实染色体；
- `is_primary=N` 的序列是否应保留；
- 是否出现未挂载 scaffold；
- `below_min_genes_per_chromosome` 是否误删真实小染色体。

### 6. 人工裁决

复核人确认数量后执行：

```powershell
python -m pgci_preprocess adjudicate `
  --project-root D:\plant_project `
  --species Mun `
  --confirmed-n 10 `
  --decision accept `
  --reviewer "<实际复核人>" `
  --rationale "<核型或文献依据>"
```

裁决记录：

```text
resources/adjudications/Mun.json
```

`confirmed_n` 必须等于程序生成的 `n_wgdi_sequences`。如果复核认为需要调整，使用 `--decision adjust`，修正输入、阈值或选择规则后重跑，再重新裁决。

### 7. 严格验收

```powershell
python -m pgci_preprocess validate `
  --project-root D:\plant_project `
  --refresh-links `
  --require-adjudication

python -m pgci_preprocess check-names --project-root D:\plant_project
python -m pgci_preprocess audit --project-root D:\plant_project
```

严格验收通过后，`analysis_v2/<Abbr>/` 中的四件套可作为 WGDI 输入。

## 染色体选择与 100 基因下限

程序按以下顺序生成候选 WGDI 序列：

1. 有主染色体标记时，将全部非细胞器主序列作为候选；
2. 无主染色体标记时，保留代表基因数不少于 50 的序列；
3. 阈值方案不足时，按基因数从高到低累计选择，目标覆盖 90% 基因，最多 200 条序列；
4. 在候选集合上应用染色体基因数下限。

默认参数：

| 参数 | CLI 默认值 | 含义 |
| --- | ---: | --- |
| `--min-genes-per-sequence` | 50 | 阈值筛选的每条序列最少代表基因数 |
| `--coverage` | 0.90 | 碎片化组装的累计基因覆盖率 |
| `--max-sequences` | 200 | WGDI 序列数上限 |
| `--min-genes-per-chromosome` | 100 | 最终保留集合中每条序列的最少代表基因数 |

`--min-genes-per-chromosome 100` 表示：

```text
99 个代表基因  删除
100 个代表基因 保留
101 个代表基因 保留
```

该过滤在序列选择完成后执行，并按过滤后的集合重新计算保留率和 `wgdi_ready`。

CLI 默认 100；直接调用 Python API 的 `PipelineOptions` 时默认 0，需要显式传入 100。传入 0 可关闭该过滤。

## `wgdi_ready`

| 状态 | 条件 | 使用方式 |
| --- | --- | --- |
| `yes` | 基因保留率 ≥ 70%，序列数不超过上限 | 进入标准下游分析 |
| `marginal` | 基因保留率 ≥ 50% | 结合组装和筛选记录复核 |
| `no` | 基因保留率 < 50% | 默认不进入下游；归档时显式使用 `--allow-not-ready` |

`wgdi_ready` 反映当前筛选集的代表基因保留率。论文分析时建议同时记录 BUSCO、组装连续性和物种倍性信息。

## 人工裁决门禁

裁决状态：

| 状态 | 条件 |
| --- | --- |
| `ok` | 记录存在，QC 版本和 SHA-256 匹配，decision 为 accept/accept_with_note，reviewer 非空，且确认数等于 `n_wgdi_sequences` |
| `missing` | 没有裁决记录 |
| `stale` | QC 版本或 SHA-256 已变化 |
| `count_mismatch` | 确认数与程序结果不一致 |
| `rejected` | decision 为 reject |
| `pending_adjustment` | decision 为 adjust |
| `no_qc` | 找不到 QC |

严格门禁：

```powershell
python -m pgci_preprocess validate `
  --project-root D:\plant_project `
  --require-adjudication
```

下游：

```powershell
$env:HOMOLOGY_REQUIRE_ADJUDICATION = "true"
```

主程序的这两个开关默认关闭；Skill 执行规则要求启用。实验室交付流程以严格门禁通过为准。

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

`status` 使用 `complete`。执行：

```powershell
python -m pgci_preprocess batch `
  --project-root D:\plant_project `
  --source-status D:\plant_project\resources\source_status.tsv `
  --expected-count 20
```

`batch` 顺序：

```text
prepare → 可选 GFF 规范化 → run → validate
```

路由记录：

```text
resources/preprocess_batches/<batch_tag>_routing.tsv
```

批量处理完成后，仍需逐物种查看染色体清单并完成人工裁决，再运行严格验收。

## 数据检查

程序检查：

- CDS 长度是否为 3 的倍数；
- CDS 与蛋白翻译是否一致；
- 是否存在内部终止密码子和非法字符；
- 坐标、链、region 和序列长度是否一致；
- gene、mRNA 和 CDS 层级是否匹配；
- gene、protein 和 FASTA ID 是否重复或冲突；
- 同一物种是否混入多个核基因组 taxon；
- GFF、CDS、PEP 和 lens 的 ID/序列集合是否闭合；
- 压缩格式、输入 SHA-256 和转换链是否记录。

被剔除的转录本记录在：

```text
clean_v2/<Abbr>/<Abbr>.removed_list.tsv
```

`validate` 检查文件存在性、名称、七列 GFF、三列 lens、坐标、顺序、ID 集合、字符集、哈希、兼容链接、`wgdi_ready` 和可选裁决记录。

## 命名标准

物种缩写从 `code/resources/species_abbreviations.tsv` 读取；新物种按拉丁名生成全局唯一缩写。

基因 ID：

```text
<Abbr><sequence_number>g<至少四位顺序号>
例如 Mun1g0001
```

旧 `s<show_id>` 项目可以迁移：

```powershell
python -m pgci_preprocess migrate-artifacts `
  --project-root D:\plant_project `
  --dry-run

python -m pgci_preprocess migrate-artifacts `
  --project-root D:\plant_project
```

迁移后运行 `check-names`，核对缩写表、下载目录、库存表、索引和产物目录。

## 可选下游分析

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
$env:HOMOLOGY_REQUIRE_ADJUDICATION = "true"
bash 0check_data.sh --preprocess --refresh-links
powershell -ExecutionPolicy Bypass -File .\run_all_pipeline.ps1
```

下游根据 `wgdi_ready` 和人工裁决状态选择物种。跳过记录写入 `skipped_species.txt`，报告中应同时列出成功和跳过物种。

**Ks 分支默认关闭**：上述命令只产出点图与共线性块（`.png`/`.icl`）。需要
Ks 值表、block-Ks 与分布峰值图时，显式启用：

```bash
WGDI_KS=true bash 6run_wgdi.sh   # 配置生成与验收同样读取 WGDI_KS
```

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all_pipeline.ps1 -WithKs
```

原因：上游 WGDI 的 blockinfo/correspondence 依赖 Ks 值，缺 Ks 时
correspondence 会过滤掉全部区块，因此这些阶段只在启用 Ks 时运行。

## 科研使用范围

- 多倍体和高杂合数据应选择能够区分等位 contig、homeolog 和亚基因组的注释来源，并在代表转录本筛选后人工抽查。
- 复杂细胞器应同时检查 `organelle_inference.inferred` 和 `review_candidates`。
- 100 基因下限按数量机械过滤，可能移除基因数少的真实染色体，需要结合 `chr_map.tsv` 复核。
- `wgdi_ready=yes` 表示当前筛选集的基因保留率达到标准，建议与 BUSCO 和组装评估联合使用。
- 大规模 Ks 分析先在代表性小样本验证参数和运行时间。
- WGD 定年、Ks 峰解释和物种演化结论由研究团队结合物种背景判断。

## 测试与验证

```powershell
cd code
python -m pytest -q
python -m ruff check src tests scripts tools
python -m mypy src --ignore-missing-imports
```

3.7.0 版本包含 212 项 Python 自动化测试，覆盖多格式输入、转换、命名、QC、染色体清单、100 基因下限、人工裁决门禁、严格染色体模式、CDS 优先规则、Ks 可选分支、下游门禁和 Windows 控制台回归。

验证回执：

```text
validation/smoke-report-3.7.0.json        # 完整 Ks 模式（4 组比较 × 8 产物）
validation/smoke-report-3.7.0-core.json   # 默认核心模式（4 组比较 × 2 产物）
```

详细环境、回执和发布门禁见 [VALIDATION.md](VALIDATION.md)。

## 仓库关系

- `plant-genome-preprocess`：程序实现、测试、发布材料和主文档；
- `plant-genome-preprocess-skill`：Agent 调用入口和同步的 `pipeline/` 运行副本。

Skill 仓库通过镜像检查保持与主程序仓一致。

## 许可

项目代码使用 MIT License。WGDI、DIAMOND、BLAST+、MUSCLE、PAML、Perl 和 PAL2NAL 按其上游许可使用，详见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
