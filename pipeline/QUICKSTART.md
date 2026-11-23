# 植物基因组 WGDI 与 DIAMOND 流程教师使用说明

本文面向需要使用本工具的老师和课题组同学，说明如何把本地植物基因组数据整理为可审计的 WGDI 输入。命令以 Windows PowerShell 为例，Linux 和 Git Bash 用户按相同参数执行即可。

## 1. 使用前确认

在仓库根目录执行：

```powershell
cd code
python -m pgci_preprocess --version
```

输出 `pgci_preprocess 3.8.0` 表示本地环境已就绪。

### 换电脑后第一步：扫描本机路径

新电脑上只有 skill 文件夹和原始数据时，不要手工改路径，也不要手工抄物种清单：

```powershell
# 数据在 <项目>\data\<show_id>_<Order>_<Family>(<Species>)/ 下：
python <skill>\pipeline\discover.py

# 数据在其他磁盘或目录：
python <skill>\pipeline\discover.py E:\raw_genomes --project-root D:\plant_project
```

这条命令不需要安装任何依赖，用系统里的 Python ≥ 3.11 直接运行。它会：

1. 定位代码根（按 `pipeline/discover.py` 自身位置，符号链接/复制都不影响）；
2. 扫描数据目录，按**文件内容**识别 gff/genome/cds/pep，并配对同一套组装；
3. 判断每个物种有无染色体证据，写 `metadata/sequence_mode.tsv`；
4. 写 `metadata/machine_paths.json` 与 `download_inventory.tsv`，并记住项目根。

屏幕上会逐物种列出选中的四个文件，`metadata/discovery_report.tsv` 列出全部候选
和得分，供人工复核。之后所有命令都可以省略 `--project-root`：

```powershell
python -m pgci_preprocess doctor
python -m pgci_preprocess run --sequence-mode from-map
```

创建项目后，可检查 Python、外部程序和输入文件：

```powershell
python -m pgci_preprocess doctor
```

新项目还没有物种数据时加 `--allow-empty`：

```powershell
python -m pgci_preprocess doctor --allow-empty
```

完整比较分析还可以执行：

```powershell
python -m pgci_preprocess doctor --workflow full
```

## 2. 准备输入文件

常规单物种输入包括：

```text
GFF/GFF3 注释
CDS FASTA
蛋白 FASTA
```

genome FASTA 作为常用输入，可帮助核对 region、坐标和序列长度。缺少 genome 时，程序能够根据 CDS 和转录本坐标完成长度计算与 region 回填。

`<Abbr>.lens` 是预处理输出，不需要提前准备。

### 支持的输入

| 已有文件 | 处理方式 |
| --- | --- |
| GFF/GFF3 | 直接进入标准预处理，支持多种注释结构 |
| GTF | 转换为 GFF3 后进入流程 |
| GBFF/GenBank | 提取编码区 GFF3、CDS、蛋白，可同时导出 genome |
| genome、CDS、蛋白 FASTA | 按文件名和内容识别角色 |
| `.gz`、`.bz2`、`.xz`、`.zip` | 自动识别压缩类型 |
| 功能注释 CSV/TSV | 规范化为辅助元数据，与坐标和序列数据一起归档 |

数据准备时建议做到：

1. 一个物种使用同一套组装版本；
2. GFF、CDS 和蛋白文件来自同一套注释；
3. 每个输入文件使用清晰、彼此不同的文件名；
4. 记录物种拉丁名、组装编号和下载来源；
5. 保留原始文件，程序会把源文件和哈希写入项目记录。

## 3. 单物种流程

### 3.1 创建项目

```powershell
cd code
python -m pgci_preprocess init-project --project-root D:\plant_project
```

程序会创建：

```text
D:\plant_project\
├── data\
│   └── {show_id}_{Order}_{Family}({Species})\   # 原始数据直接放这一层
├── resources\
│   ├── download_inventory.tsv
│   ├── species_abbreviations.tsv
│   └── format_exceptions.tsv
├── wgdi_input\
```

旧项目的 `data\downloads\<物种目录>` 仍能直接读：`data/` 里有物种目录时以 `data/` 为准，
否则回退到 `data\downloads`。把物种目录上提一层就是全部迁移工作。

### 3.2 注册物种（GFF/GFF3 输入）

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

参数说明：

| 参数 | 用途 |
| --- | --- |
| `--show-id` | 课题组内部物种编号 |
| `--scientific-name` | 物种拉丁名，用于生成唯一缩写 |
| `--order`、`--family` | 分类信息 |
| `--genome` | genome FASTA，可省略 |
| `--gff` | GFF/GFF3 注释 |
| `--cds` | CDS FASTA |
| `--pep` | 蛋白 FASTA |

`init-project` 将显式 GFF/GFF3、CDS 和蛋白文件镜像到项目目录，并写入库存表、装配信息表和来源清单。

### 3.3 使用 GTF 或 GBFF

GTF 和 GBFF 用户先完成格式转换，再用生成的文件按 3.2 注册。

GTF：

```powershell
python -m pgci_preprocess convert-annotation `
  --input D:\raw\species.gtf.gz `
  --output-dir D:\plant_project\data\derived `
  --prefix species
```

GBFF/GenBank：

```powershell
python -m pgci_preprocess convert-annotation `
  --input D:\raw\species.gbff.gz `
  --output-dir D:\plant_project\data\derived `
  --prefix species
```

GBFF 默认输出 GFF3、CDS、蛋白和 genome。转换完成后，将生成的 GFF3、CDS 和蛋白文件注册到项目。

批量处理时，`prepare` 和 `batch` 会自动完成 GTF/GBFF 转换，以及 genome + GFF 条件下缺少 CDS/蛋白时的提取。

### 3.4 运行

```powershell
python -m pgci_preprocess run --project-root D:\plant_project
```

序列筛选默认 **严格染色体模式**（`--sequence-mode chromosomes`）：只保留有染色体
挂载证据的序列；无染色体证据的物种报 `no_chromosome_sequences` 并记录原因，不会
用 scaffold 充数。碎片化组装确需保留 scaffold 时显式使用逃生门：

```powershell
python -m pgci_preprocess run --project-root D:\plant_project --sequence-mode adaptive
```

代表转录本按 SOP 规则选择：**CDS 最长优先，等长按稳定转录本 ID 字典序**，
不以蛋白长度或出现顺序为依据。

处理单个物种：

```powershell
python -m pgci_preprocess run `
  --project-root D:\plant_project `
  --species 2848
```

需要替换已有结果时：

```powershell
python -m pgci_preprocess run `
  --project-root D:\plant_project `
  --species 2848 `
  --overwrite
```

`--overwrite` 会替换相关物种的现有产物。使用前先备份 `wgdi_input` 和 `resources`。

### 3.5 验证

```powershell
python -m pgci_preprocess validate --project-root D:\plant_project

python -m pgci_preprocess check-names --project-root D:\plant_project
python -m pgci_preprocess audit --project-root D:\plant_project
```

新布局不需要兼容目录，下游脚本直接读取 `wgdi_input/<Abbr>/`。旧 `analysis_v2`
项目可先用 `migrate-layout` 迁移；只有暂不迁移的旧项目才需要 `--refresh-links`。

三个步骤全部通过后，将以下四个文件交给 WGDI：

```text
wgdi_input\<Abbr>\<Abbr>.pep
wgdi_input\<Abbr>\<Abbr>.cds
wgdi_input\<Abbr>\<Abbr>.new.gff
wgdi_input\<Abbr>\<Abbr>.lens
```

## 4. 批量流程

批量模式从状态表读取待处理物种。状态表使用 TSV 格式，必需列如下：

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

填写规则：

- `status` 使用 `complete`；
- `folder` 写物种原始文件夹的完整路径；
- `current_gff_or_gtf` 写 GFF/GFF3 或 GTF 文件名；
- 其余 `current_*` 写对应文件名；
- 同一角色有多个候选时使用 `A|B`；
- 可选列 `candidate_accession` 用于记录组装编号。

示例：

```text
show_id	Order	Family	Species	status	folder	current_genome	current_gff_or_gtf	current_cds	current_pep
1	Rosales	Moraceae	Morus notabilis	complete	D:\raw\species1	genome.fa.gz	annotation.gff3.gz	cds.fa.gz	pep.fa.gz
```

执行批量流程：

```powershell
python -m pgci_preprocess batch `
  --project-root D:\plant_project `
  --source-status D:\plant_project\resources\source_status.tsv `
  --expected-count 20
```

`batch` 包含：

```text
prepare → 可选 GFF 规范化 → run → validate
```

运行后重点查看：

```text
resources\preprocess_batches\
resources\wgdi_v2_index.tsv
resources\naming_audit.tsv
data-validation-report.tsv
```

已生成且版本一致的有效物种会标记为 `SKIPPED`。修改数据或筛选参数后，对相关物种使用 `run --species ... --overwrite` 重算。

## 5. 结果检查

### 5.1 先看总体结果

```powershell
python -m pgci_preprocess audit --project-root D:\plant_project
python -m pgci_preprocess readiness --project-root D:\plant_project --species Mun
```

重点记录：

- 处理物种数；
- `wgdi_ready` 中 `yes`、`marginal`、`no` 的数量；
- FAIL 物种及原因；
- 被剔除基因数和主要剔除原因。

### 5.2 再看单物种 QC

打开：

```text
wgdi_input\<Abbr>\<Abbr>.qc.json
wgdi_input\<Abbr>\<Abbr>.qc_record.md   # 人类可读 QC 摘要（由 qc.json 确定性生成）
```

关注：

| 字段 | 含义 |
| --- | --- |
| `n_regions` | 输入序列区域数 |
| `n_representative_genes` | 代表基因数 |
| `rejection_counts` | 各类剔除原因计数 |
| `sequence_selection` | 主序列选择模式和保留率 |
| `wgdi_bundle` | 最终 WGDI 基因数和序列数 |
| `organelle_inference` | 细胞器推断证据和复核候选 |
| `anomalies` | 格式和内容异常 |

### 5.3 查看细胞器复核项

在 `<Abbr>.qc.json` 中检查：

```text
organelle_inference.inferred
organelle_inference.review_candidates
```

`inferred` 是程序已识别并排除的序列。`review_candidates` 是建议人工复核的序列。

确认某个合成序列属于细胞器后，可指定正则重跑：

```powershell
python -m pgci_preprocess run `
  --project-root D:\plant_project `
  --species 2848 `
  --organelle-pattern "^(cp|mt)_scaffold" `
  --overwrite
```

### 5.4 处理 `wgdi_ready=no`

`wgdi_ready=no` 表示当前筛选集保留率低于 50%。若报 `no_chromosome_sequences`，
说明该物种没有染色体挂载证据（严格模式默认）：先确认组装是否染色体级别，
再决定是否用 `--sequence-mode adaptive`。其余情况建议依次检查：

1. 组装是否高度碎片化；
2. GFF 与 genome 的序列名是否一致；
3. 是否存在大量重复或等位 contig；
4. 是否需要调整 `--coverage`、`--min-genes-per-sequence` 或 `--max-sequences`；
5. 当前组装是否适合继续做 WGDI。

仅做数据完整性归档和复核时，可使用：

```powershell
python -m pgci_preprocess validate `
  --project-root D:\plant_project `
  --allow-not-ready
```

### 5.5 染色体清单与人工裁决（3.6.0 起）

`run` 后每个物种生成 `<Abbr>.chromosome_map.tsv`，列出保留了几条染色体、每条基因数/长度/来源，并在 run 输出里直接打印 `chromosomes N (id:genes) ...` 给调用方（ASCII 标签，避免非 UTF-8 控制台编码失败）。`--min-genes-per-chromosome`（CLI 默认 100）会丢弃基因数不足 100 的染色体（记入 chr_map 与 QC）。查看清单：

```powershell
python -m pgci_preprocess chromosomes --project-root D:\plant_project --species Hbr
```

老师复核后签字确认染色体数（记录绑定 QC 版本与哈希，重跑后旧签字失效）：

```powershell
python -m pgci_preprocess adjudicate --project-root D:\plant_project --species Hbr `
  --confirmed-n 21 --decision accept --reviewer "<复核人姓名>" --rationale "与核型一致"
python -m pgci_preprocess adjudication-status --project-root D:\plant_project --species Hbr
```

启用硬门禁后，未签字/过期/数量不符的物种 validate 判 FAIL、下游跳过（Skill 默认启用）：

```powershell
python -m pgci_preprocess validate --project-root D:\plant_project --require-adjudication
```

下游设 `HOMOLOGY_REQUIRE_ADJUDICATION=true` 即跳过未裁决物种。decision 可选 accept / accept_with_note / adjust / reject。

## 6. 可选下游分析

预处理验收完成后，可以继续执行 DIAMOND/BLAST 和 WGDI。

先检查完整工作流：

```powershell
python -m pgci_preprocess doctor `
  --project-root D:\plant_project `
  --workflow full
```

逐步执行：

```bash
bash 0check_data.sh --preprocess
bash run_homology.sh
bash 2gen_wgdi_conf.sh all
bash 3run_wgdi.sh
```

Windows 总控：

```powershell
bash 0check_data.sh --preprocess
powershell -ExecutionPolicy Bypass -File .\run_all_pipeline.ps1
```

`run_homology.sh` 和 Windows 总控会读取 `wgdi_ready`，并跳过未就绪物种。跳过名单写入：

```text
metadata/skipped_species.txt
```

结果汇报时同时列出成功物种、跳过物种和实施参数。

**Ks 分支默认关闭**：默认只产出点图与共线性块（`.png`/`.icl`）。需要 Ks 值表、
block-Ks 与分布峰值图时显式启用（配置生成、执行、验收三处都读同一开关）：

```bash
WGDI_KS=true bash 2gen_wgdi_conf.sh all
WGDI_KS=true bash 3run_wgdi.sh
```

Windows 总控用 `-WithKs`。原因：上游 blockinfo/correspondence 依赖 Ks 值，
缺 Ks 时 correspondence 会过滤掉全部区块，故这些阶段仅在启用时运行。

## 7. 常见问题

| 现象 | 处理方式 |
| --- | --- |
| 命令提示找不到 `pgci_preprocess` | 切换到仓库 `code` 目录，并确认当前 `python` 是已部署环境 |
| 产物名仍为 `s1234` | 提供 `--scientific-name`，或检查清单中的 `Species` 列 |
| 旧 `analysis_v2` 项目报 `compat_missing_*` | 运行 `migrate-layout`；暂不迁移时再用 `validate --refresh-links` |
| `validate` 报 `wgdi_not_ready_low_coverage` | 查看 QC 保留率；按第 5.4 节复核组装和筛选参数 |
| `run` 报 `no_chromosome_sequences` | 物种无染色体挂载证据；确认组装级别，或显式 `--sequence-mode adaptive` |
| 旧版本 QC 更新后需要重算 | 备份结果，对指定物种使用 `run --overwrite` |
| 旧项目的 `compat_missing_*` 与文件名迁移同时出现 | 先运行 `migrate-artifacts`，再运行 `migrate-layout` |
| 找不到 WGDI Python | 设置 `WGDI_PYTHON`，或使用本机已有的 `.wgdi-env` |
| `PAL2NAL` 未找到 | 设置 `WGDI_PAL2NAL` 指向本机已有脚本 |
| Windows 路径过长 | 将项目和临时目录放在较短路径下，建议总长度控制在 250 字符以内 |
| GBFF 翻译严格检查失败 | 检查 `/translation` 与 CDS 是否一致；将该记录列入数据复核 |

## 8. 目录瘦身（v3.8.0）

默认不再写筛选前的 `.all.*` 快照（实测占 `wgdi_input` 的 51%）。需要回溯全量候选时：

```powershell
python -m pgci_preprocess run --project-root D:\plant_project --keep-all
```

已有项目瘦身（同时更新索引与 SHA-256 清单，**瘦身后 `validate` 仍 PASS**）：

```powershell
# 先看能省多少
python -m pgci_preprocess prune --project-root D:\plant_project --all-snapshots --dry-run
# 确认后执行
python -m pgci_preprocess prune --project-root D:\plant_project --all-snapshots
```

根目录的 `<Abbr>/` 兼容目录是硬链接，不占额外空间；若只想让目录列表干净：

```powershell
python -m pgci_preprocess prune --project-root D:\plant_project --compat-links
```

不要手动删产物文件：手工删除不会同步索引与清单，`validate` 会报缺失。

## 9. 交付清单

项目交付给下游分析人员前，建议确认：

- [ ] `validate` 显示 PASS；
- [ ] `check-names` 显示 PASS；
- [ ] `wgdi_ready` 状态已记录；
- [ ] `review_candidates` 已完成人工复核；
- [ ] `wgdi_input/<Abbr>/` 四件套齐全；
- [ ] `<Abbr>.gene_map.tsv` 和 `<Abbr>.chr_map.tsv` 已核对；
- [ ] `<Abbr>.chromosome_map.tsv` 染色体数已人工裁决（`adjudication-status` 为 ok）；
- [ ] `<Abbr>.removed_list.tsv` 剔除原因已抽查；
- [ ] `<Abbr>.qc.json` 与 `<Abbr>.qc_record.md` 已归档；
- [ ] 使用的代表转录本、序列筛选参数和程序版本已记录；
- [ ] 发布包、验证回执和第三方工具版本已记录。

详细验证材料见 `VALIDATION.md` 和 `validation/`。
