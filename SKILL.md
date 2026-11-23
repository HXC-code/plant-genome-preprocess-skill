---
name: plant-wgdi-preprocess
description: 植物比较基因组 WGDI、DIAMOND 与共线性分析工作流（实验室整体流程第 1–4 步，Ks 可选）。完成本地数据整理、质控、代表转录本筛选、染色体清单与人工裁决、DIAMOND 同源比对和 WGDI 点图/共线性块；交付结果供系统发育、文献核对、多重共线性和圈图分析继续使用。
---

# 植物基因组 WGDI、DIAMOND 与共线性分析

按实验室《整体流程》（权威范围文档：`pipeline/code/docs/WORKFLOW.md`）实现
**第 1–4 步**。`pipeline/` 是与程序仓同步的完整可运行工程，不是一层调用壳。

## 能力边界（分析开始前必须向用户声明）

### 已覆盖（第 1–4 步）

1. 数据收集与整理：用户按网址自行下载（优先测序最完整、发表最新的来源），
   本工具注册/路由/校验本地文件；不可用物种在清单与失败记录中标注原因，不强制补全；
   染色体数目由染色体清单 + 人工裁决确认。
2. 数据预处理：产出 WGDI 四件套 `pep/cds/new.gff/lens`，CDS 与 PEP 一一对应，
   lens 由注释生成；含基因清洗、命名统一、QC 与剔除留痕。
3. DIAMOND 同源比对：物种内部 + 物种间两两比对（默认 e-value 1e-5），
   点图用于初步判断加倍事件。
4. WGDI 共线性块鉴定：默认产出点图与共线性块（`.png`/`.icl`）；
   **可选 Ks 分支**：仅用户明确要求时设 `WGDI_KS=true`，追加 Ks 值表、
   block-Ks 与分布峰值图。

### 后续分析交接（第 5–8 步）

本技能交付第 4 步的点图、共线性块和统计表；以下工作由后续分析流程承接：

- 系统发育树构建（第 5 步）
- 已发表树比对与加倍节点文献核对（第 6 步，人工分析）
- 多重共线性分析与多重比对列表（第 7 步）
- 共线性圈图绘制（第 8 步）

## 强制执行规则

1. 代表转录本筛选严格遵循「**CDS 长度最长优先 + 稳定转录本 ID 字典序消歧**」，
   不以蛋白长度或序列出现顺序为依据（引擎已实现，见 `selection.py`）。
2. 所有剔除操作必须留痕：`removed_list.tsv` 逐条记录 id 与原因，禁止无记录过滤。
3. Ks 分析默认关闭，仅用户明确要求时启用（`WGDI_KS=true` / Windows `-WithKs`）。
4. 序列筛选 CLI 默认严格染色体模式：无染色体挂载证据的物种报
   `no_chromosome_sequences` 并记录原因，不用 scaffold 充数；碎片化组装经用户
   确认后才用 `--sequence-mode adaptive`。
5. 审计链按用户决策保留关键环：输入指纹（路径+字节+SHA-256）记录在每物种
   `qc.json`，人工裁决绑定 QC 版本+指纹（重跑作废旧签字）；项目级哈希清单与
   `source-manifest` 命令已按用户决策移除，不属于必保范围。
6. 分析开始前明确告知用户上述能力边界（第 5–8 步不支持）。
7. 换电脑后先运行 `discover` 重建路径地图，禁止把某一台机器的绝对路径、符号链接
   或人工写死的物种清单当成流程前提；哪些物种可用、是否有染色体证据，都由扫描
   结果决定（`--sequence-mode from-map`）。

## 核心工作流

### 步骤 0：换电脑 / 新环境先扫描（discover）

数据放哪儿由扫描确定，**不要手工改路径、不要依赖符号链接**。新电脑上只需要
skill 目录和原始数据，任意 Python ≥ 3.11 直接运行：

```powershell
# 数据在 <项目>\data\<show_id>_<Order>_<Family>(<Species>)/ 时，在项目根执行：
python <skill>\pipeline\discover.py

# 数据在别处时，显式给出数据目录和项目位置：
python <skill>\pipeline\discover.py E:\raw_genomes --project-root D:\my_project
```

扫描会完成四件事，并在屏幕上逐物种打印选中的注释与 genome/cds/pep（`annot=`）：

| 产物 | 作用 |
|---|---|
| `resources/machine_paths.json` | 本机路径地图：代码根、项目根、数据根、解释器、外部工具 |
| `resources/download_inventory.tsv` | 由文件内容识别角色并自动配对同一套组装；已有物种行时默认不覆盖（`--force` 才重写） |
| `resources/sequence_mode.tsv` | 逐物种判定有无染色体证据，给出 `chromosomes` / `adaptive` |
| `resources/discovery_report.tsv` | 全部候选与得分，供人工复核 |

同时把项目根写进指针（`pipeline/code/.pgci-machine.json`，只读安装时回落
`~/.pgci-machine.json`），因此之后**所有命令都不用再传 `--project-root`**：

```powershell
python -m pgci_preprocess doctor
python -m pgci_preprocess run --sequence-mode from-map
```

`from-map` 表示逐物种套用 `sequence_mode.tsv`：有染色体证据的走严格模式，
scaffold-only 组装的物种才走 adaptive，不再需要人工记哪几个物种特殊。

新 Linux 机器的解释器与 WGDI 环境也必须现场确认：激活目标 conda/venv 后运行
`discover`，或在 `WGDI_PYTHON` 中显式指定；启动器只检查显式变量、当前激活环境和
项目内环境，不扫描其他用户目录或把某台机器的安装路径写死。

项目根只保留 `data/`（原始数据）、skill、`README.md`、以及两个数据目录：
`metadata/`（扫描清单、路径地图、序列模式表、`folder_list.txt`、`species_map.tsv`
等索引；`resources/` 是指向它的兼容联接）和 `wgdi_input/`（预处理四件套、DIAMOND
blast、扁平配置）。所有下游结果统一进 `wgdi_output/results/`，日志进 `logs/`、
断点标记进 `state/`。旧 `analysis_v2/` 布局可用 `migrate-layout` 迁移。

识别规则（与电脑和物种无关，按内容判断）：染色体证据只看 genome header /
GFF seqid 的 `Chr1`、`Chr01`、`Ajul.Chr01`、`Chromosome 1` 等写法，`NW_…`
scaffold 不算；`chloroplastic`/`mitochondrial` 只出现在基因产物描述里不算细胞器，
只有整套序列都标注细胞器才排除；SwissProt/InterPro/单拷贝比对等参考文件不参与配对。

### 步骤 1：数据收集与基础校验

```powershell
cd pipeline\code
python -m pip install -e ".[gbff]"        # Python ≥ 3.11；GBFF 需要 biopython
python -m pgci_preprocess init-project --project-root D:\my_project
# 单物种注册（四个显式路径，角色不靠后缀猜测）：
python -m pgci_preprocess init-project --project-root D:\my_project `
  --show-id 2848 --scientific-name "Morus notabilis" --order Rosales --family Moraceae `
  --genome genome.fasta.gz --gff annotation.gff3.gz --cds cds.fa.gz --pep pep.fa.gz
python -m pgci_preprocess doctor --project-root D:\my_project
```

批量：物种目录直接放入 `data/<show_id>_<Order>_<Family>(<Species>)/`（旧布局
`data/downloads/` 仍可读），编辑 `resources/download_inventory.tsv`
（gff/genome/cds/pep 四列，多候选写 `A|B`），然后 `prepare`/`batch`。
完整性校验（序列 ID 一致性、GFF 匹配、压缩完整性）由 `doctor`/`prepare`/`validate`
完成；不可用数据源标注原因，不强制补全。

### 步骤 2：预处理与代表序列构建

```powershell
python -m pgci_preprocess run      --project-root D:\my_project
python -m pgci_preprocess validate --project-root D:\my_project
python -m pgci_preprocess check-names --project-root D:\my_project
```

- 主染色体筛选：默认只保留染色体挂载序列，剔除未挂载 scaffold/contig 与细胞器
  （叶绿体/线粒体推断 + 复核候选）；`--min-genes-per-chromosome`（CLI 默认 100）
  丢弃基因数不足的染色体，全部记入 `chr_map.tsv` 与 QC。
- 代表转录本：CDS 最长优先，等长按稳定 ID 字典序（强制规则 1）。
- 染色体数人工裁决（SOP）：`chromosomes` 查看清单 → 复核人 `adjudicate` 签字
  （绑定 QC 版本+SHA-256，重跑作废）→ `validate --require-adjudication` 硬门禁。
  Skill 调用默认启用该门禁。

### 步骤 3：DIAMOND 同源比对（必做）

```bash
bash run_homology.sh        # 或 1run_diamond.sh / 1run_blast.sh
```

物种内 + 物种间两两比对，默认 very-sensitive、e-value 1e-5；按 `wgdi_ready`
与裁决状态跳过未就绪物种（`metadata/skipped_species.txt` 留痕）。

### 步骤 4：WGDI 共线性分析（必做）+ 可选 Ks

```bash
bash 2gen_wgdi_conf.sh all    # 生成扁平比较目录与相对路径配置（零复制）
bash 3run_wgdi.sh             # 核心：.png + .icl（阶段级断点续跑）
```

**可选 Ks 分支（仅用户明确要求时）**：

```bash
WGDI_KS=true bash 2gen_wgdi_conf.sh all
WGDI_KS=true bash 3run_wgdi.sh    # 追加 .ks/.bi.csv/.bi.c.csv/.bk.png/.kp.png/.pf.png
```

Windows 总控：`run_all_pipeline.ps1`（1–3 步），加 `-WithKs` 启用 Ks。
原因：上游 blockinfo/correspondence 依赖 Ks 值，缺 Ks 时 correspondence 会
过滤掉全部区块，故六个 Ks 阶段只在启用时运行。

### 追溯报告输出（无冗余）

每物种四类核心追溯文件 + 审计附件：

1. `chromosome_map.tsv`：染色体编号、长度、基因数、原始序列 ID
2. `gene_map.tsv`：代表基因 ID 对应与来源
3. `removed_list.tsv`：全部剔除序列与原因
4. `qc_record.md`：人类可读全流程 QC 记录（由 `qc.json` 确定性生成）

审计附件（保留）：`qc.json`（含输入指纹）、`chr_map.tsv`（逐序列取舍）、
裁决记录、`resources/` 索引与台账（不再写项目级哈希清单）。

### 输出体积控制（v3.8.0）

筛选前的 `.all.*` 快照**默认不写**（实测占产物目录的 51%）：

- 需要回溯全量候选时：`run --keep-all`；
- 已有项目瘦身：`prune --all-snapshots`（先加 `--dry-run` 看能省多少）。
  `prune` 会同时清空索引中失效路径，**瘦身后 `validate` 仍 PASS**；
- 新布局不再生成根目录 `<Abbr>/` 兼容目录；旧项目可用 `migrate-layout` 一并清理；
- **不要手动删产物文件**：手工删除不会同步索引，`validate` 会报缺失。

## 支持哪些后缀

按**内容**判定，后缀只作提示；`.gz/.bz2/.xz/.zip` 按压缩魔数识别。

| 类型/后缀 | 处理 |
| --- | --- |
| `.gff` `.gff3` | standard / transcript_only / simplified5 / gene_span / gene_cds_direct 五种结构 |
| `.gtf` | gene/transcript/exon/CDS/start_codon/stop_codon，转 GFF3 后继续 |
| `.gbff` `.gbk` `.gb` `.genbank` | 提取编码 GFF3 + CDS/PEP（可导 genome）；需 biopython |
| `.annotations.csv` `.csv` `.tsv` | 七列功能表规范化为辅助元数据，不替代坐标注释 |
| `.genome.fa` `.genome.fasta` `.genome.fna` `.genome.fas` `_genomic.fna` `.fna` | genome FASTA |
| `.cds.fa` `.cds.fasta` `.cds.fna` `.CDS.fasta` `.cds.fas` | CDS FASTA |
| `.pep.fa` `.pep.fasta` `.protein.faa` `.protein.fasta` `.pro.fa` `.faa` | 蛋白 FASTA |
| `.fa` `.fasta` `.fna` `.faa` `.fas` | 通用 FASTA，角色由文件名语义 + 首记录内容判定 |
| `.gz` `.bz2` `.xz` `.zip` | 可组合到上述后缀；按压缩魔数识别 |

常见压缩组合包括 `.gff.gz`、`.gff3.gz`、`.gtf.gz`、`.gbff.gz`、
`.genome.fasta.gz`、`.CDS.fasta.gz`、`.Protein.faa.gz`、`.pep.fa.gz` 等。
ZIP 要求压缩包内只有一个非空数据文件。

非标准后缀的文件只要内容符合上述格式，也可以进入内容识别流程；`discover`
对标准后缀和带角色语义的文件名识别更稳定。`.lens` 是输出文件，不是输入。

## 极端情况（都有明确行为）

| 情况 | 行为 |
| --- | --- |
| 无染色体挂载证据 | 默认报 `no_chromosome_sequences`，记录原因；经确认用 `--sequence-mode adaptive` |
| 只有 genome、没有注释 | 报 `annotation_missing`，不凭空造 GFF |
| 只有注释/CDS/蛋白、没有 genome | 允许；长度由 CDS 计算，region 回填 |
| 只有 GBFF / 只有 GTF | 自动转换后继续 |
| 同目录多套重复 | 清单 `A|B` 打分选择，平局报歧义 |
| 细胞器混入 | 推断排除 + `review_candidates` 人工复核 |
| 旧 `s<show_id>` 产物 | `migrate-artifacts` 幂等迁移为缩写 |

## 校验与命名

- `validate`：四件套齐全、基因数一致、坐标/链/长度合法、兼容链接内容一致；
  `--require-adjudication` 拦截缺失/过期/拒绝/数量不符的裁决。
- `check-names`：缩写表 ↔ 清单 ↔ 索引 ↔ 产物一致性（`naming_audit.tsv`）。
- 物种缩写：属名首字母 + 种加词字母，3 位起、冲突换组合/扩位；表内 show_id 沿用。

## 测试

```powershell
cd pipeline\code
python -m pytest -q        # 259 个用例；临时目录受限时加 --basetemp=<short-path>
```

## 交接与复现

- 权威范围与步骤映射：`pipeline/code/docs/WORKFLOW.md`；发布门禁：`pipeline/VALIDATION.md`。
- 验证回执：`validation/smoke-report-3.8.0.json`（Ks 模式 4×8）与
  `validation/smoke-report-3.8.0-core.json`（核心模式 4×2），均由真实运行环境生成。
- WGDI 用钉版 `wgdi==0.75` 上游包；`bin/wgdi_entry.py` 只维护必要兼容补丁；
  PAL2NAL 由 `install_pal2nal.py` 从官方源安装，不随包分发。
- CI 跑 pytest/ruff/mypy + 镜像漂移检查；代码变更先在 canonical 仓库落地再同步。
- 预处理与下游比较是两个明确模式；"四件套就绪"不等于 WGDI 已完成。
- `wgdi_input/<Abbr>/` 保存四件套与追溯文件；`wgdi_input/blast/` 保存 DIAMOND
  结果；`wgdi_input/config/` 保存每个比较的扁平配置。`wgdi_output/results/`
  保存全部点图与共线性块，`wgdi_output/logs/` 与 `wgdi_output/state/` 保存日志
  和阶段级断点标记。配置全部用项目相对路径，项目可整体搬移。旧项目用
  `pgci_preprocess migrate-layout` 迁移（`--dry-run` 可预演）。
