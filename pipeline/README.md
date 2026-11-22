# 植物基因组 WGDI 流程预处理工具

把下载到本地的植物基因组注释与序列，整理成 WGDI 可直接使用的四件套，并可继续完成
DIAMOND 同源比对与 WGDI 共线性块鉴定。

本工具按实验室《[整体流程](code/docs/WORKFLOW.md)》实现 **第 1–4 步**：

```text
1 数据收集整理  →  2 预处理（pep/cds/gff/lens）  →  3 DIAMOND 比对  →  4 WGDI 共线性
                                                                        └─ Ks 分析（可选）
```

| 信息 | 值 |
| --- | --- |
| 当前版本 | **3.7.0** |
| 运行环境 | Python ≥ 3.11（预处理纯标准库；GBFF 需 biopython） |
| 文档 | [快速上手](QUICKSTART.md) · [整体流程](code/docs/WORKFLOW.md) · [验证说明](VALIDATION.md) · [版本记录](CHANGELOG.md) |
| 下载 | [v3.7.0 Release](https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.7.0) |
| Agent Skill | [plant-wgdi-preprocess](https://github.com/HXC-code/plant-genome-preprocess-skill) |

## 能力边界

| 流程步骤 | 本工具 | 说明 |
| --- | :---: | --- |
| 1 数据收集并整理 | ✅ | 用户按来源自行下载（优先测序最完整、发表最新）；工具注册、路由、校验本地文件，不可用物种标注原因，不强制补全 |
| 2 数据预处理 | ✅ | 产出 `pep`/`cds`/`new.gff`/`lens` 四件套；含代表转录本筛选、基因清洗、命名统一、QC 与剔除留痕 |
| 3 DIAMOND 同源比对 | ✅ | 物种内 + 物种间两两比对，默认 `very-sensitive`、e-value `1e-5` |
| 4 WGDI 共线性 | ✅ | 点图与共线性块（`.png`/`.icl`） |
| 4′ Ks 分析 | ⚙️ 可选 | 默认关闭；`WGDI_KS=true`（Windows `-WithKs`）启用六个 Ks 阶段 |
| 5 系统发育树 | ❌ | 不在范围内 |
| 6 文献树比对与加倍节点核对 | ❌ | 人工分析，不在范围内 |
| 7 多重共线性与多重比对列表 | ❌ | 不在范围内 |
| 8 共线性圈图 | ❌ | 不在范围内 |

第 4 步的点图与共线性块是第 5–8 步的输入，但本工具**不执行、也不声称完成**这些后续步骤。
分析开始前应先向使用者说明上表边界。

## 安装

```powershell
cd code
python -m pip install -e ".[gbff]"     # 不处理 GBFF 时可省略 [gbff]
```

外部程序按需安装：WGDI 用 `bash bootstrap_wgdi.sh`（钉版 `wgdi==0.75`）创建独立环境；
DIAMOND/BLAST、MUSCLE、PAML `yn00`、Perl 加入 PATH 即可；PAL2NAL 由
`python code/tools/install_pal2nal.py` 从官方源安装（不随包分发）。
依赖预检：

```powershell
python -m pgci_preprocess doctor --project-root D:\plant_project --workflow full
```

## 快速开始

### 1. 创建项目并注册物种

```powershell
python -m pgci_preprocess init-project --project-root D:\plant_project

python -m pgci_preprocess init-project `
  --project-root D:\plant_project `
  --show-id 2848 --scientific-name "Morus notabilis" --order Rosales --family Moraceae `
  --genome D:\raw\Morus.genome.fasta.gz `
  --gff    D:\raw\Morus.gff3.gz `
  --cds    D:\raw\Morus.CDS.fasta.gz `
  --pep    D:\raw\Morus.Protein.faa.gz
```

`--genome` 可省略。四个角色文件按**内容**识别，不靠后缀猜测；同一角色有多个候选时
清单里写 `A|B`。GTF/GBFF 可先转换，或交给 `prepare`/`batch` 自动路由：

```powershell
python -m pgci_preprocess convert-annotation --input D:\raw\sp.gtf.gz --output-dir D:\plant_project\data\derived --prefix sp
```

### 2. 预处理 → 四件套

```powershell
python -m pgci_preprocess run      --project-root D:\plant_project
python -m pgci_preprocess validate --project-root D:\plant_project --refresh-links
python -m pgci_preprocess check-names --project-root D:\plant_project
```

`run` 会打印每个物种的染色体清单：

```text
[OK] Mun show_id=2848 Morus notabilis genes=... seqs=10 ready=yes reason=
    chromosomes 10 (id:genes) 1:4210  2:3980  3:4521 ...
    inventory: D:\plant_project\analysis_v2\Mun\Mun.chromosome_map.tsv
```

交接给 WGDI 的四个文件：

```text
analysis_v2\<Abbr>\<Abbr>.pep
analysis_v2\<Abbr>\<Abbr>.cds
analysis_v2\<Abbr>\<Abbr>.new.gff
analysis_v2\<Abbr>\<Abbr>.lens
```

### 3. 复核与人工裁决（SOP 要求）

```powershell
python -m pgci_preprocess chromosomes --project-root D:\plant_project --species Mun
python -m pgci_preprocess adjudicate  --project-root D:\plant_project --species Mun `
  --confirmed-n 10 --decision accept --reviewer "<复核人>" --rationale "<核型/文献依据>"
python -m pgci_preprocess adjudication-status --project-root D:\plant_project --species Mun
```

裁决记录绑定 QC 版本与 SHA-256，重跑后旧签字自动失效，不会被挪用到新结果上。

### 4. 严格验收

```powershell
python -m pgci_preprocess validate --project-root D:\plant_project --refresh-links --require-adjudication
python -m pgci_preprocess audit    --project-root D:\plant_project
```

### 5. 下游：DIAMOND + WGDI

```bash
bash 0check_data.sh --preprocess --refresh-links
bash run_homology.sh            # 或 1run_diamond.sh / 1run_blast.sh
bash 2run_pre_wgdi.sh
bash 3blast_to_wgdi.sh
python code/tools/generate_wgdi_configs.py --mode all
bash 6run_wgdi.sh               # 核心：.png + .icl
```

Windows 一键（1–6 步）：

```powershell
$env:HOMOLOGY_REQUIRE_ADJUDICATION = "true"
powershell -ExecutionPolicy Bypass -File .\run_all_pipeline.ps1
```

下游按 `wgdi_ready` 与裁决状态选择物种，跳过名单写入 `skipped_species.txt`。

#### 可选 Ks 分支

```bash
WGDI_KS=true python code/tools/generate_wgdi_configs.py --mode all
WGDI_KS=true bash 6run_wgdi.sh   # 追加 .ks/.bi.csv/.bi.c.csv/.bk.png/.kp.png/.pf.png
```

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all_pipeline.ps1 -WithKs
```

配置生成、`6run_wgdi.sh` 验收、Windows 总控、回执收集器读同一开关。
之所以不能"默认就跑"：上游 WGDI 的 `blockinfo`/`correspondence` 依赖 Ks 值，
缺 Ks 时 `correspondence` 会把所有区块过滤掉，产出一张空表。

## 输出与追溯

`analysis_v2/<Abbr>/`（权威输出）：

```text
<Abbr>.pep  <Abbr>.cds  <Abbr>.new.gff  <Abbr>.lens        # 四件套
<Abbr>.all.pep.fa  <Abbr>.all.cds.fa  <Abbr>.all.gff  <Abbr>.all.lens
<Abbr>.gene_map.tsv     # 代表基因 ID 与来源
<Abbr>.chr_map.tsv      # 逐序列保留/排除决定
<Abbr>.chromosome_map.tsv  # 染色体清单：编号/来源/长度/基因数
<Abbr>.qc.json          # 机器可读 QC（含输入 SHA-256）
<Abbr>.qc_record.md     # 人类可读 QC 记录（由 qc.json 确定性生成）
```

`clean_v2/<Abbr>/`（清洗记录）：

```text
<Abbr>.gene_map.tsv  <Abbr>.chr_map.tsv
<Abbr>.removed_list.tsv      # 全部剔除序列与原因
<Abbr>.chromosome_map.tsv
```

`new.gff` 固定七列：`序列号 / 基因ID / 起点 / 终点 / 链 / 序列内顺序号 / 原始基因ID`；
`lens` 固定三列：`序列号 / 长度 / 代表基因数`。

> 审计链是特性而非冗余：输入 SHA-256、`source-manifest`、绑定 QC 的裁决记录共同保证
> "删掉中间产物也能回溯到原始下载"。这是实验室 SOP 的交付依据。

## 关键规则（SOP 强制）

1. **代表转录本**：CDS 最长优先；CDS 等长时按稳定转录本 ID 字典序取第一条。
   **不以蛋白长度或序列出现顺序为依据**，保证可复现。
2. **只保留染色体挂载序列**：CLI 默认 `--sequence-mode chromosomes`，剔除未挂载
   scaffold/contig 与细胞器；无染色体证据的物种报 `no_chromosome_sequences` 并记录
   原因，不用 scaffold 充数。碎片化组装经确认后用 `--sequence-mode adaptive`。
3. **剔除必须留痕**：`removed_list.tsv` 逐条记录 id 与原因，禁止无记录过滤。
4. **Ks 默认关闭**：仅在使用者明确要求时启用。

序列筛选参数（CLI 默认）：

| 参数 | 默认 | 含义 |
| --- | ---: | --- |
| `--sequence-mode` | `chromosomes` | `adaptive` 为碎片化组装逃生门 |
| `--min-genes-per-sequence` | 50 | threshold 模式下每条序列最少代表基因数 |
| `--coverage` | 0.90 | coverage 模式的累计基因覆盖率 |
| `--max-sequences` | 200 | WGDI 序列数上限 |
| `--min-genes-per-chromosome` | 100 | 最终保留集合中每条序列的最少代表基因数 |

> Python API 的 `PipelineOptions` 默认值为中性值（`adaptive`、0），SOP 默认值在 CLI 层；
> 以命令行或 Skill 调用为准。

## `wgdi_ready` 与裁决状态

| `wgdi_ready` | 条件 | 用法 |
| --- | --- | --- |
| `yes` | 基因保留率 ≥ 70% 且序列数 ≤ 上限 | 进入标准下游 |
| `marginal` | 基因保留率 ≥ 50% | 结合组装与筛选记录复核 |
| `no` | < 50% | 默认不进下游；归档用 `--allow-not-ready` |

裁决状态：`ok` / `missing` / `stale` / `count_mismatch` / `rejected` /
`pending_adjustment` / `no_qc`。严格门禁下非 `ok` 一律 FAIL。

`wgdi_ready` 只反映当前筛选集的基因保留率，论文分析应同时记录 BUSCO、组装连续性
和倍性信息。

## 批量处理

```powershell
python -m pgci_preprocess batch `
  --project-root D:\plant_project `
  --source-status D:\plant_project\resources\source_status.tsv `
  --expected-count 20
```

顺序：`prepare → 可选 GFF 规范化 → run → validate`，路由决策写入
`resources/preprocess_batches/<tag>_routing.tsv`。批量跑完仍需逐物种复核染色体清单
并完成裁决，再走严格验收。

## 命名标准

物种缩写来自 `code/resources/species_abbreviations.tsv`（2543 条已发布物种的命名快照）；新物种按
拉丁名生成全局唯一缩写（属名首字母 + 种加词字母，3 位起，冲突换组合/扩位）。
基因 ID 形如 `Mun1g0001`。旧 `s<show_id>` 项目：

```powershell
python -m pgci_preprocess migrate-artifacts --project-root D:\plant_project --dry-run
python -m pgci_preprocess migrate-artifacts --project-root D:\plant_project
```

## 数据检查

`run`/`validate` 检查：CDS 是否为 3 的倍数、CDS 与蛋白翻译是否一致、内部终止与非法
字符、坐标/链/region/长度一致性、gene–mRNA–CDS 层级、ID 重复与冲突、混入多 taxon、
四件套 ID/序列集合闭合、压缩格式与 SHA-256。被剔除记录写入
`clean_v2/<Abbr>/<Abbr>.removed_list.tsv`。

## 科研使用范围

- 多倍体/高杂合数据应选择能区分等位 contig、homeolog 与亚基因组的注释来源，并在代表
  转录本筛选后人工抽查。
- 细胞器判定同时看 `organelle_inference.inferred` 与 `review_candidates`。
- 100 基因下限是机械过滤，可能删掉基因数少的真实小染色体，须结合 `chr_map.tsv` 复核。
- 严格染色体模式在 GWH 等"裸数字 seqid + 无 region 行"的简化注释上会判为无染色体证据，
  批量前建议先小样本实测命中率。
- 大规模 Ks 分析先在代表性小样本验证参数与耗时；WGD 定年与 Ks 峰解释由研究团队结合
  物种背景判断。

## 测试与验证

```powershell
cd code
python -m pytest -q          # 212 项
python -m ruff check src tests scripts tools
python -m mypy src --ignore-missing-imports
```

端到端冒烟回执（均由运行当时环境生成）：

```text
validation/smoke-report-3.7.0.json        # Ks 模式：4 组比较 × 8 产物
validation/smoke-report-3.7.0-core.json   # 核心模式：4 组比较 × 2 产物
```

发布门禁与环境明细见 [VALIDATION.md](VALIDATION.md)。

## 版本历史

| 版本 | 要点 |
| --- | --- |
| **v3.7.0** | 对齐整体流程 1–4 步；CDS 优先规则；默认严格染色体模式；Ks 改为可选；SOP 追溯文件名 + `qc_record.md` |
| [v3.6.2](https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.6.2) | 修复非 UTF-8 Windows 控制台打印染色体清单异常 |
| [v3.6.1](https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.6.1) | `run` 向调用方打印染色体清单 |
| [v3.6.0](https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.6.0) | 染色体清单、每染色体基因数下限、绑定 QC 的人工裁决硬门禁 |
| [v3.5.5](https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.5.5) | 线粒体推断、统一就绪门禁、产物完整性检查 |

完整变更见 [CHANGELOG.md](CHANGELOG.md)。

## 仓库关系

- `plant-genome-preprocess`：程序实现、测试、发布材料与主文档（本仓库）
- `plant-genome-preprocess-skill`：Agent 调用入口（Skill 名 `plant-wgdi-preprocess`）
  及同步的 `pipeline/` 运行副本，镜像检查保证一致

## 许可

项目代码使用 MIT License。WGDI、DIAMOND、BLAST+、MUSCLE、PAML、Perl 和 PAL2NAL
按其上游许可使用，详见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
