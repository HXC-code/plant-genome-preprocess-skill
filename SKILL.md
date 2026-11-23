---
name: plant-wgdi-preprocess
description: 植物比较基因组的本地数据整理、WGDI 输入构建、DIAMOND 同源比对和 WGDI 点图/共线性分析（实验室流程第 1–4 步，Ks 可选）。适用于跨机器重建项目路径、代表转录本筛选、染色体人工裁决、两两同源比较和可追溯结果交付。
---

# 植物基因组 WGDI 与 DIAMOND 分析

本 skill 实现实验室《整体流程》第 1–4 步：

1. 本地数据整理与配对；
2. WGDI 四件套预处理；
3. DIAMOND 物种内和物种间同源比对；
4. WGDI 点图与共线性块。

Ks 是第 4 步的可选分支。系统发育树、已发表树核对、多重共线性和圈图分别属于
第 5–8 步，不在本 skill 范围内。

分析开始前必须向使用者声明这个边界。第 4 步的 `.png` 和 `.icl` 可以作为后续分析
输入，但不能把“四件套就绪”或“点图已生成”表述为已完成后续 WGD 推断。

## 不可违反的规则

1. 代表转录本按 **CDS 长度最长优先，等长时按稳定转录本 ID 字典序**选择。
   不使用蛋白长度、文件顺序或任意排序替代。
2. 任何过滤和删除都必须写入 `removed_list.tsv` 或对应的 `chr_map.tsv` / QC；
   禁止无记录的静默过滤。
3. 默认使用严格染色体模式。没有染色体挂载证据时报告
   `no_chromosome_sequences`；只有使用者确认是碎片化组装后，才可使用
   `--sequence-mode adaptive`。
4. 人工染色体裁决绑定 QC 版本与 SHA-256。预处理重跑后，旧签字自动失效。
5. Ks 默认关闭，只有使用者明确要求时才启用。
6. 新机器必须先重新扫描路径。不要把当前机器的绝对路径、虚拟环境或其他用户的
   Linux 安装路径写死到项目、配置或代码中。
7. 不手工删除结果文件。需要瘦身时使用 `prune`，否则索引和校验会失配。

## 项目结构

推荐布局：

```text
<PROJECT_ROOT>/
├── data/                       # 原始数据，只读
├── metadata/                   # 扫描索引、映射、裁决、校验报告
├── resources/                  # metadata/ 的兼容别名（联接或回退目录）
├── wgdi_input/
│   ├── <Abbr>/                 # 四件套与追溯文件
│   ├── blast/                  # DIAMOND 数据库、结果、设置和日志
│   └── config/                 # 每个比较一份扁平 <pair>.conf
└── wgdi_output/
    ├── results/                # .png、.icl 和可选 Ks 产物
    ├── logs/                   # 每个比较的运行日志
    └── state/                  # 阶段断点与 Ks 中间文件
```

`metadata/` 是首选位置；读取代码依次回退 `metadata/`、`resources/` 和项目根。
因此即使系统不允许创建兼容联接，流程也应继续工作。

## 步骤 0：新机器先扫描

`discover` 不需要先安装 Python 包，任意 Python ≥ 3.11 即可运行。

Linux / WSL：

```bash
python3 /path/to/plant-wgdi-preprocess/pipeline/discover.py \
  /path/to/raw-data \
  --project-root /path/to/project
```

Windows PowerShell：

```powershell
python <skill>\pipeline\discover.py E:\raw_data --project-root D:\plant_project
```

如果数据已经在项目的 `data/` 下，可省略数据目录：

```bash
python3 /path/to/plant-wgdi-preprocess/pipeline/discover.py
```

扫描会：

- 按文件内容识别注释、genome、CDS 和蛋白；
- 自动配对同一套组装；
- 判断逐物种染色体证据；
- 写入 `metadata/machine_paths.json`、`download_inventory.tsv`、
  `sequence_mode.tsv` 和 `discovery_report.tsv`；
- 写入项目指针，使后续命令可以省略 `--project-root`。

指针优先位于 `pipeline/code/.pgci-machine.json`，只读安装时回落到用户主目录。
换电脑、换盘或数据搬家后必须重新扫描。

### WGDI Python 环境

WGDI 启动器只从以下来源选择解释器：

1. 显式指定的 `WGDI_PYTHON`；
2. 当前激活的 conda 环境；
3. 项目内 `.venv` / `.wgdi-env`；
4. PATH 上的 `python3` / `python`。

它不会扫描其他用户主目录。新 Linux 机器应激活目标环境后执行 `discover`，
或者由使用者明确提供：

```bash
export WGDI_PYTHON=/confirmed/path/to/python
```

启动器和环境检查可用以下命令确认：

```bash
bash /path/to/plant-wgdi-preprocess/pipeline/bin/wgdi --version
```

## 步骤 1：建立项目和登记输入

安装预处理程序：

```bash
cd /path/to/plant-wgdi-preprocess/pipeline/code
python3 -m pip install -e ".[gbff]"
```

普通预处理不需要 Biopython；GTF、GBFF/GenBank 转换和对应测试需要。

建立项目：

```bash
python3 -m pgci_preprocess init-project --project-root /path/to/project
```

登记一个拥有显式路径的物种：

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

批量数据应放在：

```text
data/<show_id>_<Order>_<Family>(<Species>)/
```

旧布局 `data/downloads/` 仍可读取。批量清单为
`metadata/download_inventory.tsv`；同一角色有多个候选时写 `A|B`，由程序校验并
选择唯一兼容文件。可使用 `prepare` 和 `batch` 进行批量路由，具体参数以
`python -m pgci_preprocess prepare --help` / `batch --help` 为准。

## 步骤 2：预处理、人工裁决和校验

按扫描结果逐物种选择模式：

```bash
python3 -m pgci_preprocess run --sequence-mode from-map
```

也可以显式指定：

```bash
python3 -m pgci_preprocess run --sequence-mode chromosomes
python3 -m pgci_preprocess run --sequence-mode adaptive
```

`adaptive` 只用于经人工确认的碎片化组装。为替换已有有效结果，必须显式增加
`--overwrite`。

查看染色体或 WGDI 序列：

```bash
python3 -m pgci_preprocess chromosomes --species <code>
```

人工确认最终数量：

```bash
python3 -m pgci_preprocess adjudicate \
  --species <code> \
  --confirmed-n <N> \
  --decision accept \
  --reviewer "<实际复核人>" \
  --rationale "<核型或文献依据>"
```

执行质量门禁：

```bash
python3 -m pgci_preprocess validate --require-adjudication
python3 -m pgci_preprocess check-names
python3 -m pgci_preprocess audit
```

单物种产物位于 `wgdi_input/<Abbr>/`，四件套为：

```text
<Abbr>.pep
<Abbr>.cds
<Abbr>.new.gff
<Abbr>.lens
```

追溯文件包括 `chromosome_map.tsv`、`gene_map.tsv`、`removed_list.tsv`、
`chr_map.tsv`、`qc.json` 和 `qc_record.md`。

默认不生成筛选前的 `.all.*` 快照。需要保留时使用 `run --keep-all`；
已有项目可用 `prune --all-snapshots --dry-run` 先评估，再去掉 `--dry-run`。

## 步骤 3：DIAMOND 同源比对

```bash
cd /path/to/plant-wgdi-preprocess/pipeline
bash run_homology.sh
```

也可使用兼容入口：

```bash
bash 1run_diamond.sh
bash 1run_blast.sh       # 显式切到 BLAST+
```

默认参数和输出：

- DIAMOND `very-sensitive`；
- e-value `1e-5`；
- 物种内和物种间有序两两比较；
- 数据库与设置：`wgdi_input/blast/`；
- 命中表：`wgdi_input/blast/<A>_pep/<A>_<B>.blast`；
- 物种映射：`metadata/species_map.tsv`；
- 跳过目录的物种：`metadata/skipped_species.txt`。

只有 `wgdi_ready=yes` 或 `marginal`，且满足当前人工裁决门禁的物种才继续运行。
不要以“文件存在”替代 readiness 和 QC 检查。

## 步骤 4：WGDI 点图和共线性

先生成配置，再执行：

```bash
bash 2gen_wgdi_conf.sh all
bash 3run_wgdi.sh
```

可选模式：

```bash
bash 2gen_wgdi_conf.sh different
bash 3run_wgdi.sh different
```

配置由 `metadata/species_map.tsv` 生成，使用项目相对路径，不复制四件套或 BLAST。
核心模式每个比较产出：

```text
wgdi_output/results/<pair>.png
wgdi_output/results/<pair>.icl
```

`3run_wgdi.sh` 按阶段记录指纹并断点续跑。配置内容、BLAST、GFF/lens、
CDS/PEP 或 WGDI 启动器发生变化时，相应阶段会自动重建。

### Ks 可选分支

只有使用者明确要求时才执行：

```bash
WGDI_KS=true bash 2gen_wgdi_conf.sh all
WGDI_KS=true bash 3run_wgdi.sh
```

追加产物：

```text
.ks
.bi.csv
.bi.c.csv
.bk.png
.kp.png
.pf.png
```

Ks 还要求 MUSCLE、PAML `yn00`、Perl 和 PAL2NAL v14。缺少这些工具时必须
明确报告缺失项，不能把核心模式的结果冒充为 Ks 结果。

## Windows 总控

Windows 用户可让总控依次执行第 1–3 步：

```powershell
powershell -ExecutionPolicy Bypass -File <skill>\pipeline\run_all_pipeline.ps1
```

启用 Ks：

```powershell
powershell -ExecutionPolicy Bypass -File <skill>\pipeline\run_all_pipeline.ps1 -WithKs
```

常用只读检查：

```powershell
powershell -ExecutionPolicy Bypass -File <skill>\pipeline\run_all_pipeline.ps1 -CheckOnly
powershell -ExecutionPolicy Bypass -File <skill>\pipeline\run_all_pipeline.ps1 -DryRun
```

## 支持的输入

文件角色以内容为主，后缀仅为提示。`.gz`、`.bz2`、`.xz` 和单成员 `.zip`
按压缩魔数识别。

| 类型 | 支持内容 |
| --- | --- |
| 注释 | GFF/GFF3、GTF、GBFF/GenBank、七列功能注释 CSV/TSV |
| GFF 结构 | standard、transcript_only、simplified5、gene_span、gene_cds_direct |
| 序列 | genome、CDS、蛋白 FASTA |
| 非标准后缀 | 内容符合格式时仍可识别；`.lens` 仅作为输出，不作为输入 |

ZIP 只能包含一个非空数据文件。同一目录存在多套候选时，清单中的 `A|B` 由评分和
兼容性校验处理；平局或无法配对时必须报错，不能任意猜测。

## 关键边界情况

| 情况 | 行为 |
| --- | --- |
| 无染色体挂载证据 | `no_chromosome_sequences`；人工确认后才可用 adaptive |
| 只有 genome，没有注释 | `annotation_missing`，不凭空生成 GFF |
| 有注释/CDS/PEP，没有 genome | 允许；由 CDS 和转录本坐标计算长度、回填 region |
| 只有 GTF 或 GBFF | 转换后继续 |
| 多套重复候选 | 评分选择；歧义时停止并要求复核 |
| 细胞器混入 | 自动推断并记录 `review_candidates` |
| 旧 `analysis_v2/` | 用 `migrate-layout --dry-run` 预演后迁移 |
| 旧 `s<show_id>` 名称 | 用 `migrate-artifacts` 幂等迁移 |

## 测试与发布检查

```bash
cd /path/to/plant-wgdi-preprocess/pipeline/code
python3 -m pytest -q
python3 -m ruff check src tests tools ../discover.py
python3 -m mypy src
```

当前测试数为 259。Shell 入口还应通过 `bash -n`；Windows 总控应通过
PowerShell parser 或 `-DryRun`。完整发布门禁见
[`pipeline/VALIDATION.md`](pipeline/VALIDATION.md)。

## 交接

- 流程范围与步骤映射：[`pipeline/code/docs/WORKFLOW.md`](pipeline/code/docs/WORKFLOW.md)
- 快速上手：[`pipeline/QUICKSTART.md`](pipeline/QUICKSTART.md)
- 用户主页：[`README.md`](README.md)
- WGDI 启动器：[`pipeline/bin/wgdi`](pipeline/bin/wgdi)
- 配置生成器：[`pipeline/code/tools/generate_wgdi_configs.py`](pipeline/code/tools/generate_wgdi_configs.py)
- WGDI 执行器：[`pipeline/3run_wgdi.sh`](pipeline/3run_wgdi.sh)
