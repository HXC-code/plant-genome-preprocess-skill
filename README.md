# Plant Genome Preprocess Skill

`plant-genome-preprocess` 是 HXC 植物基因组规范化与质控工具的 Agent Skill。它把数据识别、预处理命令、质量检查、结果复核和可选下游分析组织成一套可执行工作流，适合在 Codex、WorkBuddy 等支持 Agent Skill 的环境中直接调用。

当前版本：**3.6.2**

- [程序仓](https://github.com/HXC-code/plant-genome-preprocess)
- [发布包](https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.6.2)
- [完整使用说明](pipeline/QUICKSTART.md)
- [验证说明](pipeline/VALIDATION.md)

## 仓库组成

```text
plant-genome-preprocess-skill/
├── SKILL.md
├── README.md
├── pipeline/
│   ├── code/                    预处理源码、测试和资源表
│   ├── validation/              Windows、Linux 和真实物种验证回执
│   ├── QUICKSTART.md            教师使用说明
│   ├── VALIDATION.md            发布门禁
│   ├── run_homology.sh          DIAMOND/BLAST 执行入口
│   ├── 0check_data.sh … 6run_wgdi.sh
│   ├── run_all_pipeline.ps1     Windows 分析总控
│   ├── bootstrap_wgdi.ps1/.sh   WGDI 环境准备
│   └── bin/wgdi                 可移植 WGDI 启动器
└── .github/                     CI 和镜像一致性检查
```

`SKILL.md` 是 Agent 调用入口。`pipeline/` 是与程序仓同步的完整运行副本，包含同一套源码、资源表、脚本和验证材料。两处内容由同步检查保持一致。

## Skill 能完成的工作

| 能力 | 结果 |
| --- | --- |
| 项目初始化 | 创建标准目录、库存表、物种缩写表和资源文件 |
| 输入识别 | 识别 GFF/GFF3、GTF、GBFF、genome、CDS、蛋白和压缩格式 |
| 单物种处理 | 注册物种并生成 WGDI 四件套 |
| 批量处理 | 读取状态表，完成文件路由、转换、预处理和验收 |
| 注释修复 | 补充 region、mRNA 和 gene/CDS 层级，记录异常 |
| 代表转录本选择 | 每个基因选择一个最长蛋白候选，并按 CDS 长度和稳定 ID 消歧 |
| 主序列选择 | 根据主染色体、基因数、覆盖率和序列上限生成 WGDI 序列集 |
| 染色体清单 | 输出最终序列数、每条序列代表基因数、长度和原始序列 |
| 人工裁决 | 记录复核人确认的染色体数，并绑定 QC 版本与 SHA-256 |
| 数据验收 | 检查 ID、坐标、链、长度、翻译、哈希、链接和命名一致性 |
| 追溯报告 | 生成 gene map、chromosome map、rejected、QC、index 和 manifest |
| 旧项目迁移 | 将 `s<show_id>` 命名迁移为物种缩写 |
| 来源冻结 | 构建并验证 role-level source manifest |
| 可选比较分析 | 执行 DIAMOND/BLAST、WGDI 共线性、Ks、block-Ks 和 Ks peaks |

## 两种工作模式

### 预处理模式

```text
init-project / prepare
        │
        ▼
run → 染色体清单 → 人工裁决 → 严格 validate
        │
        ▼
check-names → audit
```

主要结果：

```text
analysis_v2/<Abbr>/
├── <Abbr>.pep
├── <Abbr>.cds
├── <Abbr>.new.gff
└── <Abbr>.lens
```

同时生成 `<Abbr>.all.*`、gene map、chromosome map、QC 和 `clean_v2/<Abbr>/.rejected.tsv`。

### 比较模式

```text
预处理四件套
      │
      ▼
DIAMOND / BLAST
      │
      ▼
WGDI dotplot、共线性、Ks
      │
      ▼
block-Ks、Ks peaks
```

## 本地使用

仓库已在本地后，将 Skill 根目录放入 Codex 或 WorkBuddy 的 skills 目录，即可按自然语言调用。

Codex 常用目录：

```text
%USERPROFILE%\.codex\skills\plant-genome-preprocess
```

WorkBuddy 常用目录：

```text
%USERPROFILE%\.workbuddy\skills\plant-genome-preprocess
```

Skill 根目录包含 `SKILL.md`，完整运行工程位于 `pipeline/`。

## 调用示例

### 预处理与验收

```text
使用 $plant-genome-preprocess 将这些植物基因组数据清洗并验证为 WGDI 输入。
```

### 批量处理

```text
使用 $plant-genome-preprocess 读取 source_status.tsv，批量处理完整物种，逐项报告 QC、wgdi_ready 和失败原因。
```

### 完整比较

```text
使用 $plant-genome-preprocess 执行预处理、DIAMOND/BLAST、WGDI 共线性和 Ks 分析。
```

### 现有项目验收

```text
使用 $plant-genome-preprocess 检查这个项目的 WGDI 四件套、命名、QC、哈希和剔除记录。
```

### 旧项目迁移

```text
使用 $plant-genome-preprocess 将这个旧项目的 s<show_id> 目录、文件名和基因 ID 迁移为物种缩写。
```

## 直接运行工程

也可以在 `pipeline/` 中使用标准命令：

```powershell
cd pipeline\code

python -m pgci_preprocess init-project --project-root D:\plant_project
python -m pgci_preprocess run      --project-root D:\plant_project
python -m pgci_preprocess chromosomes --project-root D:\plant_project --species <code>
python -m pgci_preprocess adjudicate --project-root D:\plant_project --species <code> `
  --confirmed-n <N> --decision accept --reviewer "<复核人>" --rationale "<依据>"
python -m pgci_preprocess validate --project-root D:\plant_project `
  --refresh-links --require-adjudication
python -m pgci_preprocess check-names --project-root D:\plant_project
python -m pgci_preprocess audit    --project-root D:\plant_project
```

`run` 默认会丢弃代表基因数少于 100 的 WGDI 序列。最终序列数由人工结合 `chromosomes.tsv`、`chr_map.tsv` 和物种核型确认。直接 CLI 的门禁开关默认关闭，Skill 执行时必须显式使用 `--require-adjudication`，下游必须设置 `HOMOLOGY_REQUIRE_ADJUDICATION=true`。

详细输入格式、批量状态表、QC 检查和下游命令见 [pipeline/QUICKSTART.md](pipeline/QUICKSTART.md)。

## 运行时

| 用途 | 依赖 |
| --- | --- |
| 预处理 | Python 3.11+ |
| GBFF/GenBank | Biopython 1.83+ |
| WGDI | WGDI 0.75 |
| 同源比对 | DIAMOND 或 BLAST+ |
| WGDI 外部工具 | Bash、MUSCLE、PAML `yn00`、Perl |
| Ks 密码子比对 | PAL2NAL v14 |

本地环境可通过以下命令检查：

```powershell
python -m pgci_preprocess doctor `
  --project-root D:\plant_project `
  --workflow full
```

## 质量与验证

当前 3.6.2 版本包含：

- 202 项 Python 自动化测试；
- Ruff、Mypy 和镜像漂移检查；
- Windows 与 Linux/WSL Shell 语法检查；
- Windows 完整合成 WGDI 流程，四组比较每组检查 8 个关键产物；
- 染色体清单（染色体数 + 每条基因数）与可配染色体基因数下限（默认 100）；
- 人工裁决硬门禁（绑定 QC 版本+哈希，重跑作废）；
- 植物线粒体核心基因识别；
- 统一 readiness 读取入口；
- WGDI 产物结构完整性检查；
- 41 物种真实数据试点审查；
- NCBI RefSeq 拟南芥两物种预处理；
- 真实物种子集双向 WGDI 验证。

验证材料：

```text
pipeline/validation/smoke-report-3.6.2.json
pipeline/validation/smoke-report-3.5.5.json
pipeline/validation/smoke-report-3.5.4.json
pipeline/validation/smoke-report-3.5.3.json
pipeline/validation/smoke-report-3.5.2-linux.json
pipeline/validation/real-species-cross-3.5.2.json
pipeline/validation/real-species-3.5.2.json
```

## 科研使用建议

- **多倍体和高杂合数据**：选择能够区分等位 contig、homeolog 和亚基因组的注释来源，并在代表转录本筛选后人工抽查。
- **复杂细胞器**：同时检查 `organelle_inference.inferred` 和 `review_candidates`，对匿名、稀疏或证据不足的序列进行复核。
- **准备度判断**：将 `wgdi_ready` 与 BUSCO、组装连续性和物种倍性信息联合使用。
- **大规模 Ks**：先在代表性小样本验证参数和运行时间，再扩展到全蛋白组。
- **生物学解释**：预处理提供一致、可追溯的数据基础；WGD 定年、Ks 峰解释和物种演化结论由研究团队结合证据判断。

## 同步与发布

程序仓是代码与文档的同步来源。更新 Skill 镜像时，在程序仓执行：

```powershell
python code\tools\sync_skill_package.py `
  --pipeline-target D:\path\to\plant-genome-preprocess-skill\pipeline
```

Skill CI 在 push 和 pull request 时执行 pytest、Ruff、Mypy 和镜像一致性检查。

## 许可

Skill 与程序代码使用 MIT License。PAL2NAL、WGDI、DIAMOND、BLAST+、MUSCLE、PAML 和 Perl 按其各自上游许可使用，详见 [pipeline/THIRD_PARTY_NOTICES.md](pipeline/THIRD_PARTY_NOTICES.md)。
