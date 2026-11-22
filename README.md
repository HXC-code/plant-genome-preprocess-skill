# 植物基因组 WGDI、DIAMOND 与共线性分析技能

**Plant Genome WGDI, DIAMOND and Comparative Analysis Skill**

`plant-wgdi-preprocess` 是植物比较基因组工作流的 Agent Skill。它按实验室《整体流程》覆盖数据整理、WGDI 四件套预处理、DIAMOND 同源比对和 WGDI 点图/共线性块，Ks 分析作为可选分支。

当前版本：**3.8.0**

- [程序仓](https://github.com/HXC-code/plant-genome-preprocess)
- [v3.8.0 Release](https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.8.0)
- [完整使用说明](pipeline/QUICKSTART.md)
- [实验室整体流程](pipeline/code/docs/WORKFLOW.md)
- [验证门禁](pipeline/VALIDATION.md)

Skill 调用名为 `plant-wgdi-preprocess`。首页标题使用完整分析范围，便于使用者理解它同时覆盖预处理、DIAMOND 比对和 WGDI 共线性分析。

## 交付范围

```text
1 数据收集与整理
      ↓
2 WGDI 四件套预处理
      ↓
3 DIAMOND 物种内与物种间比对
      ↓
4 WGDI 点图与共线性块
      └── Ks 可选分支
```

| 步骤 | 技能能力 | 主要结果 |
| --- | --- | --- |
| 1 | 注册本地下载文件、识别数据角色、批量路由和失败原因记录 | inventory、assembly map、routing tables |
| 1 | 输出染色体清单并取得人工裁决 | `chromosome_map.tsv`、adjudication JSON |
| 2 | 注释修复、代表转录本、序列 QC、细胞器处理和统一命名 | `pep/cds/new.gff/lens`、maps、removed list、QC |
| 3 | 运行 DIAMOND 或 BLAST 完成物种内与物种间比对 | BLAST hit tables、点图输入 |
| 4 | 生成并运行 WGDI 配置，检查核心产物 | `.png`、`.icl` |
| 4′ | 按需启用 Ks 分支 | `.ks`、blockinfo、correspondence、block-Ks、Ks peaks |

第 4 步产出的点图、共线性块与统计表可继续交给系统发育、文献核对、多重共线性和圈图分析使用。

## 仓库组成

```text
plant-genome-preprocess-skill/
├── SKILL.md
├── README.md
├── pipeline/
│   ├── code/                    完整预处理源代码、测试和资源表
│   ├── code/docs/WORKFLOW.md    实验室整体流程与实现映射
│   ├── validation/              核心与 Ks 模式验证回执
│   ├── QUICKSTART.md            教师使用说明
│   ├── README.md                程序仓首页镜像
│   ├── VALIDATION.md            发布门禁
│   ├── run_homology.sh          DIAMOND/BLAST 执行入口
│   ├── 0check_data.sh … 6run_wgdi.sh
│   ├── run_all_pipeline.ps1     Windows 分析总控
│   └── bootstrap_wgdi.ps1/.sh   WGDI 环境准备
└── .github/                     CI 与镜像一致性检查
```

`pipeline/` 是程序仓的完整版本化镜像，包含同一套源码、资源、测试和下游脚本。Skill 的差异集中在 `SKILL.md` 的执行纪律与范围说明。

## 安装

将包含 `SKILL.md` 与 `pipeline/` 的仓库根目录放入 Agent 的 skills 目录：

```powershell
xcopy /E /I . "$env:USERPROFILE\.codex\skills\plant-wgdi-preprocess"
```

WorkBuddy：

```powershell
xcopy /E /I . "$env:USERPROFILE\.workbuddy\skills\plant-wgdi-preprocess"
```

完整运行工程位于 `pipeline/`，不需要从程序仓另外复制代码。

## 调用示例

### 完整流程

```text
使用 plant-wgdi-preprocess 把 D:\raw 下的植物基因组整理成 WGDI 四件套，
完成染色体人工裁决、DIAMOND 同源比对和 WGDI 共线性分析。
```

### 只要预处理

```text
使用 plant-wgdi-preprocess 将这些数据整理成 WGDI 四件套，
逐物种报告染色体数、每条染色体基因数、wgdi_ready 和剔除原因。
```

### 启用 Ks 分支

```text
使用 plant-wgdi-preprocess 完成 DIAMOND 和 WGDI 共线性，
并明确启用 Ks 分支，输出 Ks 表、block-Ks 和 Ks 峰图。
```

### 现有项目验收

```text
使用 plant-wgdi-preprocess 检查现有项目的四件套、染色体清单、
人工裁决、命名、QC 和剔除记录。
```

## 手动运行

```powershell
cd pipeline\code

python -m pgci_preprocess init-project --project-root D:\plant_project
python -m pgci_preprocess run      --project-root D:\plant_project
python -m pgci_preprocess chromosomes `
  --project-root D:\plant_project --species <code>

python -m pgci_preprocess adjudicate `
  --project-root D:\plant_project --species <code> `
  --confirmed-n <N> --decision accept `
  --reviewer "<复核人>" --rationale "<依据>"

python -m pgci_preprocess validate `
  --project-root D:\plant_project `
  --refresh-links --require-adjudication

python -m pgci_preprocess check-names --project-root D:\plant_project
python -m pgci_preprocess audit --project-root D:\plant_project
```

下游：

```bash
bash run_homology.sh
bash 2run_pre_wgdi.sh
bash 3blast_to_wgdi.sh
python code/tools/generate_wgdi_configs.py --mode all
bash 6run_wgdi.sh
```

Ks 分支：

```bash
WGDI_KS=true python code/tools/generate_wgdi_configs.py --mode all
WGDI_KS=true bash 6run_wgdi.sh
```

Windows 总控：

```powershell
powershell -ExecutionPolicy Bypass -File .\pipeline\run_all_pipeline.ps1
powershell -ExecutionPolicy Bypass -File .\pipeline\run_all_pipeline.ps1 -WithKs
```

## 强制执行规则

1. **代表转录本**：CDS 最长优先；CDS 等长时按稳定 transcript ID 字典序；蛋白长度和文件出现顺序不参与选择。
2. **染色体模式**：CLI 默认 `--sequence-mode chromosomes`；无染色体证据时报告 `no_chromosome_sequences`。碎片化组装经确认后用 `--sequence-mode adaptive`。
3. **100 基因下限**：CLI 默认移除代表基因数少于 100 的 WGDI 序列；移除记录写入 `chr_map.tsv` 与 QC。
4. **剔除留痕**：全部剔除记录写入 `removed_list.tsv`。
5. **人工裁决**：复核人填写最终染色体数、reviewer 和 rationale；记录绑定 QC 版本与 SHA-256。
6. **严格门禁**：验收使用 `validate --require-adjudication`，下游设置 `HOMOLOGY_REQUIRE_ADJUDICATION=true`。
7. **Ks 按需开启**：默认完成点图与共线性；明确要求时使用 `WGDI_KS=true` 或 `-WithKs`。

## 关键输出

```text
analysis_v2/<Abbr>/          # 唯一输出层（不再写 clean_v2）
├── <Abbr>.pep
├── <Abbr>.cds
├── <Abbr>.new.gff
├── <Abbr>.lens
├── <Abbr>.gene_map.tsv
├── <Abbr>.chr_map.tsv
├── <Abbr>.chromosome_map.tsv
├── <Abbr>.removed_list.tsv
├── <Abbr>.qc.json
└── <Abbr>.qc_record.md
```

可选（`run --keep-all` 时才写，筛选前全量快照）：

```text
<Abbr>.all.pep.fa  <Abbr>.all.cds.fa  <Abbr>.all.gff  <Abbr>.all.lens
```

旧项目的 `clean_v2/` 仍可读；重跑后新产物只写 `analysis_v2`，旧层可整目录删除。

`new.gff` 使用七列：

```text
序列号  唯一基因ID  起点  终点  链方向  序列内顺序号  原始基因ID
```

`lens` 使用三列：

```text
序列号  序列长度  代表基因数
```

## 运行时

| 用途 | 依赖 |
| --- | --- |
| 预处理 | Python 3.11+ |
| GBFF/GenBank | Biopython 1.83+ |
| WGDI | WGDI 0.75 |
| 同源比对 | DIAMOND 或 BLAST+ |
| Ks 分支 | MUSCLE、PAML `yn00`、PAL2NAL v14、Perl |
| Shell 脚本 | Bash 或 Git Bash |

## 验证

```powershell
cd pipeline\code
python -m pytest -q
python -m ruff check src tests scripts tools
python -m mypy src --ignore-missing-imports
```

当前版本包含 213 项 Python 测试。验证回执：

```text
pipeline/validation/smoke-report-3.8.0.json        Ks 模式：4×8
pipeline/validation/smoke-report-3.8.0-core.json   核心模式：4×2
```

## 同步与版本

程序仓是代码、测试和文档的同步源。更新 `pipeline/`：

```powershell
python code\tools\sync_skill_package.py `
  --pipeline-target D:\path\to\plant-genome-preprocess-skill\pipeline
```

当前版：v3.8.0。v3.7.0 及更早版本保留为历史 Release，不改写原 tag。

## 许可

Skill 与程序代码使用 MIT License。第三方工具按其上游许可使用，详见 [pipeline/THIRD_PARTY_NOTICES.md](pipeline/THIRD_PARTY_NOTICES.md)。
