# Plant Genome Preprocess Skill

本仓库提供 `plant-genome-preprocess` Agent Skill 和与程序仓同步的完整
`pipeline/` 运行包。Skill 面向植物基因组数据清洗、WGDI 输入生成、数据验收、
命名迁移，以及可选的 DIAMOND/BLAST + WGDI 比较工作流。

程序仓：

<https://github.com/HXC-code/plant-genome-preprocess>

发布包：

<https://github.com/HXC-code/plant-genome-preprocess/releases/tag/v3.5.4>

## 3.5.4 修复的问题

| # | 原问题 | Skill 中的修复 |
|---|---|---|
| 1 | 简化注释的叶绿体/线粒体序列可能进入核基因组 WGDI 集 | 对合成 region 做可审计的细胞器推断，QC 写入 `organelle_inference`；支持 `--organelle-pattern` 与 `--no-organelle-inference` |
| 2 | `wgdi_ready=no` 没有硬门禁 | `validate` 默认 FAIL，`--allow-not-ready` 仅用于归档；下游跳过未就绪或无记录物种 |
| 3 | 只有 GFF + CDS + PEP 时无法注册项目 | `init-project` 的 genome 和 assembly ID 改为可选 |
| 4 | 同源扫描误把 `bin/`、`code/` 等目录当物种 | `run_homology.sh` 与 Windows 总控排除基础设施目录 |
| 5 | 生成的 WGDI 脚本包含 Windows 反斜杠 | 自动转换为 Git Bash/Linux 可用的 `/d/...` 路径 |
| 6 | 临时文件清理失败导致成功任务返回非零 | 清理容错，并使用按进程区分的临时结果文件 |
| 7 | Windows `TMPDIR` 和 WGDI 解释器错误缺少预检 | 启动器规范化临时路径；`doctor` 通过 `wgdi_probe` 实际启动 WGDI 入口 |

## Skill 可以完成的工作

| 能力 | 内容 |
| --- | --- |
| 项目初始化 | 创建标准目录、inventory、资源表和项目结构 |
| 输入识别 | 从名称和内容识别 GFF/GFF3、GTF、GBFF、genome、CDS、蛋白和压缩格式 |
| 单物种处理 | 注册一个显式物种 bundle，并生成 WGDI 四件套 |
| 批量处理 | 发现物种、解析候选文件、路由输入、逐物种运行和校验 |
| 注释修复 | 补 region、补 mRNA、修复 gene/mRNA/CDS 层级和拆分基因 |
| 代表转录本选择 | 每个基因保留一个最长蛋白候选，并按 CDS 长度和稳定 ID 消歧 |
| 主序列选择 | 使用主染色体、基因数阈值、累计覆盖率和序列数上限生成 WGDI 序列集 |
| 数据验收 | 检查 ID、坐标、链、长度、翻译、哈希、兼容链接和命名一致性；`wgdi_ready=no` 默认拦截 |
| 追溯报告 | 生成 gene map、chromosome map、rejected、QC、index 和 manifest |
| 旧项目迁移 | 将 `s<show_id>` 命名迁移为物种缩写并刷新兼容链接 |
| 冻结来源 | 构建和验证 role-level source manifest |
| 可选完整链 | 运行 DIAMOND/BLAST、WGDI collinearity、Ks、block-Ks 和 Ks peaks |

## 两种运行模式

| 模式 | 执行内容 | 主要结果 |
| --- | --- | --- |
| 预处理模式 | `init-project` / `prepare` → `run` → `validate` → `check-names` | WGDI 四件套、QC、映射、剔除记录、index 和 naming audit |
| 比较模式 | 预处理结果 → DIAMOND/BLAST → WGDI 配置 → WGDI 分析 | 同源命中、dotplot、共线性、Ks、block-Ks 和 peaks |

预处理输出：

```text
analysis_v2/<Abbr>/
├── <Abbr>.pep
├── <Abbr>.cds
├── <Abbr>.new.gff
└── <Abbr>.lens
```

同时生成 `<Abbr>.all.*` 全部有效代表基因集、`.gene_map.tsv`、
`.chr_map.tsv`、`.qc.json` 和 `clean_v2/<Abbr>/.rejected.tsv`。

## Skill 支持的项目

| 项目形态 | 行为 |
| --- | --- |
| GFF/GFF3 + genome + CDS + 蛋白 | 执行完整清洗、校验和 WGDI 输入生成 |
| GTF + genome/CDS/蛋白 | 转换 GTF 后进入标准流程 |
| GBFF/GenBank | 提取 GFF3、CDS、蛋白，可选导出 genome |
| GFF + CDS + 蛋白，无 genome | 从 CDS 计算长度，并在需要时由转录本回填 region；合成 region 做细胞器推断，`init-project` 可无 genome 注册 |
| 多个同角色候选 | 在 inventory 中用 `A|B` 表示候选，按内容兼容性选择 |
| 旧 `s<show_id>` 项目 | 支持 dry-run、bundle 重命名、基因 ID 迁移和兼容链接刷新 |
| 现有完整项目 | 支持 `validate`、`check-names`、`audit` 和 source manifest 验收 |

注释支持 standard、transcript-only、simplified 5 列、gene-span 和
gene-CDS-direct 五种 GFF 结构。`.gz`、`.bz2`、`.xz` 和 `.zip` 按压缩魔数
识别。

## 仓库结构

```text
plant-genome-preprocess-skill/
├── SKILL.md
├── README.md
├── pipeline/
│   ├── code/                    # 与程序仓同步的完整代码树
│   ├── validation/              # Windows、Linux、真实物种回执
│   ├── QUICKSTART.md
│   ├── VALIDATION.md
│   ├── THIRD_PARTY_NOTICES.md
│   ├── run_homology.sh
│   ├── 0check_data.sh … 6run_wgdi.sh
│   ├── run_all_pipeline.ps1
│   └── bootstrap_wgdi.ps1 / .sh
└── .github/                     # Skill CI 和镜像检查
```

`pipeline/` 是程序仓的版本化镜像，包含预处理源码、测试、资源表、可选下游脚本、
WGDI 启动器和验证材料。Skill 不仅能生成命令，也能直接运行镜像中的完整工程。

## 安装

Codex：

```powershell
git clone https://github.com/HXC-code/plant-genome-preprocess-skill.git `
  "$env:USERPROFILE\.codex\skills\plant-genome-preprocess"
```

WorkBuddy：

```powershell
git clone https://github.com/HXC-code/plant-genome-preprocess-skill.git `
  "$env:USERPROFILE\.workbuddy\skills\plant-genome-preprocess"
```

Skill 根目录包含 `SKILL.md`，运行工程位于 `pipeline/`。

## 调用方式

预处理与验收：

```text
使用 $plant-genome-preprocess 将这些植物基因组数据清洗并验证为 WGDI 输入。
```

完整比较：

```text
使用 $plant-genome-preprocess 执行预处理、DIAMOND/BLAST、WGDI 共线性和 Ks 分析。
```

现有项目验收：

```text
使用 $plant-genome-preprocess 检查这个项目的 WGDI 四件套、命名、QC、哈希和剔除记录。
```

旧项目迁移：

```text
使用 $plant-genome-preprocess 将这个旧项目的 s<show_id> 目录、文件名和基因 ID 迁移为物种缩写。
```

## 运行时

| 用途 | 依赖 |
| --- | --- |
| 预处理 | Python 3.11+ |
| GBFF/GenBank | Biopython 1.83+ |
| WGDI | WGDI 0.75 |
| 同源比对 | DIAMOND 或 BLAST+ |
| WGDI 外部工具 | Bash、MUSCLE、PAML `yn00`、Perl |
| Ks 密码子比对 | PAL2NAL v14，从上游安装并校验 SHA-256 |

WGDI 环境：

```powershell
powershell -ExecutionPolicy Bypass -File .\pipeline\bootstrap_wgdi.ps1
```

```bash
bash pipeline/bootstrap_wgdi.sh
```

PAL2NAL：

```powershell
python pipeline\code\tools\install_pal2nal.py
```

Python 包安装：

```powershell
cd pipeline\code
python -m pip install -e ".[gbff]"
```

## 验证

当前 Skill 镜像版本为 **3.5.4**，包含：

- 172 项 Python 自动化测试；
- Ruff、Mypy 和镜像漂移检查；
- Windows 与 Linux/WSL 严格 LF shell 检查；
- Windows 完整合成 WGDI smoke（3.5.4，含 readiness 门禁正路径），四组比较每组 8 个关键产物；
- 41 物种真实数据试点审查修复验收：橡胶树/银白杨叶绿体 contig 被推断排除；
  严格验收 36 PASS / 5 FAIL（拦截 5 个 `wgdi_ready=no` 物种），`--allow-not-ready` 41 PASS；
- Linux/WSL 完整合成 WGDI smoke（3.5.2），四组比较每组 8 个关键产物；
- NCBI RefSeq `Arabidopsis thaliana` TAIR10.1 与
  `Arabidopsis lyrata` v1.0 的真实物种预处理；
- 真实物种子集 `Ath` chromosome 1 对 `Aly` chromosome 364 的双方向 WGDI
  执行回执。

验证材料：

```text
pipeline/validation/smoke-report-3.5.4.json
pipeline/validation/smoke-report-3.5.3.json
pipeline/validation/smoke-report-3.5.2-linux.json
pipeline/validation/real-species-cross-3.5.2.json
pipeline/validation/real-species-3.5.2.json
```

发布门禁见 [pipeline/VALIDATION.md](pipeline/VALIDATION.md)，完整使用说明见
[pipeline/QUICKSTART.md](pipeline/QUICKSTART.md)。

## 同步与发布

程序仓是代码与文档的权威来源。Skill 镜像由以下命令刷新和核对：

```powershell
python pipeline\code\tools\sync_skill_package.py `
  --pipeline-target pipeline
```

镜像同步忽略测试缓存、虚拟环境、editable-install 元数据和 `.wgdi-env`。
Skill CI 在每次 push / pull request 中执行 pytest、Ruff、Mypy 和镜像检查。

## 许可

Skill 与程序代码使用 [MIT License](LICENSE)。PAL2NAL、WGDI、DIAMOND、
BLAST+、MUSCLE、PAML 和 Perl 作为独立上游工具使用，其许可与代码仓许可分别适用。
第三方说明见
[pipeline/THIRD_PARTY_NOTICES.md](pipeline/THIRD_PARTY_NOTICES.md)。
