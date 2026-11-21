---
name: plant-genome-preprocess
description: 植物基因组预处理与可选 WGDI 比较工作流（完整工程版）。把 GFF/GFF3/GTF/GBFF/功能注释 CSV 与 genome/CDS/蛋白 FASTA 清洗成 WGDI 四件套，含命名校验、断点续跑、基因清洗与 QC；也可继续执行 DIAMOND/BLAST、共线性区块与 block-Ks/Ks-peak 分析。
---

# 植物基因组预处理（完整工程）

这个 Skill 是植物基因组规范化与质控工具的 **Agent 调用与分发层**。权威实现、
测试和发布流程位于程序仓；`pipeline/` 是其版本化完整镜像，不是另一套独立产品。
`pipeline/` 下可以直接运行源码、测试、资源表、可选下游脚本和固定版本 WGDI
启动入口。

预处理核心到 `pgci_preprocess run` 为止，只产出 WGDI 输入四件套。下游层是可选
模式；`run` 不会调用它。只有需要 DIAMOND/BLAST、共线性和 Ks 分析时才继续。

## 何时用

- 拿到植物基因组的注释（GFF/GFF3/GTF/GBFF）与序列（genome/CDS/蛋白），要转成 WGDI 的 `pep/cds/new.gff/lens`；
- 已有 WGDI 输入想复核：命名、长度、坐标、重复基因、兼容链接；
- 批量物种（几十到几百）需要断点续跑、逐物种 QC、失败原因可追溯；
- 需要把 `s<show_id>` 旧命名迁移成物种缩写（如 `Mun`）。

## 目录

```
pipeline/
├── QUICKSTART.md               # 面向使用者的完整说明
├── run_homology.sh             # DIAMOND/BLAST 共用实现
├── 0check_data.sh … 6run_wgdi.sh
├── run_all_pipeline.ps1        # Windows 分析总控（1–6）
├── bootstrap_wgdi.ps1 / .sh    # 建 WGDI 环境，并从官方站点安装 PAL2NAL
├── bin/wgdi                    # 可移植启动器
├── requirements-wgdi-validated.txt
├── VALIDATION.md               # 发布门禁
├── validation/                 # PASS 的 WGDI smoke 回执
└── code/                       # 预处理源码（唯一代码树）
    ├── src/pgci_preprocess/    # Python 包
    ├── scripts/ tools/ tests/ docs/
    ├── resources/species_abbreviations.tsv  # 2543 行物种缩写表
    └── pyproject.toml
```

## 物种缩写规则（当前策略）

- **3 位起**：属名首字母 + 种加词里任取 2 个字母。默认取种加词前两个（`Cercis chinensis` → `Cch`）。
- **冲突换组合**：`Cch` 已被别的物种占用时，依次尝试种加词内其他两字母组合（第 1+3 位 `ci` → `Cci`、第 1+4 位 `cn` → `Ccn` …）。
- **3 位占满再扩位**：所有两字母组合都被占用，就升到 4 位（属首 + 种加词 3 个字母，如 `Cchi`），仍冲突继续加长；极端情况才追加数字后缀。
- **表是快照**：`resources/species_abbreviations.tsv`（2543 行）是已发布物种的既有命名，表里有的 show_id 一律沿用；表里没有的新物种按上面的规则现算，并自动避开表内所有 code。
- 已注册的物种 code 不会改变，历史产物目录与 gene ID 不受影响。

## 3 分钟上手

```powershell
cd pipeline\code
# Python ≥ 3.11；处理 GBFF 时安装 biopython
python -m pip install -e ".[gbff]"

# 1) 建项目（自动带上物种缩写表与格式例外表，无需手工拷贝）
python -m pgci_preprocess init-project --project-root D:\my_project
# 2) 把物种目录放进 D:\my_project\data\downloads\<show_id>_<Order>_<Family>(<Species>)
#    编辑 resources\download_inventory.tsv 四列 gff/genome/cds/pep（留空表示没有）

# 3) 预处理 + 校验
python -m pgci_preprocess run      --project-root D:\my_project
python -m pgci_preprocess validate --project-root D:\my_project --refresh-links
python -m pgci_preprocess check-names --project-root D:\my_project
python -m pgci_preprocess doctor   --project-root D:\my_project
```

单物种可用一条命令注册（文件自动镜像/转换/校验；`--genome` 可省，只有 GFF+CDS+PEP 三件套也能注册）：

```powershell
python -m pgci_preprocess init-project --project-root D:\my_project `
  --show-id 2848 --scientific-name "Morus notabilis" --order Rosales --family Moraceae `
  --genome GWH...genome.fasta.gz --gff GWH...gff.gz `
  --cds GWH...CDS.fasta.gz --pep GWH...Protein.faa.gz --species-names resources\species_abbreviations.tsv
```

## 支持哪些后缀

注释（按**内容**判定，后缀只作提示，不是白名单；`.gz/.bz2/.xz/.zip` 按压缩魔数识别，后缀写错也能读）：

| 后缀 | 处理 |
| --- | --- |
| `.gff` `.gff3`（含 `.gz`） | 直接解析，支持 standard / transcript_only / simplified5 / gene_span / gene_cds_direct 五种结构 |
| `.gtf`（含 `.gz`） | 转 GFF3（gene/transcript/exon/CDS/start_codon/stop_codon） |
| `.gbff` `.gbk` `.gb` `.genbank`（含 `.gz`） | 提取编码区 GFF3 + CDS/PEP（可顺带导 genome）；需要 biopython |
| `.annotations.csv.gz` | 识别为无表头 7 列功能注释（GO/NR/SP/TR），规范化为元数据表，**不当坐标输入** |
| `.csv` / `.csv.gz` | 仅用于功能注释规范化，缺坐标表头会被拒绝 |

序列（FASTA 家族，含 `.gz/.bz2/.xz/.zip`）：`.fa` `.fasta` `.fna` `.faa` `.fas` 及其压缩形式，以及常见组合 `.genome.fasta.gz`、`.genome.fa.gz`、`_genomic.fna.gz`、`.CDS.fasta.gz`、`.cds.fa(.gz)`、`.Protein.faa.gz`、`.protein.fa.gz`、`.pep.fa`、`.pro.fa`、`.RNA.fasta.gz`。其它扩展名只要内容是 FASTA 也会被读取，角色由文件名语义与首条记录内容共同判定。

命名不是硬编码：代码用“后缀 + FASTA 首记录预览 + 角色打分”决定 genome/cds/pep，避免把 CDS 当 genome、把 RNA 当蛋白。

## 基因清洗做了什么

一次 `run` 会完成：

1. **注释解析与层级修复**：缺 region 时用 genome 回填；缺 mRNA 时由 CDS 合成转录本；`gene`/`mRNA`/`CDS` 层级不一致时合并或拆分并记账。无 genome、无显式 region 的简化注释还会对合成 region 做**细胞器推断**（标识符信号或高密度紧凑基因集），证据写入 QC `organelle_inference`；`--organelle-pattern <正则>` 可强制排除指定合成序列，`--no-organelle-inference` 可关闭。
2. **基因模型筛选**：CDS 长度非 3 倍数、含内部终止、翻译与 `/translation` 不一致、缺 `protein_id`、坐标越界或链方向错误的记录被剔除。
3. **序列选择**：按 `--min-genes-per-sequence`（默认 50）、`--coverage`（默认 0.90）、`--max-sequences`（默认 200）挑主序列，保证 WGDI 输入既不太碎也不丢太多基因。
4. **去重与改名**：同一 `gene_id` 跨物种重复、同一蛋白对应多基因、切分基因同名等都会消歧。
5. **产出可审计**：`<Abbr>.rejected.tsv` 逐条记录被剔除的 id 与原因（如 `protein_id_missing`、`cds_length_not_multiple_of_3`、`internal_stop`），`<Abbr>.qc.json` 记录计数与阈值。

## 极端情况（都有明确行为）

| 情况 | 行为 |
| --- | --- |
| 目录为空 | 该物种不进入索引，`run` 汇总里不出现；`doctor` 计数为 0 |
| 只有 genome、没有注释 | `run` 报 `annotation_missing`，不会凭空造 GFF |
| 只有注释/CDS/蛋白、没有 genome | 允许；`init-project` 可无 `--genome` 注册；长度由 CDS 计算，region 缺失时用转录本回填，合成 region 会做细胞器推断 |
| 简化 GFF 混入叶绿体/线粒体 contig | 推断排除并计入 `n_organelle_genes`；拿不准的高密度 region 列入 `review_candidates` 供人工复核，不静默排除 |
| `wgdi_ready=no` 的碎片化物种 | `validate` 默认 FAIL（`wgdi_not_ready_low_coverage`），下游跳过；仅数据完整性复核时显式 `--allow-not-ready` |
| 只有 GBFF | 自动转 GFF3/CDS/PEP；无 CDS 的“基因组版 GBFF”会明确报错 |
| 只有 GTF | 转 GFF3 后继续；序列长度需 genome 回填 |
| 同目录多套重复（GCF + GWH） | 清单里用 `A|B` 选择器，或直接写唯一文件名；打分平局会报歧义而不是乱选 |
| 未压缩 + 压缩混放 | 两者都可读；优先内容判定 |
| 文件名无 GCA/GCF 号（GWH/figshare/手写） | accession 记 `NA`，有 `assembly_accessions.tsv` 时以其为准；只有多个不同 accession 冲突才失败 |
| 旧 `s<show_id>` 产物 | 再次 `run` 自动改名成缩写，索引与兼容链接同步 |
| GFF 写错列/坐标/链 | `validate` 逐物种报 `invalid_*`，`--refresh-links` 可修复链接 |

## 校验与命名

- `validate`：逐物种核对 pep/cds/gff/lens 是否存在、基因数是否一致、坐标/链/长度是否合法、兼容链接是否指向同一内容；`--skip-hashes` 可跳过内容哈希（大批量提速）。QC 中 `wgdi_ready=no` 默认判 FAIL；仅做数据完整性归档时显式加 `--allow-not-ready`。
- `check-names`：缩写表 ↔ 清单 ↔ 索引 ↔ 产物目录四者一致性；输出 `resources/naming_audit.tsv`。
- `doctor`：Python 版本、biopython、WGDI/BLAST/DIAMOND、清单表头、缩写表覆盖、输入文件存在性；`--workflow wgdi|full` 还会实际启动 WGDI 入口（`wgdi_probe`），暴露“存在但跑不起来”的环境问题。

## 调用规则（v3.5.4）

代理执行本 skill 时必须遵守：

1. **不得静默绕过准备度门禁**：不要为了凑 PASS 默认加 `--allow-not-ready` 或 `HOMOLOGY_ALLOW_NOT_READY=true`。先向用户报告哪些物种 `wgdi_ready=no`（保留率多少），由用户决定是否仅做完整性复核。
2. **报告细胞器推断结果**：`run` 后检查 QC 的 `organelle_inference.inferred` 与 `review_candidates`，向用户列出被排除的序列和待人工复核的高密度 region；不得把 `review_candidates` 当作已排除。
3. **报告下游跳过**：`run_homology.sh` / `run_all_pipeline.ps1` 跳过的物种记录在 `skipped_species.txt`，汇总时必须一并报告，不得只报成功比较。
4. **WGDI 环境先探针**：比较模式开始前用 `doctor --workflow full` 确认 `wgdi_probe` PASS；WGDI 必须用 `bootstrap_wgdi` 建的专用环境或 `WGDI_PYTHON` 指定的解释器（Python 3.13 + wgdi 0.75 已验证），不要复用依赖不明的 base/其他 conda 环境。
5. **旧产物重跑**：3.5.4 的细胞器推断会改变结果；QC 版本不是当前版本的旧 bundle 会被 `run` 自动重算，不要手工改 QC 凑版本。

## 测试

```powershell
cd pipeline\code
python -m pytest -q                     # 172 个用例；临时目录受限时加 --basetemp=<short-path>
```

测试覆盖：多格式项目、GBFF 往返、GTF 转换、命名迁移、兼容链接、校验失败路径、
版本一致性、共享脚本契约、发布包和 Skill 镜像。

## 交接与复现

- 发布包结构和两种使用模式见 `pipeline/QUICKSTART.md`。
- CI 在 push / PR 时跑 `pytest`、`ruff`、`mypy`，并检查 `pipeline/` 与 canonical 仓库一致；代码变更先在 canonical 仓库落地再同步。
- WGDI 使用固定 `wgdi==0.75` 上游包；`bin/wgdi_entry.py` 只维护必要的兼容补丁。
- PAL2NAL v14 不随 Skill 或发布 ZIP 分发。`bootstrap_wgdi` 会运行
  `code/tools/install_pal2nal.py`，从官方站点下载并校验 SHA-256。
- `validation/smoke-report-3.5.4.json` 是 Windows 四组比较 PASS 回执（含 readiness 门禁正路径）；
  `smoke-report-3.5.3.json` 是上一版 Windows 四组比较 PASS 回执；
  `smoke-report-3.5.2-linux.json` 是 WSL/Linux 四组比较 PASS 回执。
- `validation/real-species-3.5.2.json` 记录真实 RefSeq 拟南芥两物种的清洗、
  校验、命名、QC 和真实种间子集 WGDI 验证。
- 预处理与下游比较是两个明确模式，不应把“四件套就绪”解释成 WGDI 已完成。
- 所有中间产物可删：`analysis_v2/<Abbr>/` 是唯一权威输出，`clean_v2` 是清洗记录，根下 `<Abbr>/` 只是兼容链接。
