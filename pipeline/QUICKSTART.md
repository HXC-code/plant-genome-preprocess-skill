# 快速上手：从自己的原始下载数据到 WGDI 输入

本包只包含**代码**：预处理程序、可选的下游比对/WGDI 执行层和 WGDI 启动器。
不含任何原始数据、也不含已生成的结果目录。
拿到包的人只要有 Python 和少量外部程序，就能从自己的数据生成 WGDI 能直接跑的
`<Abbr>.pep / .cds / .new.gff / .lens`，并按缩写命名。

---

## 0. 包里有什么

```text
wgdi_pipeline/
├── QUICKSTART.md                  # 本页
├── bootstrap_wgdi.ps1 / .sh        # 一键创建 WGDI 运行环境（.wgdi-env）
├── bin/wgdi                        # 可移植启动器（自动找 .wgdi-env / PATH）
├── run_homology.sh                 # DIAMOND/BLAST 共用实现
├── 0check_data.sh … 6run_wgdi.sh   # 兼容入口
├── run_all_pipeline.ps1            # Windows 分析总控（1–6）
├── requirements-wgdi-validated.txt # 可选：WGDI 依赖钉版清单
├── VALIDATION.md                   # 发布和端到端验证门禁
├── THIRD_PARTY_NOTICES.md          # 外部工具许可说明
└── code/                           # 预处理源码（唯一代码树）
    ├── src/pgci_preprocess/        # Python 包
    ├── scripts/ tools/ tests/ docs/
    ├── resources/species_abbreviations.tsv   # 物种缩写表（关键！）
    └── pyproject.toml
```

`code/resources/species_abbreviations.tsv` 决定所有产物名，必须随包携带。
命名规则（3 位起）：属名首字母 + 种加词 2 个字母（默认前两位，如 `Cercis chinensis` → `Cch`）；被占用则换种加词内其他两字母组合（1+3 位 → `Cci`、1+4 位 → `Ccn`…）；3 位全占满扩到 4 位（`Cchi`），仍冲突继续加长。表里已有的 show_id 一律沿用表内 code；表外新物种按拉丁名现算并避开表内全部已有 code。只有既无表内记录又完全没有拉丁名时，才会保留旧式 `s<show_id>` 兼容命名，新数据不会出现这种情况。

---

## 1. 准备 Python 与外部程序

| 用途 | 需要的东西 | 说明 |
| --- | --- | --- |
| 预处理 | Python ≥ 3.11 | 纯标准库即可；要转 GBFF 才需要 `biopython>=1.83` |
| WGDI | Python + `wgdi==0.75` | 用 `bootstrap_wgdi` 或 requirements 安装；只有固定上游来源 |
| PAL2NAL | 官方 v14 下载 | `bootstrap_wgdi` 默认安装；也可单独运行 `python code/tools/install_pal2nal.py`；发布包不捆绑 |
| 同源比对 | `diamond` 或 `blastp` + `makeblastdb` | 装好后在 PATH 里即可 |
| WGDI 依赖 | `muscle`、`yn00`（PAML）、`perl` | 同上，PATH 或 `-BioinfoBin` 指定目录 |
| Windows 执行 shell 脚本 | Git Bash | `run_all_pipeline.ps1` 会自动找 `git.exe` 同级的 `bash.exe` |

安装预处理包（在仓库根目录执行，只需一次）：

```powershell
cd code
python -m venv .venv
# 激活 .venv 后：
python -m pip install -e ".[gbff]"
```

不处理 GBFF/GenBank 时可以省略 `[gbff]`。外部程序按上表安装并加入 PATH，或用 `doctor`/总控的参数显式指定。

---

## 2. 建 WGDI 运行环境（只做一次）

```powershell
powershell -ExecutionPolicy Bypass -File .\bootstrap_wgdi.ps1
# 指定解释器 / 环境目录：
# powershell -ExecutionPolicy Bypass -File .\bootstrap_wgdi.ps1 -PythonExe C:\Python311\python.exe -EnvDir D:\wgdi-venv
```

Linux / macOS / Git Bash：

```bash
bash bootstrap_wgdi.sh
```

`bin/wgdi` 找解释器的顺序：`$WGDI_PYTHON` → `.wgdi-env` → `.wgdi_python`（本机提示文件）→ PATH。
**换机器不需要改代码。**

---

## 3. 建一个自己的项目目录（数据与结果都放这里）

代码和项目分开，一个项目一个目录：

```text
my_project/
├── data/downloads/{show_id}_{Order}_{Family}({Species})/   # 每个物种一个目录
├── resources/                                              # 清单与索引（程序生成）
├── analysis_v2/<Abbr>/                                     # 预处理产物
└── clean_v2/<Abbr>/
```

每个物种最终要有四个文件（原始名不限，程序按内容识别角色）：

- 基因组蛋白 FASTA（pep）
- 编码序列 FASTA（cds）
- 注释 GFF3 / GTF
- 长度表（lens，可由预处理生成）

你手里的原始下载可以是这些形式；**后缀只是提示，真正决定处理方式的是文件内容**：

| 你有的文件 | 处理方式 |
| --- | --- |
| `.gff` / `.gff3`（含 `.gz`） | 直接解析；自动识别 standard、transcript_only、simplified5、gene_span、gene_cds_direct 五种结构 |
| `.gtf`（含 `.gz`） | 转成 GFF3 后继续 |
| `.gbff` / `.gbk` / `.gb` / `.genbank`（含 `.gz`） | 提取编码 GFF3 + CDS + PEP，可顺带导出 genome（需 `biopython>=1.83`） |
| `.gz` / `.bz2` / `.xz` / `.zip` | 按压缩魔数识别，后缀写错或缺失也能读 |
| `.fa` / `.fasta` / `.fna` / `.faa` / `.fas` 及 `.genome.fasta.gz`、`.CDS.fasta.gz`、`.Protein.faa.gz` 等组合名 | 用文件名语义 + FASTA 首条记录内容判定 genome / CDS / protein 角色 |
| `*.annotations.csv.gz`、制表符分隔的功能注释表 | 只规范化成功能元数据，不能替代 GFF/CDS/PEP，也不能当坐标文件 |

同一角色有多个候选时，在清单里写 `A|B`；程序按角色兼容性打分，平局会明确报歧义，不会乱选。

注册一个物种（会复制四个源文件到规范目录并写 inventory）：

```powershell
cd code
python -m pgci_preprocess init-project `
  --project-root D:\my_project `
  --show-id 2848 --scientific-name "Macrotyloma uniflorum" `
  --order Fabales --family Fabaceae `
  --genome  D:\raw\GWHERCC00000000.genome.fasta.gz `
  --gff     D:\raw\GWHERCC00000000.gff.gz `
  --cds     D:\raw\GWHERCC00000000.CDS.fasta.gz `
  --pep     D:\raw\GWHERCC00000000.Protein.faa.gz `
  --species-names .\resources\species_abbreviations.tsv
```

三个参数可以省：

- `--genome`：只有 GFF + CDS + PEP 三件套时可不传；长度由 CDS 计算，region 缺失时由转录本坐标回填。无 genome 的简化注释还会对合成 region 做细胞器推断（见第 4 节）。
- `--assembly-id`：不写时从 genome 文件名推断（`GWHERCC00000000.genome.fasta.gz` → `GWHERCC00000000`，NCBI 的 `GCF_xxx.x_..._genomic.fna.gz` 保留版本号）；无 genome 时可留空。
- `--species-names`：不写或指向项目内尚未生成的路径时，自动使用包内 `resources/species_abbreviations.tsv`；表里有的 show_id 沿用表内 code，表里没有的按命名规则现算并避开表内全部 code。表外物种请提供 `--scientific-name`（批量时填清单 `Species` 列），以便生成缩写。

批量模式：把 `data/downloads` 放好，写 `resources/download_inventory.tsv`，然后 `prepare`。

---

## 4. 生成 WGDI 输入（本包的核心）

```powershell
cd code
python -m pgci_preprocess run --project-root D:\my_project
python -m pgci_preprocess validate --project-root D:\my_project --refresh-links
python -m pgci_preprocess check-names --project-root D:\my_project
```

跑完在 `analysis_v2/<Abbr>/`、`clean_v2/<Abbr>/` 下得到：

```text
<Abbr>.pep  <Abbr>.cds  <Abbr>.new.gff  <Abbr>.lens
<Abbr>.all.pep.fa  <Abbr>.all.cds.fa  <Abbr>.all.gff  <Abbr>.all.lens
<Abbr>.gene_map.tsv  <Abbr>.chr_map.tsv  <Abbr>.rejected.tsv  <Abbr>.qc.json
```

`validate --refresh-links` 会在项目根建 `<Abbr>/` 兼容目录（下游脚本按这个目录取文件）。
也可以 `--no-compat-links` 跳过，只用 `analysis_v2` 里的规范文件。

### 细胞器推断与准备度门禁（3.5.4 起）

- **细胞器推断**：无 genome、无显式 region 的简化注释，`run` 会对合成 region
  识别叶绿体/线粒体（序列或基因标识符信号，或高密度、紧凑、基本无内含子的
  基因集），排除后计入 QC 的 `n_organelle_genes`，证据在 `organelle_inference`。
  拿不准的高密度 region 只列入 `review_candidates` 供人工复核，不会静默排除。
  显式 region 和 genome FASTA 记录永不被启发式改写。可用
  `--organelle-pattern "^(cp|mt)_scaffold"` 之类的正则强制排除指定合成序列，
  或 `--no-organelle-inference` 完全关闭推断。
- **植物线粒体（3.5.5 新增）**：检测一个 region 中不同的
  `nad/cob/cox/atp/matR/mttB/ccm` 核心基因家族，不套用叶绿体的紧凑形态假设。
  `--organelle-mt-min-genes` 可调整自动排除所需家族数；证据不足的候选只进入
  `review_candidates`。ID 完全匿名且没有参考相似度的稀疏线粒体仍需用
  `--organelle-pattern` 或人工检查；QC 的 `organelle_inference.limits`
  会明确记录这个边界。
- **准备度门禁**：`validate` 对 QC 中 `wgdi_ready=no` 的物种默认判 FAIL
  （`wgdi_not_ready_low_coverage`），避免碎片化物种因“四件套齐全”被误用。
  只做数据完整性归档/复核时显式加 `--allow-not-ready`。
  下游 `run_homology.sh` / `run_all_pipeline.ps1` 同样跳过 `ready=no` 物种。
- **readiness 单一权威源（3.5.5）**：`python -m pgci_preprocess readiness`
  按物种结构化读取 QC，缺失时回退 index；Python、Shell 和 PowerShell 不再各自
  实现不同的 JSON 文本匹配规则。
- **下游产物完整性（3.5.5）**：`6run_wgdi.sh` 调用
  `python -m pgci_preprocess check-wgdi-outputs` 检查 8 个产物的 PNG 魔数、
  Ks 表头、共线性块和 CSV 表头；截断但非空的文件不会标记为完成。

### 命名规则（务必理解）

- `show_id`（如 `2848`）只是内部主键；**目录名、文件名、基因 ID 一律用缩写**（如 `Mun`）。
- 基因 ID 形如 `Mun1g0001` = `{缩写}{序列号}g{四位顺序号}`；`.cds/.pep` 头行 `>Mun1g0001 source_gene=...`。
- 缩写来自 `code/resources/species_abbreviations.tsv`（由物种缩写表导入，全局唯一）。
- 旧项目里若还是 `s2848/s2848.cds`，用迁移命令改成缩写命名（可先 dry-run）：

```powershell
python -m pgci_preprocess migrate-artifacts --project-root D:\my_project --dry-run
python -m pgci_preprocess migrate-artifacts --project-root D:\my_project --prune-legacy-dirs
```

`check-names` 会把逐物种核对结果写到 `resources/naming_audit.tsv`（表↔清单↔下载目录↔索引↔产物五方对齐）。

---

## 5. 可选：下游比对和 WGDI（不属于预处理）

前面第 3–4 步已经完成预处理，产出 WGDI 可直接读取的四件套。下面的下游层是
**可选的下游执行层**：`pgci_preprocess run` 不会调用它们。只有当你还想用本包
完成 DIAMOND/BLAST 同源比对和 WGDI 分析时才需要；纯预处理用户跳到第 7 节即可。

Windows 分析总控（步骤 1–6；第 0 步负责预处理和验收）：

```powershell
bash 0check_data.sh --preprocess --refresh-links
powershell -ExecutionPolicy Bypass -File .\run_all_pipeline.ps1
# 外部程序不在 PATH 时：
# -BioinfoBin D:\bioinfo\bin -DiamondPath D:\bioinfo\bin\diamond.exe
```

逐步（Git Bash）：

```bash
bash 0check_data.sh --preprocess --refresh-links
bash run_homology.sh           # 默认 DIAMOND；设 HOMOLOGY_ENGINE=blast 改变引擎
bash 2run_pre_wgdi.sh
bash 3blast_to_wgdi.sh
python code/tools/generate_wgdi_configs.py --mode all
bash 6run_wgdi.sh
```

`1run_diamond.sh`、`1run_blast.sh`、`4run_wgdi_different.sh` 和
`5run_wgdi_same.sh` 继续可用，但只是共享实现的兼容入口。

`run_homology.sh` 会按 QC/index 的 `wgdi_ready` 跳过未就绪物种（记录在
`skipped_species.txt`），也不把 `bin/`、`resources/` 等目录当物种扫描；
确需复核时设 `HOMOLOGY_ALLOW_NOT_READY=true`。生成的 WGDI `run.sh` 使用
Bash 安全的 POSIX 路径（`/d/...`），Git Bash 与 Linux 均可直接执行。

`bin/wgdi` 会设置 `MPLBACKEND=Agg` 和临时 `MPLCONFIGDIR`，无图形界面的服务器也能跑；
Windows 反斜杠 `TMPDIR` 会被规范化，系统临时目录不可写时回退到项目内 `.wgdi-tmp`。

---

## 6. 常见问题

| 现象 | 原因 / 处理 |
| --- | --- |
| 产物名还是 `s1234` | 既没有缩写表记录，也没有拉丁名；补 `--scientific-name` 或清单 `Species` 列。若表丢失，`doctor` 会报 `species_abbreviations FAIL` |
| `validate` 报 `compat_missing_*` | 加 `--refresh-links` 建兼容目录 |
| `validate` 报 `wgdi_not_ready_low_coverage` | 该物种 `wgdi_ready=no`（碎片化、保留率过低），不适合 WGDI；仅做数据完整性复核时加 `--allow-not-ready` |
| 怀疑某条合成序列是细胞器但未被排除 | 看 QC 的 `organelle_inference.review_candidates`；确认后 `run --organelle-pattern <正则>` 重跑该物种 |
| WGDI 报 `TMPDIR` 相关错误或 `mkdir /tmp` 失败 | 新版 `bin/wgdi` 已规范化 Windows `TMPDIR` 并自动回退；仍失败时检查项目目录写权限 |
| `wgdi: not found python with wgdi` | 先跑 `bootstrap_wgdi.ps1/.sh`，或设 `WGDI_PYTHON` |
| WGDI 仍使用旧机器上的 MUSCLE/yn00 路径 | 新版启动入口会依次读取 `WGDI_MUSCLE`、`MUSCLE`、`WGDI_YN00`、`YN00` 和 PATH |
| `PAL2NAL 未安装` | 运行 `python code/tools/install_pal2nal.py`，或用 `WGDI_PAL2NAL` 指向已有脚本 |
| `doctor` 报 diamond/blastp FAIL | 装工具并加入 PATH，或用 `-BioinfoBin` / `DIAMOND` / `BLASTP` 环境变量 |
| Windows 路径过长报错 | 项目与临时目录尽量短（<250 字符），或把项目放到盘符根下 |
| 迁移中断后目录里混着新旧文件名 | 重跑 `migrate-artifacts`（幂等，会补齐 `>` 头、GFF 列与基因 ID） |
| `pytest` 报 `PermissionError: ...\Temp\pytest-of-*` | 测试会自动检测默认临时目录是否可写并重定向；也可显式加 `--basetemp=<short-path>` |
| GBFF 转换报 `extracted CDS and translation qualifier are inconsistent` | 真实 NCBI 文件个别 CDS 与 `/translation` 不一致；默认严格中止，Python API 传 `convert_gbff_to_bundle(..., strict_translation=False)` 可跳过并计数 |

---

## 7. 交付前自检（建议）

```powershell
cd code
python -m pgci_preprocess validate --project-root D:\my_project --refresh-links
python -m pgci_preprocess check-names --project-root D:\my_project
python -m pytest -q                      # 无系统临时目录权限时加 --basetemp=D:\pt
```

`validate` 只校验数据完整性：缺少 `resources/species_abbreviations.tsv` 时跳过命名审计（命名问题由 `check-names` / `doctor` 报告），不会因此判 FAIL。两项都通过，才把 `analysis_v2/<Abbr>/` 交给 WGDI 使用。

需要验证完整下游执行时，再在隔离目录运行自包含 smoke：

```powershell
python code/tools/build_synthetic_wgdi_smoke.py `
  --output-root $env:TEMP\plant-genome-smoke `
  --genes-per-species 120 --chromosomes 8 --protein-length 80
bash code/tools/run_wgdi_smoke.sh $env:TEMP\plant-genome-smoke
```

发布门禁和已通过的 3.5.5 回执见 `VALIDATION.md` 与 `validation/`。
可直接交接的压缩包为 `dist/plant-genome-preprocess-v3.5.5.zip`。
