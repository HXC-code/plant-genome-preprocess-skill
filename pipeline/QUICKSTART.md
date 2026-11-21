# 快速上手：从自己的原始下载数据到 WGDI 输入

本包只包含**代码**：预处理程序、0–6 步执行脚本、WGDI 启动器。不含任何原始数据、也不含已生成的结果目录。
拿到包的人只要有 Python 和少量外部程序，就能从自己的数据生成 WGDI 能直接跑的
`<Abbr>.pep / .cds / .new.gff / .lens`，并按缩写命名。

---

## 0. 包里有什么

```text
wgdi_pipeline/
├── QUICKSTART.md                  # 本页
├── bootstrap_wgdi.ps1 / .sh        # 一键创建 WGDI 运行环境（.wgdi-env）
├── bin/wgdi                        # 可移植启动器（自动找 .wgdi-env / PATH）
├── 0check_data.sh … 6run_wgdi.sh   # 0–6 步流程
├── run_all_pipeline.ps1            # Windows 一键串联 0–6
├── requirements-wgdi-validated.txt # WGDI 依赖钉版清单
├── vendor/                         # 内置 WGDI 兼容 shim
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
| WGDI | Python + `wgdi` 包 | 用 `bootstrap_wgdi` 自动装，不要手改脚本里的解释器路径 |
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

两个参数可以省：

- `--assembly-id`：不写时从 genome 文件名推断（`GWHERCC00000000.genome.fasta.gz` → `GWHERCC00000000`，NCBI 的 `GCF_xxx.x_..._genomic.fna.gz` 保留版本号）。
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

`validate --refresh-links` 会在项目根建 `<Abbr>/` 兼容目录（0–6 步脚本按这个目录取文件）。
也可以 `--no-compat-links` 跳过，只用 `analysis_v2` 里的规范文件。

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

## 5. 跑 0–6 步

Windows 一键：

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all_pipeline.ps1
# 外部程序不在 PATH 时：
# -BioinfoBin D:\bioinfo\bin -DiamondPath D:\bioinfo\bin\diamond.exe
```

逐步（Git Bash）：

```bash
bash 0check_data.sh --preprocess --refresh-links
bash 1run_diamond.sh          # 或 1run_blast.sh
bash 2run_pre_wgdi.sh
bash 3blast_to_wgdi.sh
bash 4run_wgdi_different.sh
bash 5run_wgdi_same.sh
bash 6run_wgdi.sh
```

`bin/wgdi` 会设置 `MPLBACKEND=Agg` 和临时 `MPLCONFIGDIR`，无图形界面的服务器也能跑。

---

## 6. 常见问题

| 现象 | 原因 / 处理 |
| --- | --- |
| 产物名还是 `s1234` | 既没有缩写表记录，也没有拉丁名；补 `--scientific-name` 或清单 `Species` 列。若表丢失，`doctor` 会报 `species_abbreviations FAIL` |
| `validate` 报 `compat_missing_*` | 加 `--refresh-links` 建兼容目录 |
| `wgdi: not found python with wgdi` | 先跑 `bootstrap_wgdi.ps1/.sh`，或设 `WGDI_PYTHON` |
| `doctor` 报 diamond/blastp FAIL | 装工具并加入 PATH，或用 `-BioinfoBin` / `DIAMOND` / `BLASTP` 环境变量 |
| Windows 路径过长报错 | 项目与临时目录尽量短（<250 字符），或把项目放到盘符根下 |
| 迁移中断后目录里混着新旧文件名 | 重跑 `migrate-artifacts`（幂等，会补齐 `>` 头、GFF 列与基因 ID） |
| `pytest` 报 `PermissionError: ...\Temp\pytest-of-*` | 系统临时目录不可写：加 `--basetemp=D:\pt`（目录要短） |
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
