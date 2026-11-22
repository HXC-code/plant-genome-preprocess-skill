# HXC 植物基因组预处理核心

本目录使用统一入口 `python -m pgci_preprocess`，将 GFF3/GTF/GBFF 注释与
CDS/蛋白 FASTA 整理为带全局唯一 ID、可直接进入 BLAST/WGDI 的一基因一代表转录本数据；
在 genome 和带 CDS 的 GFF 可用时，程序可从坐标自动提取缺失的 CDS 和蛋白序列。
`scripts/` 下的四个旧入口保留为薄兼容包装器。

## 物种代号

物种的表型/图上/目录代号来自实验室维护的缩写表（例如 `species_abbreviations.xlsx`），规范化副本是
`resources/species_abbreviations.tsv`。`show_id` 仍是内部稳定主键，分析代号使用表中的
全局唯一缩写，例如 `show_id=2720` 对应 `Lcu`，`show_id=2725` 对应 `Lsb`。

旧 `s<show_id>`（如 `s2720`）只作为兼容别名；`--species` 同时接受 `2720`、`s2720` 和
`Lcu`。Excel 缩写表更新后执行：

```powershell
python tools/import_species_names.py --input <path/to/species_abbreviations.xlsx>
python -m pgci_preprocess migrate-species-names --project-root .
```

迁移命令先备份 metadata，再更新 inventory/registry/summary/index；既有生物结果不重算，
兼容目录会生成 `Lcu.pep`、`Lcu.cds`、`Lcu.new.gff`、`Lcu.lens` 等缩写文件名。

原始下载目录不会被修改：

```text
data/downloads/
resources/download_inventory.tsv
```

## 新项目初始化与依赖预检

预处理项目根保存 `data/downloads` 和 `resources`，分析执行根保存
`0check_data.sh` 至 `6run_wgdi.sh` 以及 `<Abbr>/` 四件套兼容目录。两个目录分别管理数据与执行参数，便于长期归档和复现。

```powershell
# 自动创建嵌套项目根、data/downloads、resources 和两张标准 TSV 表头
python -m pgci_preprocess init-project --project-root <new-project>

# 空项目结构预检；正式预处理预检不加 --allow-empty
python -m pgci_preprocess doctor --project-root <new-project> --allow-empty

# 检查完整执行链；所有路径均可显式覆盖，也可用 PATH/同名环境变量
python -m pgci_preprocess doctor --project-root <project> --workflow full `
  --diamond <diamond.exe> --wgdi <bin/wgdi> --bash <bash.exe> `
  --muscle <muscle.exe> --yn00 <yn00.exe> --perl <perl.exe>
```

`doctor` 只读，不安装软件、不创建目录、不修改项目。工具检查表示“路径可解析”；
`--workflow wgdi|full` 额外运行 `wgdi_probe`：实际启动 WGDI 入口读版本，
暴露“存在但跑不起来”的解释器/环境问题（例如用错 Python 环境导致异常退出）。
完整执行链仍由独立 smoke 验收。GBFF 转换还需 `biopython>=1.83`，建议在新
虚拟环境中安装本项目（Python ≥ 3.11），不要复用依赖冲突不明的 base 环境。Windows 建议使用较短的项目和临时目录；
过长的物种目录与 pytest 临时路径可能触发传统路径长度限制。

```powershell
python -m venv .venv
# 激活新环境后：预处理 + GBFF 支持
python -m pip install -e ".[gbff]"
# 完整 WGDI：使用父目录中本次真实 smoke 验证过的精确版本
python -m pip install -r ..\requirements-wgdi-validated.txt
```


`init-project` 的一次性显式登记要求 GFF、CDS、PEP 三个不同源文件；genome 可选
（只有三件套时长度由 CDS 计算、region 由转录本回填，`--assembly-id` 也可留空）。批量
`prepare` 则支持 GTF/GBFF 自动转换，以及 `genome + GFF` 自动提取缺失的 CDS/PEP；
角色先由列名确定，后缀只作提示，最终以文件内容和角色兼容性复核。

```powershell
python -m pgci_preprocess init-project --project-root <new-project> `
  --show-id 42 --scientific-name "Example species" --order <order> --family <family> `
  --assembly-id <assembly> `
  --species-names <project-root>/resources/species_abbreviations.tsv `
  --genome <genome.fna.gz> --gff <annotation.gff3.gz> `
  --cds <cds.fa.gz> --pep <protein.faa.gz>
```

初始化源文件只读，目标优先硬链接、失败时复制并核对 SHA-256；已有物种目录、重复 show ID、
缺失/空文件和不安全 assembly ID 默认失败关闭，不覆盖既有项目。

## 运行

```powershell
cd <project-root>

# 测试在仓库的 code/ 目录执行
# cd <repo>/code
python -m pytest -q

# 正式子命令（convert 为历史兼容入口）
python -m pgci_preprocess prepare --source-status <status.tsv> --expected-count 20
python -m pgci_preprocess normalize-gff --species 4200,5630
python -m pgci_preprocess convert-annotation --input <annotation.gtf.gz> --output-dir <derived> --prefix <code>
python -m pgci_preprocess run --species 996
python -m pgci_preprocess validate --expected-count 109 --no-compat-links
python -m pgci_preprocess batch --source-status <status.tsv> --expected-count 20
python -m pgci_preprocess audit --json

# 旧入口仍可用（内部转发到统一 CLI）
python scripts/preprocess_wgdi_v2.py --species 5,225,996
python scripts/validate_wgdi_dataset.py --expected-count 109 --no-compat-links
```

`run` 默认不覆盖已有输出；确实需要重跑时显式加 `--overwrite`。旧参数
`--no-overwrite` 仍被接受，语义等同默认行为。

完整的输入路由按以下顺序执行，并写入
`resources/preprocess_batches/<tag>_routing.tsv`：

```text
后缀提示 -> 内容确认 -> 选择路径 -> 转换/提取 -> 强校验 -> 注册 inventory -> run -> validate
```

| 输入内容 | 路由 | 注册结果 |
| --- | --- | --- |
| GFF/GFF3 五种结构 | 直接 run | 原 GFF 路径 |
| GTF | `convert-annotation` 生成 GFF3 | 新 GFF3；原始 GTF 写入 `gtf` 辅助列 |
| GBFF/GenBank | 生成 GFF3 + CDS + PEP（可含 genome） | 派生四件套；原始 GBFF 写入 `gbff` 辅助列 |
| genome + 带 CDS 的 GFF，但缺 CDS/PEP | 按 CDS 坐标拼接、反向互补并翻译 | 派生 CDS/PEP 及 extraction manifest |
| inventory 的 `A\|B` 候选 | 按 gff/genome/cds/pep 角色和内容筛唯一候选 | run 可直接消费；prepare 会固化选择 |
| `.gz/.bz2/.xz/.zip` 或错误后缀压缩文件 | 按压缩魔数识别 | 后缀不再是解压前提 |

可选参数：

```text
--min-genes-per-sequence 50   # 碎片组装中一个序列至少保留的基因数
--coverage 0.90               # 碎片组装兜底保留的基因覆盖率
--max-sequences 200           # WGDI 分析集最大序列数
--overwrite                   # 已存在输出时允许替换
--gff-format auto             # auto/standard/transcript_only/simplified5/gene_span/gene_cds_direct
--protein-normalization auto  # auto/strict/evm-dot
--accession-policy strict     # strict/legacy
--organelle-pattern <regex>   # 强制排除匹配的合成序列 ID（细胞器）
--no-organelle-inference      # 关闭合成 region 的细胞器推断
--organelle-mt-min-genes 3    # 线粒体核心基因家族自动排除阈值
--min-genes-per-chromosome 100 # 丢弃基因数不足 N 的染色体（CLI 默认 100，0 关）
```

无 genome、无显式 region 的简化注释，`run` 会对合成 region 做细胞器推断：
序列/基因标识符信号（chloroplast/plastid/mitochondrion、cp/chrc/mt、plastid 基因名）
或“高密度 + 紧凑 + 基本无内含子”的基因集会被标记为细胞器并排除；拿不准的
高密度 region 只列入 QC `organelle_inference.review_candidates` 供人工复核。
显式 region 行和 genome FASTA 记录永不被启发式改写。排除证据与计数写入
`organelle_inference` 和 `n_organelle_genes`。

`validate` 默认将 QC 中 `wgdi_ready=no` 的物种判 FAIL
（`wgdi_not_ready_low_coverage`）；仅做数据完整性归档时显式加
`--allow-not-ready`（`0check_data.sh` 同名参数透传）。

### 染色体清单与人工裁决（3.6.0）

`run` 生成 `<code>.chromosome_map.tsv`（analysis_v2 与 clean_v2），列出保留染色体数、
每条基因数/长度/来源/是否 primary；QC 增 `chromosome_summary` 与顶层
`n_wgdi_sequences`。`--min-genes-per-chromosome`（CLI 默认 100）在选择后丢弃
基因数不足的染色体，排除原因 `below_min_genes_per_chromosome` 记入 chr_map。

```powershell
python -m pgci_preprocess chromosomes --project-root <root> --species <code>
python -m pgci_preprocess adjudicate --project-root <root> --species <code> `
  --confirmed-n <N> --decision accept --reviewer <name> --rationale <text>
python -m pgci_preprocess adjudication-status --project-root <root> --species <code>
python -m pgci_preprocess validate --project-root <root> --require-adjudication
```

裁决记录绑定 QC 版本 + SHA-256，重跑后旧签字失效（stale）。decision：
accept / accept_with_note（confirmed_n 须等于程序染色体数）、adjust / reject（门禁拦截）。
`validate --require-adjudication` 与下游 `HOMOLOGY_REQUIRE_ADJUDICATION=true`
构成硬门禁；裸 CLI 默认关闭以便归档，Skill 调用默认开启。

## 同源搜索与 WGDI 执行

本章脚本构成可选下游层，用于在预处理完成后执行同源比对和 WGDI 分析。
只进行数据标准化时，上一章的 `run` / `validate` / `check-names` 完成后即可交付四件套。
需要继续执行同源比对和 WGDI 分析时，按下面的步骤运行。

预处理通过只代表“输入就绪”。实际执行在 `wgdi_pipeline` 根目录进行，四件套必须先放在
`s<code>/` 中并成套使用：

```bash
export PATH="$PWD/bin:$PATH"  # 提供项目内 WGDI wrapper
# 非本机路径时设置 DIAMOND、WGDI_PYTHON、WGDI_PERL_BIN，或把工具加入 PATH
bash 0check_data.sh             # 可选：预处理项目根存在时执行完整预处理和验证
bash run_homology.sh            # 默认 DIAMOND, very-sensitive + 1e-5
bash 2run_pre_wgdi.sh
bash 3blast_to_wgdi.sh          # 复制比对结果，保留 Auto_pep 原始文件
python tools/generate_wgdi_configs.py --mode all
bash 6run_wgdi.sh               # 验收：默认 2 个核心结果（点图+共线性）；WGDI_KS=true 时 8 个
```

`run_homology.sh` 按 QC/index 的 `wgdi_ready` 跳过 `ready=no` 和没有可读准备度
记录的物种（写入 `skipped_species.txt`），并把 `bin/`、`code/`、`resources/` 等
基础设施目录排除在物种扫描之外；确需复核未就绪物种时设
`HOMOLOGY_ALLOW_NOT_READY=true`。生成的 `run.sh` 使用 Bash 安全的 POSIX 路径
（Windows `D:\...` → `/d/...`）。

`1run_diamond.sh`、`1run_blast.sh`、`4run_wgdi_different.sh` 和
`5run_wgdi_same.sh` 是兼容入口；新脚本只需要维护共享实现。

WGDI 的串联重复过滤由 `tandem=false` 控制，本项目默认开启。上游 `-r <config>` 用于
保留率分析；项目 wrapper 使用末尾单独 `-r` 或 `--remove-tandem` 控制串联重复过滤。
恢复串联重复时设置 `WGDI_TANDEM=true` 重新生成配置。

**Ks 分支默认关闭**（SOP 第 4 步可选分支）：默认只运行 dotplot + collinearity，
产出 `.png` 与 `.icl`。设置 `WGDI_KS=true` 后追加 ks/blockinfo/correspondence/
blockks/kspeaks/peaksfit 六个阶段（配置生成、`6run_wgdi.sh` 验收、Windows 总控
`-WithKs` 均读同一开关）。原因：上游 blockinfo/correspondence 依赖 Ks 值，
缺 Ks 时 correspondence 会过滤掉全部区块。

DIAMOND 可用 `DIAMOND_EVALUE`、`DIAMOND_SENSITIVITY`、`DIAMOND_THREADS`、
`DIAMOND_MAX_TARGETS` 配置；本次参数写入 `Auto_pep/diamond_run_settings.tsv`。`1e-5` 是同源
搜索阈值，由下游脚本配置。`very-sensitive` 是速度—灵敏度折中；`ultra-sensitive` 提供更高灵敏度并增加运行时间，可按数据集和验证结果选用。

Windows 总控默认选择 DIAMOND，不需要指定历史 `bioinfo_win` 环境：

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all_pipeline.ps1 `
  -HomologyEngine diamond -DiamondSensitivity very-sensitive -HomologyEvalue 1e-5
# 改用 BLAST：-HomologyEngine blast
# 启用可选 Ks 分支：-WithKs
# 非默认安装位置：-DiamondPath <path> 或 -BlastpPath <path> -MakeBlastDbPath <path>
```

总控默认不自动重跑已有完成状态；确需处理失效或未完成步骤时，检查日志后显式使用
`-RetryIncomplete`。WGDI wrapper 可用 `WGDI_PYTHON`、`WGDI_PERL_BIN`、`WGDI_PERL5LIB`
覆盖本机路径；MUSCLE、yn00、PAL2NAL 也可分别用 `WGDI_MUSCLE`、`WGDI_YN00`、
`WGDI_PAL2NAL` 覆盖。默认 Ks 列为 `ks_NG86`，可通过 `WGDI_KS_COL=ks_YN00`
在短序列冒烟或特定数据上切换。

PAL2NAL v14 不随发布包分发。`bootstrap_wgdi` 默认会从官方站点安装；也可以单独执行：

```powershell
python tools/install_pal2nal.py
```

安装器从官方站点下载并校验 SHA-256，文件写入被忽略的 `external/` 目录。

## 独立真实数据 smoke

`tools/build_wgdi_smoke_fixture.py` 可从既有四件套和历史种内 `.icl` 选择一个小型共线性子集，
复制分析脚本到全新目录；历史比对只用于选基因，DIAMOND 必须重新计算：

```powershell
python tools/build_wgdi_smoke_fixture.py --workflow-root <wgdi_pipeline> `
  --species-dir <wgdi_pipeline/Ggl> --self-icl <old/Ggl.icl> `
  --output-root <wgdi_pipeline/smoke_wgdi>
bash tools/run_wgdi_smoke.sh <wgdi_pipeline/smoke_wgdi>
```

启动器依次运行 DIAMOND、2–6 号脚本，日志在 `logs/`，最终状态在 `smoke_status.tsv`。
这一步验证的是派生小样本的 DIAMOND→WGDI 执行，不等同于重新运行全部物种；第 0 步的
预处理结构/哈希验收由项目测试和 manifest 验证单独覆盖。


## 结构兼容与格式例外

`gff_formats.py` 明确识别五种结构，不依赖文件名：

```text
standard          gene -> mRNA -> CDS
transcript_only   只有 mRNA/CDS，自动补 gene 层
simplified5       五列简化 GFF
gene_span         只有 gene 行，gene 即代表转录本
gene_cds_direct   CDS 直接挂在 gene 上
```

无法从结构安全推断的历史特例（例如 NCBI 分段 gene 需要合并）记录在
`resources/format_exceptions.tsv`，并写入每个物种的 `qc.json`。GWH `Accession`、
`cds.`/`cds-` 前缀、无 region 时的 transcript 坐标回填、EVM 蛋白尾部 `.` 规范化
都作为显式适配步骤执行，不再使用 import-time monkeypatch。

发行包副本由 canonical 源同步，不允许手工维护：

```powershell
# 旧式 skill-packages/scripts 布局
python tools/sync_skill_package.py --check
python tools/sync_skill_package.py

# 发布用的 Skill 仓库镜像（整棵 pipeline/ 树）
python tools/sync_skill_package.py --pipeline-target <skill-repo>/pipeline --check
python tools/sync_skill_package.py --pipeline-target <skill-repo>/pipeline
```

两个仓库的 GitHub Actions 会在 push / PR 时跑 `pytest`、`ruff`、`mypy`，Skill 仓库还会做上面的镜像漂移检查。

## 额外格式转换（v2.2.0）

后缀只用于提示，转换器按文件内容和数据角色处理，不把改后缀当作格式转换。
压缩层支持 gzip、bzip2、xz 和单成员 ZIP，并通过魔数识别错误后缀或缺失后缀。
原始下载文件只读，输出默认写入独立的派生目录，默认拒绝覆盖，并生成带源文件及输出
SHA-256 的 `<prefix>.conversion.json`。

```powershell
# GTF -> GFF3；随后把该 GFF3 与同一 assembly 的 genome/CDS/PEP 一起交给 run
python -m pgci_preprocess convert-annotation `
  --input <species.gtf.gz> --output-dir <derived> --prefix <code>

# GBFF -> 编码 GFF3 + CDS/PEP；默认同时提取 genome FASTA
python -m pgci_preprocess convert-annotation `
  --input <species.gbff.gz> --output-dir <derived> --prefix <code>

# 功能注释表（项目中的 *.annotations.csv.gz）-> 标准七列 TSV
python -m pgci_preprocess convert-annotation `
  --input <species.annotations.csv.gz> --output-dir <derived> --prefix <code>
```

支持边界：

| 输入内容 | 输出 | 进入 BLAST/WGDI 的方式 |
| --- | --- | --- |
| GTF | `<prefix>.gff3.gz` | `prepare/batch` 自动转换并注册；若 CDS/PEP 缺失且 genome 可用，会继续自动提取 |
| GenBank/GBFF | `<prefix>.gff3.gz`、`.cds.fa.gz`、`.pep.fa.gz`，默认另有 `.genome.fa.gz` | 四者作为同一派生 bundle 交给 `run`；导出内容覆盖编码分析所需结构 |
| 功能注释表 | `<prefix>.functional.tsv.gz` | 规范化功能元数据；WGDI 坐标分析使用 GFF/CDS/PEP |

GBFF 解析使用惰性加载的 Biopython（>=1.83），只在使用 GBFF 转换时需要：

```powershell
python -m pip install "biopython>=1.83"
```

GBFF 的 CDS/PEP 以 `protein_id` 为键，并在写出后强制回读检查 GFF/CDS/PEP 的 ID 集合一致。
NCBI Gnomon 多 isoform 会先按显式 transcript、坐标结构和 `isoform Xn`/`transcript variant Xn`
消歧；需要时从 CDS 合成 transcript，并在报告中写入警告和对应关系。
重复 protein_id 或显式 transcript 对应多个 protein 会失败关闭。`--max-records N` 仅用于
GBFF 冒烟测试，manifest 会写成 `complete=false`，便于区分测试子集和完整数据。
转换时会用 `transl_table` 重新翻译提取的 CDS，并与 `translation` qualifier 强制比对；
不一致直接拒绝。GFF 中的 `transl_table` 也会进入 run 的翻译 QC 和 `gene_map.tsv`。

功能表的实测格式为无表头、tab 分隔的七列混合功能注释（GO、NR、SP、TR）。输出保留所有来源并分别计数；描述中的换行会转为可读的 `\\n` 转义文本。
CSV/TSV/JSONL 转换依据明确 schema 执行。

## 核心规则

1. GFF3 中的 gene、mRNA、CDS 和 protein_id 建立映射。
2. 一个基因只选一个代表转录本（SOP 强制规则）：
   - 第一优先级：CDS 序列最长；
   - CDS 等长时按稳定转录本 ID 字典序升序取第一条；
   - 不以蛋白长度或序列出现顺序作为依据，保证可复现。
3. 只接受 CDS 与 PEP 能实际翻译对齐的候选。
4. 拒绝以下数据：
   - CDS 长度不是 3 的倍数；
   - 内部终止密码子；
   - 氨基酸非法字符；
   - 翻译长度或内容不一致；
   - 没有 CDS/PEP 的非编码转录本；
   - 细胞器基因；
   - 跨核基因组 taxon 污染；
   - 同一蛋白 ID 被不同基因复用；
   - 重复且序列冲突的 FASTA ID。
5. 原始 ID、是否代表转录本、翻译状态和拒绝原因都保留在映射/报告文件中。

## 输出

```text
clean_v2/{code}/
├── {code}.gene_map.tsv
├── {code}.chr_map.tsv
└── {code}.removed_list.tsv

analysis_v2/{code}/
├── {code}.cds                 # WGDI 筛选集 CDS
├── {code}.pep                 # WGDI 筛选集蛋白
├── {code}.new.gff             # WGDI 七列位置文件
├── {code}.lens                # WGDI 染色体长度/基因数
├── {code}.all.cds.fa          # 全部合格代表基因
├── {code}.all.pep.fa
├── {code}.all.gff
├── {code}.all.lens
├── {code}.qc.json
├── {code}.gene_map.tsv
└── {code}.chr_map.tsv
```

`{code}.new.gff` 为七列：

```text
简化序列号  唯一基因ID  起点  终点  链方向  基因序号  原始基因ID
```

WGDI 使用前六列，第七列用于回溯 NCBI 原始基因名。

## 汇总文件

```text
resources/stage_preprocess_v2_summary.tsv
resources/wgdi_v2_index.tsv
resources/species_registry_v2.tsv
resources/wgdi_v2_manifest_sha256.tsv
```

## 项目目录

```text
plant-genome-resources/
├── data/downloads/      # 原始下载，只读
├── resources/           # manifest、QC、哈希
├── clean_v2/            # 映射和拒绝原因
├── analysis_v2/         # WGDI 四件套
├── results/
│   ├── blast/
│   ├── wgdi/
│   ├── synteny/
│   ├── tree/
│   └── wgd/
├── logs/
├── runs/
├── scripts/
├── src/
└── tests/
```

目录决策见 `docs/adr/0001-data-layout-and-validation-contract.md`。

## WGDI 环境

不同机器上 WGDI 的安装位置不同，因此不再写死任何解释器路径。仓库根目录
`bin/wgdi` 会按 `$WGDI_PYTHON` → `.wgdi-env` → `.wgdi_python` → PATH 的顺序
自动定位解释器；首次使用直接运行：

```powershell
powershell -ExecutionPolicy Bypass -File .\bootstrap_wgdi.ps1
```

```bash
bash bootstrap_wgdi.sh
```

装好后等价于：

```text
.wgdi-env/Scripts/python.exe -m wgdi ...
```

若本机已有可用环境（例如 `conda` 的 `pced`），把该解释器路径写一行到
`.wgdi_python`（仅本机使用，不要放进发布包），或设置 `WGDI_PYTHON`。

WGDI 配置中的 CDS、PEP、GFF 和 lens 应使用同一套筛选集：

```text
{code}.cds
{code}.pep
{code}.new.gff
{code}.lens
```

## 环境注意事项（2026-11-20 实测）

- **pytest 临时目录**：`tests/conftest.py` 在 basetemp 超过 90 字符或默认临时目录不可写时自动重定向到短临时根并打印 `PytestWarning`（Windows 传统 260 字符限制与过期 ACL 会造成大批量假失败）。仍可显式用 `--basetemp=<short-temp-path>`。
- **Windows 长路径守卫**：`run` 遇到物种目录/输入/输出路径超过 250 字符时返回 `FAIL path_too_long_for_windows`（detail 给出整改建议），不再抛出裸 `FileNotFoundError`。
- **GFF 结构错误语义**：`GffFormatError` 统一转为 `FAIL gff_format_error`（不再以 `ERROR` 跳过统计）；原始 GTF/GBFF 直接进 run 则返回 `FAIL conversion_required` 并附可执行 `convert-annotation` 命令。
- **Biopython**：仅 GBFF 转换需要 `biopython>=1.83`；`doctor` 会扫描当前解释器、PATH 上的 `python/python3/py` 以及 `PGCI_PYTHON_CANDIDATES`（`;` 分隔）并列出可用的解释器版本。
- **脚本行尾**：`0check_data.sh` 至 `6run_wgdi.sh` 已统一为 LF（备份在 `99_ARCHIVE\shell_lineendings_20261120`），Git Bash 与 WSL/Linux 下 `bash -n` 均通过。

