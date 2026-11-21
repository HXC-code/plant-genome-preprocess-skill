# HXC Plant Genome Preprocess

Prepare **your own downloaded plant genome data** for WGDI. Hand the tool the
annotation and sequence files you have; it identifies each file's format and
biological role by content, normalizes gene models and IDs, and writes the four
files WGDI expects:

```text
<Abbr>.pep   <Abbr>.cds   <Abbr>.new.gff   <Abbr>.lens
```

It also writes QC reports, gene/chromosome maps, a rejection log, a SHA-256
manifest and a naming audit, so every retained gene can be traced back to the
original file.

The Python package stops before BLAST and WGDI execution. This repository also
contains the 0–6 execution scripts for the downstream comparison and WGDI
analysis. `wgdi_ready=yes` means the inputs are ready, **not** that WGDI has
already completed.

For a complete end-to-end handoff guide, start from **QUICKSTART.md**.

## Input formats

Suffixes are only hints. Format and biological role are decided from file
content, so uncommon or misleading suffixes are handled as long as the content
is recognizable. Ambiguous candidates are reported as errors rather than
silently guessed.

| Input | Handling |
| --- | --- |
| `.gff` / `.gff3` | Parsed directly; standard, transcript-only, simplified 5-column, gene-span and gene-CDS-direct structures are recognized by content |
| `.gtf` | Converted to GFF3 (gene/transcript/exon/CDS/start/stop codons) |
| `.gbff` / `.gbk` / `.gb` / `.genbank` | Coding GFF3 plus CDS/PEP extracted; genome FASTA can be exported too (requires `biopython>=1.83`) |
| `.annotations.csv.gz` and tab-separated function tables | Normalized to metadata only; **cannot** replace GFF/CDS/PEP or serve as WGDI coordinate input |
| FASTA family: `.fa`, `.fasta`, `.fna`, `.faa`, `.fas` and compound names such as `.genome.fasta.gz`, `_genomic.fna.gz`, `.CDS.fasta.gz`, `.Protein.faa.gz` | Role (genome / CDS / protein) resolved from filename semantics plus a bounded FASTA content preview |
| `.gz`, `.bz2`, `.xz`, `.zip` | Detected by compression magic; a wrong or missing suffix does not block reading |

A raw download directory may contain more than one candidate per role. List
them as `A|B` in the inventory; the tool scores each candidate and selects a
unique one. For single-species registration, `init-project` takes four explicit
paths (`--genome`, `--gff`, `--cds`, `--pep`) and never infers a biological role
from a suffix alone.

## Repository layout

```text
.
├── 0check_data.sh
├── 1run_blast.sh ... 6run_wgdi.sh
├── run_all_pipeline.ps1
├── requirements-wgdi-validated.txt
├── bin/wgdi
├── bin/wgdi_entry.py
├── vendor/sitecustomize.py
├── code/                       # single source tree (this is what gets packaged)
│   ├── src/pgci_preprocess/
│   ├── scripts/ tools/ tests/ docs/
│   ├── resources/species_abbreviations.tsv
│   └── pyproject.toml
└── plant-genome-resources/     # example project data root (NOT code)
    ├── analysis_v2/<Abbr>/
    ├── clean_v2/<Abbr>/
    └── resources/              # index, manifests, audits
```

Raw genome files and generated analysis results are intentionally not tracked.
Keep code and project data in separate directories.

## Expected data layout

Place raw downloads and the inventory under your project root:

```text
my_project/
├── data/downloads/{show_id}_{Order}_{Family}({Species})/
└── resources/download_inventory.tsv
```

A validated species ends with one GFF, one CDS FASTA and one protein FASTA.
`prepare` can derive those roles from GTF, GBFF, or a genome-plus-GFF source
before registration, and records every routing decision in
`resources/preprocess_batches/_routing.tsv`.

## Initialize, inspect dependencies, and execute

Python 3.11 or newer. From the repository root:

```powershell
cd code
python -m pip install -e ".[gbff]"
python -m pgci_preprocess init-project --project-root <new-project>
python -m pgci_preprocess doctor --project-root <new-project> --allow-empty
python -m pgci_preprocess doctor --project-root <project> --workflow full
```

The `gbff` extra (`biopython>=1.83`) is only needed for GenBank/GBFF input.
Python dependencies for preprocessing are otherwise standard-library only.
External programs (Bash, DIAMOND/BLAST, MUSCLE, yn00, Perl) and the WGDI
environment are described in QUICKSTART.md; they are not Python packages and
`doctor` only checks their availability.

`init-project` creates missing nested directories and standard TSV headers.
Registration requires four explicit inputs: `--genome`, `--gff`, `--cds` and
`--pep`. It never overwrites an existing species directory. `doctor` is
read-only: it neither installs tools nor claims that an external program has
successfully completed an analysis.

From the workflow root, the full controller defaults to DIAMOND with
`very-sensitive` and `1e-5`:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all_pipeline.ps1
# Optional overrides:
# -HomologyEngine blast
# -DiamondPath <path> -DiamondSensitivity ultra-sensitive -HomologyEvalue 1e-5
# -CondaEnv <env> -BioinfoBin <external-tool-directory>
```

Generated comparison configs default to `tandem = false`, so tandem-repeat
blocks are removed from correspondence, block-Ks and Ks-peak analyses. Set
`WGDI_TANDEM=true` before regenerating configs if tandem blocks must be
retained. The WGDI wrapper accepts `WGDI_PYTHON`, `WGDI_PERL_BIN` and
`WGDI_PERL5LIB` for non-local installations.

## Standardize and validate

From the repository root:

```powershell
& 'bash' '.\0check_data.sh' `
  --preprocess `
  --refresh-links `
  --expected-count 20
```

On a project that does not need root compatibility links:

```powershell
& 'bash' '.\0check_data.sh' `
  --no-compat-links `
  --expected-count 20
```

The project-specific guide is in `code/docs/data-validation.md`.

## End-to-end smoke receipt

Build an isolated subset from a previously validated bundle and its historical
self-collinearity file. Historical alignments are used only for gene selection;
DIAMOND is rerun with the configured settings.

```powershell
python code/tools/build_wgdi_smoke_fixture.py `
  --workflow-root <workflow-root> --species-dir <workflow-root/Ggl> `
  --self-icl <old/Ggl.icl> --output-root <workflow-root/smoke_wgdi>
bash code/tools/run_wgdi_smoke.sh <workflow-root/smoke_wgdi>
```

The runner executes DIAMOND and scripts 2–6, then requires all eight key WGDI
outputs per comparison. `smoke_report.json` records parameters, file sizes and
SHA-256 hashes. This receipt applies only to the derived smoke fixture, not to
a full rerun of every species; step 0 remains covered by preprocessing
validation.

## Test

```powershell
cd code
python -m pytest -q
ruff check src tests scripts
mypy src --ignore-missing-imports
```

## Skill package

The reusable skill package is published separately at
<https://github.com/HXC-code/plant-genome-preprocess-skill>. Copy the directory
that contains `SKILL.md` into `%USERPROFILE%\.codex\skills\` (Codex) or
`%USERPROFILE%\.workbuddy\skills\` (WorkBuddy).

```text
Use $plant-genome-preprocess to standardize and validate this plant genome dataset.
```
