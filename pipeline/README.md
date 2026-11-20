# HXC Plant Genome Preprocess

Reproducible preprocessing and validation for plant genome annotations. The
pipeline converts NCBI-style GFF3, CDS FASTA, and protein FASTA into a
one-gene-one-transcript WGDI-ready dataset.

Species analysis codes come from the project's canonical abbreviation table:
`show_id` remains the internal key, while directories, BLAST files, WGDI panels,
tree labels and figures use the globally unique abbreviation (for example
`2720 -> Lcu`). Gene IDs follow `{abbr}{sequence}g{order:04d}`, e.g.
`Lcu3g0007`. Legacy `s2720` selectors remain accepted as compatibility aliases,
but newly generated artifacts are always abbreviation-named. **Handing the
package to another lab: start from `QUICKSTART.md`.**

## Scope

This repository covers:

- raw file inventory and validation;
- GFF3, CDS, and protein parsing;
- canonical cross-species gene IDs;
- representative transcript selection;
- CDS/protein translation consistency checks;
- organelle and taxon contamination filtering;
- WGDI GFF, lens, CDS, and protein generation;
- manifest and SHA-256 validation;
- explicit project initialization and read-only dependency/input doctor;
- configurable DIAMOND/BLAST and the 0–6 WGDI execution scripts;
- isolated real-data smoke receipts;
- a reusable Codex skill package.

The `pgci_preprocess` Python package ends before BLAST and WGDI execution. This
workflow repository additionally contains the downstream execution scripts;
`wgdi_ready=yes` never means that WGDI has already completed.

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
│   ├── scripts/preprocess_wgdi_v2.py
│   ├── scripts/validate_wgdi_dataset.py
│   ├── src/pgci_preprocess/
│   ├── tools/ tests/ docs/
│   └── resources/species_abbreviations.tsv
└── plant-genome-resources/     # project data root (NOT code)
    ├── analysis_v2/<Abbr>/
    ├── clean_v2/<Abbr>/
    └── resources/              # index, manifests, audits
```

Raw genome files and generated analysis results are intentionally ignored.

## Expected data layout

Place raw downloads and the inventory under:

```text
plant-genome-resources/
├── data/downloads/{show_id}_{Order}_{Family}({Species})/
└── resources/download_inventory.tsv
```

Each validated species ends with one GFF, one CDS FASTA, and one protein FASTA.
The batch preparation path may derive those roles from GTF, GBFF, or a
genome-plus-GFF source before registration.

## Initialize, inspect dependencies, and execute

Use a new virtual environment. The preprocessing package has no mandatory
non-standard-library dependency unless GBFF conversion is requested; GBFF needs
`biopython>=1.83`. The complete WGDI environment used for the real smoke is
pinned in `requirements-wgdi-validated.txt`. Bash, DIAMOND/BLAST, MUSCLE, yn00,
and Perl are external programs, not Python dependencies.

```powershell
cd plant-genome-resources
python -m pip install -e ".[gbff]"
python -m pgci_preprocess init-project --project-root <new-project>
python -m pgci_preprocess doctor --project-root <new-project> --allow-empty
python -m pgci_preprocess doctor --project-root <project> --workflow full
```

`init-project` creates missing nested directories and standard TSV headers. A
species registration requires four explicit inputs: `--genome`, `--gff`, `--cds`,
and `--pep`. It never guesses biological roles from an arbitrary suffix and does
not overwrite an existing species directory. `doctor` is read-only and does not
install tools or claim that an executable has successfully completed an analysis.

From the workflow root, the full controller defaults to DIAMOND with
`very-sensitive` and `1e-5`:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all_pipeline.ps1
# Optional overrides:
# -HomologyEngine blast
# -DiamondPath <path> -DiamondSensitivity ultra-sensitive -HomologyEvalue 1e-5
# -CondaEnv <env> -BioinfoBin <external-tool-directory>
```

The generated `Auto_pep/diamond_run_settings.tsv` records the actual homology
settings. BLAST results are copied, not moved, into WGDI comparison directories.
The WGDI wrapper accepts `WGDI_PYTHON`, `WGDI_PERL_BIN`, and `WGDI_PERL5LIB` for
non-local installations. Windows users should keep project and pytest temporary
paths short.

Generated comparison configs default to `tandem = false`, so tandem-repeat blocks
are removed from correspondence, block-Ks and Ks-peak analyses. Upstream
`wgdi -r <config>` remains the retention-analysis option. In this workflow's
wrapper, a final standalone `-r` (or any `--remove-tandem`) is the convenience
flag that applies the same tandem filter at runtime. Set `WGDI_TANDEM=true`
before regenerating configs when tandem-repeat blocks must be retained.

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

The project-specific guide is in
`plant-genome-resources/docs/data-validation.md`.

## End-to-end smoke receipt

Build an isolated subset from a previously validated bundle and its historical
self-collinearity file. Historical alignments are used only for gene selection;
DIAMOND is rerun with the configured settings.

```powershell
python plant-genome-resources/tools/build_wgdi_smoke_fixture.py `
  --workflow-root <workflow-root> --species-dir <workflow-root/Ggl> `
  --self-icl <old/Ggl.icl> --output-root <workflow-root/smoke_wgdi>
bash plant-genome-resources/tools/run_wgdi_smoke.sh <workflow-root/smoke_wgdi>
```

The runner executes DIAMOND and scripts 2–6, then requires all eight key WGDI
outputs per comparison. `smoke_report.json` records parameters, file sizes, and
SHA-256 hashes. This receipt applies only to the derived smoke fixture, not to a
full rerun of every species. Step 0 remains covered by preprocessing validation.


## Test

```powershell
cd plant-genome-resources
python -m pytest -q
ruff check src tests scripts
mypy src --ignore-missing-imports
```

## Codex skill

The reusable skill package is under:

```text
skill-packages/plant-genome-wgdi-preprocess/
```

Install it by copying that directory into:

```text
%USERPROFILE%\.codex\skills\
```

Example prompt:

```text
Use $plant-genome-wgdi-preprocess to standardize and validate this plant genome dataset.
```

## Environment notes (2026-11-20)

- `tests/conftest.py` redirects an over-deep `--basetemp` (>90 chars) to a short temp root and emits a `PytestWarning`; deep Windows paths otherwise hit the legacy 260-char limit and produce misleading "directory exists but file cannot be written" failures.
- Shell scripts `0check_data.sh`–`6run_wgdi.sh` now use LF line endings (backup: `99_ARCHIVE\shell_lineendings_20261120`); `bash -n` passes under both Git Bash and WSL/Linux.
- GBFF conversion needs `biopython>=1.83`; `doctor` lists local interpreters that have Biopython when the current one does not.
- `run` fails fast with `path_too_long_for_windows` when a species/input/output path exceeds ~250 characters, instead of surfacing a bare `FileNotFoundError`.

