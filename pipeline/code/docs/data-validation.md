# Data Validation Guide

## Prerequisites

- Git Bash at `D:\Git\bin\bash.exe`
- Python available as `python`, or set `PYTHON`
- Canonical data root containing species folders directly under `data/`
  (legacy `data/downloads/` is still read) and
  `resources/download_inventory.tsv`

Do not run the `.sh` file directly from Windows PowerShell. Invoke Git Bash.

## Fresh preprocessing

When derived outputs are absent:

```powershell
cd C:\path\to\HXC-Project

& 'D:\Git\bin\bash.exe' '.\0check_data.sh' `
  --preprocess `
  --refresh-links `
  --expected-count 20
```

This runs:

```text
0check_data.sh
→ Python preprocessing
→ manifest validation and SHA-256 checks
→ compatibility-link checks
→ data-validation-report.tsv
```

## Existing output validation

```powershell
& 'D:\Git\bin\bash.exe' '.\0check_data.sh' --expected-count 20
```

Without legacy `s<show_id>` directories:

```powershell
& 'D:\Git\bin\bash.exe' '.\0check_data.sh' `
  --no-compat-links `
  --expected-count 20
```

## Incremental preprocessing

```powershell
& 'D:\Git\bin\bash.exe' '.\0check_data.sh' `
  --preprocess `
  --species 996 `
  --refresh-links `
  --expected-count 20
```

Incremental runs preserve manifest and summary rows for unselected species.

## Output layers

```text
data/              immutable raw input; one folder per species
resources/         inventory, index, registry, QC, adjudications
analysis_v2/       canonical CDS/PEP/GFF/lens plus all traceability files
                   (gene_map, chr_map, chromosome_map, removed_list, qc)
results/           future BLAST/WGDI analysis output
```

Single-layer layout: no `clean_v2` is written; legacy projects' `clean_v2`
remains readable and can be deleted after a re-run.

No per-file SHA-256 manifest is written any more; input provenance stays in
each species' `*.qc.json` (`inputs` section).

Standardized data shares the project root with raw data but is not written
inside individual raw species directories.

