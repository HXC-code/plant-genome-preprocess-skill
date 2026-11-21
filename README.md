# Plant Genome Preprocess Skill

This repository is the installable Skill package. The complete runnable
program and release archive live in
[HXC-code/plant-genome-preprocess](https://github.com/HXC-code/plant-genome-preprocess).

The Skill provides two explicit modes:

| Mode | Purpose |
| --- | --- |
| Preprocess only | Convert genome/CDS/protein and GFF3/GTF/GBFF into WGDI pep, CDS, new.gff and lens bundles with QC, naming and hashing. |
| Full comparison | Optionally continue through DIAMOND/BLAST, WGDI collinearity, block-Ks and Ks-peak analysis. |

`pipeline/` is a versioned mirror of the canonical program repository. Edit the
canonical repository first, then synchronize the Skill:

```powershell
python pipeline/code/tools/sync_skill_package.py `
  --pipeline-target pipeline
```

## Install

Copy this repository directory into the skills directory for your Codex host:

```text
%USERPROFILE%\.codex\skills\plant-genome-preprocess
```

The folder must contain `SKILL.md` at its root. Restart the client if the new
Skill is not discovered automatically.

## Invoke

```text
Use $plant-genome-preprocess to standardize and validate this plant genome dataset.
```

For the optional comparison mode, explicitly request:

```text
Use $plant-genome-preprocess end-to-end: preprocess, DIAMOND/BLAST, WGDI
collinearity and Ks analysis.
```

## Runtime

- Python 3.11 or newer.
- WGDI 0.75 and external alignment tools are installed separately.
- PAL2NAL v14 is not redistributed. Run
  `python pipeline/code/tools/install_pal2nal.py`, or let `bootstrap_wgdi`
  download it from the official source with SHA-256 verification.
- Linux/WSL and Windows Git Bash are supported.

## Validation

The mirrored pipeline contains:

- Windows and Linux full-chain smoke receipts;
- real RefSeq Arabidopsis thaliana and Arabidopsis lyrata validation;
- real-derived cross-species WGDI subset receipt;
- 161 automated tests and CI mirror checks.

Read `SKILL.md` for the operational instructions and `pipeline/VALIDATION.md`
for release gates and evidence.
