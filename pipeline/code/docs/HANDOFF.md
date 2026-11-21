# Handoff Guide

## Product boundary

The project contains two explicitly separated modes:

| Mode | Entry | Output |
| --- | --- | --- |
| Preprocess only | `python -m pgci_preprocess run` then `validate` | WGDI-ready pep/CDS/GFF/lens plus QC |
| Full comparison | `run_homology.sh`, scripts 2-3, config wrappers 4-5, `6run_wgdi.sh` | Pairwise homology, collinearity blocks and Ks evidence |

The name `plant-genome-preprocess` describes the authoring package. The
workflow root additionally provides an optional WGDI execution layer.

## Supported entry points

Use these commands for new work:

```powershell
cd code
python -m pgci_preprocess init-project --project-root D:\my_project
python -m pgci_preprocess run --project-root D:\my_project
python -m pgci_preprocess validate --project-root D:\my_project --refresh-links
python -m pgci_preprocess check-names --project-root D:\my_project
```

The files under `code/scripts/` are deprecated compatibility wrappers. They
forward to the canonical CLI and should not be used in new documentation.

For the optional comparison layer:

```bash
bash run_homology.sh                 # HOMOLOGY_ENGINE=diamond by default
bash 2run_pre_wgdi.sh
bash 3blast_to_wgdi.sh
python code/tools/generate_wgdi_configs.py --mode all
bash 6run_wgdi.sh
```

The numbered wrappers remain supported:

- `1run_diamond.sh` and `1run_blast.sh` select the shared homology runner.
- `4run_wgdi_different.sh` and `5run_wgdi_same.sh` select the shared config
  generator.

## Runtime ownership

- Preprocessing uses the Python package in `code/src/pgci_preprocess`.
- WGDI comes from the pinned `wgdi==0.75` requirement.
- `bin/wgdi_entry.py` owns only the matplotlib and tandem-boolean compatibility
  patches.
- `bin/wgdi` resolves the runtime from `WGDI_PYTHON`, `.wgdi-env`, `.wgdi_python`
  or an importable Python on `PATH`.

## Validation

Read `VALIDATION.md` before publishing. A source tree with passing unit tests is
not automatically an end-to-end validated WGDI release. Release archives ship
LF shell scripts (`.gitattributes` plus builder-side normalization), so the
0–6 layer runs unchanged on Linux/WSL and Windows Git Bash.
