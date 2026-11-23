# Handoff Guide

## Product boundary

Scope follows the lab's authoritative workflow (`code/docs/WORKFLOW.md`,
整体流程): steps 1–4 are implemented; steps 5–8 (phylogeny, literature
cross-check, multiple collinearity, circos) are out of scope and must be
disclosed before analysis.

| Mode | Entry | Output |
| --- | --- | --- |
| Preprocess only | `python -m pgci_preprocess run` then `validate` | WGDI-ready pep/CDS/GFF/lens plus QC, chromosome_map/removed_list/qc_record |
| Full comparison | `run_homology.sh`, `2gen_wgdi_conf.sh`, `3run_wgdi.sh` | Pairwise homology, dotplots and collinearity blocks; Ks evidence only with `WGDI_KS=true` |

SOP rules enforced by the engine: representative transcript = longest CDS,
ties by stable transcript ID (never protein length or file order); CLI sequence
selection defaults to `chromosomes` (fail-closed `no_chromosome_sequences`
without chromosome evidence; `--sequence-mode adaptive` is the escape hatch);
Ks stages are opt-in because upstream blockinfo/correspondence consume Ks
values. Raw species folders sit directly under ``data/`` (legacy
``data/downloads`` is still read). Per-file SHA-256 manifests and the
``source-manifest`` command were removed by explicit decision; input
provenance stays in each species' ``*.qc.json`` and adjudication records
remain bound to the judged QC version and fingerprint, so a re-run still
invalidates stale sign-offs.

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
bash 2gen_wgdi_conf.sh all           # flat configs in wgdi_input/config, relative paths
bash 3run_wgdi.sh                    # core: .png + .icl (stage-level resume)
# Optional Ks branch (both steps read the same switch):
WGDI_KS=true bash 2gen_wgdi_conf.sh all
WGDI_KS=true bash 3run_wgdi.sh       # adds .ks/.bi.csv/.bi.c.csv/.bk/.kp/.pf
```

The homology runner gates on the preprocessing QC: species with
`wgdi_ready=no` (or no readable readiness record) are skipped and listed in
`metadata/skipped_species.txt`; `HOMOLOGY_ALLOW_NOT_READY=true` overrides explicitly.
`validate` applies the same gate (`--allow-not-ready` for archival checks).

The numbered wrappers remain supported:

- `1run_diamond.sh` and `1run_blast.sh` select the shared homology runner.
- `2gen_wgdi_conf.sh` writes flat `<pair>.conf` files to `wgdi_input/config/`
  for `same`, `different` or `all` pairs; `3run_wgdi.sh` is the single shared
  executor. Results go to `wgdi_output/results/` with logs in `logs/` and
  resume markers in `state/`; no per-comparison folder or `run.sh` is emitted.

## Runtime ownership

- Preprocessing uses the Python package in `code/src/pgci_preprocess`.
- WGDI comes from the pinned `wgdi==0.75` requirement.
- `bin/wgdi_entry.py` owns only the matplotlib and tandem-boolean compatibility
  patches.
- `bin/wgdi` resolves the runtime from `WGDI_PYTHON`, `.wgdi-env`, `.wgdi_python`
  or an importable Python on `PATH`.
- PAL2NAL is installed separately with
  `python code/tools/install_pal2nal.py`; it is not included in release
  archives because the historical v14 distribution has no explicit license
  text. `bootstrap_wgdi` invokes the installer by default.

## Validation

Read `VALIDATION.md` before publishing. A source tree with passing unit tests is
not automatically an end-to-end validated WGDI release. Release archives ship
LF shell scripts (`.gitattributes` plus builder-side normalization), so the
0–6 layer runs unchanged on Linux/WSL and Windows Git Bash.
