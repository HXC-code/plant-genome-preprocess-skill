# Changelog

## 3.5.4 — 2026-11-21

Fixes from the 41-species real-data pilot review (WGDI full-chain execution
and program audit report). Two scientific-correctness/data-reliability issues,
three medium and two minor robustness issues are addressed.

- Organelle inference for genome-less annotations: simplified/transcript-only
  GFFs without a genome FASTA previously could not identify chloroplast or
  mitochondrial contigs, so short organelle sequences entered the nuclear WGDI
  set. `run` now infers organelles among *synthetic* regions from direct
  identifier signals (chloroplast/plastid/mitochondrion names, `cp`/`chrc`/`mt`
  accessions, plastid gene ids) or a dense, compact, mostly intronless gene set.
  Ambiguous dense regions are listed as `review_candidates` instead of being
  silently excluded. Authoritative `region` rows and genome FASTA records are
  never reclassified. New flags: `--organelle-pattern <regex>` for explicit
  synthetic-sequence exclusion and `--no-organelle-inference` to disable it.
  QC records the evidence under `organelle_inference` and counts excluded genes
  in `n_organelle_genes`. Verified on the pilot: Hevea brasiliensis
  `NC_015308.1` (83 genes) and Populus alba `NC_008235.1` (83 genes) are now
  excluded from the WGDI set.
- Readiness hard gate: `validate` now fails bundles whose QC reports
  `wgdi_ready=no` (`wgdi_not_ready_low_coverage`); `--allow-not-ready` (also
  `0check_data.sh --allow-not-ready`) restores archival/data-integrity checks.
  `run_homology.sh` and `run_all_pipeline.ps1` skip `ready=no` species and
  species without a readable readiness record; `HOMOLOGY_ALLOW_NOT_READY=true`
  overrides explicitly. Pilot acceptance: 36 PASS / 5 FAIL strict, 41 PASS with
  `--allow-not-ready`.
- `init-project` registers three-piece bundles (GFF + CDS + PEP) without a
  genome FASTA; `--genome` and `--assembly-id` are now optional.
- Generated WGDI `run.sh` entries use Bash-safe POSIX paths (`/d/...`) instead
  of Windows backslashes, so Git Bash / Linux no longer mis-escape the entry.
- `bin/wgdi` normalizes Windows `TMPDIR` values and falls back to a
  project-local temp root when the system temp dir is not writable.
  `doctor --workflow wgdi|full` now actually launches the WGDI entry
  (`wgdi_probe`) instead of only checking existence.
- `run_homology.sh` no longer scans infrastructure directories (`bin/`,
  `code/`, `resources/`, ...) as species, tolerates temp-file cleanup failure
  without polluting the exit code, and writes per-process `.part.<pid>` outputs.
- The synthetic smoke fixture records `wgdi_ready=yes` QC so the readiness gate
  is exercised on its pass path; a new functional test executes
  `run_homology.sh` against ready / not-ready / unrecorded species.
- Engineering: 172 tests, Ruff and Mypy pass; Windows synthetic full-chain smoke
  receipt `validation/smoke-report-3.5.4.json` (4 comparisons × 8 products).

## 3.5.3 — 2026-11-21

- Pinned Ruff 0.12.0 and Mypy 1.17.1 in both program and Skill CI so newer
  default lint rules cannot make previously passing commits fail without a
  deliberate repository change.
- Declared the accepted Python lint rule set explicitly in `pyproject.toml`.
- Forced UTF-8 mode in every shell wrapper that forwards to Python, fixing
  Windows cp1252 failures in `0check_data.sh` and related entry points.

## 3.5.2 — 2026-11-21

- Added a strict-LF Linux/WSL full-chain smoke and versioned Linux receipt.
- Added real-species acceptance using NCBI RefSeq Arabidopsis thaliana
  TAIR10.1 and Arabidopsis lyrata v1.0, including a real cross-species WGDI
  subset receipt.
- Stopped bundling PAL2NAL v14. Added a verified official-site installer
  (`code/tools/install_pal2nal.py`) so release archives no longer redistribute
  a historical archive whose license text is absent.
- Improved smoke receipts to resolve external tools from their `WGDI_*`
  environment paths and to record PAL2NAL provenance.
- Added a real-bundle subset fixture builder for repeatable cross-species
  validation.

## 3.5.1 — 2026-11-21

- Release archives now guarantee LF line endings for every shell script and
  `bin/wgdi`: the builder normalizes CRLF at packaging time, and `.gitattributes`
  keeps LF in every checkout. Fixes Linux/WSL execution of steps 0/2/3/6 from
  archives built on a Windows checkout (Git Bash tolerates CRLF; Linux bash
  rejects it). Third-party bytes (PAL2NAL) are never rewritten.
- Release ZIP entries now carry Unix executable bits for shell scripts.
- The smoke runner captures a public receipt (`smoke_report.public.json`) in the
  same environment as the run. The shipped 3.5.0 receipt's perl entry did not
  match the archived run; it is replaced by `validation/smoke-report-3.5.1.json`,
  regenerated from the actual validation run with `--public` (all 32 output
  hashes verified identical to the run report).
- Third-party notices: added the bundled PAL2NAL v14 script with attribution,
  citation, upstream link and an explicit license caveat (upstream ships no
  license file; confirm terms before public redistribution).
- Release tests derive the version and receipt name from the package version,
  assert that archived shell entries contain no CRLF, and unit-test the
  builder's normalization.
- Mirror drift check and sync now ignore local tool caches and environments
  (`.pytest_cache`, `.mypy_cache`, `.ruff_cache`, `__pycache__`, `.venv`,
  `.wgdi-env`); previously the skill CI drift check failed right after running
  the test suite, and a sync could delete a local WGDI environment.

## 3.5.0 — 2026-11-21

- Consolidated DIAMOND and BLAST into one shared homology runner while
  preserving `1run_diamond.sh` and `1run_blast.sh` as compatibility wrappers.
- Consolidated same-species and cross-species WGDI configuration generation
  into `code/tools/generate_wgdi_configs.py`.
- Removed the vendored WGDI source tree. The workflow now uses the pinned
  `wgdi==0.75` package and keeps only the two required launch patches.
- Replaced the obsolete v1 pipeline guide and stale review README with a current
  handoff guide and explicit validation/release gate.
- Added root README, license, changelog, third-party notices and validation
  instructions to release archives.
- Extended mirror checks to cover workflow shell/PowerShell scripts, not only
  Python files.
- Added release archive and clean extraction tests.
- Pinned pandas 2.3.3 after a real WGDI 0.75 smoke exposed a pandas 3
  string-dtype incompatibility.
- Made upstream WGDI tool paths portable by resolving MUSCLE, yn00, PAL2NAL and
  related programs from environment variables or `PATH` at launch time.
- Replaced the ambiguous generated `-r` flag with `--remove-tandem`.
- Added a self-contained synthetic WGDI smoke fixture and committed the passing
  public receipt as `validation/smoke-report-3.5.0.json`.

## 3.4.4 — 2026-11-21

- Species analysis codes: genus initial plus two specific-epithet letters, widening to three/four or more letters on collision, with a numeric suffix only as a last resort. Table rows keep their published codes; new species are derived from the Latin name and avoid every existing code.
- Documentation: README and QUICKSTART rewritten around the actual job — preprocessing the user's own raw downloads into WGDI inputs. Input-format coverage now reflects content-based detection, compression magic and `A|B` candidate routing. Install paths, Python requirement (>= 3.11), smoke-test paths and skill-package location corrected.
- Packaging: MIT `LICENSE` added to both repositories.
- Engineering: `pytest` (151 tests), `ruff` and `mypy` all pass. GitHub Actions added for push/PR. `tests/conftest.py` now redirects pytest when the default temp root is unwritable or too deep. `sync_skill_package.py` gained `--pipeline-target` to check/sync the released `pipeline/` mirror.
- Documentation: clarified that the bundled 0–6 scripts and WGDI shim are an optional downstream layer; the preprocessing core ends at `run` / `validate` / `check-names`.
- Minor: removed an unused `csv` import in `naming_audit.py` and fixed three type errors in `artifact_migration.py` / `cli.py`.
