# Changelog

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
