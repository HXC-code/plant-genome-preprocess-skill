# Changelog

## 3.4.4 — 2026-11-21

- Skill package now contains the full runnable preprocessing pipeline under `pipeline/`, including source, tests, resources, 0–6 scripts and the WGDI compatibility shim.
- Species analysis codes: genus initial plus two specific-epithet letters, widening on collision; table rows keep their published codes, new species are derived from the Latin name and avoid all existing codes.
- Documentation: SKILL.md and the mirrored README/QUICKSTART now describe content-based input routing, the full set of readable formats and how to install from `pipeline/code`. Test count corrected to 151.
- Packaging: MIT `LICENSE` added.
- Documentation: clarified that the 0–6 scripts and WGDI shim are an optional downstream layer; preprocessing ends at `run` / `validate` / `check-names`.
- Engineering: GitHub Actions runs pytest, ruff and mypy and checks that `pipeline/` has not drifted from the canonical repository. The shared pytest temp-root fallback and type-error fixes are mirrored from the canonical source.
