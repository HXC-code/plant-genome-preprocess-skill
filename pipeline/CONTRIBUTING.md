# Contributing

## Setup

Python 3.11 or newer. Work from the repository root:

```powershell
cd code
python -m venv .venv
# activate .venv
python -m pip install -e ".[test]" ruff mypy
```

## Before opening a pull request

```powershell
python -m pytest -q
python -m ruff check src tests scripts tools
python -m mypy src --ignore-missing-imports
```

All three must pass. CI runs the same commands on push and pull requests. If the
default pytest temp root is unwritable or too deep, `conftest.py` redirects it
automatically; an explicit `--basetemp=<short-path>` is still accepted.

## Where things live

- Canonical code: `code/src/pgci_preprocess/`
- Tests: `code/tests/`
- Published species table: `code/resources/species_abbreviations.tsv`
- Skill package: a separate repository, mirrored under `pipeline/`

Do not hand-edit the skill mirror. After changing canonical code, docs or tests:

```powershell
cd code
python tools/sync_skill_package.py --pipeline-target <skill-repo>/pipeline --check
python tools/sync_skill_package.py --pipeline-target <skill-repo>/pipeline
```

The skill repository's CI checks for drift between its `pipeline/` tree and this
repository.

## Rules of thumb

- Never commit raw genome data, generated results or local interpreter paths.
- Format and biological role are decided from content; suffixes are hints only.
  Ambiguous candidates must fail loudly, not be guessed.
- `show_id` is an internal key. Table rows keep their published abbreviation;
  new species derive a code from the Latin name without colliding with any
  existing code.
- Keep user-facing docs focused on the job. Environment details and
  troubleshooting belong in `QUICKSTART.md` and `code/README.md`.
- Use short, imperative commit messages; describe the user-visible effect.
