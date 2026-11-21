# Contributing

This repository is a **published mirror**. The canonical implementation lives in
[HXC-code/plant-genome-preprocess](https://github.com/HXC-code/plant-genome-preprocess).

## Setup

Python 3.11 or newer:

```powershell
cd pipeline/code
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

CI runs the same commands and additionally checks that `pipeline/` matches the
canonical repository. If that check fails, do not hand-patch individual mirrored
files: update the canonical repository first, then run:

```powershell
cd <canonical-repo>/code
python tools/sync_skill_package.py --pipeline-target <this-repo>/pipeline --check
python tools/sync_skill_package.py --pipeline-target <this-repo>/pipeline
```

`SKILL.md` is owned by this repository; everything under `pipeline/` is
generated from the canonical source.

## Rules of thumb

- Never commit raw genome data, generated results or local interpreter paths.
- Keep `SKILL.md` short and task-oriented: when to use the skill, what formats
  are accepted, and what the outputs are.
- Any change to Python code, tests or the species table must land in the
  canonical repository first and then be mirrored here.
