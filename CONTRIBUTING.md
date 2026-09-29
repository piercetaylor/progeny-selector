# Contributing

## Setup

Python 3.11 or later. The package installs shiny and pandas; the extras are `dev` (pytest, pytest-cov, ruff, mypy), `e2e` (pytest-playwright and axe-playwright-python) and `export` (shinylive).

- `dev`: `pip install -e ".[dev]"` is enough for the per-commit gates below. A pull request must pass them.
- `e2e`: the browser tests under `tests/e2e/`. Install with `pip install -e ".[dev,e2e]"`, run `playwright install chromium` once, then `pytest -m e2e`. Plain `pytest` excludes these tests.
- `export`: shinylive, for the static Shinylive export. No per-commit gate needs it.

Local venv recipe (Windows): `py -3.12 -m venv .venv`, then `.venv\Scripts\python -m pip install -e ".[dev,export,e2e]"` and `.venv\Scripts\python -m playwright install chromium`.

After changing `src/progeny_selector/_version.py`, re-run `pip install -e .` so `tests/test_version.py` sees the new metadata.

## Gates

All must pass before a pull request:

```
ruff check .
ruff format --check .
mypy
pytest
python scripts/check_contract.py
python scripts/make_fixture.py
git diff --exit-code -- tests/fixtures src/progeny_selector/examples
```

The generator is deterministic, so regenerating the fixture must leave `tests/fixtures` and `src/progeny_selector/examples` unchanged. Edit the generator, never the generated files.

## Conventions

Commits follow Conventional Commits 1.0.0 (`feat(core): ...`, `fix(io): ...`, `docs: ...`, `test: ...`, `chore: ...`; `!` or a `BREAKING CHANGE:` footer for changes to the data contract or the results columns). Versions follow SemVer 2.0.0; docs/data-formats.md is part of the public interface. CHANGELOG.md follows Keep a Changelog 1.1.0; add a line under Unreleased with every user-visible change. Decisions are MADR files in docs/adr/, numbered sequentially; supersede rather than edit an accepted record.

## Code rules

Everything in src/progeny_selector/core is pure: numpy arrays and dataclasses in, arrays and dataclasses out; no file I/O, no Shiny, no pandas. Validation happens once at the boundary (src/progeny_selector/io) and raises DataContractError or CriteriaError with the file, line and column. The UI (src/progeny_selector/app) calls `run_analysis` and renders rows; it computes nothing. Colours come only from `constants.STATE_COLORS` and `STATUS_COLORS` (Okabe-Ito). Every module starts with a docstring stating its responsibility and interface.

## Data

Never commit genotype data other than the synthetic fixture; .gitignore excludes *.vcf, *.vcf.gz and *.hmp.txt outside tests/fixtures/.

## Releasing

Releases are made by the maintainer only: set the version in `src/progeny_selector/_version.py`, CHANGELOG.md and CITATION.cff, push a `v<version>` tag, and `.github/workflows/release.yml` builds, publishes to PyPI by trusted publishing and creates the GitHub release with the CHANGELOG section as notes.
