# Contributing

## Setup

Python 3.11 or later. The package has four extras:

- `app`: shiny and pandas, to run the UI.
- `dev`: the `app` dependencies plus pytest, pytest-cov, ruff and mypy. `pip install -e ".[dev]"` is enough for the per-commit gates: `ruff check .`, `ruff format --check .`, `mypy` and `pytest` must pass before a pull request.
- `e2e`: pytest-playwright, for the browser tests under `tests/e2e/`. Install with `pip install -e ".[dev,e2e]"`, run `playwright install chromium` once, then `pytest -m e2e`. Plain `pytest` excludes these tests.
- `export`: shinylive and build, for the static Shinylive export. No test gate needs it.

Regenerate the fixture with `python scripts/make_fixture.py`; CI fails if the committed fixture differs from the generator's output.

## Conventions

Commits follow Conventional Commits 1.0.0 (`feat(core): ...`, `fix(io): ...`, `docs: ...`, `test: ...`, `chore: ...`; `!` or a `BREAKING CHANGE:` footer for changes to the data contract or the results columns). Versions follow SemVer 2.0.0; docs/data-formats.md is part of the public interface. CHANGELOG.md follows Keep a Changelog 1.1.0; add a line under Unreleased with every user-visible change. Decisions are MADR files in docs/adr/, numbered sequentially; supersede rather than edit an accepted record.

## Code rules

Everything in src/progeny_selector/core is pure: numpy arrays and dataclasses in, arrays and dataclasses out; no file I/O, no Shiny, no pandas. Validation happens once at the boundary (src/progeny_selector/io) and raises DataContractError or CriteriaError with the file, line and column. The UI (src/progeny_selector/app) calls `run_analysis` and renders rows; it computes nothing. Colours come only from `constants.STATE_COLORS` and `STATUS_COLORS` (Okabe-Ito). Every module starts with a docstring stating its responsibility and interface.

## Data

Never commit genotype data other than the synthetic fixture; .gitignore excludes *.vcf, *.vcf.gz and *.hmp.txt outside tests/fixtures/.
