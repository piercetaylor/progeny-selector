# Screens are tested as pure helpers plus browser tests; shinylive and Playwright live in their own extras

Status: accepted. Date: 2026-09-12.

## Context and Problem Statement

M0 tested the compute core and the parsers but only served the Shiny shell; no test exercised a screen. M1 adds screens whose behaviour needs checking, and some of it (uploads, reactive wiring, rendered text) exists only in a browser. At the same time the `dev` extra carried shinylive, whose dependency lzstring 1.0.4 ships only an sdist that failed to build in the scaffold container, so a packaging problem unrelated to the code could break the per-commit gates. How are screens tested, and which dependencies does each gate need?

## Considered Options

1. Unit-test `@module.server` functions with a mocked session, keep one `dev` extra.
2. Put logic the screens need into pure helpers in `core/`, `io/` and `app/present.py` with unit tests; test server wiring end to end through `shiny.pytest.create_app_fixture` and the `shiny.playwright.controller` classes; mark browser tests `e2e` and exclude them by default; split the extras into `dev`, `export` and `e2e`.
3. Browser tests only, run on every commit.

## Decision Outcome

Option 2. Helpers are where the logic worth checking lives, and they test without a browser. The wiring is checked where it runs: `create_app_fixture` starts `shiny run` in a subprocess with the parent environment merged under any `env` it is given, and the controllers drive real inputs. Every item under `tests/e2e/` gets the `e2e` marker in `tests/e2e/conftest.py`, `addopts = "-q -m 'not e2e'"` keeps them out of plain `pytest`, and `pytest -m e2e` runs them per milestone and in the CI `e2e` job. `dev` no longer installs shinylive; `export = ["shinylive>=0.8", "build"]` and `e2e = ["pytest-playwright>=0.5"]` are installed only by the jobs that need them.

### Consequences

Good: `pytest` stays fast and needs no browser; lint, type and unit gates cannot fail because lzstring does not build; each screen has at least one test that uploads real files and reads real output. Bad: browser tests need a Chromium matching the revision pinned by the installed Playwright, which is a separate download of about 115 MiB; mocked-session tests of render functions are not written, so a regression in wiring is caught only when `pytest -m e2e` runs. Neutral: without Playwright installed, `tests/e2e/` is skipped at collection rather than failing, so the `check` job reports it as skipped.
