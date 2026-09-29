# First public release packaging: single-sourced version, app command, package data, sdist, trusted publishing, CI matrix

Status: accepted. Date: 2026-09-29. Decided by the Fable planner's recommendations from the audit of 2026-09-29, accepted by the main session under the maintainer's delegation. The joint decisions with backcross (version 0.1.0, provenance columns, citation, same-day release) are in docs/adr/0035.

## Context and Problem Statement

Milestone M4 is the first public release. The audit of 2026-09-29 found that the package could not yet be installed and used by a breeder who has no checkout: the version lived in two files, `shiny run` needed a path into a source tree, the example data was reachable only from the repository, the sdist shipped agent working documents, nothing built the wheel in CI, there was no release workflow, and the citation and community files were missing. The criterion was the one used throughout: reproducible with standard tools, working for academic users on Windows, macOS and Linux, and never silently wrong. Open questions Q1, Q2, Q4, Q5 and Q7 to Q15 are decided here; Q3 and Q6 are decided in docs/adr/0035.

## Considered Options

Each decision below lists its alternative in the Decision Outcome.

## Decision Outcome

- **Version.** One source, `src/progeny_selector/_version.py`, read by hatchling for the package metadata and imported by `__init__.py`. CHANGELOG.md and CITATION.cff carry the same number by hand, and two tests plus one workflow guard fail if the three disagree. Q28: `importlib.metadata.version` reflects the install-time metadata, so after a bump the version test fails until `pip install -e .` is run again. That is accepted and stated in CONTRIBUTING; skipping when metadata is absent would be a vacuous pass.
- **Q1, shiny and pandas.** Core dependencies. The `app` extra is deleted and `dev` lists only tools. `pipx install progeny-selector` then `progeny-selector app` works in one line on all three platforms. Where the import still fails (a `--no-deps` install), `app` prints the exact install line. Alternative not taken: keep the extra and use a self-reference in `dev`.
- **Q2, command shape.** A subcommand of the existing entry point, `progeny-selector app [--host 127.0.0.1] [--port 8000] [--no-browser]`, calling `shiny.run_app` and opening a browser by default. Alternatives not taken: a second console script; `python -m progeny_selector.app`, which does not work for a pipx install.
- **Q4, Load example under Shinylive.** The same Load-example control works in the exported site. `scripts/build_shinylive.py` already stages the whole package, so the example files travel with it and nothing in the build script changes. A browser test proves the requests stay on the page's origin.
- **Q5, documentation.** Plain Markdown on GitHub, `docs/README.md` as the index, absolute links in the README so the PyPI page works. No MkDocs or Sphinx: Pages already serves the app, and a docs site would add a second deploy path, a theme dependency and a build job. Revisit past about a dozen pages.
- **Q7, sdist.** Excludes `CLAUDE.md`, `PLAN.md`, `docs/*-phases.md`, `.github`, `.env.example`, `site`, `build` and `data`. Keeps `tests/`, `docs/`, `contract/`, `scripts/`, `LICENSE`, `CHANGELOG.md` and `CITATION.cff`, so downstream packagers can run the tests. The contract cases with gzip, bgzip and CRLF stay byte-exact.
- **Q8, example command.** `progeny-selector example`, default `--out progeny-selector-example`, refusing to overwrite without `--force`. The examples ship as package data written by `scripts/make_fixture.py`, and a test proves the copies byte-identical to `tests/fixtures/`. The BC3F1 example is made self-contained at write time. `.gitignore` carries an exception for `src/progeny_selector/examples/**`, without which hatchling would drop the ignored VCF from the wheel.
- **Q9, `--top N`.** The behaviour stays (`rank_in_family <= N`, up to N per family; `--overall` for N in total). The help text and README are corrected to say so. Renaming the flag was rejected as a breaking change for no gain.
- **Q10, contacts.** No personal email in the repository. The conduct contact is the maintainer through a GitHub issue or the private reporting form; SECURITY.md uses GitHub private vulnerability reporting.
- **Q11, fixture regeneration.** Windows and macOS run unit tests only. The byte-for-byte regeneration check stays on Ubuntu.
- **Q12, `app/requirements.txt`.** Deleted, because no code path reads it; docs/adr/0009 gets a dated line.
- **Q13, TestPyPI.** No rehearsal. The `package` job installs the built wheel and sdist into fresh venvs on every push, and a failed first publish leaves nothing uploaded.
- **Q14, `unit-os` Python.** 3.12 only, one job each for Windows and macOS.
- **Q15, cffconvert.** `cffconvert --validate` runs in the `package` job. CITATION.cff carries no `date-released` until the maintainer tags, because the CFF 1.2.0 schema requires a date and rejects a placeholder.
- **Release.** `release.yml` publishes to PyPI by trusted publishing (OIDC, no token) through the `pypi` environment. Its `build` job refuses a tag that differs from the package version or a CHANGELOG section without a date (`scripts/changelog_section.py`), and it extracts the release notes from CHANGELOG.md. Nothing in the repository tags or publishes; the maintainer does.
- **CI matrix.** Python 3.11, 3.12 and 3.13 on Ubuntu, plus one Windows and one macOS unit job. Browser tests run on Linux only.

### Consequences

Good: a breeder installs with one line and reaches a working interface and an example without a checkout; the published build is the one CI tested; the citation and the release notes come from files that tests keep in step. Bad: the wheel grows by the example files; `shiny` and `pandas` install for users who only run the command line; a release needs three files edited by hand (`_version.py`, the CHANGELOG heading, CITATION.cff), checked by two tests and one workflow guard; the Windows and macOS jobs do not run the browser tests, which stays a documented limitation.
