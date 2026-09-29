# M4: first public release

This document decomposes milestone M4, the first public release of progeny-selector, into phases, one commit each. Scope is packaging and first-run usability only, as decided from the audit of 2026-09-29 (`C:/Users/pierc/AppData/Local/Temp/claude/C--Users-pierc-Projects-Plant-Breeding-Projects/c7f725a0-f419-45f6-a1fb-8200abe9211c/scratchpad/audit-2026-09-29.md`). Every commit passes the gates of CLAUDE.md ("Gates"), run with `.venv/Scripts/python.exe`. Dates are 2026-09-29. Facts that could not be verified this session are marked "not verifiable here" with the test or job that settles them. Doers receive one phase's line range, make no design decisions, never run git, never tag, never publish; each phase names its files, exact names and values, its tests with a discriminating input (one that fails if the behaviour is removed), its model, and its acceptance criterion. Every CHANGELOG line except phase 8's restructuring is written by the main session.

State at spec time: HEAD `429b0a8` plus the uncommitted results_schema 1.2.0 change (docs/adr/0033, docs/adr/0034, untracked), treated as landed before phase 1. Gate baseline from PLAN.md Handoff: 617 unit tests, `pytest -m e2e` 43 passed 1 skipped, contract mirror 1.12.0. `constants.RESULTS_SCHEMA` is `"1.2.0"`; `io/export.py` `FIXED_COLUMNS` has 38 fixed names before `token_profile` (`tests/test_results_schema.py:35-75`).

Invariants carried from M1-M3 (`docs/m3-phases.md:5`): `core/` pure; `io/` validates; `app/` computes nothing; `pytest` browser-free; screen interface `view(id)` / `server(id, state)`; fixture files authored only by `scripts/make_fixture.py`; reviewer on every `core/`, `io/`, `model/` and `contract/` diff plus one end-of-milestone pass; at most two concurrent implementers on disjoint files; no AI attribution anywhere; ADRs numbered from `0036`. **M4 introduces no contract version** and, apart from the conditional phase 0b, no `results_schema` change.

Everything below is a documented limitation in the release notes and README (phases 5 and 8), never a phase: Q14 duplicate semantics (family-blind IBS), Question B (no real two-generation dataset), Shinylive BrAPI and CORS, intercross ranking, multi-donor, non-soybean chromosome lengths, the VCF 4.4 phase indicator, the accessibility gaps in `docs/accessibility.md`.

## Decisions of 2026-09-29 (supersede the text below where they differ)

Taken by the main session after the spec was written, from the joint release research run with backcross (docs/adr/0035 here, backcross docs/adr/0030) and the Fable planner's own recommendations. Doers apply these; where a phase below says otherwise, this block wins.

- **Phase 0b is superseded.** The results_schema 1.2.0 change lands with three provenance columns, `tool`, `tool_version`, `tool_commit` (41 fixed names, 42 in an empty file's header; selected.csv carries them after `sample_db_id`), not the single `tool_version` of Q3. `tool_commit` is `g` plus seven hex, suffixed `-dirty` when tracked files differ from HEAD, and `NA` with no git; it is looked up at run time by `src/progeny_selector/provenance.py`, with `scripts/build_shinylive.py` writing `_build_commit.txt` into the staged copy. Nothing fails to import without a generated file, which answers the fragility Q3 raised. The cap column is `background_max_marker_coverage` (the 09-27 draft called it `background_max_coverage`); this spec already uses the new name. Phase 1 still creates `_version.py` as the single version source and must keep the provenance import chain free of cycles.
- **ADR numbering.** `docs/adr/0035` is the release decisions and provenance record, written with 1.2.0. Phase 8's ADR is `docs/adr/0036-first-public-release-packaging.md`.
- **Placeholders resolved.** `{{VERSION}}` = `0.1.0` (both tools; 1.0.0 waits for frozen columns and flags and a second cross or crop validated on real data). `{{SAME_DAY}}` = yes: both tools tag v0.1.0 the same day and each release note names input data contract 1.12.0 and links the sibling's release. `{{BACKCROSS_REFERENCE}}` = included, not commented out. `{{CITATION_AUTHORS}}` = `family-names: Taylor`, `given-names: Pierce`, no ORCID or affiliation unless the maintainer adds one. The CITATION.cff content is in docs/adr/0035 (title, url, keywords, the `references` entry for backcross); no `.zenodo.json`, because Zenodo ignores CITATION.cff when both exist. `{{RELEASE_DATE}}` stays the maintainer's.
- **Q6:** `Development Status :: 3 - Alpha`, the joint decision, not `4 - Beta`.
- **Q1, Q2, Q4, Q5, Q7 to Q15, items 16 to 26, 28:** accepted as recommended.
- **Item 27:** the main session restructures CHANGELOG.md itself in phase 8; no doer edits CHANGELOG.md.

## Placeholders (defined once; every phase refers to these names)

- `{{VERSION}}`: the release version. The repository declares `0.1.0` today (`pyproject.toml:7`, `src/progeny_selector/__init__.py:10`); the joint decision with backcross on `0.1.0` vs `1.0.0` is pending. Doers never write a version literal taken from this spec: after phase 1 the only code source is `src/progeny_selector/_version.py`, and CHANGELOG.md and CITATION.cff copy whatever it holds when phase 8 runs. The maintainer changes it in one edit (Maintainer steps, step 1), and the tests of phases 1 and 8 fail if the three files then disagree.
- `{{RELEASE_DATE}}`: `YYYY-MM-DD` of the tag. Written as the literal text `unreleased` by phase 8 in CHANGELOG.md and CITATION.cff; `scripts/changelog_section.py` refuses that literal, so the release workflow cannot run until the maintainer sets it.
- `{{CITATION_AUTHORS}}`: the `authors:` block of CITATION.cff (given name, family name, ORCID, affiliation). Phase 8 writes `- family-names: Taylor` / `given-names: Pierce` only; the maintainer completes it.
- `{{BACKCROSS_REFERENCE}}`: whether CITATION.cff carries a `references:` entry for backcross, and its `repository-code` and version. Phase 8 writes the entry commented out.
- `{{SAME_DAY}}`: whether backcross releases the same day. Affects only the README "Related tool" sentence (phase 8): "released together" or "see its own release notes".

## Order, dependencies and parallelism

| phase | scope | model | reviewer | depends on | files (collision set) |
| --- | --- | --- | --- | --- | --- |
| 0 | Precondition: 1.2.0 committed, tree clean, venv `[dev,export,e2e]`, site built | none | no | nothing | none |
| 0b | (Q3, conditional) `tool_version` column folded into the uncommitted 1.2.0 diff | opus | yes | Q3 answered yes, before 0 | core/pipeline.py, io/export.py, scripts/read_results.R, tests/test_results_schema.py, docs/data-formats.md, docs/adr/0033 |
| 1 | Version single source, `--version`, `--top` help text, pyproject metadata, sdist contents | sonnet | no | 0 | pyproject.toml, src/progeny_selector/_version.py, __init__.py, cli.py, tests/test_version.py, tests/test_cli_select_top.py |
| 2 | Examples as package data, `example` command, generator writes both copies | opus | no (gates) | 1 | scripts/make_fixture.py, src/progeny_selector/examples/**, cli.py, .gitignore, .github/workflows/ci.yml (one pathspec line), tests/test_examples.py |
| 3 | `progeny-selector app`, Load-example control, Shinylive proof | sonnet | no | 2 | cli.py, app/screens/load.py, tests/test_cli_app.py, tests/e2e/test_load_screen.py, tests/e2e/test_shinylive_export.py, tests/e2e/helpers.py, docs/keyboard-walkthrough.md |
| 4 | CI: 3.13, Windows and macOS unit jobs, `package` job with the fresh-venv smoke; `release.yml`; `scripts/changelog_section.py` | sonnet | no | 3 | .github/workflows/ci.yml, .github/workflows/release.yml, scripts/changelog_section.py, tests/test_changelog_section.py |
| 5 | README for breeders, absolute links, `--top` semantics, Known limitations | sonnet | no | 3 | README.md, tests/test_readme_links.py |
| 6 | docs/tutorial.md, docs/glossary.md, docs/README.md index, tutorial-commands test | sonnet | no | 3 | docs/tutorial.md, docs/glossary.md, docs/README.md, tests/test_docs_commands.py |
| 7 | Community files, CONTRIBUTING.md, SECURITY.md, CODE_OF_CONDUCT.md, templates, dependabot | sonnet | no | 1 | .github/ISSUE_TEMPLATE/**, .github/PULL_REQUEST_TEMPLATE.md, .github/dependabot.yml, SECURITY.md, CODE_OF_CONDUCT.md, CONTRIBUTING.md, tests/test_github_templates.py |
| 8 | CHANGELOG release section, CITATION.cff, cross-links, ADR 0036, docs consistency, verification block | sonnet, then main session | end-of-milestone pass | 4, 5, 6, 7 | CHANGELOG.md, CITATION.cff, README.md (Citing, Related), .github/workflows/ci.yml (cffconvert step), tests/test_changelog_section.py, docs/adr/0035, docs/adr/0003 amendment, PLAN.md, CLAUDE.md |

Serial chain on `cli.py`: 1 -> 2 -> 3. Concurrent pairs on disjoint files: 5 with 6; 5 with 7; 6 with 7; 4 with 7. Phase 4 must follow 3 because the smoke runs `example` and `app`. Phase 8 is last because it copies the version, ties CITATION.cff to it, and edits the shared CHANGELOG.md.

## Phase 0: precondition

No work. The main session commits the uncommitted results_schema 1.2.0 change (with phase 0b folded in if Q3 is answered yes) and confirms: `pytest` green, `ruff check .`, `ruff format --check .`, `mypy` clean, `python scripts/check_contract.py ../backcross` identical, `python scripts/make_fixture.py` then `git diff --exit-code -- tests/fixtures` clean, `python scripts/build_shinylive.py` has produced `site/`, `.venv` carries `[dev,export,e2e]` and Chromium. `docs/adr/` ends at `0034`.

## Phase 0b (SUPERSEDED 2026-09-29; see Decisions) (conditional on Q3 = yes): `tool_version` in results_schema 1.2.0

Runs only if Q3 is answered yes, and only inside the uncommitted 1.2.0 diff, before phase 0 commits it; otherwise the column would be schema 1.3.0 and is not part of M4. Model opus; reviewer yes (io, core).

Files: `src/progeny_selector/core/pipeline.py`, `src/progeny_selector/io/export.py`, `scripts/read_results.R`, `tests/test_results_schema.py`, `docs/data-formats.md` (results.csv and selected.csv tables), `docs/adr/0033-results-schema-1.2-cap-and-end-rule.md` (a dated paragraph under Decision Outcome), `src/progeny_selector/_version.py` (created here rather than in phase 1 if 0b runs; phase 1 then only re-points `__init__.py` and pyproject).

- `src/progeny_selector/_version.py` (leaf module, no imports): `__version__ = "0.1.0"` with a module docstring "Single source of the package version. Read by hatchling (pyproject `[tool.hatch.version]`), re-exported as `progeny_selector.__version__`, and written into results.csv and selected.csv as `tool_version`. No imports: `core/pipeline.py` imports it, and `progeny_selector/__init__.py` imports `core.pipeline`."
- `core/pipeline.py`: `from progeny_selector._version import __version__`; every row dict gains `"tool_version": __version__` immediately after `"results_schema"`.
- `io/export.py`: `FIXED_COLUMNS` gains `"tool_version"` after `"background_max_marker_coverage"` and before `"token_profile"` (39 fixed names; the empty header has 40); `SELECTION_COLUMNS` gains `"tool_version"` after `"sample_db_id"` and before `"token_profile"`; `selection_csv_text` writes `r.get("tool_version", "")`. The module docstring lines 23-24 name the column. Placement follows docs/adr/0016 line 29.
- `scripts/read_results.R`: `tool_version = col_character()` after `background_max_marker_coverage`; the comment at line 7 names it as the last fixed column.
- Tests (`tests/test_results_schema.py`): `FIXED_PREFIX` gains `"tool_version"` last (the file's comment says to edit it twice on purpose); `test_fixed_prefix_is_38_names` becomes `..._is_39_names`; `test_metadata_columns_on_every_row` asserts `row["tool_version"] == progeny_selector.__version__`; `test_selection_csv_header_and_na` asserts `cells[SELECTION_COLUMNS.index("tool_version")] == __version__`. Discriminating input: `test_tool_version_is_never_empty_or_na`: `results_csv_text(run_analysis(make_dataset(make_matrix(["A","H","A","A"])), criteria).rows)` parsed by `csv.DictReader` has `tool_version` equal to `__version__` on the row and not `""` or `"NA"`; removing the pipeline key makes the writer emit `""` and the test fails.
- No commit hash. Grounds recorded in ADR 0033's paragraph: a pip- or pipx-installed package has no repository; deriving a hash at build time (hatch-vcs) makes a fresh clone or the Shinylive staged copy (`scripts/build_shinylive.py` copies `src/progeny_selector` uninstalled) fail to import without a generated file; the reproducible provenance of a release is `tool_version` plus `results_schema`; a run from a checkout is identified by its commit in the user's own notes. The `r-reader` job reads the new column unchanged.

Acceptance: `tests/test_results_schema.py` green with 39 fixed names; `Rscript scripts/read_results.R` reads a results.csv carrying `tool_version`; `constants.RESULTS_SCHEMA` still `"1.2.0"`.

## Phase 1: version single source, `--version`, `--top` help text, pyproject metadata

Model sonnet; gates only. Files: `pyproject.toml`, `src/progeny_selector/_version.py` (new unless 0b ran), `src/progeny_selector/__init__.py`, `src/progeny_selector/cli.py`, `tests/test_version.py` (new), `tests/test_cli_select_top.py` (new).

1. `src/progeny_selector/_version.py` as in phase 0b (no `tool_version` wiring). `__init__.py`: replace `__version__ = "0.1.0"` with `from progeny_selector._version import __version__`; `__all__` unchanged.
2. `pyproject.toml` `[project]`: delete `version = "0.1.0"`; add `dynamic = ["version"]`; add
   ```toml
   [tool.hatch.version]
   path = "src/progeny_selector/_version.py"
   ```
   Classifiers become exactly:
   ```
   "Development Status :: 3 - Alpha",       # Q6, decided 2026-09-29
   "Intended Audience :: Science/Research",
   "Operating System :: OS Independent",
   "Programming Language :: Python :: 3 :: Only",
   "Programming Language :: Python :: 3.11",
   "Programming Language :: Python :: 3.12",
   "Programming Language :: Python :: 3.13",
   "Topic :: Scientific/Engineering :: Bio-Informatics",
   ```
   (`License :: OSI Approved :: MIT License` removed: PEP 639, `license = "MIT"` is the SPDX field; hatchling 1.25+ emits `License-Expression`.) `[project.urls]` becomes:
   ```toml
   Homepage = "https://github.com/piercetaylor/progeny-selector"
   Documentation = "https://github.com/piercetaylor/progeny-selector/blob/main/docs/README.md"
   Repository = "https://github.com/piercetaylor/progeny-selector"
   Issues = "https://github.com/piercetaylor/progeny-selector/issues"
   Changelog = "https://github.com/piercetaylor/progeny-selector/blob/main/CHANGELOG.md"
   "Browser app" = "https://piercetaylor.github.io/progeny-selector/"
   ```
   Dependencies per Q1 (recommended branch): `dependencies = ["numpy>=1.26", "pyyaml>=6.0", "shiny>=1.4", "pandas>=2.2"]`; `app` extra removed; `dev = ["pytest>=8", "pytest-cov>=5", "ruff>=0.6", "mypy>=1.10"]`. If Q1 is answered "keep the extra": `app` stays, `dev = ["progeny-selector[app]", "pytest>=8", ...]` (PEP 508 self-reference) so the two lists cannot drift. Sdist per Q7 (recommended branch):
   ```toml
   [tool.hatch.build.targets.sdist]
   exclude = ["/CLAUDE.md", "/PLAN.md", "/docs/*-phases.md", "/.github", "/.env.example", "/site", "/build", "/data"]
   ```
   `tests/`, `docs/`, `contract/`, `scripts/`, `LICENSE`, `README.md`, `CHANGELOG.md`, `CITATION.cff` stay in the sdist (downstream packagers run the tests; `tests/test_check_contract.py` needs `contract/` and `scripts/`).
3. `cli.py`: `from progeny_selector._version import __version__`; in `build_parser()` after `ArgumentParser(...)`: `parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")`. The `select` `--top` help becomes exactly: `"select every passing individual ranked N or better within its family, so up to N per family (ranks are dense, so ties can return more); with --overall, the N best overall"`. Module docstring interface block: add `progeny-selector --version`, and `[--coding auto|nucleotide|abh]` on validate/rank and `[--project backcross|self]` on select, which it omits today.
4. `tests/test_version.py`:
   - `test_pyproject_reads_the_version_from_the_package`: `tomllib.load` of `pyproject.toml`; `"version" not in data["project"]`, `data["project"]["dynamic"] == ["version"]`, `data["tool"]["hatch"]["version"]["path"] == "src/progeny_selector/_version.py"`.
   - `test_version_is_a_plain_semver`: `re.fullmatch(r"\d+\.\d+\.\d+", __version__)`.
   - `test_installed_metadata_matches_the_package`: `importlib.metadata.version("progeny-selector") == progeny_selector.__version__`. A `PackageNotFoundError` is a failure, not a skip: the venv recipe in CLAUDE.md is an editable install, and after any bump `pip install -e .` must be re-run (record this in CONTRIBUTING, phase 7). The fresh-venv comparison in phase 4 is the release-grade form of this test.
   - `test_version_flag(capsys)`: `with pytest.raises(SystemExit) as e: main(["--version"])`; `e.value.code == 0`; `capsys.readouterr().out == f"progeny-selector {__version__}\n"`.
   - `test_no_other_version_literal`: `grep`-style scan of `src/progeny_selector/**/*.py` for `__version__ = "`; exactly one file, `_version.py`. Discriminating: restoring the literal in `__init__.py` fails it.
5. `tests/test_cli_select_top.py` (uses the `fixture_dir` fixture from `tests/conftest.py`):
   - `test_top_is_per_family_unless_overall(tmp_path, fixture_dir)`: `main(["rank", ...fixture args..., "--out", results])` returns 0; `main(["select","--results",results,"--top","3","--out",sel])` writes 6 data rows; with `--overall` 3 rows; `--top 1` writes 2 rows (one per family, F1 and F2). Discriminating: 6 vs 3 separates per-family from overall; 2 for `--top 1` separates "N per family" from "N total".
   - `test_top_help_names_per_family`: from `build_parser()` find the `select` subparser (`parser._subparsers._group_actions[0].choices["select"]`) and its `--top` action; `"up to N per family" in action.help` and `"--overall" in action.help`.

Acceptance: `progeny-selector --version` prints `progeny-selector <_version.py value>`; `pyproject.toml` carries no version literal; `python -m build` succeeds (run locally by the main session, not a gate); `select --top 3` on the fixture still writes 6 rows and the help text says why.

## Phase 2: examples as package data and `progeny-selector example`

Model opus (package data crosses the Shinylive staging boundary and the generator is the fixture authority); gates only (no core/io/model diff). Depends on 1. Files: `scripts/make_fixture.py`, `src/progeny_selector/examples/__init__.py` (new), `src/progeny_selector/examples/synthetic_bc2f1/{genotypes.vcf,samples.csv,markers.csv,criteria.yaml}` and `src/progeny_selector/examples/synthetic_bc3f1/{genotypes.vcf,samples.csv}` (generated, never hand-edited), `src/progeny_selector/cli.py`, `.gitignore`, `.github/workflows/ci.yml` (one line), `tests/test_examples.py` (new).

1. Generator. `scripts/make_fixture.py` gains `EXAMPLES = Path(__file__).resolve().parents[1] / "src" / "progeny_selector" / "examples"` and a final step in `main()`, `copy_examples()`, that `shutil.copyfile`s the six files above from `OUT` and `OUT_BC3F1` into `EXAMPLES/<name>/` (creating directories), and prints `wrote package examples to <EXAMPLES>: 6 files`. Not copied: `expected_results.csv`, `README.md`, anything under `tests/fixtures/brapi`. The generator remains the only author of both copies; byte copy means Windows line endings cannot differ between them. Module docstring gains one line saying so.
2. `.gitignore`: add `!src/progeny_selector/examples/**` after `!contract/**`. Reason, as a comment: `*.vcf` is ignored globally, and hatchling excludes VCS-ignored files from the sdist and wheel by default, so without this line the package data would be missing from both git and the wheel (proved by the phase 4 package job's `unzip -l`).
3. `src/progeny_selector/examples/__init__.py`:
   ```python
   """Shipped example inputs: the synthetic BC2F1 and BC3F1 fixtures as package data.

   Responsibility: locate the example files inside the installed package (and inside the
   Shinylive staged copy, which is the same directory tree) and write them to a directory.
   The files are byte copies of tests/fixtures/ written by scripts/make_fixture.py, the only
   author of both; tests/test_examples.py proves the copies identical. Synthetic data: the
   rankings are test cases, not breeding recommendations.

   Interface:
       EXAMPLES: dict[str, tuple[str, ...]]          # name -> files shipped under examples/<name>/
       SHARED_FROM_BC2F1: tuple[str, ...]            # ("markers.csv", "criteria.yaml"): bc3f1 reuses bc2f1's
       EXAMPLE_FILES: tuple[str, ...]                # ("genotypes.vcf", "samples.csv", "markers.csv", "criteria.yaml")
       example_files(name) -> dict[str, Traversable] # the four inputs of an example, shared ones resolved from bc2f1
       example_paths(name) -> AbstractContextManager[dict[str, Path]]   # real paths (importlib.resources.as_file)
       write_example(name, out_dir: Path, *, force: bool = False) -> list[Path]
                                                     # writes out_dir/<name>/<four files>; FileExistsError naming the
                                                     # first existing target unless force
   """
   ```
   `EXAMPLES = {"synthetic_bc2f1": EXAMPLE_FILES, "synthetic_bc3f1": ("genotypes.vcf", "samples.csv")}`. Resources via `importlib.resources.files("progeny_selector.examples")`; an unknown name raises `KeyError` naming the valid names.
4. CLI subcommand `example`: `progeny-selector example [--out DIR] [--name {synthetic_bc2f1,synthetic_bc3f1,all}] [--force]`; `--out` default `progeny-selector-example`, `--name` default `all`. Help: `"write the shipped synthetic example inputs (BC2F1 and its BC3F1 next generation) to a directory"`. `cmd_example(args) -> int`: for each name, `write_example`; on `FileExistsError` print `error: {path} exists; pass --force to overwrite` to stderr and return 1 (nothing partially written before the check: `write_example` checks all four targets first). Stdout, exactly:
   ```
   wrote {out}/synthetic_bc2f1: genotypes.vcf, samples.csv, markers.csv, criteria.yaml
   wrote {out}/synthetic_bc3f1: genotypes.vcf, samples.csv, markers.csv, criteria.yaml
   synthetic data: the rankings are test cases, not breeding recommendations
   next: progeny-selector rank --genotypes {out}/synthetic_bc2f1/genotypes.vcf --samples {out}/synthetic_bc2f1/samples.csv --markers {out}/synthetic_bc2f1/markers.csv --criteria {out}/synthetic_bc2f1/criteria.yaml --out results.csv
   ```
   `{out}` is the path as given. Module docstring interface block gains the line.
5. `ci.yml` `check` job step "fixture is reproducible": `git diff --exit-code -- tests/fixtures src/progeny_selector/examples`. (Phase 4 edits the rest of the file later.) CLAUDE.md's gate sentence (line 24) gets the same pathspec (main session).
6. `tests/test_examples.py`:
   - `test_package_copy_is_byte_identical_to_the_fixture` parametrized over the six `(name, file)` pairs plus the two shared bc3f1 files resolved to bc2f1: `example_files(name)[file].read_bytes() == (FIXTURES / src_name / file).read_bytes()`. Discriminating: any byte edit to either copy fails exactly one case.
   - `test_example_command_writes_self_contained_directories(tmp_path)`: `main(["example", "--out", str(tmp_path / "ex")]) == 0`; `sorted(os.listdir(ex/"synthetic_bc3f1")) == ["criteria.yaml", "genotypes.vcf", "markers.csv", "samples.csv"]`; then `main(["rank", "--genotypes", ex/synthetic_bc3f1/genotypes.vcf, "--samples", ..., "--markers", ex/synthetic_bc3f1/markers.csv, "--criteria", ex/synthetic_bc3f1/criteria.yaml, "--out", ...]) == 0` and stdout contains `40 individuals` (4 selected parents x 10 progeny; proves the shared copies are the right files, since bc3f1's genotypes need bc2f1's markers.csv for the cM map).
   - `test_example_refuses_to_overwrite_unless_forced(tmp_path, capsys)`: second `main(["example","--out",...])` returns 1, stderr contains `exists; pass --force`; a third with `--force` returns 0. Discriminating: removing the check makes the second call return 0.
   - `test_unknown_example_name_is_a_usage_error`: `main(["example","--name","nope"])` raises `SystemExit` with code 2 (argparse choices).
   - `test_example_output_is_ignored_by_gitignore_but_package_data_is_not`: no git call (doers never run git); instead parse `.gitignore` text and assert the line `!src/progeny_selector/examples/**` is present after the `*.vcf` line. Weak but cheap; the wheel proof is phase 4.

Not verifiable here: whether `shinylive export` stores `.vcf` as text or base64 in `app.json`; either lands on Pyodide's file system, and phase 3's export test settles it. Acceptance: regenerating with `python scripts/make_fixture.py` leaves `git diff --exit-code -- tests/fixtures src/progeny_selector/examples` clean; `progeny-selector example` from a temp directory produces two runnable example directories.

## Phase 3: `progeny-selector app` and the Load-example control

Model sonnet; gates only, plus `pytest -m e2e` and the export smoke run by the doer (`python scripts/build_shinylive.py` then `PS_SITE_DIR=site pytest -m e2e tests/e2e/test_shinylive_export.py`). Depends on 2. Files: `src/progeny_selector/cli.py`, `src/progeny_selector/app/screens/load.py`, `tests/test_cli_app.py` (new), `tests/e2e/test_load_screen.py`, `tests/e2e/test_shinylive_export.py`, `tests/e2e/helpers.py`, `docs/keyboard-walkthrough.md` (Load section).

1. CLI subcommand `app`: `progeny-selector app [--host HOST] [--port PORT] [--no-browser]`; `--host` default `127.0.0.1`, `--port` type int default `8000`, `--no-browser` store_true ("do not open a browser tab"). Help: `"open the Shiny interface in your browser (Ctrl+C stops the server)"`. `cmd_app(args) -> int`:
   ```python
   try:
       import pandas  # noqa: F401  (the screens import it; fail here, not inside uvicorn)
       from shiny import run_app
   except ModuleNotFoundError as exc:
       print(APP_MISSING_MESSAGE.format(missing=exc.name), file=sys.stderr)
       return 2
   print(f"progeny-selector {__version__}: serving the interface on http://{args.host}:{args.port}/ (Ctrl+C stops it)")
   run_app("progeny_selector.app.app:app", host=args.host, port=args.port, launch_browser=not args.no_browser)
   return 0
   ```
   `shiny.run_app` signature verified in the installed shiny 1.7 (`.venv/Lib/site-packages/shiny/_main/_run.py:213-230`): `run_app(app, host="127.0.0.1", port=8000, *, ..., launch_browser: bool = False, dev_mode: bool = True, ...)`. `APP_MISSING_MESSAGE` under Q1 = core dependencies (recommended):
   ```
   error: the Shiny interface needs the package "{missing}", which is not installed.
   Reinstall the tool with:  python -m pip install --upgrade --force-reinstall progeny-selector
   With pipx:                pipx install --force progeny-selector
   ```
   Under Q1 = keep the extra: line 2 `Install it with:  python -m pip install "progeny-selector[app]"`, line 3 `With pipx:        pipx install "progeny-selector[app]"`.
2. Load screen (`app/screens/load.py`). Under the Run button in the Input files card: `ui.input_action_button("example", "Load example (synthetic BC2F1)", class_="btn-outline-secondary")` and `ui.p("Synthetic data shipped with the package; its rankings are test cases, not breeding recommendations.", class_="form-text")`. Element id in the page: `load-example`. Server: factor the AppState writes of `_run` (lines 198-210) into `_publish(dataset, criteria, *, first_line: str | None) -> None`, which does exactly what those lines do and sets `message` to `first_line` (when given) followed by the existing summary line and warnings. `_run` calls `_publish(dataset, criteria, first_line=None)`. New:
   ```python
   @reactive.effect
   @reactive.event(input.example)
   def _example() -> None:
       try:
           with example_paths("synthetic_bc2f1") as p:
               dataset = load_dataset(p["genotypes.vcf"], p["samples.csv"], p["markers.csv"], crop=DEFAULT_CROP_ID)
               criteria = read_criteria(p["criteria.yaml"])
           result = run_analysis(dataset, criteria)   # inside _publish or here, matching _run
           _publish(dataset, criteria, first_line=EXAMPLE_FIRST_LINE)
       except (DataContractError, CriteriaError) as exc: message.set(f"error: {exc}")
       except Exception as exc: (as _run)
   ```
   `EXAMPLE_FIRST_LINE = "loaded the shipped example: synthetic BC2F1 (test data; its rankings are not breeding recommendations)"`. The summary line stays verbatim `500 markers, 40 progeny; 11 pass hard filters`, which `tests/e2e/helpers.py::LOAD_STATUS` and the export test match with `to_contain_text`. The token-profile and crop selects are not read for the example (contract default, soybean). Module docstring: add the control and that the example is read from `progeny_selector.examples`, which under Shinylive is the staged copy of the package. `app` computes nothing new.
3. `tests/e2e/helpers.py`: `EXAMPLE_FIRST_LINE` re-exported (import from `progeny_selector.app.screens.load`).
4. Tests:
   - `tests/test_cli_app.py::test_app_without_shiny_prints_the_install_line(monkeypatch, capsys)`: `monkeypatch.setitem(sys.modules, "shiny", None)` (a `None` entry makes `import shiny` raise `ModuleNotFoundError`); `main(["app"]) == 2`; stderr contains `is not installed` and the exact second line of `APP_MISSING_MESSAGE`. Discriminating: removing the handler propagates `ModuleNotFoundError` out of `main`.
   - `test_app_passes_host_port_and_browser_flag(monkeypatch)`: `import shiny`; `monkeypatch.setattr(shiny, "run_app", recorder)`; `main(["app", "--host", "0.0.0.0", "--port", "8123", "--no-browser"]) == 0` and `recorder.calls == [(("progeny_selector.app.app:app",), {"host": "0.0.0.0", "port": 8123, "launch_browser": False})]`; `main(["app"])` records `host="127.0.0.1", port=8000, launch_browser=True`. Discriminating: any default change fails.
   - `tests/e2e/test_load_screen.py::test_load_example_button_runs_the_shipped_example(page, app)`: `page.goto(app.url)`; no uploads; `controller.InputActionButton(page, "load-example").click()`; `expect(page.locator("#load-status")).to_contain_text(EXAMPLE_FIRST_LINE)` and `.to_contain_text(LOAD_STATUS, timeout=60_000)`; click the Rank tab (`a.nav-link[data-value='rank']`) and `expect(page.locator("#rank-table.html-fill-item div.shiny-data-grid table tbody tr")).to_have_count(11)`. Discriminating: without the control the status reads "Genotypes, samples.csv and criteria.yaml are required." and Rank is empty.
   - `tests/e2e/test_shinylive_export.py::test_shinylive_load_example_stays_on_origin(site_server, page)`: same preamble as the existing test to `#load-run` attached, then `frame.locator("#load-example").click()`, `expect(frame.locator("#load-status")).to_contain_text(LOAD_STATUS, timeout=STATUS_TIMEOUT_MS)`, every recorded request starts with `ORIGIN`. Proves the package data reached Pyodide through the staged copy (docs/adr/0009) without any off-origin fetch. Not verifiable here until the doer runs it.
   - `tests/e2e/test_focus_order.py` and `test_a11y.py`: run unchanged; if the focus-order test pins the Load tab order, the doer adds `load-example` after `load-run` and records it in `docs/keyboard-walkthrough.md` Load section ("Load example (synthetic BC2F1)" after "Load and analyse").

Acceptance: from a directory that is not the checkout, `progeny-selector app --no-browser --port 8123` serves the page (curl 200, title `progeny-selector`), which phase 4 automates; the Load-example button reaches `11 pass hard filters` under `shiny run` and under the Shinylive export with no off-origin request.

## Phase 4: CI matrix, OS jobs, `package` job, release workflow

Model sonnet; gates only. YAML is outside ruff and mypy; not verifiable here, the first push settles it (push before calling the phase finished, PLAN.md Handoff rule). Depends on 3. Files: `.github/workflows/ci.yml`, `.github/workflows/release.yml` (new), `scripts/changelog_section.py` (new), `tests/test_changelog_section.py` (new).

1. `ci.yml` `check` job: `python-version: ["3.11", "3.12", "3.13"]`; mypy `if: matrix.python-version == '3.13'` (CLAUDE.md line 61: mypy on the newest supported Python; update the comment). `pip install -e ".[dev]"` unchanged (under Q1 = core, dev no longer lists shiny/pandas and they come from the core dependencies).
2. New job `unit-os` (needs nothing, runs in parallel with `check`):
   ```yaml
   unit-os:
     strategy:
       fail-fast: false
       matrix:
         os: [windows-latest, macos-latest]
     runs-on: ${{ matrix.os }}
     steps:
       - uses: actions/checkout@v4
       - uses: actions/setup-python@v5
         with: { python-version: "3.12", cache: pip }
       - run: python -m pip install --upgrade pip
       - run: pip install -e ".[dev]"
       - run: pytest
   ```
   No coverage, no fixture regeneration (Q11: the generator's CSV writers open with `newline=""` and `.gitattributes` normalises, so it probably passes on Windows, but the gate is proved on ubuntu and the maintainer's laptop; add it only if Q11 says so). Not verifiable here: macOS has never run this suite.
3. New job `package` (needs `check`), ubuntu, Python 3.12:
   ```yaml
   - run: python -m pip install --upgrade pip build twine
   - run: python -m build
   - run: twine check --strict dist/*
   - name: wheel carries the examples and the sdist excludes the working documents
     run: |
       unzip -l dist/*.whl | grep -E 'progeny_selector/examples/synthetic_bc2f1/genotypes\.vcf$'
       unzip -l dist/*.whl | grep -E 'progeny_selector/examples/synthetic_bc3f1/samples\.csv$'
       ! tar tzf dist/*.tar.gz | grep -E '/(CLAUDE|PLAN)\.md$'
       tar tzf dist/*.tar.gz | grep -E '/tests/fixtures/synthetic_bc2f1/genotypes\.vcf$'
       tar tzf dist/*.tar.gz | grep -E '/LICENSE$'
   - name: fresh venv from the wheel
     run: |
       python -m venv "$RUNNER_TEMP/fresh"
       "$RUNNER_TEMP/fresh/bin/pip" install --upgrade pip
       "$RUNNER_TEMP/fresh/bin/pip" install dist/*.whl        # Q1 extra branch: "$(ls dist/*.whl)[app]"
   - name: smoke from outside the checkout
     working-directory: ${{ runner.temp }}
     run: |
       set -euo pipefail
       PS=fresh/bin/progeny-selector; PY=fresh/bin/python
       $PS --version | tee version.txt
       v=$(cut -d' ' -f2 version.txt)
       ls "$GITHUB_WORKSPACE"/dist/progeny_selector-"$v"-py3-none-any.whl
       $PY -c "import progeny_selector, importlib.metadata as m; assert m.version('progeny-selector') == progeny_selector.__version__"
       $PS example --out ex
       $PS validate --genotypes ex/synthetic_bc2f1/genotypes.vcf --samples ex/synthetic_bc2f1/samples.csv --markers ex/synthetic_bc2f1/markers.csv --criteria ex/synthetic_bc2f1/criteria.yaml
       $PS rank --genotypes ex/synthetic_bc2f1/genotypes.vcf --samples ex/synthetic_bc2f1/samples.csv --markers ex/synthetic_bc2f1/markers.csv --criteria ex/synthetic_bc2f1/criteria.yaml --out results.csv | grep -F "40 individuals, 11 pass hard filters"
       $PS select --results results.csv --top 3 --out selected.csv --next-manifest next_samples.csv --next-generation BC3F1 --samples ex/synthetic_bc2f1/samples.csv
       test "$(tail -n +2 selected.csv | wc -l)" -eq 6
       $PS rank --genotypes ex/synthetic_bc3f1/genotypes.vcf --samples ex/synthetic_bc3f1/samples.csv --markers ex/synthetic_bc3f1/markers.csv --criteria ex/synthetic_bc3f1/criteria.yaml --out results_bc3f1.csv | grep -F "40 individuals"
       $PY -c "import progeny_selector.app.app as a; assert a.app is not None"
       $PS app --help
       $PS app --no-browser --port 8765 & pid=$!
       for i in $(seq 1 60); do curl -fsS http://127.0.0.1:8765/ -o page.html && break; sleep 1; done
       grep -F "<title>progeny-selector</title>" page.html
       kill $pid
   - name: sdist installs too
     run: |
       python -m venv "$RUNNER_TEMP/fromsdist"
       "$RUNNER_TEMP/fromsdist/bin/pip" install dist/*.tar.gz
       "$RUNNER_TEMP/fromsdist/bin/progeny-selector" --version
   - uses: actions/upload-artifact@v4
     with: { name: dist, path: dist/ }
   ```
   `working-directory: ${{ runner.temp }}` is load-bearing: run from the checkout, `src/` or `tests/fixtures` could mask a broken wheel. The wheel filename check ties the tag-time version to the built artefact (`hatchling` names the wheel from the dynamic version). This job is the milestone's usability acceptance.
4. `scripts/changelog_section.py`: `python scripts/changelog_section.py VERSION [--changelog CHANGELOG.md] [--out FILE] [--check-citation CITATION.cff]`. Finds `^## \[VERSION\] - (?P<date>.+)$`, requires `date` to match `\d{4}-\d{2}-\d{2}`, prints (or writes) the body up to the next `^## \[` or the `^\[` link-reference block, whichever comes first, trimmed. Exit 1 with `CHANGELOG.md: no section "## [VERSION] - YYYY-MM-DD"` or `CHANGELOG.md: section [VERSION] has date "unreleased", not YYYY-MM-DD`. With `--check-citation`, also requires CITATION.cff's `version:` to equal VERSION and `date-released:` to equal the section date (plain line parsing, no YAML dependency beyond pyyaml which is a core dependency). Pure function `section(text: str, version: str) -> tuple[str, str]` (date, body) for the tests.
5. `.github/workflows/release.yml`, full text:
   ```yaml
   name: release

   on:
     push:
       tags: ["v*"]

   permissions:
     contents: read

   jobs:
     build:
       runs-on: ubuntu-latest
       outputs:
         version: ${{ steps.version.outputs.v }}
       steps:
         - uses: actions/checkout@v4
         - uses: actions/setup-python@v5
           with:
             python-version: "3.12"
         - run: python -m pip install --upgrade pip build twine
         - run: python -m build
         - run: twine check --strict dist/*
         - name: the tag names the package version
           id: version
           run: |
             python -m venv "$RUNNER_TEMP/fresh"
             "$RUNNER_TEMP/fresh/bin/pip" install dist/*.whl
             v=$("$RUNNER_TEMP/fresh/bin/progeny-selector" --version | cut -d' ' -f2)
             test "v$v" = "$GITHUB_REF_NAME" || { echo "tag $GITHUB_REF_NAME does not match package version $v"; exit 1; }
             echo "v=$v" >> "$GITHUB_OUTPUT"
         - name: release notes from CHANGELOG.md
           run: python scripts/changelog_section.py "${{ steps.version.outputs.v }}" --out notes.md --check-citation CITATION.cff
         - uses: actions/upload-artifact@v4
           with:
             name: release
             path: |
               dist/
               notes.md

     publish-pypi:
       needs: build
       runs-on: ubuntu-latest
       environment:
         name: pypi
         url: https://pypi.org/project/progeny-selector/${{ needs.build.outputs.version }}/
       permissions:
         id-token: write
       steps:
         - uses: actions/download-artifact@v4
           with:
             name: release
         - uses: pypa/gh-action-pypi-publish@release/v1
           with:
             packages-dir: dist/

     github-release:
       needs: [build, publish-pypi]
       runs-on: ubuntu-latest
       permissions:
         contents: write
       steps:
         - uses: actions/download-artifact@v4
           with:
             name: release
         - run: gh release create "$GITHUB_REF_NAME" dist/* --repo "$GITHUB_REPOSITORY" --title "progeny-selector ${{ needs.build.outputs.version }}" --notes-file notes.md
           env:
             GH_TOKEN: ${{ github.token }}
   ```
   No API token anywhere: PyPI trusted publishing (OIDC) through `id-token: write` on the `pypi` environment; the GitHub release uses the workflow's own token. `--check-citation` fails the build until phase 8's CITATION.cff exists, so this file is committed in phase 4 but no tag is pushed before phase 8 (Maintainer steps).
6. `tests/test_changelog_section.py`: a three-section text (`[Unreleased]`, `[1.1.0] - 2026-10-02`, `[1.0.0] - 2026-09-30`, then the link block); `section(text, "1.0.0")` returns `("2026-09-30", body)` with the 1.0.0 lines only; `section(text, "1.1.0")` excludes every 1.0.0 line and the link block; `[0.9.0] - unreleased` raises `SystemExit(1)` with the "not YYYY-MM-DD" message; a missing version raises `SystemExit(1)`. Discriminating: the body of 1.1.0 must not contain a 1.0.0 line.

Acceptance: CI green on `check (3.11, 3.12, 3.13)`, `unit-os (windows-latest, macos-latest)`, `package`, plus the existing `e2e`, `export`, `r-reader`, `deploy-pages`; the `package` smoke passes from `runner.temp`.

## Phase 5: README for breeders

Model sonnet; gates only. Depends on 3 (the commands must exist). Files: `README.md`, `tests/test_readme_links.py` (new). Every link absolute with base `https://github.com/piercetaylor/progeny-selector/blob/main/` (directories with `/tree/main/`), so the PyPI page renders them. Headings and content, in order:

1. `# progeny-selector`: two sentences from the current README line 3, then: "Status: version {{VERSION}} (copy from `_version.py`); see [CHANGELOG.md](abs) and [Known limitations](#known-limitations)."
2. `## Install` (Python 3.11 or newer):
   ```sh
   python -m pip install progeny-selector
   ```
   ```sh
   pipx install progeny-selector
   ```
   Windows note: `py -3 -m pip install progeny-selector`; if `progeny-selector` is not found afterwards, use `py -3 -m progeny_selector ...` or install with pipx (`py -3 -m pip install --user pipx && py -3 -m pipx ensurepath`). Under Q1 = keep the extra, each line reads `"progeny-selector[app]"` and a sentence says the CLI-only install is `progeny-selector`.
3. `## Open the interface`: `progeny-selector app` (opens a browser tab at http://127.0.0.1:8000/; `--port`, `--host`, `--no-browser`; Ctrl+C stops it). "Or use the [browser site](https://piercetaylor.github.io/progeny-selector/): files are processed in the tab and never uploaded. Both offer **Load example (synthetic BC2F1)** on the Load screen."
4. `## Run the example` (the fenced `sh` block is executed by phase 6's test, so paths are exactly these):
   ```sh
   progeny-selector example --out example
   progeny-selector validate --genotypes example/synthetic_bc2f1/genotypes.vcf --samples example/synthetic_bc2f1/samples.csv --markers example/synthetic_bc2f1/markers.csv --criteria example/synthetic_bc2f1/criteria.yaml
   progeny-selector rank --genotypes example/synthetic_bc2f1/genotypes.vcf --samples example/synthetic_bc2f1/samples.csv --markers example/synthetic_bc2f1/markers.csv --criteria example/synthetic_bc2f1/criteria.yaml --out results.csv
   progeny-selector select --results results.csv --top 3 --out selected.csv --next-manifest next_samples.csv --next-generation BC3F1 --samples example/synthetic_bc2f1/samples.csv
   ```
   Then one paragraph: 500 markers, 2 parents, 40 BC2F1 progeny in 2 families, 11 pass the hard filters; the example is synthetic and its rankings are test cases.
5. `## What the outputs mean`: `results.csv` (one row per individual: `rank_overall`, `rank_in_family`, `passes_filters`, `exclusion_reason`, `composite_score`, `rpp_total`/`rpp_carrier`/`rpp_noncarrier`, `drag_total_est`/`drag_total_max` with `drag_unit`, `missing_rate`, `qc_flags`, one status column per target and avoid locus, one `rpp_<chromosome>` per chromosome; `results_schema` names the column layout, `background_model`, `background_unit`, `background_max_marker_coverage`, `assembly` and `crop` say how it was computed); `selected.csv` (the chosen individuals with `notes`); `next_samples.csv` (a samples.csv for the next genotyping round). Link to the tutorial, the glossary and `docs/data-formats.md`.
6. `### What `--top N` selects`: exactly the phase 1 help text in prose, with the example: `--top 3` on the example writes 6 rows (3 in family F1, 3 in F2); `--top 3 --overall` writes 3; ranks are dense, so ties can return more than N; the Selection screen's "add top N per family" does the same.
7. `## Your own data`: formats (VCF, VCF.gz, HapMap, wide CSV; BrAPI v2.1 on the command line), `samples.csv` (one recurrent and one donor parent), `criteria.yaml`, `markers.csv` (optional map), crops (the 13 schemes), assembly, `--profile` token profiles; links to `contract/data-contract.md`, `docs/data-formats.md`, `docs/soybean-inputs.md` (not linked today: audit line 25).
8. `## Known limitations`: eleven bullets, verbatim in CHANGELOG (phase 8): (1) `possible_duplicate` is raw pairwise IBS >= 0.995 over informative markers, blind to family; advisory only (docs/adr/0017). (2) The two-generation round trip is verified on synthetic data only; no real linked BC(n) -> BC(n+1) dataset has been run. (3) The browser site has no BrAPI source; BrAPI is command-line only until a server sends CORS headers for the site's origin (docs/adr/0024). (4) One donor per analysis; an intercross or pyramiding population is not ranked (contract line 65, docs/adr/0002). (5) Chromosome-length tables exist for soybean only (Wm82.a2, a4); under any other crop a chromosome ends at its last marker (docs/adr/0020). (6) The bp coverage cap default of 2 Mb is a soybean translation of 10 cM; set `background.max_marker_coverage` for other crops (docs/adr/0033). (7) A VCF 4.4 leading phase indicator in GT is rejected as `genotypes.invalid_gt` (contract 1.10.0). (8) Accessibility: multi-row selection in the Rank grid needs a mouse; three axe rules are recorded exceptions; no screen-reader, Firefox, Safari or 320 px reflow testing (docs/accessibility.md). (9) Purdy labels: single crosses only, and no heterozygosity checks because the label carries no filial generation (docs/adr/0034). (10) The long-format KASP header names are unconfirmed against a real LGC export (docs/adr/0019). (11) CI runs the browser tests on Linux only; Windows and macOS run the unit tests.
9. `## Verification`: the gate commands; browser and export checks per PLAN.md; `docs/limits.md`; `docs/accessibility.md`.
10. `## Related tool`: one sentence on backcross with its repository link (phase 8 finalises the wording per `{{SAME_DAY}}`).
11. `## Citing`: "See [CITATION.cff](abs) (GitHub shows a Cite this repository button)." Phase 8 completes it.
12. `## Licence`: MIT.

`tests/test_readme_links.py::test_readme_links_are_absolute_and_resolve`: every `[text](target)` in README.md has `target` starting with `https://` or `#`; every `https://github.com/piercetaylor/progeny-selector/(blob|tree)/main/<path>` has `<path>` existing relative to the repo root (`#fragment` stripped); every `#anchor` matches a heading slug of README.md (lowercase, spaces to `-`, characters outside `[a-z0-9-]` removed). Discriminating: today's nine relative links fail it, and a typo in any path fails the existence check. Acceptance: the test is green and `twine check` (phase 4) reports no rendering warning.

## Phase 6: tutorial, glossary, docs index, commands test

Model sonnet; gates only. Depends on 3. Files: `docs/tutorial.md` (new), `docs/glossary.md` (new), `docs/README.md` (new), `tests/test_docs_commands.py` (new). Plain language; each definition that comes from literature names the source already cited in PLAN.md or an ADR, with the ADR number.

1. `docs/tutorial.md`, "One run on the shipped example". Sections: `## Get the example` (`progeny-selector example --out example`, or the Load-example button); `## Command line` (the four README commands in one `sh` block, then `select` again with `--per-selected 10`, then `rank` on `example/synthetic_bc3f1` as the next generation, all in `sh` blocks that the test runs); `## What the example contains` (2 parents, 40 BC2F1 progeny in families F1 and F2, 500 markers on Gm01-Gm20, 475 informative, target `syn_Gm06_13` donor allele required, avoid `syn_Gm13_10` must be recurrent, flank windows 6 cM, count-model background; planted individuals from `tests/fixtures/synthetic_bc2f1/README.md`: BC2F1-F1-001 ranks first; F1-002 fails the target; F2-001 fails the avoid locus; F2-002 shows homozygous-donor calls and is flagged `possible_self_or_outcross`; F1-003 has 30 % missing calls and is `high_missing`); `## The screens` with one subsection each for Load, Validate and QC, Navigate, Rank (what each column and chip means, the hide-excluded switch, ctrl-click and shift-click to select rows), Compare (chips, per-chromosome RPP, drag bounds, the 20-row strip in Okabe-Ito colours), Selection list (notes, "add top N per family" and why 3 gives 6, the projected next generation: backcross halves donor content, expected RPP (1 + RPP)/2), Export; `## Reading results.csv` as a table of the fixed columns (the 38 names of `tests/test_results_schema.py:35-74`, plus `tool_version` if 0b ran) with one-line meanings and the dynamic column families `target_<id>_status`, `avoid_<id>_status`, `rpp_<chrom>`, `drag_<id>_*`, `recombinant_<id>_*` as named in `docs/data-formats.md` (the doer copies the exact prefixes from that file, never invents them); `## Next generation` (edit `next_samples.csv` ids after planting; `generation` BC3F1; the BC3F1 example is that file's shape). Every screen control name is copied from the screen modules (`app/screens/*.py`), never guessed.
2. `docs/glossary.md`, alphabetical, each entry two to five sentences: avoid locus (PLAN.md 5; Flapjack's QTL Source RP, docs/adr/0006); background selection and RPP, count vs weighted (Flapjack MABC, Hospital and Charcosset 1997, docs/adr/0006, 0033; the weight formula and the chromosome-end rule); carrier and non-carrier chromosome (Hospital and Charcosset 1997; Frisch, Bohn and Melchinger 1999; PLAN.md 3); composite score and hard filters (docs/adr/0007); crop scheme and assembly (docs/adr/0015, 0020); expected RPP by generation, 1 - (1/2)^(n+1) (PLAN.md 3, the Iowa State chapter cited there); `family_donor_outlier` (docs/adr/0012); flanking window and recombinant flags (Flapjack "first recombination on each side", docs/adr/0006); foreground selection and `required_state` (PLAN.md 2, docs/adr/0011 for `rule: run`); generation labels BCnFm and Purdy (docs/adr/0034); IBS and `possible_rp_sample`/`possible_donor_sample` (SNPRelate, PLAN.md 6); linkage drag, min and max bounds, estimate (PLAN.md 4); `possible_duplicate` (docs/adr/0017, with the family-blind caveat and the BC3F1 fixture note); `possible_self_or_outcross`, `possible_outcross`, `het_rate_deviates`, `high_missing`, `generation_unparsed` (PLAN.md 8); `results_schema` (docs/adr/0016); the six genotype states A, H, B, X, N, U with the backcross labels rp_hom, donor_hom, het, nonparental, missing, uninformative and the ABHgenotypeR coding (PLAN.md 1, docs/adr/0006); staged vs weighted ranking (docs/adr/0007 amendment; Frisch, Bohn and Melchinger 1999 four-stage order); token profile (docs/adr/0014); uninformative marker and its reasons (PLAN.md 1).
3. `docs/README.md`: an index with one line per document (tutorial, glossary, data-formats, soybean-inputs, keyboard-walkthrough, accessibility, limits, reference-repos, adr/, the phase specs as history), absolute links so it renders from the Documentation URL.
4. `tests/test_docs_commands.py::test_fenced_sh_commands_run(doc, tmp_path)` parametrized over `README.md` and `docs/tutorial.md`: extract every ```` ```sh ```` block; for each line starting with `progeny-selector `, run `subprocess.run([sys.executable, "-m", "progeny_selector", *shlex.split(rest)], cwd=tmp_path, capture_output=True, text=True)` in document order in one `tmp_path` per document, and assert `returncode == 0` with the command and stderr in the message; skip lines starting with `progeny-selector app` and any line containing `--brapi-url`; lines not starting with `progeny-selector ` (pip, pipx, py -3) are ignored. Discriminating: a wrong flag, a wrong path or a command placed before `example --out` returns 2 or 1.

Acceptance: both docs exist, the commands test is green for both files, `docs/README.md` links resolve (extend `tests/test_readme_links.py` to include `docs/README.md` in its parametrization, one-line edit permitted here despite the phase-5 file list, because phases 5 and 6 may run concurrently and the second to land adds the parameter).

## Phase 7: community files and CONTRIBUTING

Model sonnet; gates only. Depends on 1. Files: `.github/ISSUE_TEMPLATE/bug_report.yml`, `.github/ISSUE_TEMPLATE/dataset_or_format_question.yml`, `.github/ISSUE_TEMPLATE/config.yml`, `.github/PULL_REQUEST_TEMPLATE.md`, `.github/dependabot.yml`, `SECURITY.md`, `CODE_OF_CONDUCT.md`, `CONTRIBUTING.md`, `tests/test_github_templates.py` (new).

1. `bug_report.yml` (GitHub issue form): `name: Bug report`, `description: Something is wrong or fails`, `labels: [bug]`; body: markdown "Never attach real genotype data. Reproduce with the shipped example (`progeny-selector example --out example`) or a synthetic file, and paste headers only."; inputs, all `required: true` unless noted: `version` ("Output of `progeny-selector --version`, or the browser site's date"), `how_installed` dropdown (`pip`, `pipx`, `browser site`, `checkout`), `os` dropdown (Windows, macOS, Linux, other), `python` ("Output of `python --version`", not required for the browser site), `command` textarea ("The exact command or the screen and control"), `error` textarea ("The exact error text or the wrong value"), `expected` textarea.
2. `dataset_or_format_question.yml`: `name: Dataset or format question`, `description: My file does not load or my chromosome names are not recognised`, `labels: [question, data-format]`; inputs: `crop` dropdown (the 13 scheme names plus other), `platform` (KASP, array, GBS or sequencing VCF, other), `format` dropdown (VCF, VCF.gz, HapMap, wide CSV, BrAPI), `header` textarea ("The first five lines of the file with sample names replaced; never real calls beyond a header line"), `chromosomes` ("How chromosomes are named in the file, for example `Gm01`, `chr1`, `1`"), `error` textarea (the exact message, not required).
3. `config.yml`: `blank_issues_enabled: true`; a contact link to `docs/README.md` titled "Documentation and tutorial".
4. `PULL_REQUEST_TEMPLATE.md`: checklist: gates run (`ruff check .`, `ruff format --check .`, `mypy`, `pytest`, `python scripts/check_contract.py`), fixture regenerated if the generator changed, CHANGELOG line under `[Unreleased]`, an ADR for a decision, no real genotype data, no AI attribution trailer, contract changes start in backcross.
5. `dependabot.yml`: `package-ecosystem: github-actions`, `directory: /`, `schedule: interval: monthly`. No pip ecosystem (the pins are lower bounds; churn without benefit).
6. `SECURITY.md`: supported version (the latest release); what the tool does with data (files are read locally or in the browser tab; the only network use is the opt-in BrAPI loader on the command line, which sends the token only in the `Authorization` header, docs/adr/0024); how to report: GitHub private vulnerability reporting on the repository (the maintainer enables it: Maintainer steps, step 6), response within 14 days; no bug bounty.
7. `CODE_OF_CONDUCT.md`: the Contributor Covenant 2.1 text (CC BY 4.0; keep its attribution paragraph and link `https://www.contributor-covenant.org/version/2/1/code_of_conduct/`), enforcement contact `{{CONDUCT_CONTACT}}` written as the literal text "the maintainer, through a GitHub issue or the repository's private reporting form" (Q10).
8. `CONTRIBUTING.md`: Setup rewritten against pyproject: under Q1 = core, "The package installs shiny and pandas; the extras are `dev` (pytest, pytest-cov, ruff, mypy), `e2e` (pytest-playwright and axe-playwright-python), `export` (shinylive)"; the stale line 10 ("shinylive and build") corrected; the M3 venv recipe from CLAUDE.md line 62; "after changing `_version.py`, re-run `pip install -e .` so `tests/test_version.py` sees the new metadata"; the gates including `git diff --exit-code -- tests/fixtures src/progeny_selector/examples`; a `## Releasing` section: "Releases are made by the maintainer only: set the version in `src/progeny_selector/_version.py`, CHANGELOG.md and CITATION.cff, push a `v<version>` tag, and `.github/workflows/release.yml` builds, publishes to PyPI by trusted publishing and creates the GitHub release with the CHANGELOG section as notes." Conventions and Code rules unchanged.
9. `tests/test_github_templates.py`: each `.github/ISSUE_TEMPLATE/*.yml` loads with `yaml.safe_load` and has `name`, `description`, a non-empty `body`; the bug template's body text contains `progeny-selector --version` and `Never attach real genotype data`; `config.yml` has `blank_issues_enabled`; `dependabot.yml` has one update with `package-ecosystem: github-actions`. Discriminating: removing the data warning or breaking the YAML fails.

Acceptance: GitHub's issue chooser shows the two forms (verified by the maintainer after push); `tests/test_github_templates.py` green.

## Phase 8: CHANGELOG release section, CITATION.cff, cross-links, ADR 0036, docs consistency

Model sonnet for the files, then the main session for PLAN.md, CLAUDE.md and the verification block; end-of-milestone reviewer pass over the whole M4 diff. Depends on 4, 5, 6, 7. Files: `CHANGELOG.md`, `CITATION.cff` (new), `README.md` (Citing, Related tool), `.github/workflows/ci.yml` (cffconvert step), `tests/test_changelog_section.py` (repo-level tests), `docs/adr/0036-first-public-release-packaging.md` (new), `docs/adr/0003-hosting-local-shiny-and-shinylive-pages.md` (dated amendment), `PLAN.md`, `CLAUDE.md`.

1. `CHANGELOG.md`: keep `## [Unreleased]` empty at the top (one line: "Nothing yet."); merge the current `[Unreleased]` entries (lines 5-106) and the `## [0.1.0] - unreleased` section (line 107 onward) into one `## [{{VERSION}}] - unreleased` section (the version copied from `_version.py`), grouped Added / Changed / Fixed / Documentation with the duplicate `### Added` and `### Changed` headings merged and entry order preserved; add `### Known limitations` with the eleven README bullets verbatim; add the M4 lines the main session has collected (the `app` and `example` commands, the Load-example control, `--version`, `--top` help text, packaging metadata, CI jobs, release workflow, community files, docs). Link block at the end:
   ```
   [Unreleased]: https://github.com/piercetaylor/progeny-selector/compare/v{{VERSION}}...HEAD
   [{{VERSION}}]: https://github.com/piercetaylor/progeny-selector/releases/tag/v{{VERSION}}
   ```
   (`{{VERSION}}` written as the `_version.py` value.) This is the one CHANGELOG edit a doer makes in M4; the main session reviews the diff.
2. `CITATION.cff`:
   ```yaml
   cff-version: 1.2.0
   message: "If you use this software, please cite it as below."
   type: software
   title: progeny-selector
   version: <copy of _version.py>
   date-released: unreleased
   license: MIT
   repository-code: https://github.com/piercetaylor/progeny-selector
   url: https://piercetaylor.github.io/progeny-selector/
   abstract: "Ranks and selects progeny in marker-assisted backcross programs: foreground and avoid loci, recurrent-parent genome recovery on carrier and non-carrier chromosomes, linkage-drag bounds, QC flags and a weighted composite score, from VCF, HapMap or wide CSV genotypes."
   keywords: [plant breeding, marker-assisted backcrossing, recurrent parent genome, linkage drag, soybean, genotype]
   authors:
     - family-names: Taylor
       given-names: Pierce
   # references:            # {{BACKCROSS_REFERENCE}}: uncomment and complete if the maintainer wants the sibling cited
   #   - type: software
   #     title: backcross
   #     repository-code: https://github.com/piercetaylor/backcross
   #     authors:
   #       - family-names: Taylor
   #         given-names: Pierce
   ```
   `{{CITATION_AUTHORS}}` (ORCID, affiliation) and `date-released` are the maintainer's. `doi` is added after the first Zenodo archive (Maintainer steps, step 12).
3. `ci.yml` `package` job: `pip install cffconvert` and `cffconvert --validate` before the build (CFF 1.2.0 schema validation; reproducible, standard tool).
4. README: `## Citing` becomes: the CFF file, "GitHub shows a Cite this repository button; a Zenodo DOI is added after the first release", and a plain-text citation line `Taylor, P. ({{year}}). progeny-selector (version {{VERSION}}) [Computer software]. https://github.com/piercetaylor/progeny-selector` with the values copied from CITATION.cff at write time. `## Related tool`: "backcross characterises finished near-isogenic lines from the same input files (shared contract 1.12.0): https://github.com/piercetaylor/backcross" plus the `{{SAME_DAY}}` sentence.
5. `tests/test_changelog_section.py` additions: `test_changelog_has_a_section_for_the_current_version` (`section(CHANGELOG.md text, __version__)` finds a heading, whatever its date; the date is not asserted here because it is `unreleased` until the tag); `test_citation_version_matches_the_package` (parse `CITATION.cff` with `yaml.safe_load`; `data["version"] == __version__`; `data["cff-version"] == "1.2.0"`; `data["authors"]` non-empty). Discriminating: bumping `_version.py` alone fails both.
6. `docs/adr/0036-first-public-release-packaging.md` (MADR): context (audit 2026-09-29); decisions: version single-sourced in `_version.py` read by hatchling; shiny and pandas as core dependencies (or the extra, per Q1) with the exact-install-line message; `progeny-selector app` and `example` commands; examples shipped as package data written by the generator, proven byte-identical; sdist contents (Q7); trusted publishing with OIDC and the `pypi` environment, the tag-equals-version and CHANGELOG-date guards; plain-markdown docs (Q5); Windows and macOS unit jobs, browser tests on Linux only; consequences (wheel grows by the six example files; a release needs three files edited by hand, checked by two tests and one workflow guard).
7. `docs/adr/0003` amendment, dated: "The CI workflow builds the wheel" became true with the `package` job on this date; the `release.yml` publishes it.
8. Main session, `PLAN.md`: "Milestones" gains `M4 first public release (complete)` with its acceptance; "Deployment and cost" line 121: `pip install progeny-selector` then `progeny-selector app`; "Testing and CI" line 117: the 3.11-3.13 matrix, `unit-os`, `package`, `release.yml`; "Technology decisions" line 84: the matrix; "Core algorithms" line 64: the 2 Mb default (docs/adr/0033); "QC" line 74: Purdy labels and `filial_known` (docs/adr/0034); "UI walkthrough" line 80: the token-profile and crop selects and the Load-example control; "Repository layout" lines 88-113: `_version.py`, `examples/`, the bc3f1 and brapi fixtures, the new workflows and community files; a dated "M4 verification" block with the CI run ids and the fresh-venv smoke output; the Handoff block. `CLAUDE.md`: State paragraph, the gate pathspec (line 24), the extras sentence (line 62).

Acceptance: `tests/test_changelog_section.py` green; `cffconvert --validate` passes in `package`; the end-of-milestone reviewer pass is recorded; nothing is tagged.

## M4 acceptance checklist

- [ ] `pip install` of the built wheel into a fresh venv, from outside the checkout, runs `progeny-selector --version`, `example`, `validate`, `rank` (11 pass), `select --top 3` (6 rows), `rank` on the BC3F1 example, imports `progeny_selector.app.app`, and `app --no-browser` serves a page titled `progeny-selector` (CI `package` job).
- [ ] The Load-example control reaches `500 markers, 40 progeny; 11 pass hard filters` under `shiny run` and under the Shinylive export with every request on the page's origin.
- [ ] `tests/fixtures/` and `src/progeny_selector/examples/` regenerate byte-identical from `scripts/make_fixture.py`; `tests/test_examples.py` compares the copies.
- [ ] One version source: `pyproject.toml` has no literal, `--version` prints it, `tests/test_version.py` and the CITATION test fail on disagreement.
- [ ] CI green on 3.11, 3.12, 3.13, Windows, macOS, `package`, `e2e`, `export`, `r-reader`; `release.yml` committed and never triggered by a doer.
- [ ] README links absolute and resolving; README and tutorial `sh` blocks execute; glossary and tutorial exist; `docs/README.md` index exists.
- [ ] CITATION.cff validates; issue forms, PR template, SECURITY.md, CODE_OF_CONDUCT.md, dependabot present; CONTRIBUTING.md correct against pyproject.
- [ ] CHANGELOG has one release section with compare links and Known limitations; ADR 0036 written; PLAN.md and CLAUDE.md updated; `constants.RESULTS_SCHEMA` is `"1.2.0"`; `contract/` untouched.

# Maintainer steps

Executed only by Pierce, after phase 8 is merged and CI is green on `main`. Nothing in any phase depends on these having happened.

1. Decide `{{VERSION}}` with backcross. Set it in exactly three places: `src/progeny_selector/_version.py`, the `## [VERSION] - ...` heading and link block of `CHANGELOG.md`, and `version:` in `CITATION.cff`. Run `pip install -e .` and `pytest` (`tests/test_version.py`, `tests/test_changelog_section.py` fail if the three disagree).
2. Set `{{RELEASE_DATE}}`: replace `unreleased` with `YYYY-MM-DD` in the CHANGELOG heading and in `CITATION.cff` `date-released`. Run `python scripts/changelog_section.py <VERSION> --check-citation CITATION.cff` and read the printed notes.
3. Complete `{{CITATION_AUTHORS}}` (ORCID, affiliation) and decide `{{BACKCROSS_REFERENCE}}` (uncomment or delete the `references:` block); `cffconvert --validate`. Decide `{{SAME_DAY}}` and fix the README "Related tool" sentence.
4. Commit and push these edits to `main`; wait for CI green, including `package`.
5. Register the PyPI trusted publisher (no token). At https://pypi.org/manage/account/publishing/ under "Add a new pending publisher", GitHub tab: PyPI Project Name `progeny-selector`; Owner `piercetaylor`; Repository name `progeny-selector`; Workflow name `release.yml`; Environment name `pypi`. The project name is free as of the audit; a pending publisher reserves it on first publish.
6. On GitHub, Settings: create the environment `pypi` (Settings, Environments, New environment; optionally require your own review); under Code security, enable "Private vulnerability reporting" (SECURITY.md points there). No repository variable is needed for the release; `ENABLE_PAGES` stays `true`.
7. Optional, Zenodo: at https://zenodo.org/account/settings/github/ switch on `piercetaylor/progeny-selector` before tagging; Zenodo archives each GitHub release and reads CITATION.cff for metadata.
8. Tag and push: `git tag -a v<VERSION> -m "progeny-selector <VERSION>"` on the CI-green commit, then `git push origin v<VERSION>`. Watch the `release` workflow: `build` (tag equals version, notes extracted), `publish-pypi`, `github-release`.
9. If `build` or `publish-pypi` fails, nothing was uploaded: fix on `main`, delete the tag locally and remotely (`git tag -d v<VERSION>`; `git push origin :refs/tags/v<VERSION>`), re-tag the fixed commit. A version that did reach PyPI cannot be re-uploaded; the next attempt must be a new version.
10. Verify the PyPI page https://pypi.org/project/progeny-selector/: README renders with working links, classifiers, project URLs, licence expression.
11. Verify the fresh install from PyPI on this machine, from a directory that is not the checkout: `py -3.12 -m venv %TEMP%\ps-fresh`, `%TEMP%\ps-fresh\Scripts\pip install progeny-selector`, then `progeny-selector --version`, `progeny-selector example --out ex`, the `rank` and `select` lines from the README, and `progeny-selector app` (a browser tab opens; press Load example; Ctrl+C). Record the output in the PLAN.md M4 verification block.
12. After the GitHub release exists: if Zenodo was enabled, copy the concept DOI into `CITATION.cff` (`identifiers: - type: doi`) and the README Citing section in a follow-up commit; that commit is not a release.
13. Update the PLAN.md Handoff block and the CLAUDE.md State line with the tag, the PyPI URL and the release run id.

# Flagged questions

Each item is a decision the repository does not settle, or an inconsistency between PLAN.md, the ADRs and the code. Nothing is resolved silently; the recommendation is what the phases assume unless overruled. Criterion throughout: reproducible with standard tools, works for academic users on Windows, macOS and Linux, never silently wrong.

1. **Q1, shiny and pandas: `app` extra or core dependencies.** Recommendation: core dependencies; delete the `app` extra; `dev` lists only the tools. Grounds: the interface is the breeder's primary surface; `pipx install progeny-selector` then `progeny-selector app` works with one line on all three platforms; a forgotten `[app]` is the most likely first-run failure for a non-developer; the cost is about fifteen pure-Python packages for an HPC user who only runs `rank`, which changes nothing in behaviour. Either way `cmd_app` prints the exact install line when the import fails (a `--no-deps` install), so a bare install never fails confusingly. Alternative: keep the extra and use the PEP 508 self-reference `progeny-selector[app]` in `dev` so the lists cannot drift.
2. **Q2, console command name and shape.** Recommendation: a subcommand of the existing entry point, `progeny-selector app [--host 127.0.0.1] [--port 8000] [--no-browser]`, calling `shiny.run_app("progeny_selector.app.app:app", ...)` with `launch_browser=True` by default. One console script to find on PATH, one `--help`. Alternatives: a second script `progeny-selector-app`; `python -m progeny_selector.app` (fails the pipx case, where only console scripts are on PATH).
3. **Q3, a provenance column in the uncommitted results_schema 1.2.0.** Recommendation: yes, one column `tool_version` (the value of `__version__`, for example `0.1.0`) in results.csv after `background_max_marker_coverage` and in selected.csv after `sample_db_id`, both before `token_profile`, folded into the 1.2.0 diff before it is committed (phase 0b, opus, reviewer, ADR 0033 paragraph). No commit hash: a pip or pipx install has no repository; a build-time hash (hatch-vcs) makes a fresh clone and the Shinylive staged copy fail to import without a generated file, which is exactly the fragility a first release should not carry; the reproducible provenance of a release is `tool_version` plus `results_schema`; a run from a checkout is identified by its commit in the user's notes. This keeps the column count change at one bump. backcross adds app-version plus commit; the two files then differ in that its column carries a hash and this one does not, which should be recorded in both ADRs. If answered no, `tool_version` waits for 1.3.0 and is not in M4.
4. **Q4, Load example under Shinylive.** Recommendation: yes, the same control. The example files reach Pyodide because `scripts/build_shinylive.py` copies the whole `src/progeny_selector` tree (`examples/` included) into the staging directory and `shinylive export` bundles every file there into `app.json`; `importlib.resources.files` resolves to real paths on Pyodide's file system, so `load_dataset` reads them as on CPython. Cost: about six files, roughly 300 KB (not measured; the export test prints nothing about size, `scripts/build_shinylive.py` prints the site size), against a 44.5 MB site. Not verifiable here whether `.vcf` is stored as text or base64; either works, and phase 3's `test_shinylive_load_example_stays_on_origin` is the proof. Nothing in `build_shinylive.py` changes.
5. **Q5, documentation site.** Recommendation: plain Markdown on GitHub with `docs/README.md` as the index and the `Documentation` project URL pointing at it; absolute links in README so the PyPI page works. No MkDocs or Sphinx: Pages already serves the app at the site root, a docs site would need a second deploy path, a theme dependency and a build job, and the audience reads GitHub. Revisit if the docs exceed a dozen pages.
6. **Q6, Development Status classifier.** Recommendation: `4 - Beta` whichever version is chosen: the pipeline is verified on real SoySNP50K data (PLAN.md M2 verification), measured (docs/limits.md) and reviewed for accessibility. Alternative `3 - Alpha` if the version stays 0.1.0 and the maintainer wants the classifier to say so; `5 - Production/Stable` only with 1.0.0.
7. **Q7, sdist contents.** Recommendation: exclude `CLAUDE.md`, `PLAN.md`, `docs/*-phases.md`, `.github`, `.env.example`, `site`, `build`, `data`; keep `tests/`, `docs/`, `contract/`, `scripts/`, `LICENSE`, `CHANGELOG.md`, `CITATION.cff`. The working documents describe an agent workflow, not the software, and the audit lists them as shipped; downstream packagers want the tests. `contract/cases/` includes gzip and bgzip cases and a CRLF case by design; the sdist keeps them byte-exact (hatchling does not normalise). Alternative: exclude `tests/` too and ship a 100 KB smaller sdist.
8. **Q8, `example` command name and default.** Recommendation: `example`, default `--out progeny-selector-example`, refuses to overwrite without `--force`. Alternatives: `demo` (backcross's `?demo=synthetic` vocabulary), `init`. The BC3F1 directory is made self-contained by copying bc2f1's `markers.csv` and `criteria.yaml` at write time, so package data stays exactly the fixture files.
9. **Q9, `--top N` behaviour.** Recommendation: keep the behaviour (`rank_in_family <= N`, up to N per family; `--overall` for N total) and fix the help text and README; it matches `core/selection.py::select_top_n` and the Selection screen's "add top N per family". Alternative: rename to `--per-family N` with `--overall N`, a breaking CLI change one release before anyone depends on it.
10. **Q10, conduct and security contacts.** Recommendation: no personal email in the repository; the CODE_OF_CONDUCT enforcement contact is "the maintainer, through a GitHub issue or the repository's private reporting form", and SECURITY.md uses GitHub private vulnerability reporting, which the maintainer enables (step 6). Alternative: a dedicated address.
11. **Q11, fixture regeneration on the Windows and macOS jobs.** Recommendation: unit tests only on those runners; the reproducibility diff stays on ubuntu (proved there and on the maintainer's Windows laptop, where `.gitattributes` normalises). Adding it costs one step and may reveal a `csv.writer` CRLF difference on the runner's git configuration; if wanted, it is a one-line addition to `unit-os`.
12. **Q12, `src/progeny_selector/app/requirements.txt`.** The file says it is unused by both `shiny run` and `scripts/build_shinylive.py`. Recommendation: delete it in phase 3 and drop the sentence from docs/adr/0009's Neutral clause with a dated line; it ships in the wheel today for no reader. Alternative: keep it as the documented fallback for a direct `shinylive export`.
13. **Q13, a TestPyPI rehearsal.** Recommendation: none. The `package` job installs the built wheel and the sdist into fresh venvs on every push, which is the rehearsal that matters; a TestPyPI publisher would need a second pending-publisher registration and a workflow input. If the first real publish fails, step 9 covers it.
14. **Q14, Windows and macOS Python version in `unit-os`.** Recommendation: 3.12 only, one job each, to keep the matrix at seven jobs. Alternative: 3.11 and 3.13 too, four more jobs per push.
15. **Q15, cffconvert in CI.** Recommendation: yes, `cffconvert --validate` in the `package` job (phase 8); it is the standard CFF validator and takes seconds. Alternative: rely on GitHub's own rendering of the Cite button.
16. **Inconsistency, PLAN.md:64 vs docs/adr/0033.** PLAN.md says the coverage cap default is "10 cM or 4 Mb"; ADR 0033 and `constants.DEFAULT_MAX_COVERAGE["bp"]` say 2,000,000 once 1.2.0 lands. Phase 8 corrects PLAN.md.
17. **Inconsistency, PLAN.md:84 and :117 vs the CI matrix.** "GitHub Actions on Python 3.11 and 3.12" becomes 3.11-3.13 plus Windows and macOS in phase 4. Phase 8 corrects PLAN.md; CLAUDE.md line 60 ("the 3.11 half of the CI matrix") becomes "the 3.11 third".
18. **Inconsistency, PLAN.md:121 and README:24 vs the code.** "`pip install progeny-selector` and `shiny run`" cannot work: `shiny run` needs a path into a checkout. Phase 3 adds the command; phase 8 rewrites the line.
19. **Inconsistency, docs/adr/0003 Decision Outcome vs CI.** "The CI workflow builds the wheel" has been false since M0 (the audit: no build job). Phase 4 makes it true; phase 8 adds a dated amendment rather than editing the accepted text.
20. **Inconsistency, CONTRIBUTING.md:8-10 vs pyproject.** Line 10 lists `build` in `export` (removed after the M2 low in PLAN.md Handoff; `export = ["shinylive>=0.8"]`); line 9 omits `axe-playwright-python` from `e2e`; line 8's description of `dev` changes under Q1. docs/adr/0008 records `export = ["shinylive>=0.8", "build"]` as history and is left alone. Phase 7 fixes CONTRIBUTING.
21. **Inconsistency, cli.py docstring vs parser.** The interface block omits `--coding` (validate, rank) and `--project` (select). Phase 1 completes it.
22. **Inconsistency, PLAN.md:74 vs docs/adr/0034.** The QC paragraph predates Purdy labels and `filial_known`. Phase 8 adds one sentence.
23. **Inconsistency, PLAN.md:80 vs app/screens/load.py.** Screen 1's description lacks the token-profile select, the custom-profile file input and the crop select (M2, S6), and will lack the Load-example control. Phase 8 updates it.
24. **Inconsistency, PLAN.md:88-113 vs the tree.** The layout lists one fixture directory and four test files; the tree has `synthetic_bc3f1/`, `brapi/`, 56 test files, ten scripts. Phase 8 updates the layout to directories and counts rather than file lists.
25. **Inconsistency, README:36 vs CITATION.cff.** "cite this repository with the commit or version used" is replaced by the CFF file in phase 8; the same sentence in `../backcross/README.md:34` is the sibling's to change (main session notes it in that repository's handoff; nothing here edits `../backcross`).
26. **Inconsistency, `.gitignore` vs package data.** `*.vcf` is ignored globally with exceptions only for `tests/fixtures/**` and `contract/**`; hatchling excludes VCS-ignored files, so the shipped example would be missing from git and the wheel without phase 2's `!src/progeny_selector/examples/**`. Recorded here so the exception is understood as load-bearing.
27. **Inconsistency, PLAN.md Handoff rule vs phase 8.** "CHANGELOG under `[Unreleased]`, written only by the main session, never by a doer" (docs/m3-phases.md:5); phase 8 has a doer restructure the file into the release section. Recommendation: allow it for this one phase with the main session reviewing the diff, or the main session performs step 1 of phase 8 itself.
28. **Inconsistency, `tests/test_version.py` and editable installs.** `importlib.metadata.version` reflects the metadata at install time; after a version bump the test fails until `pip install -e .` is re-run. Recommendation: accept (CONTRIBUTING says so); the alternative, skipping when metadata is absent, is the vacuous-pass pattern PLAN.md warns against.
29. **Open placeholders, not mine to decide:** `{{VERSION}}` (0.1.0 vs 1.0.0, joint with backcross), `{{RELEASE_DATE}}`, `{{CITATION_AUTHORS}}`, `{{BACKCROSS_REFERENCE}}`, `{{SAME_DAY}}`. Every phase runs with the repository's current `0.1.0`; only Maintainer steps 1-3 change them.