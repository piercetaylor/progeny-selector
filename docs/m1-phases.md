# M1 vertical slice: implementation phases

This document decomposes milestone M1 (PLAN.md, "Milestones") into six phases, one commit each, in this repository only. Every commit passes the gates of CLAUDE.md. It adds only what an implementer needs to make no further design decision; where a fact could not be verified this session it says so instead of asserting it. Dates are 2026-09-12; every `[web]` citation is a page fetched today. The M0 verification block in PLAN.md ("Not verified": `shinylive export`, the screens beyond the shell, the workflow itself, real-data behaviour) is M1's inheritance; each item is allocated to a phase below.

Status: questions 1, 2, 3, 5, 6, 7 and 9 were settled by the maintainer on 2026-09-12. Questions 4 and 8 are open and must be answered before phase 6. Browser tests run locally as well as in CI; the maintainer approved installing pytest-playwright and Chromium if not already present.

## Invariants that every phase respects

1. `core/` stays pure: numpy and dataclasses in and out; no I/O, no pandas, no Shiny. `io/` is the only place that validates and serialises. `app/` renders rows and calls `run_analysis`; it computes no genetics. New logic that M1 needs (tree derivation, row filtering, strip segments, QC table rows, criteria serialisation) goes in `core/` or `io/` with unit tests. The one exception is `app/present.py`: Shiny-free colour and geometry helpers (status to style dict, strip rectangles) that import only `constants` and `core`, are unit-tested without a browser, and touch no metric (Q3).
2. No parser or auto-detection behaviour changes; `src/progeny_selector/io/calls.py`, `vcf.py`, `hapmap.py`, `wide_csv.py`, `manifest.py` are not edited. The shared input sections of `docs/data-formats.md` (Chromosome names through markers.csv) are not edited; S1 in `../isoline-browser/docs/m3-phases.md` owns them and is not scheduled here. M1 may edit only the criteria.yaml and Outputs sections, and only additively.
3. No documented column, parameter name or output format changes. `AnalysisResult.rows` keys and the three writers are unchanged. A new file the tool writes (criteria.yaml download) uses the existing criteria contract.
4. Never commit real genotype data. `tests/fixtures/` is not touched in M1: no metric is added, so `scripts/make_fixture.py` is not edited and the fixture gate passes trivially. Hand-built matrices from `tests/conftest.py::make_matrix` carry the new unit tests, which is what CLAUDE.md asks for when independence of interpretation matters.
5. Colours come only from `constants.STATE_COLORS` and `STATUS_COLORS`. Every chip and strip carries a text label or `<title>`, so colour is never the only cue (PLAN.md, "UI walkthrough").
6. `pytest` with no arguments stays the per-commit gate and needs no browser: browser tests live under `tests/e2e/`, carry the `e2e` marker, are excluded by `addopts`, and run through `pytest -m e2e` after `playwright install chromium`. CI runs them; the laptop runs them per milestone.
7. Screen modules keep the M0 interface `view(id) -> Tag`, `server(id, state) -> None`; `AppState` keeps its six fields and gains none. `state.breadcrumb` keeps the shape `{"cross", "family", "generation"}`.
8. No AI attribution in any commit. Conventional Commits; phases 3, 4, 5 and 6 add `CHANGELOG.md` entries under `[Unreleased]`; phases 1 and 2 add none (infrastructure and internal helpers).
9. `docs/adr/` gains 0008 (phase 1) and 0009 (phase 6); ADR 0003 gets a dated amendment rather than an edit (CONTRIBUTING.md: supersede, do not edit).
10. Python 3.11 stays the floor (`ruff target-version py311`); nothing in M1 uses 3.12 syntax. `shiny>=1.4` stays the floor because every API used below exists at 1.4.0 (resolution 1).
11. Implementers never run git. The main session reviews, re-runs the gates and commits.

## Order, and why

| phase | scope                                                                                         | depends on |
| ----- | --------------------------------------------------------------------------------------------- | ---------- |
| 0     | Genotype extension handling: accept `.vcf.bgz`, reject `.bcf`                                 | nothing    |
| 1     | Extras split, e2e test infrastructure, first screen test                                      | nothing    |
| 2     | Pure helpers in core/io/present with unit tests; no UI change                                 | nothing    |
| 3     | Navigate tree with breadcrumbs; Rank grid chips, per-locus columns, view-aware selection      | 1, 2       |
| 4     | Validate/QC table; Compare screen with side-by-side cards and chromosome strips               | 1, 2       |
| 5     | Criteria editor and criteria.yaml download; download handlers migrated off `render.download`  | 1, 2       |
| 6     | Shinylive export staged with the package, smoke-tested in CI, deployed to Pages               | 1, 3–5     |

Phases 3–6 assert their acceptance through browser tests, so the harness comes first. Phase 2 has no UI and can run in parallel with phase 1; it is listed second so that phases 3–5 are thin view code over tested functions. Phase 6 is last because its smoke test drives the screens phases 3–5 build, and because it is the phase most likely to need a fallback (resolution 5), which should not block the screens.

### Allocation of the M0 "Not verified" list

| item                                   | closed by | test or record                                                  |
| -------------------------------------- | --------- | --------------------------------------------------------------- |
| `shinylive export`                     | 6         | `tests/e2e/test_shinylive_export.py`; CI `export` job           |
| Shiny screens beyond serving the shell | 1, 3–5    | `tests/e2e/test_*.py` through `shiny.pytest`                    |
| GitHub Actions workflow                | 1, 6      | first green run recorded in the PLAN.md M1 verification block   |
| Behaviour on real data                 | 6, partly | maintainer protocol, Q7; not testable in CI                     |

## Environment facts found today

- `.venv` was created from Python 3.12.4 (`C:\Users\pierc\AppData\Local\Programs\Python\Python312\python.exe`) and `pip install -e ".[dev]"` completed: shiny 1.7.0, shinylive 0.8.11. On this machine lzstring installed without the wheel-build failure PLAN.md records for the scaffold container. In that venv on 2026-09-12: ruff check and format pass (61 files), mypy passes (39 files), pytest 24 passed, the regenerated fixture has no git diff.
- `shinylive>=0.8` requires `lzstring>=1.0.4` [web] https://raw.githubusercontent.com/posit-dev/py-shinylive/main/setup.cfg, and lzstring 1.0.4 (2018) ships only an sdist [web] https://pypi.org/pypi/lzstring/json. That it builds here does not show it builds on every CI image; the scaffold container is the counter-example.
- `scripts/make_fixture.py` writes CRLF line endings on Windows. `.gitattributes` normalises them, so the fixture gate passes, but the working-copy files differ byte for byte from the committed ones. Not fixed in M1 (Q9).
- The system `python3` is 3.13 without shiny; `python` is 2.7 (CLAUDE.md). Playwright and pytest-playwright are not installed in any interpreter, so every browser-facing claim below is verified against documentation only.
- Latest shiny is 1.7.0 (2026-07-29) [web] https://pypi.org/pypi/shiny/json; latest shinylive is 0.8.11 (2026-07-29), `python_requires >=3.10` [web] https://pypi.org/project/shinylive/.
- `git remote -v` and `gh repo view` both give `piercetaylor/progeny-selector`, public; `pyproject.toml` says `piercetaylor2020`.

## Resolutions

1. **Download handlers use `@render.download_button`.** `@render.download` is deprecated since shiny 1.4.0 and `@render.download_button` / `@render.download_link` "pair 1:1 with ui.download_button() and ui.download_link()" from the same version [web] https://raw.githubusercontent.com/posit-dev/py-shiny/main/CHANGELOG.md. Signature `render.download_button(fn=None, *, filename=None, media_type=None, encoding='utf-8', label='Download', width=None)`; the decorated function returns a path or yields strings or bytes [web] https://shiny.posit.co/py/api/core/render.download_button.html. Generated files are yielded as text, so nothing is written to a temp directory in the browser tab.
2. **Chips are DataGrid cell styles, not HTML.** `render.DataGrid(data, *, width='fit-content', height=None, summary=True, filters=False, editable=False, selection_mode='none', styles=None)`; `styles` is a style-info dict, a list of them, or a callable returning a list; each has `location: "body"`, `rows` (row numbers or None), `cols` (column numbers or None), `style` (CSS dict) and `class`; `selection_mode` is `"none" | "row" | "rows"` [web] https://shiny.posit.co/py/api/core/render.DataGrid.html; `styles` exists since 1.4.0 (CHANGELOG above). Row and column indices are integers; `app/present.py` computes them from the frame's column list.
3. **Selected rows come from `data_view(selected=True)`**, "the data how it is viewed within the browser" with sort and filter applied [web] https://shiny.posit.co/py/api/core/render.data_frame.html. The M0 code indexes the unfiltered list with `cell_selection()["rows"]`, whose frame of reference after a header filter is not established, so it is replaced.
4. **The tree is a Bootstrap accordion of families with a radio group of generations per family.** py-shiny has no tree component; `ui.accordion(*args, id=None, open=None, multiple=True, ...)` reports the open panels through `input.<id>()` [web] https://shiny.posit.co/py/api/core/ui.accordion.html and `ui.accordion_panel(title, *args, value=MISSING, icon=None, **kwargs)` identifies a panel by `value` [web] https://shiny.posit.co/py/api/core/ui.accordion_panel.html; `ui.input_radio_buttons(id, label, choices, *, selected=None, inline=False, width=None)` takes a value-to-label dict [web] https://shiny.posit.co/py/api/core/ui.input_radio_buttons.html. Accordion headers are buttons and radios take arrow keys, which covers the keyboard requirement without JavaScript, and a Playwright `Accordion` controller exists [web] https://shiny.posit.co/py/api/testing/playwright.controller.Accordion.html. Breadcrumbs are Bootstrap 5 markup (`nav[aria-label=breadcrumb] > ol.breadcrumb > li.breadcrumb-item`) with `ui.input_action_link(id, label, *, icon=None, **kwargs)` [web] https://shiny.posit.co/py/api/core/ui.input_action_link.html on every ancestor crumb. Q1 offers the alternative.
5. **The Shinylive export stages a copy of the package next to a stub `app.py`.** `shinylive export` bundles every non-hidden file under the app directory except `__pycache__`, `venv` and `.venv` into `app.json` [web] https://raw.githubusercontent.com/posit-dev/py-shinylive/main/shinylive/_app_json.py; it finds Pyodide packages by scanning imports and `requirements.txt`, ignores `http(s)://` lines ("If it's a URL, then it must be a wheel file. Ignore it"), and a module it cannot map is dropped with the note "Assuming it is in base Pyodide or in requirements.txt" [web] https://raw.githubusercontent.com/posit-dev/py-shinylive/main/shinylive/_deps.py. Today's export therefore ships neither `progeny_selector` nor numpy (only `app/app.py` imports are scanned, and `app/requirements.txt` lists only pyyaml), so the "shinylive export builds" CI step proves only that the command exits 0. Staging the package makes its `import numpy` visible to the scanner and puts the code in the bundle; whether the staged directory is on `sys.path` in Pyodide is **not verified here** and is exactly what the phase 6 smoke test checks. Fallback, if it is not: CI copies the built wheel into `site/` and `requirements.txt` gains the line `https://piercetaylor.github.io/progeny-selector/progeny_selector-0.1.0-py3-none-any.whl`, which the exporter ignores and micropip installs at run time. ADR 0003's sentence "the app's requirements.txt names the package wheel" is amended either way (Q4).
6. **Screens are tested at two levels.** Pure helpers in `core/`, `io/` and `app/present.py` have pytest unit tests. Server behaviour is tested end to end with `shiny.pytest.create_app_fixture(app, scope='module', timeout_secs=30, env=None)`, which yields a `ShinyAppProc` with `.url`, resolving the path relative to the collecting test file [web] https://shiny.posit.co/py/api/testing/pytest.create_app_fixture.html, and the controllers `OutputDataFrame`, `InputFile`, `InputSelect`, `InputSwitch`, `InputActionButton`, `InputTextArea`, `PageNavbar`, `Accordion`, `DownloadButton` [web] https://shiny.posit.co/py/api/testing/. Requires `pip install pytest pytest-playwright` and a browser install [web] https://shiny.posit.co/py/docs/end-to-end-testing.html; `playwright install chromium` locally (Windows cache `%USERPROFILE%\AppData\Local\ms-playwright`), `playwright install --with-deps chromium` in CI. No render-function unit testing of `@module.server` code is attempted: the helpers are what is worth unit-testing, and the wiring is what the browser test checks.
7. **Criteria editing is a YAML text editor, not a form**, in M1: `ui.input_text_area(id, label, value='', *, width=None, height=None, cols=None, rows=None, placeholder=None, resize=None, autoresize=False, autocomplete=None, spellcheck=None, update_on='change')` [web] https://shiny.posit.co/py/api/core/ui.input_text_area.html, parsed by the strict reader of ADR 0005 so errors are the same verbatim messages the CLI prints. The applied `Criteria` in `AppState` is the source of truth; the text area is a view over it, so a field-by-field form can be added in M2 as a second view without rework. Settled 2026-09-12 (Q2).
8. **The strip is inline SVG built with htmltools tags**: no matplotlib, no JavaScript, works under Pyodide. One row per chromosome in `chrom_sort_key` order; each marker owns the interval between the midpoints to its neighbours (chromosome ends: 0 and the assembly length from `chrom_length_bp`, else the last marker); adjacent same-state intervals merge. Nothing is inferred between disagreeing markers, which keeps the strip a summary rather than a browser (PLAN.md non-goals).

## Phase 0: genotype extension handling (settled Q6)

Goal: the loader accepts exactly the extensions `docs/data-formats.md` documents, as isoline-browser already does (`src/io/loaders.ts:29` strips `.gz` or `.bgz`). No documented contract changes; the code moves to match the document.

| file | change |
| ---- | ------ |
| `src/progeny_selector/io/__init__.py` | `load_genotypes` strips a trailing `.gz` **or** `.bgz` before choosing the format; `.bcf` is no longer accepted, so `x.bcf` raises the existing "cannot infer genotype format" error (binary BCF was being read as text, and neither the contract nor the sibling accepts it). |
| `src/progeny_selector/io/vcf.py` (only if needed) | The gzip branch opens a `.bgz` file exactly as it opens `.gz`. If the branch already keys on magic bytes rather than the name, no change. |
| `tests/test_io.py` | `test_vcf_bgz_extension`: the fixture VCF gzipped to `tmp_path/"g.vcf.bgz"` loads with the same marker and sample counts as the plain file; `test_bcf_extension_rejected`: a `tmp_path/"g.bcf"` raises `DataContractError` matching "cannot infer genotype format". |
| `CHANGELOG.md` | Under Fixed: "Genotype files named `.vcf.bgz` now load, as documented; `.bcf`, which was never documented and was read as text, is rejected." |

Commit: `fix(io): accept .vcf.bgz and stop accepting .bcf`. Parsers are otherwise untouched (invariant 2 still holds for M1 phases 1–6).

## Phase 1: extras split and e2e infrastructure

Goal: the per-commit gates never depend on shinylive or a browser; a browser test drives the real Load screen with the fixture.

| file | change |
| ---- | ------ |
| `pyproject.toml` | `[project.urls] Repository = "https://github.com/piercetaylor/progeny-selector"`. Extras: `app = ["shiny>=1.4", "pandas>=2.2"]` unchanged; `dev = ["pytest>=8", "pytest-cov>=5", "ruff>=0.6", "mypy>=1.10", "shiny>=1.4", "pandas>=2.2"]` (shinylive removed); new `export = ["shinylive>=0.8", "build"]`; new `e2e = ["pytest-playwright>=0.5"]`. `[tool.pytest.ini_options]` gains `markers = ["e2e: browser tests through Playwright; excluded by default"]` and `addopts = "-q -m 'not e2e'"`. |
| `tests/e2e/__init__.py`, `tests/e2e/conftest.py` | conftest: `pytest.importorskip("playwright")`; `pytest_collection_modifyitems` adds `pytest.mark.e2e` to every item under `tests/e2e/` (a module-level `pytestmark` does not cover a directory); sets `os.environ.setdefault("PYTHONPATH", <repo>/src)` before any fixture starts so `shiny run` in the subprocess imports the package even without an editable install (whether `create_app_fixture(env=)` replaces or extends the environment is unverified; the editable install CI performs is the primary path). |
| `tests/e2e/test_load_screen.py` | Below. |
| `.github/workflows/ci.yml` | `check` job: `pip install -e ".[dev]"` unchanged; the `shinylive export builds` step is removed from `check` (it moves to phase 6's `export` job). New job `e2e` (needs `check`, Python 3.12): `pip install -e ".[dev,e2e]"`, `playwright install --with-deps chromium`, `pytest -m e2e`. |
| `CLAUDE.md` | Gates: `pytest -m e2e` listed as a per-milestone gate needing `playwright install chromium`; Environment gotchas gains the extras split and the local venv recipe (`py -3.12 -m venv .venv`, `.venv\Scripts\python -m pip install -e ".[dev,export,e2e]"`). |
| `CONTRIBUTING.md` | Setup names the four extras and which gate needs which. |
| `docs/adr/0008-screen-testing-and-extras.md` | MADR: pure helpers unit-tested, server code tested through `shiny.pytest` and Playwright controllers, `e2e` marker excluded by default; shinylive isolated in `export` so the per-commit gates never depend on lzstring's sdist. |

`tests/e2e/test_load_screen.py`: `app = create_app_fixture(["../../src/progeny_selector/app/app.py"])`; `test_load_fixture(page, app)`: `page.goto(app.url)`; `controller.InputFile(page, "load-genotypes").set(FIXTURE / "genotypes.vcf")` and likewise `load-samples`, `load-markers`, `load-criteria` (module ids are `<module>-<id>`); `controller.InputActionButton(page, "load-run").click()`; a verbatim-text controller under a known name is not verified, so the test asserts `page.locator("#load-status")` `to_contain_text("500 markers, 40 progeny; 11 pass hard filters")` with a 60 s timeout. `InputFile.set(file_path, *, timeout=None, expect_complete_timeout=30000)` [web] https://shiny.posit.co/py/api/testing/playwright.controller.InputFile.html. The implementer checks the exact status text the M0 Load screen emits and uses it verbatim.

Acceptance as tests: `pytest` (26 passing after phase 0, e2e excluded, no browser needed); `pytest -m e2e` passes with Chromium installed; `ruff`, `mypy`, fixture gate pass; the CI `e2e` job is green on the first push.

Not verifiable here: `create_app_fixture(env=)` semantics; the name of the verbatim-text controller; lzstring on the CI image.

## Phase 2: pure helpers with unit tests

Goal: everything phases 3–5 display is a tested function.

| file | change |
| ---- | ------ |
| `src/progeny_selector/core/navigation.py` (new) | Interface below. |
| `src/progeny_selector/core/strip.py` (new) | Interface below. |
| `src/progeny_selector/core/qc.py` | Adds `uninformative_summary` and `qc_table_rows` (below); fixes the header docstring, which lists `sample_qc(states, gm, …)` while the function is `sample_qc(gm, dataset, classification, filters)`. |
| `src/progeny_selector/core/score.py` | `QC_EXCLUDING_FLAGS` unchanged; `qc_table_rows` imports it. |
| `src/progeny_selector/io/criteria.py` | Adds `criteria_to_dict`, `dump_criteria_yaml`, `read_criteria_text` (below). |
| `src/progeny_selector/app/present.py` (new) | Shiny-free; interface below. |
| `tests/test_navigation.py`, `tests/test_strip.py`, `tests/test_qc_table.py`, `tests/test_criteria_roundtrip.py`, `tests/test_present.py` (new) | Below. |

```python
# core/navigation.py
@dataclass(frozen=True)
class GenerationNode: generation: str | None; n: int
@dataclass(frozen=True)
class FamilyNode: family_id: str | None; n: int; generations: tuple[GenerationNode, ...]
@dataclass(frozen=True)
class NavTree: cross: str; n: int; families: tuple[FamilyNode, ...]
def build_tree(dataset: Dataset) -> NavTree            # cross = f"{rp.line_name} x {donor.line_name}"; progeny only; None sorts last
def filter_rows(rows: list[dict], family: str | None, generation: str | None) -> list[dict]   # None = no filter; order preserved
def crumb_labels(tree: NavTree, family: str | None, generation: str | None) -> list[str]     # ["cross"], ["cross","F1"], ["cross","F1","BC2F1"]

# core/strip.py
@dataclass(frozen=True)
class Segment: start: float; end: float; state: int    # fractions of the chromosome length, 0 <= start < end <= 1
@dataclass(frozen=True)
class ChromStrip: chrom: str; length_bp: float; n_markers: int; segments: tuple[Segment, ...]
@dataclass(frozen=True)
class LocusTick: chrom: str; locus_id: str; kind: str; start: float; end: float   # kind 'target' | 'avoid'
def chromosome_strips(states: np.ndarray, gm: GenotypeMatrix) -> list[ChromStrip]   # states: int8 (n_markers,), gm sorted by position
def locus_ticks(targets: dict[str, ResolvedLocus], avoid: dict[str, ResolvedLocus], strips: list[ChromStrip]) -> list[LocusTick]

# core/qc.py
def uninformative_summary(classification: Classification) -> list[tuple[str, int]]   # (reason, count) sorted by count desc, then reason; '' excluded
def qc_table_rows(qc: list[SampleQC], dataset: Dataset) -> list[dict]
    # keys: sample_id, line_name, family_id, generation, missing_rate, het_rate, expected_het, hom_donor_rate,
    # nonparental_rate, expected_rpp, ibs_rp, ibs_donor, flags ('|'-joined), excluded (bool: any flag in QC_EXCLUDING_FLAGS)

# io/criteria.py
def criteria_to_dict(criteria: Criteria) -> dict
    # key order: name (if set), targets, avoid, flank_window, flank_unit, background, weights, filters;
    # loci carry locus_id, exactly the keys of their kind, rule, min_markers, notes (if set), then required_state/flank_left/flank_right or allow_het;
    # background omits max_marker_coverage when None; weights and filters list every field explicitly
def dump_criteria_yaml(criteria: Criteria) -> str    # yaml.safe_dump(doc, sort_keys=False, allow_unicode=True), LF line endings
def read_criteria_text(text: str) -> Criteria        # yaml.safe_load then criteria_from_dict; empty text -> CriteriaError("criteria document must be a mapping")

# app/present.py
STATUS_COLUMN_RE = re.compile(r"^(target|avoid)_.+_status$")
def status_columns(columns: Sequence[str]) -> list[str]
def status_cell_styles(columns: Sequence[str], rows: list[dict]) -> list[dict]
    # one entry per (status column, status value) present:
    # {"location":"body","rows":[...],"cols":[j],"style":{"background-color":STATUS_COLORS[v],"color":"#000000","font-weight":"600"}}
def qc_row_styles(rows: list[dict]) -> list[dict]    # excluded rows: vermilion at 25% (rgba from STATUS_COLORS["fail"]); flagged-not-excluded: yellow at 35%
def chip_style(status: str) -> str                   # inline CSS for a ui.span badge; black text on STATUS_COLORS[status]
def strip_rects(strips: list[ChromStrip], ticks: list[LocusTick], width_px: int = 320, row_px: int = 8, gap_px: int = 3, label_px: int = 34) -> StripGeometry
    # StripGeometry(width, height, rows: list[StripRow(chrom, y, label_x, rects: list[Rect(x, w, fill, state_label)], ticks: list[TickShape])]);
    # bar length proportional to length_bp / max length_bp; every rect fill from STATE_COLORS
```

Acceptance as tests (with `make_matrix` / `make_dataset` from conftest, or the fixture as loaded by `test_smoke_pipeline`). Literal values below (line names, counts, sample ids) were read from the fixture by the planner; the implementer confirms each against `tests/fixtures/synthetic_bc2f1/` before hard-coding it and reports any that differ rather than adjusting silently.

- `test_navigation.py`: fixture dataset gives `cross == "Williams 82 (synthetic) x PI synthetic donor"`, two families F1, F2 with n 20 each, each one generation BC2F1; a hand-built manifest with `family_id=None` on one sample yields a trailing `FamilyNode(None, 1, …)`; `filter_rows(rows, "F1", None)` returns 20 rows in the incoming order; `filter_rows(rows, None, "BC3F1")` returns `[]`; `crumb_labels` for the three depths.
- `test_strip.py`: `make_matrix(["A","H","A"])` (1, 2, 3 Mb on Gm01, length 57,932,356) gives one strip with segments `[(0, 1.5e6/L, A), (1.5e6/L, 2.5e6/L, H), (2.5e6/L, 1, A)]` within 1e-9; states `["A","A","A"]` merge into one segment; `["N","A"]` keeps N as its own segment (missing is drawn grey, not filled); a marker beyond the assembly length extends `length_bp` to that position; two chromosomes come out in `chrom_sort_key` order; `locus_ticks` for a region locus spans `start_bp/L .. end_bp/L`, for a marker locus a zero-width tick widened to `0.004` by `strip_rects`.
- `test_qc_table.py`: on the fixture result, `qc_table_rows` has 40 rows, `BC2F1-F2-002` has `excluded is True`, `BC2F1-F1-003` has `high_missing` in `flags` and `excluded is False` (high_missing is not a QC-excluding flag; the missing-rate filter excludes it), `uninformative_summary` sums to 25 with `"parents identical (monomorphic)"` first at 20.
- `test_criteria_roundtrip.py`: `criteria_from_dict(criteria_to_dict(c)) == c` for the fixture criteria and for a hand-built `Criteria` with a region target carrying `flank_left`, a flanking avoid with `allow_het=True`, and `background.max_marker_coverage=8`; `dump_criteria_yaml` output starts with `name:` for the fixture and with `targets:` when `name` is None; `read_criteria_text(dump_criteria_yaml(c)) == c`; `read_criteria_text("")` and `read_criteria_text("- a")` raise `CriteriaError`; a dumped region appears as `chrom`/`start_bp`/`end_bp` (never the shorthand).
- `test_present.py`: `status_cell_styles(["sample_id","target_T1_status"], rows)` with statuses pass/fail/pass yields two entries with `cols == [1]` and rows `[0, 2]` and `[1]`; every colour string equals a `STATUS_COLORS` value; `strip_rects` on the `test_strip` strips yields rect widths summing to the bar length within 0.5 px and one `StripRow` per strip; importing `progeny_selector.app.present` succeeds with `shiny` absent (`monkeypatch.setitem(sys.modules, "shiny", None)` then reload).

Not verifiable here: nothing; this phase is pure Python.

## Phase 3: Navigate tree and Rank grid

Goal: the breeder picks cross, family and generation with mouse or keyboard, sees where they are, and the Rank grid shows only that node with status chips and per-locus columns.

| file | change |
| ---- | ------ |
| `src/progeny_selector/app/screens/navigate.py` | `ui_`: `ui.card(ui.card_header("Navigate: cross > family > generation"), ui.output_ui("crumb"), ui.output_ui("tree"))`. Server: `tree = reactive.calc(lambda: build_tree(state.dataset()) if state.dataset() else None)`; `@render.ui def tree()` returns `ui.accordion(*panels, id="families", multiple=False, open=False)` with, per family index `i`, `ui.accordion_panel(f"{label} ({fam.n})", ui.input_radio_buttons(session.ns(f"gen_{i}"), "Generation", choices={"": f"(all generations) ({fam.n})", **{g.generation: f"{g.generation} ({g.n})" for g in fam.generations if g.generation}}, selected=""), value=f"fam_{i}")`; label is `family_id` or `"(no family)"`. `@reactive.effect def _publish()` reads `input.families()` (open panel values; `None` or empty means no family), maps `fam_i` back to the `FamilyNode`, reads `input[f"gen_{i}"]()` (`""` means None), and sets `state.breadcrumb`. `@render.ui def crumb()` builds `ui.tags.nav({"aria-label": "breadcrumb"}, ui.tags.ol({"class": "breadcrumb"}, *items))` from `crumb_labels`; every ancestor item is `ui.input_action_link(session.ns("crumb_<depth>"), label)`; the last item is `li.breadcrumb-item.active` with `aria-current="page"`. Clicking `crumb_0` calls `ui.update_accordion("families", show=False)`; clicking `crumb_1` calls `ui.update_radio_buttons(f"gen_{i}", selected="")`. Both update functions are referenced from the accordion page but their signatures were not fetched: verify on install and record in the report. The M0 `@__import__("shiny").render.ui` hack is removed. |
| `src/progeny_selector/app/screens/rank.py` | `rows()` becomes `filter_rows(result.rows, crumb["family"], crumb["generation"])` then the `only_pass` switch. `DISPLAY_COLUMNS` unchanged; the frame's columns are `DISPLAY_COLUMNS + status_columns(keys) + [k for k in keys if k.startswith("recomb_")]`. `render.DataGrid(frame, selection_mode="rows", filters=True, height="70vh", styles=status_cell_styles(list(frame.columns), data))`. `_publish_selection` uses `table.data_view(selected=True)["sample_id"].tolist()` (resolution 3). A caption under the header reads `f"{len(rows)} individuals shown ({crumb text}); select rows, then open Compare"` via `@render.text def caption`. |
| `src/progeny_selector/app/app.py` | Header docstring status line updated (no longer "M0 scaffold"). |
| `tests/e2e/test_navigate_rank.py` (new) | Below. |
| `CHANGELOG.md` | "Navigate screen: family and generation tree with breadcrumbs that filter the Rank grid. Rank grid shows pass, fail and unknown chips and recombinant flags per locus." |

Acceptance as tests (`tests/e2e/test_navigate_rank.py`, after the phase 1 load sequence in a helper `load_fixture(page)`):

- `controller.PageNavbar(page, "screen").set("navigate")`; `controller.Accordion(page, "navigate-families").expect_panels(["fam_0", "fam_1"])`; `.set("fam_1")`; the breadcrumb `nav[aria-label=breadcrumb]` contains text `F2`; `PageNavbar.set("rank")`; with the switch `rank-only_pass` on, `controller.OutputDataFrame(page, "rank-table").expect_nrow(5)` (fixture: 11 passing, 6 in F1 and 5 in F2, per `expected_results.csv`; confirm before hard-coding); `controller.InputSwitch(page, "rank-only_pass").set(False)`; `expect_nrow(20)`.
- `OutputDataFrame.expect_column_labels` includes `target_T1_status`, `avoid_AV1_status`, `recomb_T1_left`, `recomb_T1_right`; `expect_cell("pass", row=0, col=<index of target_T1_status>)`; the cell's computed `background-color` equals `rgb(0, 158, 115)` (`#009E73`).
- Keyboard: focus the first accordion header, press `Enter`, then `Tab` and `ArrowDown`; `ArrowDown` from `(all generations)` lands on `BC2F1` and the crumb shows `BC2F1`.
- Selection: `OutputDataFrame(page, "rank-table").select_rows([0, 2])` after `set_sort` on `rpp_total` descending; `PageNavbar.set("compare")`; the Compare screen (in this phase still the M0 text panel) lists exactly the two `sample_id`s of view rows 0 and 2. This is the check that `data_view(selected=True)` follows the view.

Not verifiable here: `input.families()` value type when no panel is open (`None` versus `[]`; the effect treats both as no family); `update_accordion` and `update_radio_buttons` signatures; module-namespacing of ids created inside `render.ui` (`session.ns`), which the e2e test will expose.

## Phase 4: Validate/QC table and Compare screen

Goal: QC is a sortable, filterable table with flagged rows highlighted; Compare shows up to six individuals side by side with chips, per-chromosome RPP, drag bounds, recombinant flags and a chromosome strip.

| file | change |
| ---- | ------ |
| `src/progeny_selector/app/screens/qc.py` | `ui_`: `ui.layout_columns(ui.card(ui.card_header("Summary"), ui.output_ui("summary")), ui.card(ui.card_header("Per-individual QC"), ui.output_data_frame("table")), col_widths=(4, 8))`. `summary` renders a `ui.tags.dl`: markers, samples, progeny, informative markers, uninformative by reason (`uninformative_summary`), map unit (`result.unit`), then `ui.tags.ul` of dataset and result warnings verbatim. `table`: `render.DataGrid(pd.DataFrame(qc_table_rows(result.qc, dataset)), filters=True, selection_mode="none", height="70vh", styles=qc_row_styles(rows))`. |
| `src/progeny_selector/app/screens/compare.py` | `ui_`: `ui.card(ui.card_header("Compare selected individuals"), ui.output_ui("legend"), ui.output_ui("panels"))`. `legend`: six `ui.span` badges A/H/B/X/N/U with `STATE_COLORS` backgrounds and text labels, plus target and avoid tick glyphs. `panels`: for the first six of `state.selected_ids()` (a `ui.p` reads "showing 6 of N" beyond that), `ui.layout_columns(*cards, col_widths=[12 // n])`; each card: header `sample_id` (`line_name`), rank and score line, chips `ui.span(label, style=chip_style(status))` per target and avoid locus, a `ui.tags.table` of per-chromosome RPP (`rpp_<chrom>` keys, carrier chromosomes bold with "(carrier)"), a drag table per target (`left_max`, `right_max`, `total_est`, `total_max` with `drag_unit`, recombinant left and right as "yes"/"no"), then the SVG from `strip_rects(chromosome_strips(states[:, j], gm), locus_ticks(...))` where `j = gm.sample_index(sid)` on `result.classification.states` and `gm = state.dataset().genotypes.sorted_by_position()` (`run_analysis` sorted its own copy; the app sorts the same way so marker order matches `classification.states`). SVG: `ui.tags.svg({"viewBox": f"0 0 {g.width} {g.height}", "role": "img", "aria-label": f"chromosome strip for {sid}"}, ui.tags.title(...), ...)` with rows of `ui.tags.rect` each carrying a `ui.tags.title(f"{chrom} {state_label}")` child, `ui.tags.text` chromosome labels, and tick `polygon`s. |
| `tests/e2e/test_qc_compare.py` (new) | Below. |
| `CHANGELOG.md` | "Validate screen: QC table with flagged rows highlighted and the reasons markers were uninformative. Compare screen: side-by-side statuses, per-chromosome RPP, drag bounds and Okabe-Ito chromosome strips." |

Acceptance as tests (`tests/e2e/test_qc_compare.py`):

- After load, `PageNavbar.set("qc")`; `OutputDataFrame(page, "qc-table").expect_nrow(40)`; `expect_column_labels` contains `flags` and `excluded`; the row for `BC2F1-F2-002` has `expect_cell("True", row=r, col=<excluded>)` and its first cell's computed background is not white; the summary `dl` contains "475" and "parents identical (monomorphic): 20".
- `PageNavbar.set("rank")`, `OutputDataFrame(page, "rank-table").select_rows([0, 1])`, `PageNavbar.set("compare")`: two `.card` elements inside `#compare-panels`; the first card contains `BC2F1-F1-001`, a badge with text `pass` and computed background `rgb(0, 158, 115)`, a table cell `Gm06 (carrier)`, and an `svg[role=img]` with 20 `text` chromosome labels and at least one `rect` filled `#0072B2`; the number of `rect` elements in the Gm06 row is at least 3 (A, donor segment, A).
- Selecting seven rows gives six cards and the text "showing 6 of 7".

Not verifiable here: rendered legibility of an 8 px row on the maintainer's laptop screen; whether Bootstrap's card padding fits six columns at 1280 px, which the implementer checks with `page.set_viewport_size` and reports.

## Phase 5: criteria editor and download; download handlers migrated

Goal: criteria are editable in the browser, re-run without reloading data, and downloadable as the same strict YAML the CLI reads.

| file | change |
| ---- | ------ |
| `src/progeny_selector/app/screens/load.py` | Layout becomes two rows: the M0 `layout_columns` (files, status), then `ui.card(ui.card_header("Criteria (editable)"), ui.input_text_area("criteria_text", "criteria.yaml", rows=18, placeholder="Load a criteria.yaml above or paste one here", spellcheck="false"), ui.input_action_button("apply", "Apply criteria and re-analyse"), ui.download_button("download_criteria", "Download criteria.yaml"), ui.output_text_verbatim("criteria_status"))`. After a successful Load **and** after a successful Apply: `ui.update_text_area("criteria_text", value=dump_criteria_yaml(state.criteria()))`, signature `ui.update_text_area(id, *, label=None, value=None, placeholder=None, session=None)` [web] https://shiny.posit.co/py/api/core/ui.update_text_area.html. `apply` effect: `criteria = read_criteria_text(input.criteria_text())`, `result = run_analysis(state.dataset(), criteria)`, set `state.criteria`, `state.result`, clear `state.selected_ids`, rewrite the text area in canonical form, status "re-analysed: N pass hard filters"; `CriteriaError` and `DataContractError` shown verbatim; with no dataset, "Load data first." `@render.download_button(filename="criteria.yaml", media_type="application/yaml") def download_criteria(): yield dump_criteria_yaml(state.criteria())`. Download always serialises the applied `Criteria`, never the raw editor text, so the file a breeder keeps is what the app ran; unapplied edits are not downloaded. |
| `docs/adr/0005-criteria-file-yaml.md` | Dated amendment appended (text below). |
| `src/progeny_selector/app/screens/export.py` | The three handlers become `@render.download_button(filename=...)` and yield text. Thin wrappers `results_csv_text(rows) -> str`, `selection_csv_text(rows, notes) -> str`, `next_round_manifest_text(rows, next_generation, rp, donor) -> str` are added in `io/export.py` by refactoring each `write_*` into `_write_*(fh, ...)` plus the path-opening wrapper, so `write_*` behaviour and `tests/test_smoke_pipeline.py` are untouched. `tempfile` import removed. |
| `docs/data-formats.md` | criteria.yaml section gains one sentence: "The UI's Download criteria.yaml writes this schema with every key explicit (defaults included) and regions as `chrom`/`start_bp`/`end_bp`; the file reloads unchanged." Outputs unchanged. |
| `tests/test_export_text.py` (new) | `results_csv_text(result.rows)` equals the text `write_results_csv` writes to `tmp_path` (decoded, `\r\n` normalised, since `csv.writer` uses `\r\n` with `newline=""`); same for the other two. |
| `tests/e2e/test_criteria_editor.py` (new) | Below. |
| `CHANGELOG.md` | "Criteria can be edited on the Load screen, re-applied without reloading genotypes, and downloaded as criteria.yaml. Downloads no longer write temporary files." |

ADR 0005 amendment, dated 2026-09-12: "The in-memory `Criteria` object is the single source of truth for a run; any editor, whether the M1 YAML text area or a later field-by-field form, is a view that reads from it and writes back to it. Serialisation goes through one path, `criteria_to_dict` and `dump_criteria_yaml` in `progeny_selector.io.criteria`, the inverse of `criteria_from_dict`, and `criteria_from_dict(criteria_to_dict(c)) == c` is tested. The emitted form is canonical: every scalar key explicit with defaults included, loci written with only the keys of their kind, and the `region:` shorthand accepted on input and normalised to `chrom`/`start_bp`/`end_bp` on output. Download always serialises the applied `Criteria`, never the raw editor text. Saved presets are named criteria.yaml files in the project folder, not an in-app store."

Acceptance as tests (`tests/e2e/test_criteria_editor.py`):

- After load, `controller.InputTextArea(page, "load-criteria_text").expect_value(re.compile(r"^name: synthetic BC2F1 fixture"))`; replace `max_missing_rate: 0.2` with `max_missing_rate: 0.5` through `.set(text)`; `InputActionButton(page, "load-apply").click()`; `#load-criteria_status` contains "12 pass hard filters" (BC2F1-F1-003 at 30 % missing now passes if its target and avoid statuses pass and no QC flag excludes it; the implementer computes the number with the CLI first and uses that); `PageNavbar.set("rank")`, `expect_nrow` of the same number.
- Set the text to `targets: []` and Apply: status contains "criteria must define at least one target locus"; the grid still shows the previous rows (a failed Apply leaves state untouched), and a download at this point still yields the last applied criteria.
- Paste criteria using the `region:` shorthand and Apply: the text area now shows `chrom`/`start_bp`/`end_bp`.
- Download: `with page.expect_download() as dl: controller.DownloadButton(page, "load-download_criteria").click()`; `dl.value.suggested_filename == "criteria.yaml"`; `read_criteria_text(Path(dl.value.path()).read_text())` equals the fixture criteria with `filters.max_missing_rate == 0.5`. `DownloadButton` exposes only `click`, so Playwright's own `expect_download` is used [web] https://shiny.posit.co/py/api/testing/playwright.controller.DownloadButton.html.
- Export screen: `expect_download` on `export-results`, `export-selected` (after selecting two rows on Rank) and `export-manifest`; the results file's header starts `rank_overall,rank_in_family,sample_id`; the manifest header equals the samples.csv contract line.

Not verifiable here: `update_text_area` and `InputTextArea.set` signatures; whether `page.expect_download` works under headless Chromium with Shiny's download anchor (it is the documented Playwright path).

## Phase 6: Shinylive export staged, smoke-tested in CI, deployed

Goal: the static site on GitHub Pages runs the whole pipeline in the browser tab; CI proves it with the fixture and proves no request leaves the page's origin.

| file | change |
| ---- | ------ |
| `scripts/build_shinylive.py` (new) | `python3 scripts/build_shinylive.py [--staging build/shinylive-app] [--out site]`: deletes and recreates the staging dir; writes `app.py` containing `from progeny_selector.app.app import app  # noqa: F401` and `requirements.txt` containing `pyyaml`, `pandas`, `numpy` (explicit, so the scanner's handling of inline imports does not matter); `shutil.copytree("src/progeny_selector", staging/"progeny_selector", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))`; runs `[sys.executable, "-m", "shinylive", "export", staging, out]` (whether `python -m shinylive` works is unverified; the `shinylive` console script is the documented form [web] https://github.com/posit-dev/py-shinylive and is the fallback); asserts `out/app.json` exists and contains `"name": "progeny_selector/core/pipeline.py"`; prints the size of `out/`. `build/` is already gitignored. |
| `src/progeny_selector/app/requirements.txt` | Comment updated to say the staging script supersedes it; contents unchanged (`shiny run` does not use the file). |
| `tests/e2e/test_shinylive_export.py` (new) | Skipped unless `PS_SITE_DIR` is set. Starts `python -m http.server --directory $PS_SITE_DIR --bind 127.0.0.1 8008` as a subprocess (the documented preview command [web] https://shiny.posit.co/py/get-started/shinylive.html); uses the pytest-playwright `page` fixture directly; records every request URL through `page.on("request", ...)`; `page.goto("http://127.0.0.1:8008/")`; waits for `#load-run` attached with a 240 s timeout (Pyodide plus numpy and pandas is roughly 13 MB + 7.5 MB + 13 MB per the same page); `page.set_input_files("#load-genotypes", ...)` for the four fixture files; clicks `#load-run`; expects `#load-status` to contain the phase 1 status text within 180 s; asserts every recorded URL starts with `http://127.0.0.1:8008/` (a PyPI fetch by micropip would fail this and would mean resolution 5's primary path did not work). Then a short version of phases 3–5: select two rows on Rank, open Compare, expect two cards; `expect_download` on `export-selected` and check two data rows. |
| `.github/workflows/ci.yml` | New job `export` (needs `check`, Python 3.12): `pip install -e ".[dev,export,e2e]"`; `actions/cache@v4` on the shinylive asset cache keyed `shinylive-${{ hashFiles('pyproject.toml') }}` (the Linux cache path is what `shinylive assets info` prints; the step runs it first and the implementer pins the reported path); `python scripts/build_shinylive.py`; `playwright install --with-deps chromium`; `PS_SITE_DIR=site pytest -m e2e tests/e2e/test_shinylive_export.py`; `upload-pages-artifact` moves here with its existing guards. `deploy-pages` now `needs: export`. |
| `docs/adr/0003-…` | Dated amendment at the end: the package is bundled into the export by staging a copy, not installed from a wheel URL; the wheel-URL line is the recorded fallback (resolution 5). |
| `docs/adr/0009-shinylive-export-staging.md` | MADR: staging script, explicit requirements, CI smoke test as the definition of "export verified", origin-only request assertion as the definition of "data stays in the tab". |
| `README.md` | "What exists now" rewritten for M1; quickstart gains `python3 scripts/build_shinylive.py` in place of the bare `shinylive export`. |
| `PLAN.md` | M1 verification block (below); "UI walkthrough" sentences for screens 2, 3 and 5 lose "placeholder"; "Testing and CI" says Playwright screen tests arrived in M1. |
| `CLAUDE.md` | State: M1 complete, M2 next; gates gain `python3 scripts/build_shinylive.py && PS_SITE_DIR=site pytest -m e2e tests/e2e/test_shinylive_export.py` per milestone. |
| `CHANGELOG.md` | "The Shinylive static export includes the package and is smoke-tested in CI with the fixture; it deploys to GitHub Pages when ENABLE_PAGES is set." |

M1 verification block for PLAN.md, dated, from the maintainer's laptop: `pytest` count; `pytest -m e2e` count and Chromium version; `build_shinylive.py` output size and export wall-clock; the smoke test's time to `#load-run` and time to status; which of resolution 5's two paths is in use; the Pages URL once ENABLE_PAGES is true; and the result of the Q7 real-data protocol, or the statement that it has not been run.

Acceptance as tests: `test_shinylive_export.py` passes in the CI `export` job; the site artifact deploys when `ENABLE_PAGES == 'true'`; `pytest` (unit) is unaffected.

Not verifiable here: Pyodide `sys.path` for the staged directory (the smoke test is the check); the asset cache path on Linux; CI wall-clock for the Pyodide load (240 s is a hang guard, tightened to 3× the recorded figure afterwards, as the sibling did); real-data behaviour (Q7).

## Gates and CI after M1

Per commit: `ruff check . && ruff format --check .`, `mypy`, `pytest`, `python3 scripts/make_fixture.py && git diff --exit-code -- tests/fixtures`. Per milestone and in CI: `pytest -m e2e` (Chromium installed), `python3 scripts/build_shinylive.py && PS_SITE_DIR=site pytest -m e2e tests/e2e/test_shinylive_export.py`. CI jobs: `check` (3.11, 3.12), `e2e` (3.12), `export` (3.12, uploads the Pages artifact), `deploy-pages`.

## Inconsistencies found (recorded, not all fixed in M1)

1. `pyproject.toml` Repository URL owner is `piercetaylor2020`; the remote and `gh` say `piercetaylor`. Fixed in phase 1.
2. ADR 0003 says `requirements.txt` names the package wheel; `app/requirements.txt` lists only pyyaml and no release exists; the CI step "shinylive export builds" checks only the exit code. Fixed in phase 6 with an ADR amendment.
3. `docs/data-formats.md` accepts `.vcf.bgz`, but `io/__init__.py::load_genotypes` strips only `.gz`, so `x.vcf.bgz` raises "cannot infer genotype format"; it also accepts an undocumented `.bcf` extension and would read binary BCF as text. Not fixed: parser behaviour is frozen for S1 and Q2 of the sibling plan (Q6).
4. `core/qc.py` docstring signature disagrees with the code (phase 2 fixes the docstring only).
5. `rank.py` selection indexes the unfiltered list (resolution 3).
6. `navigate.py` uses `@__import__("shiny").render.ui` (phase 3 removes it).
7. `dev` extras include shinylive, so the per-commit gates depend on a 2018 sdist-only package that failed to build in the scaffold container, though it built on the maintainer's machine (phase 1).
8. PLAN.md "Testing and CI" places Playwright tests of the screens in a later milestone; M1's acceptance needs them now (phase 6 edits the sentence).
9. `scripts/make_fixture.py` writes CRLF on Windows (Q9).

## Questions for the maintainer

1. **Tree control.** SETTLED 2026-09-12 by the maintainer as recommended. Recommendation was: an accordion of families with a radio group of generations per panel (resolution 4). Grounds: py-shiny has no tree component; the accordion gives keyboard behaviour for free; a Playwright `Accordion` controller exists; the M0 select inputs remain a fallback. Alternative: `ui.navset_pill_list` with one pill per family, which looks more like a tree but nests a second navset inside the page navbar.
2. **Criteria editor form.** SETTLED 2026-09-12: the maintainer asked what comparable software does and delegated the choice. YAML text area in M1; a field-by-field form planned for M2. Grounds, fetched 2026-09-12: Flapjack reads loci from a QTL file and takes run options in a small dialog, but applies thresholds after the run without saving them [web] https://flapjack.hutton.ac.uk/en/latest/mabc.html; OptiMAS reads loci and favourable alleles from a map file (Valente et al. 2013, J Hered 104:586) [web] https://academic.oup.com/jhered/article/104/4/586/775686; BMS Molecular Breeding Planner is the one verified tool offering both, parameters "from an external file (.mars, .mabc, or .mas) or set by hand" [web] https://bmspro.io/222. EiB, Intertek, Genovix and Phenome documentation could not be verified. M2 form sketch: scalars through `ui.input_numeric`, `ui.input_select` and `ui.input_switch`; loci through an editable `render.DataGrid` with one row per locus; the form and the text area are two views of the same applied `Criteria`, and both validate through `Criteria.validate()` with CLI-identical errors.
3. **Where presentational helpers live.** SETTLED 2026-09-12 by the maintainer as recommended. Recommendation was: `app/present.py`, Shiny-free and unit-tested, holding only colour, style and geometry mapping; all genetics, filtering and serialisation in `core/` and `io/`. Grounds: "app computes nothing" is about metrics; colour lookups in `core` would make the pure core know about CSS. Alternative: a `progeny_selector/view/` package.
4. **Shinylive packaging.** Recommendation: stage a package copy (resolution 5); fall back to a wheel URL on the Pages site only if Pyodide does not import from the staged directory. Grounds: no release and no PyPI dependency needed, and it works from any directory; the URL fallback ties local preview to the internet and to a deploy that has already happened.
5. **Extras split and e2e policy.** SETTLED 2026-09-12 by the maintainer as recommended. Recommendation was: `dev` without shinylive; separate `export` and `e2e` extras; `pytest` stays browser-free per commit; `pytest -m e2e` per milestone and in CI. Grounds: the per-commit gate should not depend on lzstring's sdist or a browser install; this mirrors the sibling's node/browser split.
6. **`.vcf.bgz` and `.bcf` extensions.** SETTLED 2026-09-12 by the maintainer: fix now, as phase 0. The recommendation had been to wait for the contract. The fix is safe to take early because it moves the code to the documented behaviour and to what isoline-browser already does, rather than changing the contract.
7. **Real-data acceptance.** SETTLED 2026-09-12 by the maintainer: mine the SoyBase backcross-derived NILs genotyped on SoySNP50K, where both parent and progeny were genotyped (as used by PATRIOT, Frontiers in Plant Science 2021, https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2021.676269/full; SNP data https://soybase.org/snps/). Before any download: confirm the SoyBase/USDA data licence (its terms page returned 403 on 2026-09-12) and identify recurrent parent, donor and NIL sets from GRIN pedigrees. The data is downloaded by hand to the gitignored `data/` folder, never committed, and converted with a script that writes `samples.csv` from the pedigrees. Caveats recorded now: NILs are advanced inbred lines, not segregating BC2F1 progeny, so generation labels and expected heterozygosity differ from the BCnF1 QC expectations; positions are on Wm82.a2, while the chromosome-length table is Wm82.a4.v1. This is planned as separate work that runs alongside the M1 phases and does not block them. The original recommendation follows for the record: CI verifies the synthetic fixture through the exported site, including the origin-only request assertion (phase 6); the maintainer runs the same steps once on the program's BC2F1 KASP or SoySNP6K file from the Pages URL or a local `http.server`, with the browser's network panel open, and records in the PLAN.md verification block the marker and progeny counts, the number passing, wall-clock to status, and "no requests after page load left the origin". The file is never committed and never attached to an issue. Grounds: real data cannot enter the repository, and the network observation is the honest form of "without leaving the browser".
8. **Shiny version floor.** Recommendation: keep `shiny>=1.4` and record the versions the laptop and CI actually resolve in the verification block. Grounds: everything used exists at 1.4.0; raising the floor buys nothing verified here.
9. **Fixture line endings on Windows.** SETTLED and done 2026-09-12 in commit c5a7f8b; the repository URL was corrected in 28ef760. Recommendation was: a one-line fix in a separate `fix(fixture)` commit before phase 1, opening output files with `newline="\n"`, so the working copy matches the committed bytes. Grounds: the gate passes only because git normalises; any tool reading the files directly sees different bytes. Alternative: leave it and document it in CLAUDE.md.
