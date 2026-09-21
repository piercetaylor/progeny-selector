# M2 usable: implementation phases

This document decomposes milestone M2 (PLAN.md, "Milestones") plus the M1 gap (a breeder takes real segregating data to selected.csv in the browser) and the deferred reviewer findings (PLAN.md Handoff) into phases, one commit each. Every commit passes the gates of CLAUDE.md. Dates are 2026-09-16. Facts that could not be verified this session are marked "not verifiable here" with the test that settles them. Doers receive one phase's line range and make no design decisions; each phase names its files, functions, exact names and values, edge cases, tests, model, whether the reviewer runs, and its acceptance criterion.

Invariants carried from M1 (docs/m1-phases.md): `core/` pure; `io/` validates and serialises; `app/` computes nothing; colours only from `constants`; `pytest` browser-free; screen interface `view(id)`, `server(id, state)`; no AI attribution; Conventional Commits; CHANGELOG under `[Unreleased]`; decisions as MADR ADRs numbered from 0015; fixture expectations computed by the generator's independent implementation; at most two concurrent implementers on disjoint files; reviewer on core/io/model/contract diffs plus one end-of-milestone pass; doers never run git.

## Questions for the maintainer

Each item is a decision the documents do not settle or an inconsistency between PLAN.md, the ADRs and the code. Nothing below is resolved silently; the recommendation is what the phases assume unless overruled.

1. **`qc_excluded` semantics (Validate vs Rank, deferred medium finding).** Code: `core/qc.py` `qc_table_rows` sets `qc_excluded` from flags alone; `core/score.py` `hard_filters` excludes only when `filters.exclude_qc_flagged`. Recommendation: policy-aware. Add `qc_excluding_hits(flags, filters) -> list[str]` to `core/score.py`, used by both `hard_filters` and `qc_table_rows(qc, dataset, filters)`, so the rule lives once and Validate and Rank agree under `exclude_qc_flagged: false`. Tradeoff: the column then means "would be excluded under the applied criteria", so it changes when criteria change; the alternative keeps the flag-only meaning and renames the column `qc_excluding_flag` (policy-free, but Validate still cannot say whether the plant is excluded). The column is UI-only (not in docs/data-formats.md), so no contract approval is needed either way.
2. **results.csv schema start value and placement.** Recommendation: `results_schema = "1.0.0"`, written on every row; the four new columns (`background_model`, `background_unit`, `rank_mode`, `results_schema`) follow `frac_b` and precede the per-locus columns, so the fixed prefix (34 columns, listed in Phase 4) is the header written when there are no rows. Tradeoff: leading placement would make them visible first in a spreadsheet, but would move every documented leading column. This edits the results.csv section of docs/data-formats.md; no sibling file is affected (outputs are this repo's).
3. **Fifth column `assembly` and the assembly selector.** The maintainer named four columns; the assembly selector (item 7b) produces a fifth candidate. Recommendation: add `assembly` (value of the criteria key, below) now, because the schema freezes at 1.0.0 and adding it later bumps the minor version. Selector location: top-level criteria.yaml key `assembly` (default `Wm82.a4`), because chromosome lengths are a run parameter (RPP weights, drag ends, strips), criteria.yaml is this repo's file, and it round-trips through the criteria download; markers.csv is contract-owned (a change starts in backcross) and a UI-only control would not reach the CLI. Assemblies: `Wm82.a1`, `Wm82.a2`, `Wm82.a4` and `none`; `Wm82.a5` and `Wm82.a6` only if the doer obtains their lengths (Q7). Tradeoff: a criteria key is one more thing to set per run; a2 data run against the a4 default today with only a README caveat. Affects docs/data-formats.md (criteria section), no sibling file.
4. **NA scope.** Decided for results.csv. Recommendation: selected.csv also writes `NA` for a missing `family_id`, `generation`, rank or score, and carries `results_schema` as its last column; next_samples.csv keeps empty cells and no schema column because it is the input contract's samples.csv (`generation` and `family_id` "default to empty" there) and must load unchanged in backcross. Tradeoff: two output files with NA, one without; the alternative (NA everywhere) risks a backcross reader treating "NA" as a family name.
5. **Duplicate detection as a QC flag value.** `core/qc.py` `duplicate_pairs` exists but `run_analysis` never calls it, although PLAN.md algorithm 8 and the qc.py docstring say pairs are reported. Recommendation: call it in the pipeline, keep `AnalysisResult.duplicates` as pairs, and append the advisory (non-excluding) value `possible_duplicate` to `qc_flags` of both members, as ADR 0012 did for `family_donor_outlier` (an added value in an existing column, not a new column). Tradeoff: `qc_flags` content changes for datasets with duplicates; the alternative surfaces pairs in the UI only, leaving results.csv silent.
6. **criteria.yaml download: unset optional keys.** docs/data-formats.md says the download writes "every key explicit"; `criteria_to_dict` omits `name`, `notes`, `max_marker_coverage`, `flank_left`, `flank_right` when None (deferred low finding). Recommendation: write them as `null`, and change the doc sentence to "every key explicit; an optional key with no value is written as `null`". Tradeoff: `name: null` in a downloaded file reads oddly; the alternative rewrites the doc to "every key with a value". A one-sentence edit to the criteria section of docs/data-formats.md.
7. **SoyBase position table (item 7a): storage and URLs.** Recommendation: the script reads GFF3 files the maintainer downloads by hand into gitignored `data/soybase/` (an optional `--download` fetches from a URL table inside the script, never run in CI); the derived table (about 50K rows, several MB) is written to `data/soybase/` and not committed. Chromosome lengths for a1/a5/a6 are taken from the same SoyBase genome_info page as the existing constants (fetch date recorded) or from `##sequence-region` pragmas if the GFF3s carry them; an assembly whose lengths cannot be obtained is left out of the `assembly` enum. Not verifiable here: the Data Store file names and URLs (`Wm82.gnm{1,2,4,5,6}.mrk.SoySNP50K`, `.SoySNP6K`); the doer resolves them from the directory listing and records them in docs/reference-repos.md. Tradeoff: committing the table would let CI test on it but adds megabytes of third-party data to the repo.
8. **Song et al. 2016 Table S1 (item 7c): file type and smoothing.** Not verifiable here: whether Table S1 is xlsx and its column headers. Recommendation: the maintainer exports it to CSV by hand into gitignored `data/song2016/` (no openpyxl dependency); the script takes a `--columns` mapping with defaults the doer fills after inspecting the file. Monotone Marey interpolation: cM made non-decreasing along bp by a running maximum before linear interpolation, clamped to the outermost mapped cM per chromosome. Tradeoff: running maximum is the simplest monotone smoother and flattens local inversions to plateaus; pool-adjacent-violators isotonic regression spreads them but is more code without numpy support.
9. **KASP export layout (item 7d).** Assumed: LGC Kraken/SNPviewer long-format genotyping CSV, one row per sample x SNP with at least the columns `SubjectID`, `SNPID`, `Call`, calls written `A:A`, `A:G`, `G:G`, and `Uncallable`, `Missing`, `?`, `Bad`, `Dupe`, `NTC`, empty as missing; NTC/control rows dropped by `SubjectID`. Tradeoff: the grid export (SNPs as columns) is the other common shape; the script accepts only the long form in M2 and says so. Please confirm the program's export.
10. **Question B: two-generation dataset.** No public linked consecutive-generation dataset exists; the maintainer is sourcing one. Phase 8 verifies the round trip on a synthetic two-generation fixture first; the real-data step (Phase 11) is blocked on a dataset meeting: both parents genotyped; generation n progeny with samples.csv; generation n+1 progeny genotyped on the same panel whose samples.csv is this tool's next_samples.csv written from the generation n selection, unchanged; positions on one assembly. Recommendation: accept a two-generation subset of the program's own KASP or SoySNP6K data, converted by the Phase 10 KASP script, run locally from `data/`, never committed.
11. **PLAN M2 "weighted model the default" is already the code's default** (`model/criteria.py` `BackgroundOptions.model = "weighted"`, docs/data-formats.md). The fixture criteria pins `model: count`, and every e2e constant (11 pass, row orders, `15.28 cm`) depends on it. Recommendation: keep `count` in the fixture criteria and test weighted through `dataclasses.replace` against weighted columns the generator computes independently (Phase 6). Tradeoff: the fixture then does not exercise the default end to end in the browser; the real-data protocol (Phase 11) does.
12. **PLAN M2 "per-target windows" is already implemented** (`TargetSpec.flank_left/flank_right`; `core/pipeline.py` lines 118-119). Only a hand-built unit test and a docs sentence remain (Phase 6). Confirm nothing else was meant.
13. **`n_per_selected` exposure for the round trip.** `write_next_round_manifest` already accepts `n_per_selected` (default 1); the CLI and the Export screen do not expose it. Recommendation: `select --per-selected N` and an Export numeric input, so next_samples.csv can be used unchanged when one row per planted plant is requested; docs/data-formats.md next_samples.csv changes "Edit ids after planting" to say so. No column change; affects no sibling file.
14. **Contract 1.3.0 precondition is global.** CHANGELOG.md and docs/data-formats.md are in the uncommitted 1.3.0 set and nearly every phase touches CHANGELOG.md. Recommendation: no M2 phase is dispatched before that work is committed (Phase 0). Alternative: the main session stages by name, which risks committing a half-finished CHANGELOG entry.
15. **Criteria editor form** (docs/m1-phases.md Q2 said "planned for M2") is not in the maintainer's M2 scope list. Recommendation: out of M2; the YAML editor covers `ranking`, `assembly` and per-target windows.
16. **DataGrid keyboard selection** (PLAN.md "UI walkthrough": arrow keys and space/enter). Not verifiable here. Recommendation: the walkthrough documents observed behaviour, and the e2e test asserts what the walkthrough says (Phase 7).
17. **Stale text to correct at end of milestone**: `cli.py` docstring "downstream Shiny dashboards"; PLAN.md names `scripts/read_results.R` before it exists; `app/app.py` docstring "M1 in progress ... placeholders"; `select.py` docstring "Placeholder in M0"; `core/chrom.py` docstring cites contract 1.1.0 and hardcodes a4 while `constants.py` carries an unused a2 table; PLAN.md "Repository layout" lists ADRs 0001..0007; CLAUDE.md says `../backcross` while the folder is `../isoline-browser` (Handoff). All docs-only, folded into Phases 5, 7, 11.

## Decisions, 2026-09-16

The maintainer delegated questions 1-9 and 11-17 to research ("ask fable"), with the criterion: whatever serves academic and open-source plant-breeding users best (reproducible with standard tools, tolerant of spreadsheet, pandas and R exports, never silently wrong). The main session checked the key claims against the code and the URLs below. These decisions override the recommendations above where they differ; ADRs record who decided and why. Question 10 (B) stays open: the maintainer is sourcing the dataset.

1. `qc_excluded` is policy-aware, as recommended.
2. `results_schema = "1.0.0"` on every row; the five metadata columns follow `frac_b`. Contract 1.4.0 (docs/adr/0014, in progress in another session) adds `token_profile` as the last results column, so `results_schema` precedes it and the empty-file header is the 34-name prefix plus `token_profile` (35 names). `FIXED_COLUMNS`, the schema test and the R `col_types` (`token_profile = col_character()`) include it. Phases 3 and 4 wait for contract 1.4.0 to be committed.
3. `assembly` is added now, as a criteria key, default `Wm82.a4`, values `Wm82.a1, Wm82.a2, Wm82.a4, Wm82.a5, Wm82.a6, none`. a1 lengths from https://www.soybase.org/resources/genome_info/ : Gm01 55,915,596; Gm02 51,656,714; Gm03 47,781,077; Gm04 49,243,853; Gm05 41,936,505; Gm06 50,722,822; Gm07 44,683,158; Gm08 46,995,533; Gm09 46,843,751; Gm10 50,969,636; Gm11 39,172,791; Gm12 40,113,141; Gm13 44,408,972; Gm14 49,711,205; Gm15 50,939,161; Gm16 37,397,386; Gm17 41,906,775; Gm18 62,308,141; Gm19 50,589,442; Gm20 46,773,168. a5 and a6 were computed from the SoyBase Data Store `genome_main.fna.gz` files (`Wm82.gnm5.NRKG`, `Wm82.gnm6.S97D`); no published table exists, so `constants.py` cites those files and the ADR says the lengths are computed. a5 Gm01..Gm20: 58,848,252; 54,986,233; 49,126,026; 55,565,417; 44,498,915; 53,747,728; 47,217,708; 50,302,612; 51,410,279; 55,407,096; 41,206,322; 43,862,403; 46,254,787; 52,806,141; 55,385,218; 40,798,560; 44,045,024; 61,540,830; 53,925,450; 51,073,742. a6 Gm01..Gm20: 59,600,650; 54,790,513; 49,119,424; 55,533,332; 44,449,196; 53,606,008; 47,205,696; 50,288,313; 51,135,784; 55,391,814; 42,152,160; 43,682,070; 46,106,413; 52,776,670; 55,366,112; 40,198,974; 44,073,542; 60,749,638; 53,875,614; 51,026,854. The doer recomputes one chromosome per assembly from the FASTA as a spot check before committing the tables.
4. `NA` in results.csv and selected.csv; next_samples.csv keeps empty cells, because backcross's manifest reader (`src/io/manifest.ts`) would read `NA` as a family named "NA". readr (`na = c("", "NA")`) and pandas read both forms as missing.
5. `possible_duplicate`, advisory, threshold 0.995, on both members; IBS is computed over informative markers called in both (subsampled among those), not all markers, because with a small informative fraction backcross siblings reach 0.995 over all markers. The Phase 6 generator uses the same marker set. The ADR states that the threshold is stricter than PLINK/KING practice on purpose, since BC5+ and selfed siblings legitimately approach it.
6. Unset optional keys are written as `null`.
7. As recommended; `--download` also verifies the directory's `CHECKSUM.*.md5` and sends a browser User-Agent. URLs: `https://data.soybase.org/Glycine/max/markers/<dir>/glyma.<dir>.gff3.gz` with `<dir>` in `Wm82.gnm{1,2,4,5,6}.mrk.SoySNP50K` and `Wm82.gnm{1,2,4,5,6}.mrk.SoySNP6K` (licence Open in each `README.*.yml`). No `##sequence-region` pragmas.
8. Table S1 is `.xls` (https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12864-015-2344-0/MediaObjects/12864_2015_2344_MOESM1_ESM.xls, CC BY 4.0), exported to CSV by hand. Header on row 2 after a title row: `ss ID` (join key `"ss" + int`), `SNP ID`, `Glyma1.01 Chromosome`, `Glyma1.01 Coordinate`, `Wm82.a2.v1 Chromosome` (`Chr01`; scaffold cells may carry a leading newline), `Wm82.a2.v1 Coordinate`, gene columns, `WP Linkage Group`, `WP linkage position`, `EW Linkage Group ` (trailing space), `EW linkage position`. The two maps are never averaged: `--map WP|EW`, default `WP`. Cleaning replaces the running maximum: drop rows whose linkage group is not the chromosome, then drop the markers outside the longest non-decreasing subsequence of cM along bp, then `np.interp` with clamping; the log counts each step. This follows MareyMap's practice of invalidating outlier markers rather than smoothing them. `core/genetic_map.py` gains `longest_nondecreasing_mask(cm_sorted_by_bp) -> np.ndarray[bool]` in place of `monotone_cm`, and its test uses the case `cm [0, 5, 4, 10]` → mask `[T, T, F, T]`.
9. The KASP converter accepts both shapes: long, detected by the headers `SubjectID`, `SNPID`, `Call` (unverified against a primary source), and grid (SNPviewer's export), detected by matching row or column ids against `--markers` in either orientation. Calls match `^[ACGT][:\-.]?[ACGT]$`; `?` and the listed tokens are missing.
11-17. As recommended.

Spec corrections applied by the decisions: ADRs are numbered 0015-0019 (0014 is contract 1.4.0). Phase 10 takes `marker_id` from `Name=` (the `ID=` is `glyma.Wm82.gnmN.ss...`); allele attributes are `alleles=` (gnm1/2 50K and all 6K) or `ref_allele=` (gnm4/5/6 50K); gnm5 seqids are `Chr01`-`Chr20`. Phase 1: `min_markers: 2.0` and whole-float locus positions are accepted by skipping `_check_int` for whole floats at the locus call sites and converting in `_normalise_locus`. Phase 3's `chrom_length_bp` docstring cites contract 1.4.0.

## Order, dependencies and parallelism

| phase | scope | model | reviewer | depends on |
| ----- | ----- | ----- | -------- | ---------- |
| 0 | Precondition: contract 1.3.0 working tree committed | none | no | nothing |
| 1 | Deferred findings in core/io: `qc_excluded`, criteria hardening | opus | yes | 0 |
| 2 | Deferred findings in UI/CI/tests: Load accept and catch, CI cache and extras, smoke-test requests, no-family test, Compare drag assertion | sonnet | no | 0 |
| 3 | criteria `ranking.mode` (staged ranking), `assembly` key and chromosome-length tables, map warnings | opus | yes | 1 |
| 4 | results.csv freeze: NA, header on empty, metadata columns, `results_schema`, selected.csv | opus | yes | 3 |
| 5 | `scripts/read_results.R`, CI `r-reader` job | sonnet | no | 4 |
| 6 | Weighted-model and staged expectations in the generator, duplicate detection in the pipeline, per-target window test | opus | yes | 3 |
| 7 | UI: Selection notes, Validate map/model block and duplicates, Export per-selected, keyboard walkthrough | sonnet | no | 4, 6 |
| 8 | Two-generation synthetic round trip: generator extension, `--per-selected`, round-trip test | opus | no (scripts/tests/cli) | 4, 6 |
| 9 | `core/genetic_map.py` and `scripts/song2016_map.py` | opus | yes | 3 |
| 10 | `scripts/soysnp_positions.py` and `scripts/kasp_to_wide.py` | sonnet | no | 3 |
| 11 | Real-data acceptance (browser, local and Pages), Question B real-data step (blocked), end-of-milestone review, docs | main session and maintainer | end-of-milestone pass | all |

Parallel pairs on disjoint files: 1 with 2; 5 with 6; 7 with 8; 9 with 10. CHANGELOG.md is the collision point of every pair: the second phase of a pair puts its CHANGELOG line in its report and the main session inserts it. Phase 3 and Phase 4 both touch docs/data-formats.md and run sequentially.

## Phase 0: precondition

No work. The uncommitted contract 1.3.0 mirror (io/delimited.py, hapmap.py, manifest.py, vcf.py, wide_csv.py, tests/test_io.py, tests/test_contract_cases.py, CHANGELOG.md, docs/data-formats.md, contract/, docs/adr/0013) is committed by its own session before any M2 phase is dispatched (Q14). No M2 phase edits `contract/`. New criteria-hardening tests go in a new file, never in tests/test_io.py.

## Phase 1: deferred findings in core and io

Goal: Validate and Rank agree on QC exclusion under every `exclude_qc_flagged` value; criteria.yaml cannot exhaust memory or smuggle NaN, int ids or fractional positions.

| file | change |
| ---- | ------ |
| `src/progeny_selector/core/score.py` | Add `def qc_excluding_hits(flags: list[str], filters: Filters) -> list[str]`: returns `[f for f in flags if f in QC_EXCLUDING_FLAGS]` when `filters.exclude_qc_flagged` is true, else `[]`. `hard_filters` uses it (behaviour unchanged: reason `"qc:" + "|".join(hits)` when non-empty). |
| `src/progeny_selector/core/qc.py` | `qc_table_rows(qc: list[SampleQC], dataset: Dataset, filters: Filters) -> list[dict]`; `qc_excluded = bool(qc_excluding_hits(q.flags, filters))`. Docstring: "true when the applied filters exclude the individual on QC flags alone". Keys and order unchanged. |
| `src/progeny_selector/app/screens/qc.py` | `qc_table_rows(result.qc, dataset, state.criteria().filters)`; when `state.criteria()` is None, return None from `table`. |
| `src/progeny_selector/io/criteria.py` | (a) `class _NoAliasLoader(yaml.SafeLoader)` overriding `compose_node(self, parent, index)`: if `self.check_event(yaml.AliasEvent)` raise `CriteriaError("criteria.yaml: YAML anchors and aliases are not accepted")`, else `super()`. `read_criteria_text` uses `yaml.load(text, Loader=_NoAliasLoader)` and first raises `CriteriaError("criteria.yaml is larger than 1 MB")` when `len(text.encode("utf-8")) > 1_048_576`. (b) `_check_number` additionally rejects non-finite values: `math.isfinite(value)` false -> `CriteriaError(f"{where}{key} must be a finite number, got {value!r}")`. (c) `_normalise_locus` coerces `marker_id`, `left_marker`, `right_marker`, `chrom`, `locus_id` to `str` when they are int or float (`str(value)`); a float `start_bp`, `end_bp`, `anchor_bp` or `min_markers` that is not whole (`value != int(value)`) raises `CriteriaError(f"{where}.{key} must be an integer, got {value!r}")`; whole floats become int. |
| `tests/test_qc_table.py` | `qc_table_rows(result.qc, dataset, criteria.filters)`; new `test_qc_excluded_follows_filters`: with `dataclasses.replace(criteria, filters=dataclasses.replace(criteria.filters, exclude_qc_flagged=False))`, `run_analysis` gives `BC2F1-F2-002` `passes_filters is True` and `exclusion_reason == ""`, and `qc_table_rows(..., filters)` gives `qc_excluded is False` for it; with the fixture filters both say excluded. |
| `tests/test_criteria_hardening.py` (new) | `read_criteria_text("a: &x [1]\nb: *x")` raises `CriteriaError` matching "anchors and aliases"; a document of 1.1 MB raises matching "larger than 1 MB"; `weights: {rpp_noncarrier: .nan}` and `.inf` raise matching "finite"; `targets: [{locus_id: 1, marker_id: 12345}]` yields `TargetSpec.marker_id == "12345"` and `locus_id == "1"`; `start_bp: 100.5` raises matching "must be an integer"; `start_bp: 100.0` gives `100`. |
| `CHANGELOG.md` | Fixed: "The Validate screen's `qc_excluded` follows `filters.exclude_qc_flagged`, so it agrees with the Rank screen. criteria.yaml rejects YAML anchors and aliases, documents over 1 MB, NaN or infinite numbers and fractional positions; numeric marker ids are read as text." |

Edge cases: an empty `flags` list; `exclude_qc_flagged` true with only advisory flags (`qc_excluded` False); `min_markers: 2.0` accepted as 2; `min_markers: true` still rejected by `_check_int`.

Tests that must pass: the two files above, `pytest`, gates. Reviewer: yes (core, io). Acceptance: `test_qc_excluded_follows_filters` passes and `tests/test_criteria_hardening.py` passes with every case above.

## Phase 2: deferred findings in UI, CI and tests

Goal: the remaining low findings and test gaps close without touching core, io or model.

| file | change |
| ---- | ------ |
| `src/progeny_selector/app/screens/load.py` | `accept=[".vcf", ".gz", ".bgz", ".txt", ".csv", ".tsv", ".hmp", ".hapmap"]` on `genotypes`; `samples` and `markers` accept `[".csv", ".tsv", ".txt"]`. `_run` and `_apply` catch `(DataContractError, CriteriaError)` as now and add a second `except Exception as exc:` that sets the message to `f"unexpected error: {type(exc).__name__}: {exc}"` and calls `logging.getLogger(__name__).exception("load failed")`; AppState untouched on either path. |
| `pyproject.toml` | `export = ["shinylive>=0.8"]` (drop `build`). |
| `.github/workflows/ci.yml` | `export` job: a step `id: shinylive_version` runs `python -c "import importlib.metadata as m; print(m.version('shinylive'))"` into `$GITHUB_OUTPUT` as `v`; cache key `shinylive-${{ steps.shinylive_version.outputs.v }}-${{ hashFiles('pyproject.toml') }}`. |
| `tests/e2e/test_shinylive_export.py` | Record requests on the browser context: `page.context.on("request", ...)` instead of `page.on`, so service-worker fetches are seen; do not block service workers (Shinylive's app iframe depends on one). Print the count of recorded URLs. |
| `tests/e2e/test_qc_compare.py` | In `test_compare_two_cards`, after the first card's assertions, sort the Rank grid by `rpp_total` descending before selecting so view row 2 is `BC2F1-F2-015` (as `test_selection_follows_sorted_view` establishes), open Compare and assert the second card's `table[data-locus=T1]` row "total max" contains the value from `expected_results.csv` for `BC2F1-F2-015` formatted `f"{v:.2f} cm"` (Handoff: `drag_total_max_cm = 109.534`, so `109.53 cm`; the doer confirms against the file before hard-coding). |
| `tests/e2e/test_navigate_rank.py` | `test_no_family_node`: copy the fixture to `tmp_path`, blank `family_id` for `BC2F1-F1-004` in samples.csv, load through the Load screen, open Navigate, expect panels `["fam_0", "fam_1", "fam_2"]`, set `fam_2`, breadcrumb current text `(no family)`, Rank with `only_pass` off shows 1 row with `sample_id` `BC2F1-F1-004`. |
| `CHANGELOG.md` (via report) | Fixed: "Load accepts `.bgz`, `.tsv` and `.hapmap` uploads and reports unexpected errors instead of a blank status. The Shinylive asset cache key includes the shinylive version." |

Not verifiable here: whether `page.context.on("request")` sees the service worker's fetches in this Chromium; the test prints the request count so the report can say. Reviewer: no. Acceptance: `pytest -m e2e` passes including the two new tests; the `export` job is green with the new cache key.

## Phase 3: ranking mode, assembly selector, map warnings

Goal: criteria.yaml carries `ranking: {mode: weighted|staged}` (ADR 0007 amendment) and `assembly`; chromosome lengths come from a keyed table; every fallback from cM to bp and every missing or exceeded chromosome length is a warning the UI already displays.

Precondition: Q3 answered (which assemblies; `assembly` key approved as a documented parameter change to docs/data-formats.md, criteria section; no sibling file affected).

```python
# model/criteria.py
RANK_MODES: tuple[str, ...] = ("weighted", "staged")
ASSEMBLIES: tuple[str, ...] = ("Wm82.a1", "Wm82.a2", "Wm82.a4", "none")   # plus "Wm82.a5", "Wm82.a6" when their lengths are in constants
DEFAULT_ASSEMBLY = "Wm82.a4"
@dataclass
class RankingOptions:
    mode: str = "weighted"
    def validate(self) -> None: ...   # CriteriaError(f"ranking.mode must be one of {RANK_MODES}")
class Criteria:  # new fields, defaults
    ranking: RankingOptions = field(default_factory=RankingOptions)
    assembly: str = DEFAULT_ASSEMBLY      # validate: CriteriaError(f"assembly must be one of {ASSEMBLIES}")

# constants.py
SOYBEAN_CHROM_LENGTHS_BP: dict[str, dict[str, int]] = {"Wm82.a4": SOYBEAN_CHROM_LENGTHS_BP_WM82A4, "Wm82.a2": SOYBEAN_CHROM_LENGTHS_BP_WM82A2, "Wm82.a1": {...}}
# keep SOYBEAN_CHROM_LENGTHS_BP_WM82A4 and _WM82A2 by name (scripts/make_fixture.py imports the a4 table); each table's comment records source URL and fetch date

# core/chrom.py
def chrom_length_bp(name: str, fallback: int | None = None, assembly: str = DEFAULT_ASSEMBLY) -> int | None
    # "none" or an unknown assembly -> fallback; docstring cites contract 1.3.0

# core/score.py
def rank_rows_staged(rec_count, rpp_carrier, rpp_noncarrier, drag_est, missing_rate, sample_ids, passes, family_ids) -> (rank_overall, rank_in_family)
    # order: rec_count desc, rpp_carrier desc, rpp_noncarrier desc, drag_est asc, missing_rate asc, sample_id asc;
    # NaN sorts last on every key (-inf for desc keys, +inf for asc); dense ranks among passing rows, NaN for excluded

# core/strip.py
def chromosome_strips(states, gm, assembly: str = DEFAULT_ASSEMBLY) -> list[ChromStrip]   # length_bp from chrom_length_bp(chrom, None, assembly)
```

| file | change |
| ---- | ------ |
| `src/progeny_selector/io/criteria.py` | `TOP_KEYS` gains `"ranking"`, `"assembly"`; `_check_shapes` covers `ranking`; `assembly` must be a string (`CriteriaError(f"assembly must be text, got {value!r}")`); `criteria_from_dict` builds `RankingOptions` via `_build` and reads `assembly=str(doc.get("assembly", DEFAULT_ASSEMBLY))`; `criteria_to_dict` order: `name, targets, avoid, flank_window, flank_unit, assembly, background, ranking, weights, filters`; `ranking` written as `{"mode": ...}`; unset optional keys written as `null` (Q6): `name`, `notes`, `max_marker_coverage`, `flank_left`, `flank_right`. |
| `src/progeny_selector/core/pipeline.py` | `chrom_lengths` from `chrom_length_bp(c, None, criteria.assembly)`; when a chromosome has no length, use the last marker position and warn `f"chromosome {c}: no {criteria.assembly} length; the last marker is the chromosome end for RPP weights and drag bounds"` (one warning per chromosome, at most 5 listed then `"... and N more"`); when the maximum marker position exceeds the length, use the maximum position and warn `f"chromosome {c}: marker positions reach {max_pos:.0f} beyond the {criteria.assembly} length {length}; check the assembly setting"`. `_pick_unit` becomes `_pick_unit(dataset, requested, what) -> tuple[str, str | None]` returning a warning `f"{what} is cm but the map has no cM; {'RPP weights computed' if what == 'background.map_unit' else 'windows interpreted'} in bp"` when `requested == "cm"` and no cM; the existing flank warning text is preserved verbatim. When the weighted model runs in bp, append `f"weighted RPP in bp: no cM map, weights capped at {cap:.0f} bp per marker"` with `cap` the effective `max_marker_coverage`. Ranking: `rec_count = rec_frac` before division (the integer sum of flagged flanks over targets); `rank_rows_staged(...)` when `criteria.ranking.mode == "staged"`, else `rank_rows` unchanged; `composite_score` computed in both modes. `chromosome_strips` callers (`app/screens/compare.py`) pass `state.criteria().assembly`. |
| `src/progeny_selector/core/strip.py` | Signature above. |
| `docs/data-formats.md` | criteria.yaml section: new top-level keys `assembly` (values `ASSEMBLIES`, default `Wm82.a4`; "chromosome ends for RPP weights, drag bounds and strips when positions are in bp; `none` uses the last marker") and `ranking` with `mode` (`weighted` default, `staged`: "hard filters first, then recombinant flanks (count) descending, `rpp_carrier` descending, `rpp_noncarrier` descending, `drag_total_est` ascending, `missing_rate` ascending, then `sample_id`, with exact ties at each key and no bins; `composite_score` is still written"); the download sentence per Q6; the input-contract paragraph's "Wm82.a4.v1 chromosome lengths in `constants.py`" becomes "the chromosome lengths of the assembly named by `assembly`". |
| `docs/adr/0015-assembly-selector-and-map-warnings.md` (new) | MADR: `assembly` in criteria.yaml, not markers.csv (contract-owned) nor UI-only (CLI parity); keyed length tables; fallback and mismatch warnings; the a2 NIL data as the motivating case; a5/a6 included only with fetched lengths. |
| `tests/test_ranking_staged.py` (new) | Hand-built arrays for 5 samples: s1 rec 2, carrier .5, noncarrier .5; s2 rec 1, carrier .9, noncarrier .99; s3 rec 1, carrier .9, noncarrier .95, drag 5; s4 rec 1, carrier .9, noncarrier .95, drag 3; s5 rec 1, carrier NaN. Expect order s1, s2, s4, s3, s5 (ranks 1..5); with `passes[s2] = False`, s2 gets NaN and s4 rank 2; `rank_in_family` with families `["a","a","b","b","b"]` is `[1, NaN, 2, 1, 3]` in the excluded case; `sample_id` breaks a full tie. |
| `tests/test_chrom_assembly.py` (new) | `chrom_length_bp("Gm01", None, "Wm82.a2") == 56_831_625`; `"none"` returns the fallback; `"Wm82.a4"` default equals the old behaviour; `normalize_chrom` untouched. |
| `tests/test_criteria_roundtrip.py` | Update for `ranking`, `assembly` and `null` keys: dumped fixture criteria starts `name: synthetic BC2F1 fixture`; a `Criteria(name=None, ...)` dump starts `name: null`; `ranking:\n  mode: weighted` and `assembly: Wm82.a4` present; `read_criteria_text` of `ranking: {mode: lexicographic}` raises matching "ranking.mode"; `assembly: Wm82.a9` raises matching "assembly". |
| `tests/test_map_warnings.py` (new) | With `make_matrix` markers lacking cM (build `Marker(..., cm=None)`), `criteria.background.map_unit = "cm"` and `flank_unit = "cm"`: `result.warnings` contains both cm fallback warnings verbatim; a marker at `pos_bp = 60_000_000` on Gm01 with assembly `Wm82.a4` yields the "beyond the Wm82.a4 length" warning and `drag` right bounds are non-negative; chromosome `scaffold_1` yields the "no Wm82.a4 length" warning. |
| `CHANGELOG.md` | Added: "`ranking: {mode: staged}` orders survivors lexicographically (docs/adr/0007); `assembly` selects the chromosome-length table (Wm82.a1, a2, a4 or none; docs/adr/0015). Warnings now report a cM request the map cannot honour, chromosomes without an assembly length, and marker positions beyond the assembly length." |

Edge cases: `assembly: none` with a weighted model in bp (terminal weights `cap/2`, drag right end at the last marker); all targets on chromosomes without length; `staged` with zero targets is impossible (validate requires one); staged mode with every sample NaN on `rpp_carrier` ranks by the next key.

Reviewer: yes. Acceptance: the three new test files pass; `tests/test_smoke_pipeline.py` unchanged and passing (fixture criteria has no `ranking`/`assembly`, defaults apply); `pytest -m e2e` unchanged.

## Phase 4: results.csv freeze

Renamed 2026-09-19, during the phase-4 review: the column the spec called `rpp_unit` is written as `background_unit`. `rpp_<chrom>` is the only dynamic family with no suffix, so a reader selects it by the `rpp_` prefix (the phase 5 R reader does), and a fixed column sharing that prefix broke the Compare screen's chromosome table. Renaming after the freeze would be a major schema bump. The text below is updated to the shipped name.

Goal: results.csv is schema 1.0.0: `NA` for missing cells, a header even with no rows, metadata columns, no literal `nan`; selected.csv carries the schema.

Precondition: Q2, Q3 (fifth column), Q4 answered; approval to edit the Outputs section of docs/data-formats.md (no sibling file affected).

```python
# io/export.py
RESULTS_SCHEMA = "1.0.0"
NA = "NA"
FIXED_COLUMNS = LEADING_COLUMNS + ("role", "n_informative_called", "frac_a", "frac_h", "frac_b",
                                   "background_model", "background_unit", "rank_mode", "assembly", "results_schema")
# 34 columns; drop "assembly" if Q3 is declined. Header written with no rows = FIXED_COLUMNS exactly.
SELECTION_COLUMNS = ("sample_id", "line_name", "family_id", "generation", "rank_overall", "rank_in_family", "composite_score", "rpp_total", "notes", "results_schema")
def _fmt(v): None -> NA; float NaN -> NA; bool -> TRUE/FALSE; float -> as now; str -> as is (an empty string stays empty)
```

| file | change |
| ---- | ------ |
| `src/progeny_selector/core/pipeline.py` | `_num(x: float | None) -> float | None` accepts None; rows route `het_rate`, `expected_het`, `expected_rpp` through it; rows gain `"background_model": criteria.background.model`, `"background_unit": unit`, `"rank_mode": criteria.ranking.mode`, `"assembly": criteria.assembly`, `"results_schema": RESULTS_SCHEMA` (import the constant from `io.export`? No: `core` must not import `io`; define `RESULTS_SCHEMA` in `constants.py` and have `io/export.py` re-export it). Insert the five keys after `frac_b` in the row dict. |
| `src/progeny_selector/io/export.py` | `_columns([])` returns `list(FIXED_COLUMNS)`; with rows, fixed columns present in `rows[0]` first, then the rest in dict order. `_write_selection_csv` writes `SELECTION_COLUMNS`; `family_id`/`generation` None -> `NA`; `notes.get(sid, "")` stays `""`; `results_schema` = `RESULTS_SCHEMA`. `_write_next_round_manifest` unchanged (empty cells, no schema). |
| `src/progeny_selector/cli.py` | `_read_results`: every cell equal to `"NA"` becomes None before the numeric conversions; docstring drops "downstream Shiny dashboards" for "the R reader `scripts/read_results.R`". |
| `docs/data-formats.md` | results.csv: the fixed prefix listed in order (34 names), then per-target, per-avoid, per-chromosome columns; "missing values are `NA`; empty text (`exclusion_reason`, `qc_flags`, `notes`) is an empty cell, not `NA`"; `background_model` (`count`/`weighted`), `background_unit` (`bp`/`cm`, the unit RPP weights used), `rank_mode` (`weighted`/`staged`), `assembly`, `results_schema` (SemVer; 1.0.0; a column added bumps minor, a column removed or renamed bumps major); "the header row is written when there are no results and then holds exactly the fixed prefix". selected.csv: the ten columns; `NA` rule; `results_schema`. next_samples.csv: unchanged columns and empty cells; "with `--per-selected N` (CLI) or the Export count, one row per planted plant is written and the file loads unchanged as the next generation's samples.csv". |
| `docs/adr/0016-results-schema-freeze.md` (new) | MADR: NA (as backcross), header on empty, metadata columns, schema string and bump rules, selected.csv carries it, next_samples.csv does not (input contract). |
| `tests/test_results_schema.py` (new) | Pins the 34-name list literally in the test (not imported); `results_csv_text([])` == that header joined by `,` plus `\r\n`; the fixture result's header starts with the list; every row's `results_schema == "1.0.0"`, `background_model == "count"`, `background_unit == "cm"`, `rank_mode == "weighted"`, `assembly == "Wm82.a4"`; the text contains no `,nan,`/`,nan\r\n`; the excluded `BC2F1-F1-002` row has `NA` in `rank_overall`; `selection_csv_text([row])` header equals `SELECTION_COLUMNS` and a row with `family_id=None` writes `NA` in that column and `""` in `notes`. |
| `tests/test_export_text.py`, `tests/test_smoke_pipeline.py` | `test_results_csv_roundtrip` also asserts the header prefix; CLI test unchanged in behaviour. |
| `tests/e2e/test_criteria_editor.py` | `test_results_before_load_is_header_only`: on a fresh page, Export `export-results` yields exactly the 34-name header plus `\r\n`. |
| `CHANGELOG.md` | Changed (with `!`/BREAKING CHANGE footer per CONTRIBUTING): "results.csv schema 1.0.0: missing cells are `NA`, the header is written with no results, new columns `background_model`, `background_unit`, `rank_mode`, `assembly`, `results_schema`; `het_rate`, `expected_het` and `expected_rpp` never print `nan`. selected.csv gains `results_schema` and writes `NA` for missing family or generation (docs/adr/0016)." |

Edge cases: `rank_in_family` NaN for excluded rows -> `NA`; `_read_results` on a pre-freeze file (empty cells) still works; a `qc_flags` of `""` is not NA; `frac_*` are floats never None.

Reviewer: yes. Acceptance: `tests/test_results_schema.py` passes; `pytest -m e2e` passes with the new header test.

## Phase 5: R reader and CI job

Goal: `scripts/read_results.R` reads results.csv unchanged with explicit column types in CI, on the fixture and on the empty case; not part of the per-commit gates.

| file | change |
| ---- | ------ |
| `scripts/read_results.R` (new) | `#!/usr/bin/env Rscript`; `suppressPackageStartupMessages(library(readr))`; `FIXED_COLUMNS <- c(...)` the 34 names; `read_results <- function(path)` calling `readr::read_csv(path, col_types = cols(rank_overall = col_integer(), rank_in_family = col_integer(), sample_id = col_character(), line_name = col_character(), family_id = col_character(), generation = col_character(), passes_filters = col_logical(), exclusion_reason = col_character(), composite_score = col_double(), foreground_all_pass = col_logical(), avoid_all_pass = col_logical(), rpp_total = col_double(), rpp_carrier = col_double(), rpp_noncarrier = col_double(), expected_rpp = col_double(), drag_total_est = col_double(), drag_total_max = col_double(), drag_unit = col_character(), ibs_rp = col_double(), ibs_donor = col_double(), missing_rate = col_double(), het_rate = col_double(), expected_het = col_double(), qc_flags = col_character(), role = col_character(), n_informative_called = col_integer(), frac_a = col_double(), frac_h = col_double(), frac_b = col_double(), background_model = col_character(), background_unit = col_character(), rank_mode = col_character(), assembly = col_character(), results_schema = col_character(), .default = col_guess()), na = "NA", show_col_types = FALSE)`; when run as a script (`commandArgs(trailingOnly = TRUE)[1]`): `stop()` unless `names(df)[seq_along(FIXED_COLUMNS)]` equals `FIXED_COLUMNS`; when `nrow(df) > 0`, `stop()` unless every `results_schema` starts with `"1."`; print `nrow`, number with `passes_filters`, and the per-locus column names (those matching `^(target|avoid)_.*_status$`); exit 0. Dynamic columns: `target_*_status`/`avoid_*_status` character, `drag_*`/`rpp_<chrom>` double, `recomb_*` logical, by guess. |
| `.github/workflows/ci.yml` | New job `r-reader` (needs `check`, ubuntu-latest): checkout; setup-python 3.12; `pip install -e .`; `r-lib/actions/setup-r@v2` with `use-public-rspm: true`; `r-lib/actions/setup-r-dependencies@v2` with `packages: any::readr` (not verifiable here whether it accepts `packages` without a DESCRIPTION; fallback step `Rscript -e 'install.packages("readr", repos = "https://cloud.r-project.org")'`); `progeny-selector rank --genotypes tests/fixtures/synthetic_bc2f1/genotypes.vcf --samples ... --markers ... --criteria ... --out results.csv`; `python -c "from progeny_selector.io.export import write_results_csv; write_results_csv([], 'empty_results.csv')"`; `Rscript scripts/read_results.R results.csv`; `Rscript scripts/read_results.R empty_results.csv`. Not in `check`; not in the local gates. |
| `README.md`, `CLAUDE.md` | README "Quickstart" gains `Rscript scripts/read_results.R results.csv`; CLAUDE.md Gates: "CI job `r-reader` (R and readr) reads the fixture's results.csv; not a per-commit gate." |
| `CHANGELOG.md` (via report) | Added: "`scripts/read_results.R` reads results.csv with readr and explicit column types; CI verifies it on the fixture and on an empty results file." |

Reviewer: no. Acceptance: the `r-reader` job is green on both files on the first push; the Python schema test (Phase 4) and the R column list are identical (the doer diffs them by eye and says so in the report).

## Phase 6: weighted and staged expectations, duplicate detection, per-target windows

Corrected 2026-09-21, before dispatch: the table below says the duplicate IBS is computed over all 500 markers, but decision 5 says informative markers called in both members, and the decisions override this section. The shipped rule is decision 5's, and it needs one change the table does not list. `core/similarity.py` `pairwise_ibs` gains `marker_idx: np.ndarray | None = None`, defaulting to all markers so no other caller changes, and subsamples at most `max_markers` from `marker_idx` with the existing seed-0 generator; `core/qc.py` `duplicate_pairs` gains the same parameter and passes it through; `core/pipeline.py` passes `np.flatnonzero(classification.informative)`, which is already computed before `sample_qc`. The generator computes the same way over the fixture's 475 informative markers, so no subsampling happens there and the expectation is exact. `tests/test_duplicates.py` must build all 30 markers informative (RP and DONOR both called, both homozygous, different alleles) or the restricted denominator is zero and nothing is flagged at all. ADR 0017 says "at most 2,000 informative markers" and carries decision 5's sentence about the threshold being stricter than PLINK and KING practice on purpose. The ADR number stays 0017: 0018 is the sibling's number for a backcross-side decision, not a gap here.

Amended 2026-09-21, after the phase 6 review: duplicate IBS also carries a minimum-overlap floor. A pair is reported only when the two individuals share calls at `MIN_DUPLICATE_OVERLAP_FRAC` (0.5) or more of the markers used for that comparison, because without it a sample with almost no calls scores IBS 1.0 against anyone it agrees with on its handful of calls. The maintainer decided it on 2026-09-21; `docs/adr/0017-duplicate-flag.md` and its amendment section carry the reasoning and the recorded trade-off. Phase 7 shows duplicates on the Validate screen and should describe what the screen reports in those terms.

Goal: the fixture pins the weighted model and staged order from the generator's independent implementation; duplicates are detected in the pipeline; per-target windows have a hand-built test.

Precondition: Q5, Q11, Q12 answered.

| file | change |
| ---- | ------ |
| `scripts/make_fixture.py` | `expected_metrics` adds, computed independently: `rpp_total_weighted`, `rpp_carrier_weighted`, `rpp_noncarrier_weighted` (cM, cap 10 cM: per informative marker in cM order per chromosome, `w = min(gap_left/2, 5) + min(gap_right/2, 5)` with `5` on the outer side of the first and last informative marker; per individual `sum(w*s)/sum(w)` over called A/H/B informative markers with `s` 1/0.5/0; NaN when the denominator is 0); `rank_overall_staged`, `rank_in_family_staged` (passing rows sorted by `-(rec_left + rec_right)`, `-rpp_carrier`, `-rpp_noncarrier`, `drag_total_est_cm`, `missing_rate`, `sample_id`; dense; `""` for excluded); `possible_duplicate` (pairwise over all 40 progeny and all 500 markers: shared alleles per marker from the written calls, `1` if identical pair, `0.5` if one shared, else `0`, averaged over markers called in both; flag both members of any pair `>= 0.995`; assert the set is empty for this fixture). Docstring updated. |
| `src/progeny_selector/core/pipeline.py` | After `qc = sample_qc(...)`: `dups = duplicate_pairs(gm, sample_ids)`; append `"possible_duplicate"` to `qc[i].flags` for each member (once); `AnalysisResult.duplicates: list[tuple[str, str, float]]`. |
| `src/progeny_selector/core/qc.py` | Docstring lists `possible_duplicate` as advisory; `QC_EXCLUDING_FLAGS` unchanged. |
| `src/progeny_selector/cli.py` | `validate` prints `duplicate pairs: N` and up to 20 lines `  a, b: ibs 0.998`. |
| `docs/adr/0017-duplicate-flag.md` (new) | MADR: advisory `possible_duplicate` on both members, threshold 0.995 on at most 2,000 subsampled markers, skipped above 2,000 individuals; not excluding because a true duplicate sample is a bookkeeping error the breeder resolves, not a genotype fault. |
| `tests/test_smoke_pipeline.py` | `test_weighted_matches_expected`: `run_analysis(dataset, dataclasses.replace(criteria, background=BackgroundOptions(model="weighted", map_unit="cm")))` and compare `rpp_total`, `rpp_carrier`, `rpp_noncarrier` with the `*_weighted` columns (1e-6); rows carry `background_model == "weighted"`. `test_staged_ranks_match_expected`: `ranking=RankingOptions(mode="staged")`; compare ranks; `composite_score` still equal to the weighted run's. `test_no_duplicates_in_fixture`: `result.duplicates == []` and no row has `possible_duplicate`. |
| `tests/test_duplicates.py` (new) | Build a `GenotypeMatrix` with columns RP, DONOR, P1, P2, P3 over 30 markers where P2 equals P1 except one missing call and P3 differs at 10 markers; `make_dataset`-style samples; `run_analysis` flags P1 and P2 `possible_duplicate`, not P3; `result.duplicates == [("P1", "P2", ibs)]` with `ibs >= 0.995`; both still `passes_filters` when their loci pass. |
| `tests/test_status_metrics.py` | `test_per_target_windows`: two targets on the `STATES` case with `flank_left=1.0`, `flank_right=20.0` on one and defaults on the other; `recomb_<id>_left/right` differ exactly as the windows dictate against the hand-computed `left_max`/`right_max`. |
| `docs/data-formats.md` | targets table row for `flank_left, flank_right` gains "per target; overrides `flank_window` on that side only" (additive wording). |
| `CHANGELOG.md` | Added: "Advisory QC flag `possible_duplicate` and duplicate pairs in `validate` and the Validate screen (docs/adr/0017). The fixture pins weighted RPP and staged ranks computed independently." |

Edge cases: a duplicate pair where one member is high_missing (still flagged; IBS uses markers called in both); more than 2,000 individuals returns no pairs and a warning `"duplicate detection skipped above 2000 individuals"` appended to `result.warnings`; a progeny identical to a parent is already `possible_rp_sample`/`possible_donor_sample` and is not in the progeny-only pairwise set.

Reviewer: yes. Acceptance: regenerated fixture committed via the generator with no hand edits; the new tests pass; `git diff --exit-code -- tests/fixtures` clean after a second run.

## Phase 7: UI: notes, map block, duplicates, per-selected, keyboard walkthrough

Goal: notes typed on the Selection screen reach selected.csv; the Validate screen states model, units, assembly and duplicates; Export writes N placeholders per selected; the keyboard walkthrough is documented and tested.

| file | change |
| ---- | ------ |
| `src/progeny_selector/app/screens/select.py` | Layout: `ui.card(ui.card_header("Selection list"), ui.layout_columns(ui.input_numeric("top_n", ...), ui.input_action_button("apply_top_n", ...), ui.input_radio_buttons("step", ...), col_widths=(3, 3, 6)), ui.output_data_frame("table"), ui.output_text_verbatim("summary"))`. `table`: `render.DataGrid(frame, editable=True, selection_mode="none", height="50vh")` over columns `NOTES_COLUMNS = ("sample_id", "line_name", "family_id", "rank_overall", "composite_score", "rpp_total", "notes")`, rows in `state.selected_ids()` order, `notes` from `state.notes()`. `@table.set_patch_fn` `async def _(*, patch: render.CellPatch) -> render.CellValue`: if `patch["column_index"] != NOTES_COLUMNS.index("notes")`, return the original cell value from `table.data()`; else `sid = table.data().iat[patch["row_index"], 0]`, `state.notes.set({**state.notes(), sid: str(patch["value"])})`, return `patch["value"]`. Notes are keyed by `sample_id`, never by row. Docstring rewritten (no "Placeholder"). |
| `src/progeny_selector/app/screens/load.py` | `_run` resets `state.notes.set({})` on success; `_apply` keeps notes. |
| `src/progeny_selector/app/screens/export.py` | `ui.input_numeric("per_selected", "Placeholder rows per selected individual", value=1, min=1, max=999)`; `manifest` passes `n_per_selected=int(input.per_selected())`. |
| `src/progeny_selector/app/screens/qc.py` | Summary `dl` gains, after "map unit": `background model` (`state.criteria().background.model`), `RPP unit` (`result.unit`), `drag unit` (first row's `drag_unit` or `-`), `assembly` (`state.criteria().assembly`), `rank mode`, `duplicate pairs` (`ui.tags.ul` of `f"{a}, {b}: IBS {ibs:.3f}"` or `"none"`). |
| `src/progeny_selector/app/app.py` | Docstring status: "M2 in progress; every screen runs the real pipeline". |
| `docs/keyboard-walkthrough.md` (new) | Step by step from page load to selected.csv using only the keyboard: Tab order on Load (four file inputs, Load and analyse), navbar tabs (`Tab` to the tab list, arrows between tabs), Navigate (accordion `Enter`, radios with arrows), Rank (Tab into the grid; arrow keys move; the selection keys as observed), Selection (Tab into the grid, `Enter` edits a cell, `Enter` commits, `Escape` cancels), Export (Tab to a button, `Enter` downloads). Records the observed behaviour of DataGrid row selection under shiny 1.7 (Q16) and any step that needs a mouse. |
| `tests/e2e/test_selection_notes.py` (new) | Load fixture; Rank `select_rows([0, 2])`; Selection: `OutputDataFrame(page, "select-table").expect_nrow(2)`; sort by `sample_id` descending (`set_sort`) so view row 0 is `BC2F1-F1-010`; edit its `notes` cell to `"keep; vigorous"` (controller method not verifiable here: try `set_cell(text, row=0, col=NOTES_COL)`, else double-click the cell, `fill`, `Enter`; report which); Export `export-selected`: the row with `sample_id` `BC2F1-F1-010` has `notes == "keep; vigorous"` and `BC2F1-F1-001` has `""`; header equals `SELECTION_COLUMNS`. Then "Add top N per family" with N=1 keeps the note. Set `export-per_selected` to 3, download `export-manifest`: `2 + 2*3` data rows, ids ending `-BC3F1-001..003`. This is the check that `CellPatch.row_index` indexes `.data()` and not the sorted view. |
| `tests/e2e/test_keyboard_walkthrough.py` (new) | Executes the documented steps: after load, reach Rank by keyboard, Tab into `#rank-table`, press `ArrowDown` twice and the documented selection key; assert what the walkthrough states (a selected row count, or focus movement alone if selection needs a mouse). |
| `tests/e2e/test_qc_compare.py` | `test_qc_table_and_summary` also expects the summary `dl` to contain `background model`, `count`, `assembly`, `Wm82.a4`, `duplicate pairs`, `none`. |
| `CHANGELOG.md` (via report when concurrent) | Added: "Selection list: editable notes written to selected.csv. Validate: background model, units, assembly, rank mode and duplicate pairs. Export: placeholder rows per selected individual. docs/keyboard-walkthrough.md." |

Not verifiable here: `render.CellPatch` key names (`row_index`, `column_index`, `value`), `set_patch_fn` signature, the controller's cell-edit method, DataGrid keyboard selection. The e2e tests are the checks; the report names what was found.

Reviewer: no. Acceptance: both new e2e tests pass; `pytest` unaffected.

## Phase 8: two-generation synthetic round trip

Goal: next_samples.csv written from the BC2F1 selection is, unchanged, the samples.csv of a BC3F1 fixture generated independently, and that generation ranks against its own expected results.

Precondition: Q13 answered.

| file | change |
| ---- | ------ |
| `scripts/make_fixture.py` | Keep the true (pre-missing) states of every BC2F1 individual. After `expected_metrics`, pick the generator's own top 2 per family by `rank_in_family` (expected: `BC2F1-F1-001`, `BC2F1-F1-010`, `BC2F1-F2-005`, `BC2F1-F2-019`; assert). For each selected parent, 10 BC3F1 progeny: parent haplotypes (RP all 0; donor-side 1 where the true state is H), gamete by `haldane_gamete`, other gamete RP, state H where the gamete is 1; 1 % random missing never at the target or avoid marker. Ids `f"{parent}-BC3F1-{k:03d}"`, `line_name` = the parent's line_name (what the writer copies), role `progeny`, generation `BC3F1`, family_id = parent id, notes `f"derived from {parent}"`. Write `tests/fixtures/synthetic_bc3f1/`: `genotypes.vcf` (same markers, columns RP, donor, 40 progeny), `samples.csv` with rows in `select_top_n` order (sorted by family_id then rank_in_family, then k), header `sample_id,line_name,role,generation,family_id,notes`, parent rows `[RP_ID, "Williams 82 (synthetic)", "recurrent_parent", "", "", "synthetic recurrent parent"]` and the donor equivalent, LF endings; `expected_results.csv` from `expected_metrics` parametrised with `expected_het=0.125` (`het_rate_deviates`), families of 10 (`family_donor_outlier` rule applies), the same target/avoid/flank constants; `README.md` naming the design. `criteria.yaml` and `markers.csv` are not duplicated: the BC3F1 fixture reuses `synthetic_bc2f1/criteria.yaml` and `markers.csv`. Print line gains the BC3F1 counts. |
| `src/progeny_selector/cli.py` | `select --per-selected N` (int, default 1, `>= 1`), passed to `write_next_round_manifest`. |
| `tests/test_round_trip.py` (new) | (1) gen1: `select_top_n(result.rows, 2, per_family=True)`; `next_round_manifest_text(chosen, "BC3F1", rp, donor, n_per_selected=10).replace("\r\n", "\n")` equals `synthetic_bc3f1/samples.csv` text exactly. (2) gen2: `load_dataset(bc3f1/genotypes.vcf, bc3f1/samples.csv, bc2f1/markers.csv)` with the bc2f1 criteria; compare statuses, exclusion reasons, RPP, drag, recombinant flags, composite scores, ranks and advisory flags with `bc3f1/expected_results.csv` as `test_smoke_pipeline` does; every row `expected_rpp == 0.9375` and `expected_het == 0.125`; `build_tree` gives 4 families of 10, one generation `BC3F1` each. (3) CLI: `rank` gen1 -> `select --top 2 --per-selected 10 --next-manifest m.csv --samples bc2f1/samples.csv` -> `m.csv` (CRLF normalised) equals the fixture samples.csv -> `rank --samples m.csv` on gen2 genotypes prints `40 individuals`. |
| `.gitignore` | Unchanged (`!tests/fixtures/**` already un-ignores the new folder). |
| `docs/adr/0018-two-generation-round-trip.md` (new) | MADR: the manifest format is pinned by two implementations (writer vs generator); `--per-selected`; minimal requirements for a real dataset (Q10). |
| `docs/data-formats.md` | next_samples.csv sentence per Phase 4 (already edited there; this phase adds the CLI flag name). |
| `CHANGELOG.md` | Added: "`select --per-selected N`; a synthetic BC3F1 fixture generated from the BC2F1 selection proves next_samples.csv loads unchanged as the next generation (docs/adr/0018)." |

Edge cases: a selected parent with a missing call at the target in gen1 (true state used for haplotypes, so gen2 is consistent); `--per-selected 0` is a usage error (exit 2); the writer's `line_name` copy means ten progeny share a line_name (documented in the fixture README).

Reviewer: no (scripts, tests, cli). Acceptance: `test_round_trip.py` passes; fixture gate clean; CI `check` green.

## Phase 9: genetic map interpolation and the Song 2016 converter

Goal: unmapped markers on the same assembly get cM by monotone, clamped Marey-map interpolation; Song et al. 2016 Table S1 converts to markers.csv.

Precondition: Q8 answered; Table S1 exported to CSV in gitignored `data/song2016/`.

```python
# core/genetic_map.py (new, pure numpy)
def monotone_cm(bp: np.ndarray, cm: np.ndarray) -> np.ndarray
    # sort by bp; running maximum of cm in that order; returned in the input order; ties in bp keep the max
def interpolate_cm(mapped_bp: np.ndarray, mapped_cm: np.ndarray, query_bp: np.ndarray) -> np.ndarray
    # mapped_cm made monotone first; np.interp on sorted mapped_bp; queries beyond the ends are clamped to the end cM;
    # fewer than 2 mapped markers -> the single cM for every query, or NaN when none
```

| file | change |
| ---- | ------ |
| `scripts/song2016_map.py` (new) | `python3 scripts/song2016_map.py --table data/song2016/table_s1.csv --assembly a2|a1 --out markers.csv [--markers-in genotypes_markers.csv] [--columns marker=SNP,chrom=LG,pos_a1=...,pos_a2=...,cm=...]`. Reads the CSV (delimiter sniffed by `io.delimited.sniff_delimiter`), maps columns per `--columns` (defaults filled by the doer after inspecting the file and recorded in the docstring with the citation "Song et al. 2016, BMC Genomics 17:33, Table S1, CC BY 4.0"), normalises chromosomes with `core.chrom.normalize_chrom`, drops rows without a position on the chosen assembly (counted), writes `marker_id,chrom,pos_bp,cm` sorted by chromosome and position. With `--markers-in` (a markers.csv-shaped file, e.g. from Phase 10's position table), every marker absent from the table gets `cm` from `interpolate_cm` per chromosome; chromosomes with no mapped marker get an empty `cm` and a warning that cM mode will be disabled. Writes a `.log` beside `--out` with counts (mapped, interpolated, clamped, unplaced). Positions are on one assembly; mixing is refused. |
| `tests/test_genetic_map.py` (new) | Hand-computed: mapped bp `[1e6, 2e6, 3e6, 4e6]`, cm `[0, 5, 4, 10]` -> monotone `[0, 5, 5, 10]`; query `[0.5e6, 2.5e6, 3.5e6, 9e6]` -> `[0, 5, 7.5, 10]`; a single mapped marker gives its cM everywhere; none gives NaN; input order preserved. |
| `tests/test_song2016_map.py` (new) | A 12-row hand-written table with two chromosomes and one unmapped-assembly row; output rows, cM values and the log counts asserted; `load_dataset` on a tiny wide CSV plus the written markers.csv has `has_cm()` true. |
| `docs/soybean-inputs.md` (new, shared with Phase 10) | The three scripts, inputs, licences (SoyBase Open; Song 2016 CC BY 4.0), what is never committed, the assembly caveat. |
| `docs/adr/0019-soybean-input-tooling.md` (new, shared with Phase 10) | MADR: Marey interpolation by running maximum and clamping; no extrapolation; scripts, not io, because these are one-off conversions outside the contract; KASP long-format assumption. |
| `CHANGELOG.md` | Added: "`scripts/song2016_map.py` converts Song et al. 2016 Table S1 to markers.csv with monotone, clamped Marey-map interpolation for unmapped markers (`core/genetic_map.py`, docs/adr/0019)." |

Reviewer: yes (`core/genetic_map.py`). Acceptance: the two test files pass on hand-built data; the doer runs the script on the real table locally and reports counts (not committed).

## Phase 10: SoySNP position table and KASP converter

Goal: a marker position table across assemblies from SoyBase GFF3s, and a KASP export converted to the contract's wide CSV.

Precondition: Q7, Q9 answered.

| file | change |
| ---- | ------ |
| `scripts/soysnp_positions.py` (new) | `python3 scripts/soysnp_positions.py --gff3 Wm82.a2=data/soybase/glyma.Wm82.gnm2.mrk.SoySNP50K.gff3.gz [--gff3 ...] [--panel SoySNP50K|SoySNP6K per file via name] --out data/soybase/soysnp_positions.csv [--emit-markers-csv Wm82.a2 --markers-out markers.csv] [--download]`. GFF3 rows (`seqid source type start end score strand phase attributes`): `marker_id` from `ID=` (fallback `Name=`), chromosome from `seqid` with the `glyma.Wm82.gnmN.` prefix stripped as `scripts/soysnp50k_nils.py::strip_chrom` does and normalised; scaffolds kept with their names; `start` as `pos_bp`. Output wide table: `marker_id, in_SoySNP50K, in_SoySNP6K, chrom_Wm82.a1, pos_bp_Wm82.a1, chrom_Wm82.a2, pos_bp_Wm82.a2, chrom_Wm82.a4, pos_bp_Wm82.a4, chrom_Wm82.a5, pos_bp_Wm82.a5, chrom_Wm82.a6, pos_bp_Wm82.a6` (empty where absent), joined on `marker_id`. `--emit-markers-csv` writes `marker_id,chrom,pos_bp` for one assembly, Gm01..Gm20 only, sorted. `--download` fetches the URL table `SOYBASE_URLS` at the top of the script (entries resolved by the doer from the Data Store listing, recorded with date in docs/reference-repos.md; not verifiable here) into `data/soybase/`; never run in CI. Also prints per-assembly max position per chromosome so lengths can be sanity-checked against `constants.py`. |
| `scripts/kasp_to_wide.py` (new) | `python3 scripts/kasp_to_wide.py --kasp export.csv --markers markers.csv --out genotypes.csv [--drop-unplaced] [--controls NTC,H2O]`. Assumed layout per Q9: long CSV with `SubjectID`, `SNPID`, `Call` (case-insensitive header match; delimiter sniffed). Calls `X:Y` become the nucleotide pair `XY`; `Uncallable`, `Missing`, `?`, `Bad`, `Dupe`, `NTC`, empty and any call not matching `^[ACGT]:[ACGT]$` become `N` (counted per reason); rows whose `SubjectID` is in `--controls` (default `NTC`) are dropped; conflicting calls for one sample x SNP are an error naming both rows. Output: wide CSV `marker_id,chrom,pos_bp,<samples...>` with chrom/pos from `--markers`; a SNP absent from the map is an error listing ids unless `--drop-unplaced` (counted). Column order: samples in first-seen order; markers sorted by chromosome and position. |
| `tests/test_soysnp_positions.py` (new) | Two hand-written 6-line GFF3 snippets (a2 and a4) with one marker only in a2 and one scaffold; the joined table has the expected cells and blanks; `--emit-markers-csv Wm82.a4` excludes the scaffold and loads through `read_markers`. |
| `tests/test_kasp_to_wide.py` (new) | A 10-row long CSV with two samples, one NTC row, one `Uncallable`, one `A:G`; output loads with `load_dataset` (nucleotide coding detected), the het is `H` after classification with hand-built parents, the missing is `N`; a conflicting duplicate row raises; an unplaced SNP raises without `--drop-unplaced` and is dropped with it. |
| `docs/soybean-inputs.md`, `docs/adr/0019-...`, `docs/reference-repos.md` | Sections for the two scripts; URLs and fetch dates. |
| `CHANGELOG.md` (via report when concurrent) | Added: "`scripts/soysnp_positions.py` builds a SoySNP50K/6K position table across Wm82 assemblies from SoyBase GFF3s; `scripts/kasp_to_wide.py` converts an LGC long-format KASP export to the wide CSV contract." |

Reviewer: no. Acceptance: both test files pass on hand-built inputs; the doer runs `soysnp_positions.py` on the downloaded GFF3s locally and reports row counts and per-assembly maxima (nothing under `data/` committed).

## Phase 11: real-data acceptance, Question B, end of milestone

Goal: the M1 gap and M2 acceptance are exercised on real data by the maintainer; the milestone diff gets its reviewer pass; the documents are current.

| item | who | record |
| ---- | --- | ------ |
| Browser run on real segregating data, local `shiny run`: load `data/nils/` (or the program's KASP/SoySNP6K file converted by Phase 10) with `assembly: Wm82.a2`, weighted model; take it to selected.csv with notes; network panel open. | maintainer | PLAN.md M2 verification block: markers, progeny, passing, wall-clock to status, warnings shown, "no request left the origin". |
| Same on the Pages site (or `python -m http.server` on `site/`). | maintainer | same block; the Shinylive size and time to `#load-run`. |
| Question B real-data round trip: generation n -> next_samples.csv (`--per-selected` as planted) -> generation n+1 loaded unchanged -> results. **Blocked on the maintainer's dataset (Q10).** | maintainer | same block, or the statement that no dataset was available. |
| End-of-milestone reviewer pass over the whole M2 diff. | main session | findings to the responsible implementer. |
| Docs: PLAN.md ("Milestones" M2 done, "Testing and CI" names `r-reader`, "Repository layout" ADRs 0001..0019, M2 verification block, Handoff); CLAUDE.md State "M2 complete, M3 next", gates; README "What exists now"; `core/chrom.py` docstring; CLAUDE.md `../isoline-browser` note (Q17). | sonnet (docs-only, gates only) | |

Acceptance: the verification block exists with the browser figures; the reviewer pass has no open high or medium finding.

## M2 acceptance checklist

- [ ] Validate and Rank agree on QC exclusion under `exclude_qc_flagged: true` and `false` (`tests/test_qc_table.py::test_qc_excluded_follows_filters`).
- [ ] criteria.yaml rejects aliases, oversized documents, NaN/inf, fractional positions; numeric ids read as text (`tests/test_criteria_hardening.py`).
- [ ] Remaining deferred findings closed: Load accept list and error catch, `export` extra without `build`, versioned shinylive cache key, context-level request recording, "(no family)" browser test, second Compare drag assertion.
- [ ] `ranking: {mode: staged}` orders exactly as ADR 0007's amendment states; `composite_score` and `rank_mode` written (`tests/test_ranking_staged.py`, fixture staged columns).
- [ ] `assembly` selects the chromosome-length table; cM fallbacks, missing lengths and positions beyond the assembly are warnings shown on Load and Validate (`tests/test_map_warnings.py`, e2e summary).
- [ ] results.csv schema 1.0.0: `NA` for missing cells, header with no rows, `background_model`, `background_unit`, `rank_mode`, `assembly`, `results_schema`; no literal `nan` (`tests/test_results_schema.py`); selected.csv carries `results_schema`.
- [ ] `scripts/read_results.R` reads the fixture's results.csv and the empty file in the CI `r-reader` job with explicit col_types and `na = "NA"`.
- [ ] Weighted RPP and staged ranks on the fixture match the generator's independent implementation; `possible_duplicate` detected in the pipeline and shown on Validate and in `validate` (`tests/test_duplicates.py`).
- [ ] Per-target windows covered by a hand-built test and documented.
- [ ] Notes edited on the Selection screen appear in selected.csv on the right `sample_id` after sorting (`tests/e2e/test_selection_notes.py`).
- [ ] next_samples.csv written from the BC2F1 selection with `--per-selected 10` equals the BC3F1 fixture's samples.csv and loads unchanged; BC3F1 results match expected (`tests/test_round_trip.py`).
- [ ] Keyboard walkthrough in `docs/keyboard-walkthrough.md` with a browser test asserting what it states.
- [ ] Soybean input scripts with hand-built tests: position table, Song 2016 converter with Marey interpolation, KASP converter; nothing under `data/` committed.
- [ ] Real segregating data taken to selected.csv in the browser locally and on the Pages site, recorded in PLAN.md; two-generation real-data step recorded or marked blocked.
- [ ] ADRs 0015-0019 present; CHANGELOG entries under `[Unreleased]`; end-of-milestone reviewer pass done; PLAN.md Handoff updated.
