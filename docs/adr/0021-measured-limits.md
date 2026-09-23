# Measured limits, not stated

Status: accepted. Date: 2026-09-22.

## Context and Problem Statement

An earlier plan (struck) claimed "2,000 individuals under 60 s" without a machine, a Python, a
numpy version or a script anyone could re-run. `docs/m3-phases.md` (Phase 2) replaces that
claim with two scripts, `scripts/bench_pipeline.py` and `scripts/bench_shinylive.py`, that
generate a deterministic synthetic dataset at a given size and measure it, and a table,
`docs/limits.md`, that records what they measured, on which machine, with which package
versions, so a reader can tell what changed and re-run the same measurement themselves.

## Decision Outcome

- **The CPython memory figure is the `tracemalloc` peak.** `tracemalloc.get_traced_memory()[1]`
  taken around each stage separately, because numpy reports its own array allocations to
  `tracemalloc` (it has since numpy 1.13) and this project's hot arrays are numpy's. Bracketing
  each stage with its own `tracemalloc.start()`/`tracemalloc.stop()` means a stage's peak
  excludes memory already retained by an earlier stage (`docs/m3-phases.md`, Q4 line, "measure_case
  lays out"); the parser's peak does not include the loaded `Dataset` it returns, for instance.
- **The browser memory figure is `performance.measureUserAgentSpecificMemory()`.** It is the API
  that sees the Pyodide worker's heap rather than only the tab's own JS heap, which is why it is
  worth the isolation requirement it carries. It requires cross-origin isolation
  (`crossOriginIsolated`); `scripts/bench_shinylive.py` serves the exported site with
  `Cross-Origin-Opener-Policy: same-origin` and `Cross-Origin-Embedder-Policy: require-corp`,
  asserts `crossOriginIsolated` in the page, and exits 3 rather than printing a memory figure of
  0 or otherwise fabricated when isolation is not available (Q4). Whether Shinylive's own
  same-origin iframe and service worker load under `require-corp` is not verifiable in this
  session; if they do not, that returns as a question, not a doer's choice (Q4).
- **Measured per milestone, not per commit (Q5).** The two scripts are not run in CI; CI runs
  only the deterministic, fixture-scale `tracemalloc` bound tests (`tests/test_bench_generator.py`
  here; the analysis-memory bound tests of Phase 6). A wall-clock regression is caught at the next
  milestone's re-measurement, not the next commit.
- **Bench inputs live under gitignored `data/bench/`.** Nothing under `data/` is ever committed
  (`.gitignore`); the generator is deterministic (`np.random.default_rng(seed)`) so anyone can
  reproduce a case's inputs from the command line rather than from a committed file, the same
  reasoning that keeps the real fixture generated rather than hand-edited.
- **The reference machine is the maintainer's laptop**, as the M1 and M2 limits blocks were
  measured on; `docs/limits.md` records its CPU, RAM and OS alongside the Python, numpy, shiny
  and shinylive versions, so a figure can be told apart from a different machine's.
- **The generated target region's half-width is `max(1_000_000, d)`, not a flat 1 Mb**, where
  `d` is that chromosome's marker spacing at the requested marker count: an absolute 1 Mb window
  is empty whenever marker spacing exceeds it, which is every small case (`generate_case`'s own
  smoke-test size among them), so a fixed width would make the smoke test's target locus resolve
  to no markers for a reason no later reader would find obvious.

## Consequences

A performance claim in this project is now a row a reader can regenerate: `scripts/bench_pipeline.py
--markers M --progeny N --out data/bench/<case>` writes the same inputs and prints the same
Markdown row on any machine with the same seed. A regression shows up as a smaller or larger
number in the next milestone's re-measurement of the same script against the same generated
inputs (`docs/m3-phases.md`, Phase 9), not as a stale adjective.

## Revisit when

- Shinylive is confirmed (or not) to load under `require-corp`; if not, the fallback (Chromium
  process-tree RSS via `psutil`, or a figure the app itself prints) becomes this ADR's amendment,
  not a doer's silent substitution.
- The streaming parser (Phase 3) or the sample-chunked core (Phase 6) lands: Phase 9 re-runs
  `scripts/bench_pipeline.py` on the same generated inputs and records the before/after in
  `docs/limits.md`.
