# Measured limits

Date: 2026-09-22.

Machine: Intel(R) Core(TM) Ultra 5 125U; 16,573,128,704 bytes (about 15.4 GiB) RAM; Windows 11 Home (build 10.0.26200).

Python: 3.12.4 (`MSC v.1940 64 bit (AMD64)`).
numpy: 2.5.3.
shiny: 1.7.0.
shinylive: 0.8.11.
Chromium: pending (measured by `scripts/bench_shinylive.py`, not run this phase).
Commit: pending (written by the main session).

## CPython

| case | markers | progeny | VCF MB plain / gz | parse s plain / gz | parse peak MiB | resident matrix MiB | analysis s | analysis peak MiB | duplicate scan s (ran) | export s | total s | outcome |
| ---- | ------- | ------- | ------------------ | -------------------- | ---------------- | ---------------------- | ----------- | -------------------- | ------------------------ | -------- | ------- | ------- |
| 6k_x_2000 | 6000 | 2000 | 48.3 / 2.3 | 262.29 / 24.04 | 49.7 | 22.9 | 555.97 | 412.7 | 517.28 (ran) | 0.26 | 1335.80 | ok |
| 50k_x_200 | 50000 | 200 | 42.2 / 2.4 | 218.01 / 22.66 | 69.1 | 19.3 | 33.28 | 345.2 | 4.70 (ran) | 0.02 | 256.01 | ok |
| 50k_x_2000 | pending | pending | pending | pending | pending | pending | pending | pending | pending | pending | pending | pending |

`50k_x_2000` is deferred: on today's un-chunked `core`, `rpp` alone allocates roughly 3.2 GB of
float64 temporaries at this size (Phase 6, not yet landed), on top of the roughly 400 MB VCF; left
for the maintainer to run on the reference machine rather than attempted here.

Duplicate-scan check (docs/m3-phases.md edge case): a 400-marker, 2001-progeny run printed
`duplicate scan: skipped` and its row read `0.00 (skipped)`; a 400-marker, 2000-progeny run's row
read `87.91 (ran)`. `core/pipeline.py` skips `duplicate_pairs` only when `len(sample_ids) >
_MAX_DUPLICATE_SAMPLES` (2000), so exactly 2,000 progeny runs the scan and 2,001 does not — both
rows above (`6k_x_2000`, `50k_x_200`, both at or under 2,000 progeny) are marked `ran`.

## Shinylive (Chromium)

| case | gz MB | to #load-run s | uploads s | run->status s | memory before / after / delta MiB | outcome |
| ---- | ----- | --------------- | --------- | -------------- | ----------------------------------- | ------- |
| 6k_x_2000 | pending | pending | pending | pending | pending | pending |
| 50k_x_200 | pending | pending | pending | pending | pending | pending |
| 50k_x_2000 | pending | pending | pending | pending | pending | pending |

## Method

- **CPython memory** (`parse peak MiB`, `analysis peak MiB`) is the `tracemalloc` peak
  (`tracemalloc.get_traced_memory()[1]`) taken around each stage separately, with
  `tracemalloc.start()`/`tracemalloc.stop()` bracketing only that stage; numpy reports its
  allocations to `tracemalloc`, so the figure includes numpy's own arrays, but each stage's peak
  excludes whatever was already live from a previous stage (the loaded `Dataset`, for example, is
  not counted again during `analysis peak MiB`).
- **Browser memory** (`memory before / after / delta MiB`) is
  `performance.measureUserAgentSpecificMemory()`, which requires cross-origin isolation
  (`crossOriginIsolated`); it includes the Pyodide worker's heap. `scripts/bench_shinylive.py`
  serves `site/` with `Cross-Origin-Opener-Policy: same-origin` and
  `Cross-Origin-Embedder-Policy: require-corp` and asserts `crossOriginIsolated` before reading
  a memory figure, exiting 3 (Q4, docs/m3-phases.md) rather than printing a fabricated figure when
  isolation is unavailable.
- **`parse s plain / gz`**: `parse s` (plain) is timed with `tracemalloc` running (it brackets the
  same stage as `parse peak MiB`); `parse s (gz)` is a second, untimed-for-memory
  `load_dataset` call on the `.gz` input with `tracemalloc` not running, so the two wall-clock
  figures are not on an identical footing with respect to profiling overhead (docs/m3-phases.md,
  Q15: the file is decompressed a second time rather than growing a buffer).
- **`total s`**: the sum of `parse s` (plain), `analysis s`, `duplicate scan s` and `export s`; it
  does not include `parse s (gz)`, which is a second, separate parse of the same data.
- **Units**: file sizes (`VCF MB plain / gz`, `gz MB`) are decimal (1 MB = 1,000,000 bytes);
  memory figures (`... peak MiB`, `resident matrix MiB`, the browser memory columns) are binary
  (1 MiB = 1,048,576 bytes).
- **`duplicate scan s (ran)`**: `core/pipeline.py` skips `duplicate_pairs` above 2,000 progeny
  (`_MAX_DUPLICATE_SAMPLES = 2000`); at exactly 2,000 progeny the scan runs and this column reads
  `ran`, so the quadratic scan is inside the measured `total s`. Above 2,000 progeny the column
  reads `skipped` and `scripts/bench_pipeline.py` prints `duplicate scan: skipped`.
- **`resident matrix MiB`**: `dataset.genotypes.calls.nbytes` after parsing, converted to MiB —
  the size of the retained `int8` genotype array, not a peak.

### What a user can expect

From the two cases measured so far, three things the tables above make plain. The third row and
the Shinylive rows are still pending, so nothing here is extrapolated to 50,000 markers by
2,000 individuals.

- **Analysis memory is about eighteen times the genotype matrix, and parsing is about two to
  three times it.** At 6,000 x 2,000 the matrix is 22.9 MiB resident and the analysis peaks at
  412.7 MiB; at 50,000 x 200, 19.3 MiB against 345.2 MiB. The ratio is the float64 temporaries
  in `core`, not the parser: `rpp` alone builds several arrays of the full marker-by-sample
  shape at 32 bytes per call against the matrix's 2 (docs/m3-phases.md, Q1). This is the
  measurement the chunked-core phase is judged against.
- **The duplicate scan dominates whenever it runs.** At exactly 2,000 progeny it accounts for
  517 s of the 556 s analysis. It runs at 2,000 and is skipped at 2,001 (docs/adr/0017), so the
  cost appears abruptly at the boundary rather than growing into it: a run of 2,001 individuals
  is faster than a run of 2,000.
- **Parsing is the single largest cost at these sizes.** 262 s for a 48.3 MB plain VCF and 218 s
  for 42.2 MB, roughly 50,000 calls per second either way, against 24 s for the same data
  gzipped — the `.gz` figure is lower because the file is smaller to read, not because parsing
  is faster.

### Not measured

- Firefox, Safari (Chromium only, via Playwright).
- A run in which individuals pass the filters. `n_pass` is 0 in both generated cases: progeny
  states are drawn independently per marker, so a het call at every marker of a multi-marker
  target region is vanishingly unlikely. Everything these figures measure — parsing,
  classification, RPP, drag, IBS, the duplicate scan, ranking and the export — runs for every
  individual whether or not it passes, so the totals stand as written; what the rows do not
  exercise is a long selection list or a large next-round manifest.
- An unsorted VCF or a manifest in a different sample order from the genotype columns: each
  triggers an extra copy of the matrix (`GenotypeMatrix.sorted_by_position` or `select_samples`)
  that these figures do not include, since the generated bench inputs are already sorted and
  ordered.

Not verifiable here (docs/m3-phases.md): whether Shinylive loads under `require-corp` (Q4); the
wasm heap ceiling of the Pyodide build — the outcome column of the Shinylive table records what
happens.

Figures: not verifiable here — filled by the phase run.
