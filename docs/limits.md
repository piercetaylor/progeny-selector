# Measured limits

Measured on 2026-09-22 before streaming (commit `c3e54ab`, the Phase 2 rows) and on 2026-09-25 and 2026-09-26 after M3 Phases 3 and 6 (CPython at commit `674a2ff`, Shinylive at commit `829f247`).

Machine: Intel(R) Core(TM) Ultra 5 125U; 16,573,128,704 bytes (about 15.4 GiB) RAM; Windows 11 Home (build 10.0.26200).

Python: 3.12.4 (`MSC v.1940 64 bit (AMD64)`).
numpy: 2.5.3.
shiny: 1.7.0.
shinylive: 0.8.11.
Playwright: 1.62.0.
Chromium: 151.0.7922.34, the full browser in new headless mode (`channel="chromium"`).

## CPython

After streaming (Phase 3, docs/adr/0022) and sample-chunked analysis (Phase 6, docs/adr/0025), commit `674a2ff`:

| case | markers | progeny | VCF MB plain / gz | parse s plain / gz | parse peak MiB | resident matrix MiB | analysis s | analysis peak MiB | duplicate scan s (ran) | export s | total s | outcome |
| ---- | ------- | ------- | ------------------ | -------------------- | ---------------- | ---------------------- | ----------- | -------------------- | ------------------------ | -------- | ------- | ------- |
| 6k_x_2000 | 6000 | 2000 | 48.3 / 2.3 | 99.84 / 6.20 | 29.1 | 22.9 | 305.95 | 201.7 | 204.50 (ran) | 0.11 | 405.90 | ok |
| 50k_x_200 | 50000 | 200 | 42.2 / 2.4 | 65.81 / 5.75 | 69.1 | 19.3 | 16.08 | 103.0 | 2.11 (ran) | 0.01 | 81.90 | ok |
| 50k_x_2000 | 50000 | 2000 | 402.2 / 19.0 | 580.58 / 43.72 | 241.0 | 190.9 | 314.70 | 481.0 | 219.23 (ran) | 0.12 | 895.40 | ok |

Before streaming, commit `c3e54ab` (Phase 2):

| case | markers | progeny | VCF MB plain / gz | parse s plain / gz | parse peak MiB | resident matrix MiB | analysis s | analysis peak MiB | duplicate scan s (ran) | export s | total s | outcome |
| ---- | ------- | ------- | ------------------ | -------------------- | ---------------- | ---------------------- | ----------- | -------------------- | ------------------------ | -------- | ------- | ------- |
| 6k_x_2000 | 6000 | 2000 | 48.3 / 2.3 | 262.29 / 24.04 | 49.7 | 22.9 | 555.97 | 412.7 | 517.28 (ran) | 0.26 | 818.52 | ok |
| 50k_x_200 | 50000 | 200 | 42.2 / 2.4 | 218.01 / 22.66 | 69.1 | 19.3 | 33.28 | 345.2 | 4.70 (ran) | 0.02 | 251.31 | ok |

The Phase 2 rows first printed totals of 1335.80 s and 256.01 s. Those added the duplicate scan to a total that already contained it, because `run_analysis` runs the scan (`core/pipeline.py`) and `scripts/bench_pipeline.py` then timed it a second time on its own. Every total above is `parse s` (plain) + `analysis s` + `export s`, recomputed from each case's `bench.json`, and the script now computes it that way.

`50k_x_2000` was not measured before Phase 6: the un-chunked `rpp` alone would have allocated roughly 3.2 GB of float64 temporaries at this size, by arithmetic from its 32 bytes per call.

## Shinylive (Chromium)

Commit `829f247`, `site/` rebuilt from it:

| case | gz MB | to #load-run s | uploads s | run->status s | memory before / after / delta MiB | outcome |
| ---- | ----- | --------------- | --------- | -------------- | ----------------------------------- | ------- |
| 6k_x_2000 | 2.3 | 5.31 | 0.09 | 136.44 | 150.1 / 400.8 / 250.7 | ok |
| 50k_x_200 | 2.4 | 5.25 | 0.16 | 12.61 | 150.1 / 340.4 / 190.3 | ok |
| 50k_x_2000 | 19.0 | 4.83 | 0.34 | 225.84 | 150.1 / 870.3 / 720.2 | ok |

All three cases load and analyze in the browser, including 50,000 markers by 2,000 individuals. Two defects stopped the first attempts, and both were fixed before these rows were taken. `scripts/bench_shinylive.py` launched `chromium-headless-shell`, whose `measureUserAgentSpecificMemory()` throws `SecurityError` on a cross-origin-isolated page; it now launches the full browser. The Load screen could not read a `.vcf.gz` upload at all, because Shiny stores an upload under its last suffix only (`0.gz`); it now restores the upload's own name.

## Acceptance

Phase 3 (docs/adr/0022): the 50K x 2,000 parse peak is 241.0 MiB, below the matrix plus 64 MiB (190.9 + 64 = 254.9 MiB).

Phase 6 (docs/adr/0025): the 50K x 2,000 analysis peak is 481.0 MiB, 2.52 times the 190.9 MiB matrix, within the 4x threshold.

## Measurement conditions

Every run shared the machine with a WSL virtual machine that grew from 1.9 GB to 4.5 GB over the two days and was left running. The logs of available memory and page-ins, sampled every 30 s, are the basis for how far each row's timings can be trusted. The memory columns are unaffected either way, because `tracemalloc` counts allocations and not resident pages.

- CPython `50k_x_2000` (2026-09-25): available memory 1.8 to 2.4 GB; page-ins had a median of 32 per second and exceeded 500 per second in 3 of 38 samples, two of them while the gzip copy was being written before any timed stage. The timings stand.
- CPython `6k_x_2000` and `50k_x_200` (2026-09-26): available memory fell to 475 MB; page-ins had a median of 575 per second, exceeded 500 per second in 13 of 23 samples and peaked at 110,396. Part of these rows' wall-clock time is paging. Their timings show the direction of the change against the Phase 2 rows and not its size.
- Shinylive (2026-09-26): available memory 1.4 to 2.3 GB; page-ins had a median of 17 per second, with one sample above 500. The timings stand.

## Method

- **CPython memory** (`parse peak MiB`, `analysis peak MiB`) is the `tracemalloc` peak (`tracemalloc.get_traced_memory()[1]`) taken around each stage separately, with `tracemalloc.start()`/`tracemalloc.stop()` bracketing only that stage. numpy reports its allocations to `tracemalloc`, so the figure includes numpy's own arrays, but each stage's peak excludes whatever was already live from a previous stage (the loaded `Dataset` is not counted again during `analysis peak MiB`).
- **`parse peak MiB`** brackets `load_dataset`, which reads the VCF, then markers.csv, then samples.csv. It is not the peak of `read_vcf` alone. At 50,000 x 200, `read_vcf` alone peaks at 38.4 MiB and `read_markers` alone at 31.9 MiB, and the column reads 69.1 MiB because the marker table is read while the matrix is held; without markers.csv, `load_dataset` also peaks at 38.4 MiB. The column is therefore unchanged at that shape before and after streaming, although the parser's own peak fell from 68.4 to 38.4 MiB (docs/adr/0022).
- **`analysis s` and `analysis peak MiB`** bracket `run_analysis`, which includes the duplicate scan whenever it runs. At 2,000 individuals the scan alone peaks at 175.6 MiB, measured separately; at 200 individuals, 14.1 MiB.
- **`duplicate scan s (ran)`** is the scan timed a second time on its own, after `run_analysis`, so its share of `analysis s` can be read off. It is not added to `total s`. `core/pipeline.py` skips the scan above 2,000 individuals (`_MAX_DUPLICATE_SAMPLES = 2000`); at exactly 2,000 it runs and the column reads `ran`. A 400-marker, 2,001-progeny run printed `duplicate scan: skipped` and a 400-marker, 2,000-progeny run read `87.91 (ran)` (2026-09-22).
- **`parse s plain / gz`**: the plain parse is timed with `tracemalloc` running, since it brackets the same stage as `parse peak MiB`; the `.gz` parse is a second `load_dataset` call with `tracemalloc` stopped. The two are not on the same footing. Every other CPython timing also runs under `tracemalloc`.
- **`total s`**: `parse s` (plain) + `analysis s` + `export s`. It excludes the `.gz` parse and the repeated duplicate scan.
- **`resident matrix MiB`**: `dataset.genotypes.calls.nbytes` after parsing, the size of the retained `int8` genotype array and not a peak.
- **Browser memory** (`memory before / after / delta MiB`) is `performance.measureUserAgentSpecificMemory()` before the uploads and after `#load-status` reports the load. It requires cross-origin isolation; `scripts/bench_shinylive.py` serves `site/` with `Cross-Origin-Opener-Policy: same-origin` and `Cross-Origin-Embedder-Policy: require-corp`, asserts `crossOriginIsolated`, and exits 3 when isolation is unavailable. The figure covers the page, including the Pyodide worker's WebAssembly heap, which grows and is never returned, so the delta approximates the load's high-water mark.
- **Browser timings**: `to #load-run s` is page load until the Load button is ready, `uploads s` the four file uploads, and `run->status s` pressing Load until `#load-status` reports the marker and progeny count, which spans parsing, analysis and the duplicate scan. The browser uploads the `.gz` file and runs without `tracemalloc`.
- **Units**: file sizes are decimal (1 MB = 1,000,000 bytes); memory figures are binary (1 MiB = 1,048,576 bytes).

### What a user can expect

- **Analysis memory at 50,000 markers is set by the marker count and the duplicate scan, not by the call count.** At 50,000 x 2,000 the analysis peaks at 481.0 MiB against a 190.9 MiB matrix. Before Phase 6 the chunked functions alone would have needed about 3.2 GB there. A ratio to the matrix depends on shape and is not a bound: 2.52x at 50,000 x 2,000, 5.3x at 50,000 x 200 and 8.8x at 6,000 x 2,000.
- **At 2,000 individuals the duplicate scan dominates analysis memory and time.** At 6,000 x 2,000 the scan alone peaks at 175.6 MiB of the 201.7 MiB analysis peak and takes 204.50 s of the 305.95 s analysis. The scan runs at 2,000 individuals and is skipped at 2,001 (docs/adr/0017), so its cost appears at that boundary: a run of 2,001 individuals is cheaper than a run of 2,000.
- **The browser handles every measured case.** At 50,000 x 2,000, Shinylive reaches a loaded and analyzed dataset in 226 s after pressing Load, and the page grows by 720.2 MiB. At 6,000 x 2,000 the same step takes 136 s, most of it the duplicate scan.
- **The CPython wall-clock figures overstate what a user waits.** They are taken under `tracemalloc`, which charges every allocation. The same 50,000 x 2,000 records parse in 580.58 s from the plain file with `tracemalloc` running and in 43.72 s from the `.gz` file without it; the browser, uninstrumented, finishes parsing and analysis in 225.84 s. Read the CPython timings as comparisons between rows and between commits, and the browser timings as the time a user waits.

### Not measured

- Firefox and Safari (Chromium only, via Playwright).
- A run in which individuals pass the filters. `n_pass` is 0 in every generated case: progeny states are drawn independently per marker, so a het call at every marker of a multi-marker target region is vanishingly unlikely. Everything these figures measure (parsing, classification, RPP, drag, IBS, the duplicate scan, ranking and the export) runs for every individual whether or not it passes, so the totals stand as written; the rows do not exercise a long selection list or a large next-round manifest.
- An unsorted VCF or a manifest in a different sample order from the genotype columns. Either triggers an extra copy of the matrix (`GenotypeMatrix.sorted_by_position` or `select_samples`) that these figures do not include, since the generated inputs are already sorted and ordered.
- The Pyodide build's WebAssembly heap ceiling. The largest case measured grew the page to 870.3 MiB and completed.
