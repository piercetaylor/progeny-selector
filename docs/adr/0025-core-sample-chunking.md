# Core analysis in blocks of sample columns

Status: accepted. Date: 2026-09-23.

## Context and Problem Statement

`classify`, `rpp` and `ibs_to_sample` each built full `(n_markers, n_samples)` temporaries and
then reduced them along the marker axis. The genotype matrix itself is int8, two bytes per call.
A float64 temporary of the same shape is eight bytes per call, and `rpp` held four of them at
once: `rpp_contribution`'s output, `np.nan_to_num` of it, `w2`, and the product. That is 32 bytes
per call against the matrix's 2 (`docs/m3-phases.md`, Q1). At 50,000 markers by 2,000
individuals, 100 million calls, the matrix is 200 MB and `rpp` alone asks for about 3.2 GB.

ADR 0022 bounded the VCF parser to the size of the matrix it returns and named analysis memory as
the larger, separate decision. It is the binding constraint on the 50K x 2,000 row of
`docs/limits.md`, which was left as `pending` because the analysis, not the parse, is what a
16 GiB machine could not finish.

Measured on a 4,000 x 1,000 matrix (`tests/test_core_memory.py`), the pre-change `tracemalloc`
peaks were 136,017,736 bytes for `rpp` (34.0 x the 4 MB state matrix), 76,002,696 bytes for
`ibs_to_sample` (9.5 x the 8 MB call matrix) and 28,066,408 bytes for `classify` (3.5 x).

## Decision Drivers

- Analysis peak at 50K x 2,000 must be a few times the resident matrix.
- Output must be bit-identical. Every figure in `tests/fixtures/` and every written RPP, IBS and
  state in `docs/` was produced by the old code, and a change of a single ulp in RPP can reorder
  two candidates whose composite scores are close.
- The change must stay inside the three functions that hold float64 temporaries. Rewriting the
  whole core against a chunked representation would touch `pipeline.py`, `drag.py`, `qc.py` and
  `selection.py` for no measured gain.

## Considered Options

1. **Block over sample columns inside each function.** The signatures, the callers and the return
   shapes are unchanged; only the loop is new.
2. **Reduce in float32.** Halves every temporary and changes every result, so the fixture and the
   sibling's exchanged files move. Rejected on bit-identity.
3. **A chunked `GenotypeMatrix`** that yields column blocks, with every core entry point rewritten
   against it. This is option 1's benefit plus a rewrite of code that is already bounded.

## Decision Outcome

**Option 1, blocks of `CORE_SAMPLE_CHUNK = 64` sample columns** (`docs/m3-phases.md`, Phase 6,
Q1 answered yes).

- **Why 64.** A float64 `(50_000, 64)` temporary is 25.6 MB, so the three or four a block holds
  at once stay under 100 MB at 50K markers, independent of how many individuals are in the file.
  The block is wide enough that the per-block Python overhead is amortized: at 2,000 individuals
  the loop runs 32 times, and the marker axis inside each block is 50,000 long, so numpy still
  does the work in long vectorized reductions. A narrower block would lower the constant and raise
  the iteration count; 64 was chosen because it puts the bound where the milestone asked for it
  with a loop short enough to be irrelevant to run time.
- **Which functions change.** `classify` blocks the `a1`/`a2` comparison and the state assignment,
  assigning H, then B, then A, then N into a view of the output block, with
  `states[~informative, :] = STATE_U` after the loop. `rpp` folds the marker mask into the weight
  vector (`w = np.where(marker_mask, w, 0.0)`) and accumulates `numer` and `denom` per block, with
  the NaN contributions replaced once at import time by the read-only `RPP_LOOKUP = np.nan_to_num(
  RPP_CONTRIBUTION)`. `ibs_to_sample` blocks the `called`, `shared` and `denom` arrays.
- **The weight mask inside a block is `np.where`, not multiplication by the boolean.**
  `w[:, None] * counted` evaluates `inf * False` as NaN, which the sum then propagates to every
  individual, where the unchunked `np.where(counted, w[:, None], 0.0)` produced NaN only in the
  columns the infinite weight actually entered. This is reachable: `marker_weights` returns `cap`
  for a lone marker on a chromosome, and `cap` is infinite when a caller passes
  `max_coverage=inf` with no chromosome lengths. The chunked code keeps the `np.where` and keeps
  `np.nansum` for the numerator, so the infinite-weight result is the pre-change result.
  `test_rpp_with_an_infinite_weight` pins it.
- **A one-column tail is folded into the block before it.** `constants.sample_blocks` returns
  65-wide last blocks at `n_samples` 65, 129, 193 and 2049 in place of 64 plus 1. numpy reduces a
  single-column array along axis 0 by pairwise summation and a wider one by accumulating rows in
  order, so a width-1 final block gave the last individual a weighted RPP differing from the
  unchunked value by up to 3.3e-16 at 3,000 markers and 4.7e-15 at 50,000. A sweep of remainders 0
  through 8 found remainder 1 to be the only width that broke. The unweighted case was exact at
  every width because its summands are 0, 0.5 and 1. All three functions call the same helper, so
  the rule cannot drift between them.
- **Which functions do not change, and why.** `pipeline.py`'s `frac_a`, `frac_h`, `frac_b` and the
  `counted` masks are boolean, one byte per call, reduced immediately; their peak is half the
  matrix. `drag.donor_segment` already works on one chromosome's markers at a time, a twentieth of
  the matrix in soybean. `qc.sample_qc` walks one sample column per iteration. `pairwise_ibs`
  subsamples to at most 2,000 markers before it allocates anything. None of them holds a float64
  array of the full shape, so none of them was touched.
- **Results are bit-identical, and the test is the evidence.** The per-column terms are the same
  numbers in both versions: `w`, `w/2` and `0` for `rpp`, and 1.0, 0.5 or 0.0 shared-allele
  fractions for `ibs_to_sample`. numpy reduces a C-contiguous two-dimensional array of two or more
  columns along axis 0 by accumulating rows into the output buffer in row order, so for every
  block width `sample_blocks` produces the additions happen in the order the unchunked code used.
  A one-column block is the exception, and the tail rule above exists to avoid emitting one.
  `tests/test_core_memory.py` holds the pre-change bodies of all three functions verbatim as
  `_rpp_reference`, `_ibs_reference` and `_classify_reference`, and asserts
  `np.array_equal(..., equal_nan=True)` against them at `n_samples` 64, 65, 128, 129, 193, 700 and
  2049, weighted, masked and unweighted. The all-false-mask, all-zero-weight, infinite-weight,
  zero-sample, uninformative-marker and all-missing-column cases are asserted equal as well, and a
  120 x 40 matrix covers a file with fewer individuals than one block.
- **How the bound is measured.** The `tracemalloc` peak of one call, as ADR 0021 defines the
  CPython memory figure, with the inputs constructed before `tracemalloc.start()` and
  `reset_peak()` called after a warm-up so that neither the inputs nor numpy's first-call
  allocations enter the peak. The same 4,000 x 1,000 matrix gives 6,674,928 bytes for `rpp`,
  6,929,132 bytes for `ibs_to_sample` and 6,115,376 bytes for `classify`; against that matrix
  those are 1.67 x, 0.87 x and 0.76 x. The suite asserts `< 4 x`, `< 4 x` and `< 3 x`, thresholds
  the pre-change code failed at 34.0 x, 9.5 x and 3.5 x on the same inputs.
- **The bound is absolute, and the ratio is not.** The block temporaries are
  `(n_markers, 64)` float64 and do not depend on the number of individuals, so the peak is set by
  the marker count alone: 6.7 MB for `rpp` at 4,000 markers and 83.2 MB at 50,000, measured at both
  100 and 2,000 individuals. Expressed against the state matrix the same 6.7 MB reads as 16.65 x at
  4,000 x 100 and 1.67 x at 4,000 x 1,000. Quote the megabytes and the marker count; a ratio
  belongs only with the shape it was measured at.

## Consequences

Good: `rpp` holds 83.2 MB of block temporaries at 50,000 markers by 2,000 individuals, on top of
the 200 MB matrix, and the same 83.2 MB at 200. The pre-change figure was a multiple of
the call count, about 3.2 GB at that shape, which is why the 50K x 2,000 row of `docs/limits.md`
stood at `pending`. The `analysis peak MiB` figure for that case is written in Phase 9 against the
acceptance criterion of at or below 4 x `resident matrix MiB`, a ratio that holds at 2,000
individuals and not at 100.

Bad: three functions now contain a loop where they contained one expression, and the block width
is a constant a reader has to find in `constants.py` to understand the arithmetic. The bit-identity
claim rests on numpy's axis-0 reduction order, an implementation property of the library that its
documented interface does not promise. The width-1 case above is one instance of that property
already biting; if a future numpy version applied pairwise summation across the marker axis at
other widths, the equality tests would fail at the ulp.

Follow-up, not this phase: `classify.rpp_contribution` is now unused in `src/` and holds a second
copy of the 1 / 0.5 / 0 table that `RPP_LOOKUP` holds. Two tables that must agree and no caller to
keep them honest is a drift risk. Removing a public function that `__all__` exports is a separate
decision; `RPP_LOOKUP` is marked read-only (`setflags(write=False)`) so that at least one of the
two cannot be mutated in place by a caller.

Not addressed: the peak scales with the marker count, because `RPP_LOOKUP[c]` and the
`shared` array are `(n_markers, 64)`. A file with 500,000 markers would hold 832 MB of blocks. The
marker axis is chunkable by the same method if that case arrives; no run in `docs/limits.md`
needs it.

## Revisit when

- A run exceeds 100,000 markers, where the per-block temporaries pass 50 MB each and the marker
  axis becomes the one worth blocking.
- `pipeline.py` grows a float64 per-call metric, which would put it back in the class of functions
  this record covers.
