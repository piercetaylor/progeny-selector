# Streaming VCF parse: two passes over the file, one preallocated matrix

Status: accepted. Date: 2026-09-22.

## Context and Problem Statement

`read_vcf` held a list of per-record `(n_samples, 2)` int8 arrays and then `np.stack`ed them,
so at the moment of the stack both the list and the finished matrix were resident: a peak of
about twice the genotype matrix (measured on a 2,000 x 1,000 VCF: 9,417,127 bytes against a
4,000,000-byte matrix, a ratio of 2.35). The milestone's scope line (PLAN.md:131) asks the
parser to bound its memory to the size of what it returns.

## Decision Drivers

- The peak must be the matrix plus a small constant, not a multiple of it.
- Every contract rule must stay identical, and every message and physical line number must stay
  identical *for a file that can be read to the end*: the shared contract cases
  (`contract/cases/err-vcf-*`) and `tests/test_io.py` pin them. The narrow form is deliberate; the
  two places where the wider claim fails (an unreadable stream, and an empty-versus-headerless
  file) are named in the Decision Outcome and the Consequences rather than glossed.
- `core/` must not change. Every core entry point (`classify`, `rpp`, `donor_segment`,
  `ibs_to_sample`, `pairwise_ibs`) indexes whole arrays.

## Considered Options

1. **List of per-record arrays, then `np.stack`** (the code as it stood): one pass, simple,
   peak about 2 x the matrix.
2. **A growth buffer**, as the sibling's `GenotypeBuilder` (`../backcross/docs/adr/0012`):
   one pass, a doubling or 1.5 x growth policy, peak about 1.57 x the inflated bytes, with a
   retained over-allocation to trim or keep.
3. **Two passes with preallocation**: pass 1 counts data lines, pass 2 fills
   `np.empty((n_records, n_samples, 2), dtype=np.int8)`.

## Decision Outcome

**Option 3, two passes** (`docs/m3-phases.md`, Q15).

- **The bound is the matrix plus one line plus the marker tables.** `_scan` tokenises nothing
  and keeps only the sample ids and a count; `_fill` holds one line's fields, the `Marker` list
  and the per-marker allele lists (about 15 MB at 50K markers). Measured after the change on the
  same 2,000 x 1,000 VCF: 4,854,579 bytes, a ratio of 1.21. There is no growth policy and no
  retained over-allocation to reason about; the allocation is exactly the size of the result.
- **`Dataset` keeps a whole int8 `GenotypeMatrix` and `core/` is untouched.** A chunked
  representation would rewrite every core entry point for a gain the parser does not need.
  Analysis memory is a separate decision, and the larger one: `core/background.py::rpp` builds
  four float64 `(m, s)` temporaries, 32 bytes per call against the matrix's 2, so about 3.2 GB
  at 100 million calls, against the parser's 400 MB at 50K x 2,000 (`docs/m3-phases.md`, Q1);
  that is Phase 6 and ADR 0025, not this record.
- **A `.gz` is inflated twice.** The cost is one extra decompression, recorded in the
  `parse s (gz)` column of `docs/limits.md` rather than hidden.
- **Pass 1 judges nothing.** Not the header's shape either: `_scan` takes
  `line.rstrip("\r\n").split("\t")[9:]` as the sample ids and moves on. The shape check, the
  `#`-after-header check, the data-line-before-header check and the second-`#CHROM` check all sit
  in pass 2, at the physical line they always sat at, so precedence between errors *inside* the
  file is exactly what it was — including the reviewer's case of a data row followed by a `#CHROM`
  line with no FORMAT column, which must still say "VCF data line before #CHROM header" rather
  than blaming the header two lines later. The shared cases `err-vcf-hash-line-after-header` and
  `err-vcf-second-header` check the row-error half of this; `tests/test_vcf_streaming.py`
  `::test_a_malformed_header_never_outranks_an_earlier_error` checks the header half.
  `_scan`'s "header seen" state is an explicit flag, never the truthiness of `sample_ids`: a
  malformed `#CHROM` line yields an empty list, and a list-truthiness sentinel would stop counting
  every record after it and hand pass 2 a short matrix for a file that must be refused instead.
  A malformed header therefore reaches `np.empty((n, 0, 2))`, which costs nothing and which `_fill`
  refuses before indexing.
  The one thing precedence does *not* survive is an unreadable stream: see the Consequences.
- **A per-call GT cache.** `_fill` keys a `dict[str, tuple[int, int]]` on the GT sub-field alone
  (`token.split(":", 1)[0]`), so a `GT:DP` file with thousands of distinct depths cannot grow it
  beyond the number of distinct genotypes; `_parse_gt` itself is unchanged, and a half-missing
  `./1` still reads `(-1, 1)` (Q13: the BrAPI reader's `(-1, -1)` is equivalent downstream, and
  changing this would move fixture bytes for no observable gain).
- **Deviation from `docs/m3-phases.md`, Phase 3, recorded deliberately: the emptiness check runs
  after pass 2, not before the allocation.** The spec places `if n_records == 0: raise
  DataContractError("VCF contains no variant records")` "before any allocation". Taken literally,
  a VCF whose records precede its `#CHROM` header (or that has no header at all) counts zero
  records and is told it is empty, which is a false statement about the file and would send a
  breeder looking for missing data instead of a missing header; the pre-change reader said "VCF
  data line before #CHROM header" at that line. So `calls = np.empty((n_records, len(sample_ids),
  2), dtype=np.int8)` is allocated unconditionally — at `n_records == 0` that is a zero-row array
  and costs nothing — `_fill` always runs and raises the structural message at its physical line,
  and the emptiness check runs on the way out, where only a genuinely header-only file reaches it.
  The spec's intent, never allocating a full-size matrix for a bad file, is honoured exactly.
  `_fill` reaches "VCF data line before #CHROM header" before it uses `len(sample_ids)`, so the
  empty sample-id list of a headerless file is harmless; that is pinned by
  `test_fill_judges_a_missing_header_before_it_uses_the_sample_ids`.
- **`EOFError` and `OSError` become a contract error naming the file; `FileNotFoundError` does
  not.** `EOFError`, `gzip.BadGzipFile` and any other `OSError` from either pass become
  `DataContractError(f"{path}: cannot read: {exc}")`, where a truncated `.gz` previously escaped as
  a raw `EOFError` through the CLI and the Load screen. The claim is scoped to those exception
  types and no wider: a Latin-1 encoded VCF still escapes as a raw `UnicodeDecodeError` from
  `open_text`, unchanged from the single-pass reader and not addressed here.
  `FileNotFoundError` alone is re-raised unchanged, because `cli.py:217` catches it by name
  alongside `DataContractError`; wrapping it would make that catch dead code, and a later cleanup
  would then turn a missing criteria file or token profile back into a traceback. This is a
  by-name exemption, not a rule about when the failure happens: a directory passed as a path
  raises `PermissionError` at open on Windows and *is* wrapped.
- **A file that changes between the passes is refused, whichever way it changes.** `filled >=
  calls.shape[0]` before each write catches growth, which would otherwise be an `IndexError`, and
  `filled != calls.shape[0]` after the loop catches shrink, which `GenotypeMatrix.__post_init__`
  would otherwise report as the internal "calls shape (7, 4, 2) != (6, 4, 2)"; both raise
  "VCF changed while it was being read". Both guards have a test that patches `_scan` to
  miscount, because a guard with no test is a guard a refactor deletes.
- **How the bound is measured.** The `tracemalloc` peak around the parse, as ADR 0021 defines the
  CPython memory figure and as the `parse peak MiB` column of `docs/limits.md` records it; numpy
  reports its own allocations to `tracemalloc`. The check in the suite is
  `tests/test_vcf_streaming.py::test_parse_peak_is_bounded_by_the_matrix`, which asserts a peak
  below 1.5 x the matrix and fails at 2.35 x against the pre-change reader.

## Consequences

Good: a VCF parses within the size of the matrix it produces, so the parser is no longer the
binding constraint on file size, and the measurement is a test rather than a claim. A truncated
download now reports which file could not be read.

Bad: the file is read twice, so a `.gz` is inflated twice and a plain file is walked twice; a
file that changes between the passes is refused ("VCF changed while it was being read") rather
than parsed as far as it goes.

Bad, and intrinsic to reading twice: **an unreadable stream now outranks every row error in the
same file.** A 400-record `.gz` truncated to two thirds, whose line 3 also has the wrong column
count, used to report `line 3: expected 11 columns, found 10`; it now reports
`cannot read: Compressed file ended before the end-of-stream marker was reached`, because pass 1
consumes the whole stream to count records before pass 2 judges any row. This cannot be fixed
without abandoning the two-pass design, so it is a cost of option 3 rather than a defect: the
phase's "every error message and line number unchanged" holds for errors inside a readable file,
and its two exceptions are this one and the emptiness check's move above. A user with both faults
sees the truncation first and the row error on the next attempt.

Bad, and a bound on how far "bounded" goes: the marker and allele tables are per record, not per
call, so a very tall, very narrow VCF is dominated by them rather than by the matrix. Measured at
20,000 records x 2 samples: a peak of 9.06 MB against an 80 KB matrix, about 450 bytes per record
for the `Marker` object, its id string and the allele list. The bound is the matrix plus those
tables, which is what makes millions of records cheap in the matrix and not free overall.

## Revisit when

- Phase 6 chunks the analysis: with `rpp` bounded, the parser's peak stops being a rounding
  error in the total and the 50K x 2,000 row of `docs/limits.md` becomes a measurement.
- A format arrives whose record count cannot be counted cheaply (a stream without a seekable
  source, for instance a BrAPI page), where the growth buffer of option 2 is the right answer
  and this record does not apply.
