# results.csv is a versioned schema: NA for missing, a header when empty, metadata columns

Status: accepted. Date: 2026-09-19.

## Context and Problem Statement

results.csv is the file a breeder opens in a spreadsheet, reads with pandas or readr, and archives beside the genotype data. Until now a missing value was an empty cell, a run with no passing individuals produced an empty file with no header, and nothing in the file said which background model, which unit, which ranking mode or which assembly produced the numbers. A float that was `NaN` in a column the pipeline did not route through `_num` (`het_rate`, `expected_het`, `expected_rpp`) printed the literal `nan`, which pandas reads as missing, readr reads as the string `"nan"` and a spreadsheet leaves as text. Reproducing a ranking from the file alone was not possible, and a downstream reader could not tell a pre-freeze file from a later one.

## Considered Options

1. Leave the format as it is and document the empty cell.
2. Empty cells for missing values, plus the metadata columns.
3. `NA` for missing values, a header row even with no results, metadata columns, and a schema version written on every row.

## Decision Outcome

Option 3.

`NA` is what the sibling isoline-browser writes, what readr reads as missing by default (`na = c("", "NA")`), and what pandas reads as missing; an empty cell is ambiguous between "missing" and "the empty string". The two are therefore given different meanings here: a missing number, rank or identifier is `NA`, and text that is genuinely empty (`exclusion_reason`, `qc_flags`, `notes`) stays an empty cell. The literal `nan` never appears: `het_rate`, `expected_het` and `expected_rpp` now go through the same `_num` as the other floats, and the writer maps `None` and `NaN` alike to `NA`.

The header is written even with no rows, so a run that excludes everyone, or a download taken before any data is loaded, still produces a readable file rather than zero bytes. With no rows the column list cannot be derived from the data, so it is the fixed prefix: the 34 documented columns plus `token_profile` (contract 1.4.0), which keeps its place as the last column when the dynamic per-locus and per-chromosome columns are present.

Five metadata columns, written on every row rather than in a header comment, because a comment line breaks every CSV reader this file is meant to survive: `background_model`, `background_unit`, `rank_mode` (docs/adr/0007's amendment), `assembly` (docs/adr/0015) and `results_schema`. They follow `frac_b`, so no documented leading column moves. `results_schema` is SemVer and starts at `1.0.0`: adding a column bumps the minor version, removing or renaming one bumps the major. A reader that accepts `1.x` therefore keeps working when a column is added.

selected.csv carries `results_schema` too, and writes `NA` for a missing `family_id` or `generation` and for the rank or score of an individual that has none. next_samples.csv does not: it is the input contract's samples.csv (docs/adr/0010), where `generation` and `family_id` default to empty, and isoline-browser's manifest reader would take `NA` for a family named "NA". Two output files with `NA`, one without, is the price of the round trip.

### Versioning rules and column placement

A new fixed column goes after `results_schema` and before `token_profile`, which stays last; a reader that takes the first 34 columns by position therefore keeps working across a minor bump, and one that takes the last column keeps finding the token profile. Adding a fixed column, a new dynamic family (a new `<prefix>_<id>_<metric>` group) or a column to selected.csv is a minor bump, because a reader that selects by name or by the documented prefix is unaffected. A new value in an enumerated column (`background_model`, `background_unit`, `rank_mode`, `assembly`, `drag_unit`, a status or a QC flag) is not a schema change at all and does not bump the version: readers must treat these as open vocabularies. Removing or renaming any column, changing a column's type, or changing what `NA` means is a major bump.

The prefix `rpp_` is reserved: only `rpp_total`, `rpp_carrier`, `rpp_noncarrier` and the per-chromosome columns use it, which is why the unit column is `background_unit` and not `rpp_unit`. It pairs with `background_model`, and it leaves `^rpp_` minus the three summaries as a usable selector for the per-chromosome columns, which `scripts/read_results.R` relies on. Chromosome names are not controlled by this tool (`normalize_chrom` keeps an unrecognised name verbatim), so a chromosome that would produce one of the three summary names is rejected with a `DataContractError` naming the chromosome, rather than silently overwriting a metadata cell.

`NA` is meaningful only where a value can be missing. `cli.py`'s `_read_results` therefore maps `NA` to a missing value only in the numeric columns and in `family_id` and `generation`; an identifier or free-text cell is read as written, so a line genuinely named `NA` survives `rank` -> `select` -> next_samples.csv. The same ambiguity is not solvable for a free-text note whose text is exactly `NA`: quoting it would not help (pandas reads a quoted `NA` as missing, and `csv.writer` cannot quote a single field without quoting every field), so docs/data-formats.md records the caveat instead.

### Test coverage

The synthetic fixture does not pin the metadata columns, on purpose: they are configuration passed through from criteria.yaml and the resolved unit, not quantities the fixture's independent implementation computes. Their coverage is in `tests/test_results_schema.py` (the defaults on the fixture, the `bp` fallback for `background_unit`, `NA` for a NaN or absent QC number, the chromosome-name collision) and in phase 6 for the weighted model and the staged rank mode.

### Consequences

Good: the file says how it was produced; a reader can branch on the schema version; missing values are unambiguous in R, pandas and a spreadsheet; an empty run is still a valid CSV. Bad: existing scripts that test a cell for `""` to mean missing must test for `"NA"`; the columns are wider; the format change is breaking and is recorded as such in CHANGELOG.md. `cli.py`'s `_read_results` reads both forms, so a pre-freeze results.csv still feeds `select`.
