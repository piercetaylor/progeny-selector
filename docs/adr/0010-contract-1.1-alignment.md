# Contract 1.1.0: mirror the shared contract and align the loaders with it

Status: accepted. Date: 2026-09-14.

## Context and Problem Statement

`docs/data-formats.md` was said to be shared verbatim with isoline-browser, but the two had drifted and nothing ran one repository's inputs through the other's parser. isoline-browser's ADR 0013 moved the shared input specification into a `contract/` directory with hand-authored cases, canonical there and mirrored here. Its ADR 0014 records the maintainer's decisions of 2026-09-14 on every input the two loaders read differently, including one nucleotide cell vocabulary for HapMap and wide CSV that is safe for any diploid crop, with crop-specific token profiles deferred to contract 1.2.0. What does this repository change to read every case identically?

## Decision Outcome

`contract/` is mirrored byte for byte from `../isoline-browser/contract` (its README gives the mirror rule) and `tests/test_contract_cases.py` loads every case through `load_dataset`. `scripts/check_contract.py` recomputes `MANIFEST.sha256` and, given the sibling's path, byte-compares the mirror with the canonical copy; it is a per-commit gate here, and `node scripts/check-contract-mirror.mjs ../progeny-selector` is the same check from the other side. The shared sections of `docs/data-formats.md` are replaced by a pointer to the contract version.

Loader changes, each from isoline-browser docs/adr/0014: samples.csv, markers.csv and the wide CSV sniff comma or tab from the header line and accept RFC 4180 quoting (`io/delimited.py`); `line_name` is optional; genotype columns absent from samples.csv are dropped and the loaded columns follow manifest order (`build_dataset`); synthesised coded parents are recorded in `Dataset.synthetic_sample_ids`, so the contract's sample list, which excludes them, can be reported; a VCF record with ID `.` is named `<CHROM>_<POS>` from CHROM as written; the chromosome pattern gains the `LG` prefix and the space and `-` separators and loses `ch`; non-soybean names sort in natural order (`chrom_sort_key`); `io/calls.py` is the one cell vocabulary, now mirrored by isoline-browser's `src/io/calls.ts`: the missing tokens are exactly the contract's per wide-CSV mode and for HapMap (`HAPMAP_MISSING` = the nucleotide list plus `X`, `XX`; `?` leaves every set), IUPAC codes still expand, and a single character outside A, C, G, T and the IUPAC codes now raises instead of reading as missing, so `?`, a stray `B` or `H`, `0`, `+` and `A?` are errors in both formats (shared cases `err-nucleotide-*`, `err-hapmap-*`); A/B/H detection scans every row instead of the first 200; every reader strips a leading UTF-8 byte-order mark. The wide-CSV fixed column order this reader always required is now the contract's rule.

### Consequences

Good: every contract case loads here with the same markers, samples and calls as in isoline-browser, and a future divergence is a failing test in the repository that diverged; the one vocabulary module makes the next token decision a one-line edit mirrored in both repositories. Bad: files that relied on `?` as missing, on the `ch6` spelling, on a single unknown letter reading as missing, or on genotype-file column order now load differently or fail with a message naming the line and cell. Open, listed in isoline-browser docs/adr/0014 and PLAN.md: a pair of one nucleotide and one of N, `-`, `.` (`AN`, `A-`), read as missing here and there; named token profiles and a crop selector in contract 1.2.0.

## More Information

isoline-browser docs/adr/0013 and 0014 and docs/input-coding.md; contract/README.md; docs/adr/0002.
