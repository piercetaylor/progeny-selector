# Contract 1.8.0: a pair of two missing characters is read as missing

Status: accepted. Date: 2026-09-26. Mirrors backcross docs/adr/0023, the canonical record for this contract version.

## Context and Problem Statement

Contract 1.6.0 (docs/adr/0026 here, backcross docs/adr/0021) stated that a pair of one nucleotide and one of `N`, `-`, `.` is read as missing and left one cell out on purpose: a pair of two of `N`, `-`, `.` that is not itself a missing token (`N/N`, `N-`, `./N`, `..`, `-|N`). The contract said such a pair "is not defined by this version". Both tools have always read it as missing. Contract 1.8.0 states the reading. This record is the sibling half: what this repository already did, what it checks, and what did not change. The decision, the options weighed and the maintainer's criterion ("best for academic and open-source plant-breeding users: reproducible with standard tools, tolerant of spreadsheet, pandas and R exports, never silently wrong") are in backcross docs/adr/0023 and are not restated here.

## Decision Outcome

The canonical record's decision applies here unchanged: the pair is read as missing in HapMap and in a nucleotide wide CSV, in any order and in all three spellings; it is part of the pair grammar and not a missing token, so a pair with `X` that is not itself a missing token (`X/X`, `XN`, `N/X`, `X.`) stays an error while `XX` stays a HapMap missing token; `auto` detection sees such a cell that is not itself a missing token as a non-coded cell, so a wide CSV holding one reads as nucleotide (`NN`, a missing token, is skipped by detection as before); under a token profile whose `base` is `nucleotide` the rule applies below the profile's own tokens, and under `base: none` (`dart`, `axiom`, `kasp`) such a cell is `genotypes.unknown_cell` unless the profile lists that exact token, as `axiom` lists `--`.

- **No source change.** `parse_nucleotide_call` (`src/progeny_selector/io/calls.py`) already returned `None` for such a pair, for every caller: its pair grammar accepts two characters drawn from A, C, G, T, `N`, `-`, `.` and returns `None` when either is not a nucleotide, so two missing characters were covered by the same branch as the half-missing pair. Under `base: none` the profile branch raises before the grammar is reached. `detect_coding` skips only `WIDE_NUCLEOTIDE_MISSING`. Only the docstring changed, to name contract 1.8.0.
- **The rule is pinned by hand-built tests, not only by the shared cases.** `tests/test_io.py` asserts that all 27 cells (`N`, `-`, `.` paired with each other, in the two-character, slash and bar spellings) are missing under `HAPMAP_MISSING` and `WIDE_NUCLEOTIDE_MISSING`, upper and lower case, with no profile and under `tassel` and `soybase-report`; that every one raises the `base: none` message under `dart` and `kasp`, and under `axiom` all but `--`, which `axiom` lists; that `X/X`, `XN`, `N/X` and `X.` raise under both missing sets and `X/X` under `tassel`; and that a wide CSV whose only non-A/B cell is `N/N` is detected as nucleotide and then fails naming the coded letter `B`; that test pins detection only, since `B` is read before the pair, and the pair branch is pinned by the 27-cell test. Each assertion over the 23 spellings that are not missing tokens (`NN`, `--`, `./.` and `.|.` are, and return `None` before the grammar) fails if the pair branch raised, if `base: none` fell through to the grammar, or if `X` joined the pair characters.
- **The contract cases come from the sibling unchanged.** `contract/cases/hapmap-two-missing-pair`, `wide-two-missing-pair`, `profile-tassel-two-missing-pair` and `err-profile-none-two-missing-pair` are generated in backcross and run here through `tests/test_contract_cases.py`; `scripts/check_contract.py ../backcross` byte-compares the mirror.
- **`docs/data-formats.md`** moves its version reference to 1.8.0 and states the pair reading in its summary of the cell vocabulary.

## Consequences

No input changes meaning, no output changes, no error kind is added, and `results_schema` does not move. A file carrying `N/N` or `..` loads exactly as before; a future reading which rejected it would now be a major bump. As backcross docs/adr/0023 records, `-/-` (a homozygous deletion in TASSEL's reading) and `..` are both read as missing, so a deletion call in that spelling is discarded rather than carried.

## Revisit when

- A platform is met whose `-` is a real deletion allele worth carrying rather than discarding, which would be a token profile question, not a change to this rule.
- backcross docs/adr/0023 is amended or superseded.
