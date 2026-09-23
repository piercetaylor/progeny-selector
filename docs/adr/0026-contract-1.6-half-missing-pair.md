# Contract 1.6.0: a half-missing pair is read as missing

Status: accepted. Date: 2026-09-23. Mirrors backcross docs/adr/0021, the canonical record for this contract version.

## Context and Problem Statement

Contract 1.1.0 (docs/adr/0010 here, backcross docs/adr/0014) fixed the per-cell nucleotide vocabulary of HapMap and wide CSV but left one cell out: a pair of one nucleotide and one of `N`, `-`, `.`, in either order and in any of the three spellings (`AN`, `-A`, `A/N`). Both tools have always read it as missing; `docs/data-formats.md` said so and added "which the contract has not yet decided". Contract 1.6.0 states the reading. This record is the sibling half: what this repository already did, what it checks, and what did not change.

## Decision Outcome

The canonical record's decision applies here unchanged and is not restated: the pair is read as missing in HapMap and in a nucleotide wide CSV, in either order and in all three spellings; it is part of the pair grammar and not a missing token, so a pair of a nucleotide with any other character (`AX`, `A?`) stays an error and `auto` detection still sees a non-missing cell; under a token profile whose `base` is `nucleotide` the rule applies below the profile's own tokens, and under `base: none` such a cell is `genotypes.unknown_cell` unless the profile lists that exact token. A pair of two of `N`, `-`, `.` that is not itself a missing token (`N/N`, `N-`) stays "not defined by this version", as 1.1.0 left it.

- **No source change.** `parse_nucleotide_call` (`src/progeny_selector/io/calls.py:92-95`) already returned `None` for such a pair, for every caller: `read_hapmap` with `HAPMAP_MISSING`, `read_wide_csv` with `WIDE_NUCLEOTIDE_MISSING`, and both with a compiled profile. Only the docstring changed, from "undecided in contract 1.1.0" to the contract version that states it. The two implementations were compared cell by cell before the contract text was written and agreed everywhere, which is why this is a minor bump rather than a behaviour change here or in the sibling.
- **The rule is pinned by hand-built tests, not only by the shared cases.** `tests/test_io.py` asserts all four nucleotides against all three missing characters, in both orders, in all three spellings, upper and lower case, for both formats and under all five built-in profiles; that `AX` and `AU` stay errors; that `NA` is missing under `dart`, `axiom` and `kasp` because those profiles list it; and that a wide CSV whose only non-A/B cell is the pair is detected as nucleotide and then fails naming the coded letter `B`, not the pair. Both repositories' fixture and case generators are deterministic, so the hand-built tests carry the weight of the interpretation, as CLAUDE.md requires.
- **The contract cases come from the sibling unchanged.** `contract/cases/hapmap-half-missing-pair`, `wide-half-missing-pair`, `profile-tassel-half-missing-pair` and `err-profile-none-half-missing-pair` are generated in backcross and run here through `tests/test_contract_cases.py`; `scripts/check_contract.py ../backcross` byte-compares the mirror.
- **`docs/data-formats.md` loses a difference.** The clause "a pair of one nucleotide and one of N, `-`, `.` (`AN`, `A-`) is read as missing, which the contract has not yet decided" listed this reading among the things this tool adds on top of the contract. It is contract now, so the clause is gone and the version reference moves to 1.6.0.

## Consequences

No input changes meaning, no output changes, no error kind is added, and `results_schema` does not move. A file carrying `A-` loads exactly as before; what changed is that the contract now says why, and that a future reading which rejected it would be a major bump rather than a "stricter reading" of something unstated.

## Revisit when

- The maintainer decides the two-missing pair (`N/N`, `N-`), which both tools read as missing and the contract still leaves undefined.
- A platform is met whose `-` is a real deletion allele worth carrying rather than discarding, which would be a token profile question, not a change to this rule.
