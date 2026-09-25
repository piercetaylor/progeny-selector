# Contract 1.7.0: crop schemes for cowpea, pea and peanut

Status: accepted. Date: 2026-09-25. Mirrors backcross docs/adr/0022, the canonical record for this contract version.

## Context and Problem Statement

Contract 1.5.0 (docs/adr/0020 here, backcross docs/adr/0020) shipped chromosome schemes for nine crops and deferred cowpea, pea, sunflower and peanut as ambiguous. Contract 1.7.0 (backcross docs/adr/0022, the canonical record) adds three of the four as optional crop ids, a minor bump. This record is the sibling half: what this repository implements and what it checks.

## Decision Outcome

The canonical record's decisions apply here unchanged and are not restated, nor are its evidence, sources and list of unverified points: `cowpea` (canonical `Vu01`..`Vu11`) keeps accepting a bare `1`..`11` and accepts NCBI's parenthetical old-LG form only in its eleven exact pairings; `pea` (canonical `chr1LG6`..`chr7LG7`) accepts `chr`/`chromosome` with the karyotype number, or the karyotype number with its matching `LG` suffix, and refuses a bare `1`..`7`, because published pea tables use bare numbers for both the karyotype and the linkage-group numbering; `peanut` (canonical `Arahy.01`..`Arahy.20`) leaves `A01`..`B10`, `Aradu.` and `Araip.` spellings unmapped, since the scheme builds a key from a substring of the name and cannot map `B01` to `11` without an alias table; sunflower stays deferred. Both bare-number decisions were delegated by the maintainer on 2026-09-25, with the criterion "never silently wrong, tolerant of standard-tool exports".

- **Three literals in `io/crops.py`, nothing else in the source.** `BUILTIN_CROPS` gains `cowpea`, `pea` and `peanut` after `cotton`, with each `pattern` written as a raw string so `\s`, `\.` and `\(` reach `re` as the JSON's decoded value. `test_literal_equals_file` compares each literal with its file under `contract/crops/`, as for the nine. The CLI's `--crop` help and the Load screen's Crop select are built from `BUILTIN_CROPS`, so neither changed.
- **No change to `core/chrom.py`.** The cowpea and pea patterns are alternations whose other groups are undefined on a match; `_key_of` joins only the groups that are not `None`, so `Vu11(old9)` yields `11` and `4LG4` yields `4`. This was confirmed against the new patterns rather than assumed, and the pattern is compiled case-insensitively, so `chr1lg6` maps.
- **The hand-built tests carry the weight.** `tests/test_crops.py` asserts, with the same spellings as backcross `tests/crops.test.ts`, what each scheme maps and what falls through: cowpea `Vu01(old7)`, `LG4`, `VuLG4`, `Vu12`; pea `4`, `04`, `chr1LG1`, `1LG1`, `LG4`, `LG6`, `chr8`, `8`, `chr0`; peanut `A01`, `B01`, `Aradu.A09`, `Araip.B08`, `Arahy.21`, `Arahy.00`; and that under pea a bare `4` orders after all seven canonical names. The `crop-<id>-spellings` cases are generated in backcross and run here through `tests/test_contract_cases.py`.
- **The unknown-crop tests now use `sunflower`.** Three tests used `cowpea` as the example of an id that is not built in; it is built in now, so they use the crop that is still deferred.

## Consequences

A dataset loaded with any of the nine earlier crops, or with none, reads exactly as before, and `results_schema` does not move: a new crop id is a new value in the existing `crop` column. A pea file whose chromosome column holds bare `1`..`7` keeps those names as written and orders them after the canonical names, which is visible in the Compare strips and in the per-chromosome columns; the user renames the column once rather than reading five of seven chromosomes under the wrong numbering. Peanut's subgenome spellings do the same.

## Revisit when

- The scheme schema gains an alias table, which would let peanut map `A01`..`B10` and the `Aradu.`/`Araip.` spellings.
- A statement or whole-genome synteny table shows that sunflower's HA412-HO and HanXRQ chromosome numbers agree (the unblock condition in backcross docs/adr/0022).
