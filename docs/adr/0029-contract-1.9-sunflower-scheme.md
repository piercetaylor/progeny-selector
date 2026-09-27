# Contract 1.9.0: a sunflower chromosome scheme

Status: accepted. Date: 2026-09-26. Mirrors backcross docs/adr/0024, the canonical record for this contract version.

## Context and Problem Statement

Contract 1.5.0 (docs/adr/0020 here, backcross docs/adr/0020) deferred cowpea, pea, sunflower and peanut as ambiguous, and contract 1.7.0 (docs/adr/0027 here, backcross docs/adr/0022) shipped three of them and deferred sunflower until a statement or a whole-genome synteny table showed that the HA412-HO and HanXRQ chromosome numbers agree. Contract 1.9.0 (backcross docs/adr/0024, the canonical record) adds `sunflower` as an optional crop id, a minor bump. This record is the sibling half: what this repository implements and what it checks.

## Decision Outcome

The canonical record's decisions apply here unchanged and are not restated, nor are its evidence, sources and list of unverified points: `sunflower` (canonical bare `1`..`17`) reads the `HanXRQChr` and `Ha412HOChr` prefixes, `chr`, `chromosome` and a bare number, and keeps unplaced scaffolds, `0` and `18`, organelles, the `LG` prefix, NCBI accessions and HA412-HO v1.1 `Ha1`..`Ha17` names as written. The evidence chain is in backcross docs/adr/0024. Its weakest link is HA412-HOv2: no method statement says how that assembly's chromosomes were numbered, so its pairing with XRQ rests on a whole-genome MUMmer alignment figure (PNAS 2023 SI, Fig. S2) whose XRQ labels are clipped in the PDF, read from the rank-matching diagonal rather than from printed labels. The maintainer accepted this on 2026-09-26 as meeting the unblock condition of backcross docs/adr/0022, although it is a figure and not a table.

- **One literal in `io/crops.py`, nothing else in the source.** `BUILTIN_CROPS` gains `sunflower` after `peanut`, making thirteen schemes, with its `pattern` written as a raw string so `\s` reaches `re` as the JSON's decoded value. `test_literal_equals_file` compares it with `contract/crops/sunflower.json`, as for the twelve. The CLI's `--crop` help and the Load screen's Crop select are built from `BUILTIN_CROPS`, so neither changed beyond a count in a docstring.
- **No change to `core/chrom.py`** beyond the count in its module docstring. The pattern has one capture group, so `_key_of` yields the bare number (the pattern's `0*` sits outside the group, so `01` yields `1`), as for the other single-group schemes; the pattern is compiled case-insensitively, so `HANXRQCHR09` maps.
- **The hand-built tests carry the weight.** `tests/test_crops.py` asserts, with the same spellings as backcross `tests/crops.test.ts`, what the scheme maps (`Ha412HOChr01`, `HanXRQChr17`, `HANXRQCHR09`, `chr1`, `Chr_17`, `chromosome-4`, `17`, `01`) and what falls through (`HanXRQChr00c001`, `Ha412HOChr00`, `chr0`, `chr18`, `18`, `MT`, `Pltd`, `HanXRQMT`, `HanXRQCP`, `Ha412HOv2Chr01`, `LG1`, `NC_035433.2`, `CM007890.2`, `Ha1`, `Ha10`). The `crop-sunflower-spellings` case is generated in backcross and runs here through `tests/test_contract_cases.py`.
- **The unknown-crop tests now use `potato`.** Three tests used `sunflower` as the example of an id that is not built in; it is built in now, so they use `potato`, which is out of scope, as backcross does.

## Consequences

A dataset loaded with any of the twelve earlier crops, or with none, reads exactly as before, and `results_schema` does not move: a new crop id is a new value in the existing `crop` column. A sunflower export shows bare `1`..`17` rather than either assembly's prefix. A file on HA412-HO v1.1 names, or one carrying accessions, keeps those names as written and orders them after the canonical names, which is visible in the Compare strips and the per-chromosome columns.

## Revisit when

- Backcross docs/adr/0024's revisit conditions fire: a v1.1 header list confirms `Ha1`..`Ha17`, a real sunflower export uses `LG` for chromosome numbers, or evidence shows the HA412-HOv2 numbering differs from XRQ's on any chromosome.
- The scheme schema gains an alias table, which would let sunflower map NCBI accessions as it would peanut's subgenome spellings.
