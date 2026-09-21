# The manifest format is pinned by a two-generation round trip, not by the writer alone

Status: accepted. Date: 2026-09-21.

## Context and Problem Statement

The point of next_samples.csv is that a breeder ranks one generation, selects, plants the selected lines, genotypes their progeny and loads the same file back as the next generation's samples.csv. Until now nothing tested that loop. `write_next_round_manifest` was covered only by tests that read what it had just written, so any change to its columns, its empty cells or its derived ids would have stayed green: the writer was being compared with itself. The failure this invites is concrete. results.csv and selected.csv write `NA` for a missing value (docs/adr/0016); next_samples.csv must not, because it is the input contract's samples.csv and isoline-browser's manifest reader (`src/io/manifest.ts`) would take `NA` for a family named "NA". A writer that drifted to `NA` for consistency would break the sibling silently.

A second gap: the writer already accepted `n_per_selected`, but neither the CLI nor the Export screen exposed it, so a breeder planting ten seeds per selected line had to edit the file by hand — and an edited file is no longer the file the tool wrote.

The real two-generation dataset this should ultimately run on does not exist yet (M2 question B; the maintainer is sourcing one).

## Considered Options

1. Unit-test `write_next_round_manifest` harder: more columns asserted, more edge cases.
2. Generate a second fixture with the writer itself and load it back.
3. Generate a second fixture independently of the writer, and compare the writer's output with it.

## Decision Outcome

Option 3, with option 1's literal assertions kept as a backstop.

`scripts/make_fixture.py` breeds a BC3F1 generation from the BC2F1 fixture: it picks its own top two per family by its own `rank_in_family` (asserting they are `BC2F1-F1-001`, `BC2F1-F1-010`, `BC2F1-F2-005`, `BC2F1-F2-019`), takes each selected plant's true, pre-missing states as its haplotypes, draws a Haldane gamete against an all-recurrent one, and writes ten progeny per parent to `tests/fixtures/synthetic_bc3f1/`. It builds `samples.csv` in the manifest layout from its own simulation, without importing `progeny_selector.io`, exactly as it already computes `expected_results.csv` without importing `progeny_selector.core` (CLAUDE.md: fixture expectations come from a second implementation). `tests/test_round_trip.py` then asserts that `next_round_manifest_text(...)` from the BC2F1 selection equals that file character for character, that the generation loads and ranks against its own expectations, and that the same loop through the CLI produces the same bytes.

Independence is the whole value: if the generator imported the writer, or the fixture were written by the writer, the test would prove only that the code is self-consistent. The test also pins the parent rows literally and asserts the token `NA` appears nowhere in the file, so the two implementations cannot drift together.

`select --per-selected N` (default 1, an integer of at least 1; `0` is a usage error, exit 2) exposes `n_per_selected` on the CLI, and the Export screen gains the same control. A file with ten placeholder rows per selected line is what the round trip uses, so the flag is exercised rather than merely offered.

### Consequences

- The BC3F1 fixture's `criteria.yaml` and `markers.csv` are not duplicated; the tests read the BC2F1 copies. One criteria file, one map, no chance of the two drifting.
- The writer copies a selected parent's `line_name` to all of its placeholder progeny, so ten rows share a line name. That is the documented behaviour (the ids and names are placeholders to be edited after planting), and the fixture README says so rather than the generator inventing distinct names the writer would not produce.
- The generator's random draws for the BC3F1 generation all happen after the BC2F1 fixture is complete, so the BC2F1 files stay byte-identical.
- The BC3F1 fixture is synthetic. It proves the file format round-trips and that the second generation's metrics are right; it does not prove anything about real segregation. The real-data step needs a dataset with both parents genotyped, generation *n* progeny with a samples.csv, generation *n+1* progeny genotyped on the same panel whose samples.csv is this tool's next_samples.csv unchanged, and positions on one assembly. Until such a dataset exists, this fixture is the acceptance test.
