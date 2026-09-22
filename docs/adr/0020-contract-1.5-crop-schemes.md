# Contract 1.5.0: crop chromosome schemes

Status: accepted. Date: 2026-09-22. Mirrors backcross docs/adr/0020, the canonical record for this contract version.

## Context and Problem Statement

Up to contract 1.4.0 chromosome names were soybean-only: `Gm01`..`Gm20`, one alias pattern, every other name kept as written and ordered after `Gm20`. A maize `chr1` was displayed as `Gm01`, a wheat `Chr1A` was an unrecognised name, and a rice `Chr12` matched only by accident of the soybean pattern. Contract 1.5.0 (backcross docs/adr/0020, the canonical record; docs/m4-phases.md section 8) gives every crop its own scheme. This record is the sibling half, S6: what this repository implements and where it deliberately differs from the canonical record.

## Decision Outcome

The canonical record's decisions D6.1 to D6.6 apply here unchanged and are not restated: a scheme is a JSON document under `contract/crops/` with an ordered list of canonical names, a parallel list of keys and one alias regular expression whose capture groups yield the key (join the defined groups, upper-case, strip the leading zeros of each digit run, look the key up in `keys`; a miss keeps the name as written and orders it after the canonical names in natural order, exactly the 1.2.0 rule); `soybean` is the default and reproduces 1.2.0 exactly; nine crops ship (soybean, maize, rice, sorghum, wheat, barley, oat, common bean, cotton) and cowpea, pea, sunflower and peanut are deferred as ambiguous.

- **`Dataset` carries the compiled scheme and `Dataset.crop` is derived from it.** `Dataset.scheme: CompiledScheme` defaults to `SOYBEAN`, `load_dataset` sets it from `resolve_crop(crop)`, and `crop` is a property returning `scheme.scheme.id`, so the two cannot disagree. This is what lets `core` use the chosen scheme without importing `io`: `run_analysis` reads `dataset.scheme` and passes it to `sorted_by_position`, `_chrom_lengths`, `resolve_locus`, `marker_weights` and `rpp_per_chromosome`, and the Compare screen passes it to `chromosome_strips`. It matches the sibling, where `assembleDataset` takes `meta.crop` and records `meta.crop.scheme.id`.
- **The scheme lives in `core/chrom.py`, the nine literals in `io/crops.py`.** `core` imports no `io`, so `CropScheme`, `CompiledScheme`, `compile_scheme` and the `soybean.json` literal (`SOYBEAN_SCHEME`, compiled as `SOYBEAN`) are in `core/chrom.py`, which every default argument uses; `io/crops.py` holds `BUILTIN_CROPS`, `validate_scheme` and `resolve_crop`, and re-exports the types. This is the layering of backcross's `src/core/chromosomes.ts` and `src/io/crops.ts`. Shinylive stages the package and not `contract/`, so the literals are the runtime copy and `tests/test_crops.py` checks each against its JSON file, as `tests/test_profiles.py` does for token profiles.
- **`crop` is a fixed column after `results_schema` and before `token_profile`, not a trailing column.** The canonical record's D6.5 appends `crop` after `token_profile`; docs/m4-phases.md section 8.7a corrects that for this repository, because docs/adr/0016 fixes the position of a new fixed column here and `scripts/read_results.R` (commit `2b966bc`) takes the last column as the token profile. results.csv and selected.csv therefore read `..., results_schema, crop, token_profile`, and `RESULTS_SCHEMA` bumps to `1.1.0`, the minor bump docs/adr/0016 prescribes for a new fixed column.
- **The chromosome-length table stays soybean's, and `chrom_length_bp` takes the scheme.** Its signature is `chrom_length_bp(name, fallback=None, assembly=DEFAULT_ASSEMBLY, scheme=SOYBEAN)` and it returns `fallback` for every scheme but soybean, so a maize `chr1` never takes `Gm01`'s Williams 82 length: its callers (`core/pipeline.py`, `core/strip.py`) already substitute the last marker position on the chromosome, and the run warns that the assembly gives no length instead of falsely warning that markers reach beyond one. Choosing a crop therefore changes names, their order and which chromosomes have a recorded length, never a coordinate.
- **No user-supplied scheme.** `resolve_crop` accepts a built-in id or `None`; an unknown id is a `DataContractError` naming the built-ins. `validate_scheme` exists for the equality test against `contract/crops/` and for a future user-supplied file.

## Consequences

A maize `chr1` is no longer read as `Gm01`. A dataset loaded without `--crop` (CLI) or with the Load screen's default is read exactly as contract 1.2.0 read it, so no earlier output changes except for the new `crop` cell and the schema version. The nine `crop-<id>-spellings` cases under `contract/cases/` pin each scheme through `load_dataset`, and every threading test in `tests/test_crops.py` uses a non-soybean crop, since a soybean test cannot tell a threaded scheme from the default.

A target region written in the dataset's own spelling resolves: under maize, `chrom: chr7` stays `chr7` and matches the markers there, where the soybean default would have made it `Gm07` and matched nothing. Every chromosome ordering the analysis produces (marker order, per-chromosome RPP columns, Compare strips) is the chosen crop's, so oat's unplaced `Un0` follows `chr1A` and `chr7D` rather than leading them.

## Revisit when

- A crop needs chromosome lengths of its own: the tables in `constants.py` are Williams 82's, and `chrom_length_bp` returns `fallback` for every other crop rather than guessing.
- A second assembly of a shipped crop needs its own names (an `assemblies` list in the schema rather than a second id).
- A user-supplied scheme is asked for, or one of the deferred crops gets an unambiguous public convention.
