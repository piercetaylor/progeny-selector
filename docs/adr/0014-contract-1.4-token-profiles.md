# Contract 1.4.0 mirrored: named token profiles

Status: accepted. Date: 2026-09-16.

## Context and Problem Statement

Contract 1.4.0 is decided in the canonical repository, backcross docs/adr/0019 ("Contract 1.4.0: named token profiles"), and mirrored here byte for byte under `contract/`, including the five profile files under `contract/profiles/`. That record carries the decisions and their reasons: a profile is a JSON document keyed by platform with homozygous, heterozygous and missing vocabularies and a base flag; profile tokens take precedence over the format's vocabulary, which `base: "nucleotide"` keeps underneath and `base: "none"` replaces; a heterozygote token written `*` resolves to the marker's two alleles on the row, HapMap's `alleles` column included, and any other number of alleles is the error kind `genotypes.ambiguous_heterozygote`; a profile with a VCF is the error kind `genotypes.profile_format`; coded A/B/H detection runs only without a profile or under `base: "nucleotide"`, where the profile's tokens do not vote; the built-ins are `tassel`, `soybase-report`, `dart`, `axiom` and `kasp`; a user-supplied profile is recorded as `custom:<id>`, no profile as `default`. What does this repository change to read files that way and record the profile?

## Decision Outcome

- **`src/progeny_selector/io/profiles.py` (new)**: `TokenProfile` (the JSON fields; heterozygous pairs, `missing` and `sources` are lists so `dataclasses.asdict` equals the parsed file), `BUILTIN_PROFILES` as Python literals copied from the five files, `validate_profile`, `resolve_profile`, `profile_label`, `compile_profile`. Shinylive stages the package, not the repository root, so the literals are the runtime copy; `tests/test_profiles.py` compares each with its file.
- **`io/calls.py`**: `parse_nucleotide_call(text, missing, profile)` implements the resolution order and returns the sentinel `HET_OF_MARKER` for a `*` token; `encode_marker(calls, seed_alleles, cells)` resolves the sentinel with `resolve_het_of_marker` once the row's alleles are known. Before this version the HapMap reader did not read the `alleles` column at all; it now passes it (minus `N`) as the seed, used only to resolve a `*` token, so the allele table of every other row is unchanged.
- **`io/hapmap.py`, `io/wide_csv.py`**: a `profile` parameter; `detect_coding` is skipped under `base: "none"` and does not see the profile's tokens otherwise; `coding="abh"` with a profile is an error.
- **`io/__init__.py`**: `load_genotypes(path, coding, profile=None)` rejects a profile with a VCF; `load_dataset(..., profile=None)` sets the new `Dataset.token_profile` (default `"default"`).
- **CLI** `--profile ID_OR_FILE` on `validate` and `rank` (a value containing `/` or `\` or ending `.json` is read as a file, as in backcross); **Load screen** a "Token profile" select and a "Custom token profile (JSON, optional)" file input.
- **Exports**: every results row carries `token_profile` as its last key, so results.csv ends with that column; selected.csv appends `token_profile` after `notes`.

### Amendment, 2026-09-16: review findings (maintainer decisions, still 1.4.0; backcross docs/adr/0019 amendment)

- `read_wide_csv` raises `genotypes.profile_format` when a profile meets a file whose coding is `abh`, requested or detected, naming the profile and which; the claimed-token filter stays, so a profile's tokens never decide coding silently.
- `profile_label` records an id only for the `BUILTIN_PROFILES` object itself; a validated file is `custom:<id>` even with a built-in's id.
- `validate_profile` rejects unknown top-level keys, a token listed twice within a set after normalisation, and a heterozygote pair of identical symbols.
- The HapMap `alleles` seed is upper-cased; `resolve_het_of_marker` rejects a row whose alleles include `-`.
- Load screen: the Token profile select is rendered disabled while a custom profile file is set, and "Clear custom token profile" forgets the file and re-renders its input.

### Consequences

Good: SoyBase, DArT, Axiom, GenomeStudio, KASP and TASSEL-converted wide CSV files load here exactly as they load in backcross, and results.csv and selected.csv say which vocabulary produced them.

Bad: results.csv and selected.csv each gain a trailing column; a VCF with a profile is an error rather than a silent no-op.
