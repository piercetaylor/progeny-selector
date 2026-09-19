# Changelog

All notable changes to this project are documented in this file. The format follows Keep a Changelog 1.1.0 (https://keepachangelog.com/en/1.1.0/) and the project uses Semantic Versioning (https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- Contract 1.3.0 mirrored (docs/adr/0013): a blank line before the header is skipped in HapMap, wide CSV, samples.csv and markers.csv; the delimiter is sniffed from the first non-blank line; a blank line holds only spaces and tabs; a quoted field left open at the end of a delimited file, a `#` line after a VCF `#CHROM` line and an unparseable markers.csv `cm` are errors naming the line; line breaks inside quoted fields are read as LF; error kinds `genotypes.column_count` and `delimited.unterminated_quote`.
- Data contract 1.2.1 (wording only): the sibling browser tool is now called Backcross.

### Added

- Token profiles (contract 1.4.0 mirrored, docs/adr/0014; backcross docs/adr/0019): the Load screen, the CLI (`--profile`) and backcross read HapMap and wide-CSV cells under a named vocabulary, `tassel`, `soybase-report` (H heterozygous, U missing), `dart` (0/1/2/-), `axiom` (AA/AB/BB, NoCall, or 0/1/2) or `kasp` (X:X/X:Y/?), or a JSON file of the same shape. results.csv and selected.csv gain a trailing `token_profile` column; readers by position keep their columns. A profile with a VCF or with a wide CSV that is coded A/B/H is an error, as is an `H` at a marker with an indel allele; a profile file is always recorded as `custom:<id>`; the Load screen's custom profile can be cleared and disables the select while set; error kinds `genotypes.ambiguous_heterozygote` and `genotypes.profile_format`.
- Package scaffold with a pure compute core, boundary parsers, CLI and a Shiny for Python screen shell.
- Parsers for VCF 4.2+ (plain, gzip, bgzip), HapMap, wide CSV (nucleotide and A/B/H), samples.csv, markers.csv and criteria.yaml, with boundary validation.
- Parent-of-origin classification; foreground status for marker, region and flanking-marker loci; avoid-locus status; recurrent-parent proportion (count and map-weighted, per chromosome, carrier and non-carrier); donor-segment bounds and flank recombinant flags; identity-by-state to each parent; QC metrics and flags; hard filters, composite score, ranking with fixed tie-breaks; top-N selection and next-generation projection.
- results.csv, selected.csv and next-round manifest writers; `progeny-selector validate | rank | select` commands.
- Synthetic BC2F1 fixture (2 parents, 40 progeny in 2 families, 500 markers, planted target, avoid locus, contaminant, high-missing individual) with generator and expected results; unit and smoke tests.
- PLAN.md, data-format and reference-repository documents, MADR decision records, CI workflow.
- Navigate screen: family and generation tree with breadcrumbs that filter the Rank grid. Rank grid shows pass, fail and unknown chips and recombinant flags per locus.
- Validate screen: QC table with flagged rows highlighted and the reasons markers were uninformative. Compare screen: side-by-side statuses, per-chromosome RPP, drag bounds and Okabe-Ito chromosome strips.
- Criteria can be edited on the Load screen, re-applied without reloading genotypes, and downloaded as criteria.yaml. Downloads no longer write temporary files.
- The Shinylive static export includes the package and is smoke-tested in CI with the fixture; it deploys to GitHub Pages when ENABLE_PAGES is set.
- `scripts/soysnp50k_nils.py` converts one backcross-derived NIL family from the SoyBase SoySNP50K VCF and PATRIOT's pedigree file into a genotype file and samples.csv for this tool. It streams the 144 MB VCF, strips SoyBase's chromosome prefix, and drops scaffold and malformed records.
- The shared input contract, version 1.1.0, mirrored byte for byte from isoline-browser under `contract/` and checked by `tests/test_contract_cases.py`, which loads every case through `load_dataset` and verifies the manifest hashes (docs/adr/0010).
- Foreground `rule: run` for region targets, with `min_run`, `anchor_bp` and `tolerate_isolated`: a target passes only on a contiguous run of donor calls through the anchor, so scattered array calls no longer pass a wide window (docs/adr/0011).
- Advisory QC flag `family_donor_outlier`: an individual whose donor fraction (count model) lies above its family median by more than 2.5 robust scales (1.4826 x MAD, floor 0.01), in families of at least 6; it does not exclude (docs/adr/0012).
- `scripts/check_contract.py`, a per-commit gate that recomputes `contract/MANIFEST.sha256` and, given `../isoline-browser`, byte-compares the mirror with the canonical copy.

### Changed

- Contract 1.2.0, a minor version (docs/adr/0010, amended 2026-09-14): VCF POS must be decimal digits; a VCF record with ID `.` is named from the parsed POS. A position written as a whole-valued float or in exponent notation (`1000.0`, `1e3`, `1.9E+07`) is read as that integer in HapMap, wide CSV and markers.csv; a fractional `pos_bp` in wide CSV or markers.csv is now an error naming the line and value instead of being truncated, a non-integer VCF `POS` is a `DataContractError` instead of an uncaught `ValueError`, an empty `marker_id` in a wide CSV is an error, and wide-CSV and markers.csv errors name the physical line even after skipped rows.
- Loaders aligned with contract 1.1.0 (docs/adr/0010): samples.csv, markers.csv and wide CSV may be tab-delimited and quoted; `line_name` is optional; genotype columns not in samples.csv are dropped and samples load in manifest order; a VCF record with ID `.` is named from CHROM as written (`chr13_19000000`, not `Gm13_19000000`); chromosome names accept the `LG` prefix and space or `-` separators, no longer accept `ch6`, and non-soybean names order naturally (`scaffold_2` before `scaffold_10`); `?` is no longer a missing token anywhere, and in HapMap and nucleotide wide CSV a single character outside A, C, G, T and the IUPAC codes (`?`, `B`, `H`, `0`, `+`) is an error naming the line and cell rather than a missing call; the HapMap missing tokens are exactly the contract's eleven (`X`, `XX` added, `?` removed); A/B/H auto-detection reads the whole file instead of its first 200 rows; a leading byte-order mark is accepted on every input.

### Fixed

- Genotype files named `.vcf.bgz` now load, as documented; `.bcf`, which was never documented and was read as text, is rejected.

## [0.1.0] - unreleased

Initial scaffold; no release has been tagged.
