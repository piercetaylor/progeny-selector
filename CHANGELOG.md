# Changelog

All notable changes to this project are documented in this file. The format follows Keep a Changelog 1.1.0 (https://keepachangelog.com/en/1.1.0/) and the project uses Semantic Versioning (https://semver.org/spec/v2.0.0.html).

### Added

- `scripts/bench_pipeline.py` and `scripts/bench_shinylive.py` measure wall-clock and peak memory on generated datasets; `docs/limits.md` records the measured limits for CPython and Shinylive (docs/adr/0021).

## [Unreleased]

### Documentation

- Replace the milestone-heavy landing page with a concise workflow and fixture example; preserve the previous README in `docs/legacy-readme.md`.

### Fixed

- The next-round manifest refuses an empty next-generation label in `io.export`, so `select --next-generation ""` is a usage error and the Export screen shows the refusal and writes no file, where both previously wrote ids like `BC2F1-F1-001--001` and a blank `generation` column; a padded label is trimmed, and a refused write no longer leaves an empty file behind.
- results.csv and the Validate screen record the chromosome-length table the run used, and not the criteria.yaml text. A maize run with an unset `assembly` recorded `Wm82.a4` while using no length table at all, and now records `none` (docs/adr/0015, amendment 2026-09-22).
- An unset criteria.yaml `assembly` means whatever the crop implies: `Wm82.a4` under soybean, `none` under every other crop. An explicit `Wm82.*` under a non-soybean crop is refused naming the key, and an explicit `none` is accepted under any crop. A downloaded criteria.yaml writes an unset key as `assembly: null` and stays reusable across crops.

### Added

- Crop selector (contract 1.5.0, docs/adr/0020): soybean (default), maize, rice, sorghum, wheat, barley, oat, common bean and cotton chromosome schemes; a chosen crop normalizes and orders chromosome names, resolves target regions and selects chromosome lengths by that crop's convention, so a maize `chr1` is no longer read as `Gm01` and a maize target `chr7` resolves. `--crop ID` on the CLI and a Crop select on the Load screen.

### Changed

- results.csv schema 1.1.0: a fixed `crop` column after `results_schema` and before `token_profile`, which stays last (docs/adr/0016); selected.csv gains it in the same position and `scripts/read_results.R` reads it. Under any crop but soybean a chromosome ends at its last marker, because the chromosome length tables are Williams 82's; the per-chromosome "assembly gives no length" warning says so.

### Changed

- **BREAKING CHANGE:** results.csv schema 1.0.0: missing cells are `NA`, the header is written with no results, new columns `background_model`, `background_unit`, `rank_mode`, `assembly`, `results_schema`; `het_rate`, `expected_het` and `expected_rpp` never print `nan`. selected.csv gains `results_schema` and writes `NA` for missing family or generation (docs/adr/0016). next_samples.csv is unchanged: empty cells, no schema column. A chromosome whose name would collide with `rpp_total`, `rpp_carrier` or `rpp_noncarrier` is now an error naming the chromosome.
- Contract 1.3.0 mirrored (docs/adr/0013): a blank line before the header is skipped in HapMap, wide CSV, samples.csv and markers.csv; the delimiter is sniffed from the first non-blank line; a blank line holds only spaces and tabs; a quoted field left open at the end of a delimited file, a `#` line after a VCF `#CHROM` line and an unparseable markers.csv `cm` are errors naming the line; line breaks inside quoted fields are read as LF; error kinds `genotypes.column_count` and `delimited.unterminated_quote`.
- Data contract 1.2.1 (wording only): the sibling browser tool is now called Backcross.

### Fixed

- Weighted RPP: a chromosome whose markers run past the assembly length now weighs its terminal markers like a chromosome with no recorded length at all, instead of giving the terminal marker zero outer weight (docs/adr/0015, amendment 2026-09-21).
- `scripts/kasp_to_wide.py` accepts `A:G` and `G:A` for the same sample and SNP as one call instead of a conflict.
- Clearing a numeric field no longer breaks the screen: Export still downloads next_samples.csv and the Selection list's "Add top N per family" adds nothing instead of failing.
- The next-round manifest refuses a placeholder-row count below 1 in `io.export` rather than only in the CLI, so the Export screen shows the refusal and writes no file where it previously downloaded a manifest holding the two parents and no progeny.

### Added

- `scripts/soysnp_positions.py` builds a SoySNP50K/6K position table across Wm82 assemblies from SoyBase GFF3s; `scripts/kasp_to_wide.py` converts an LGC long-format KASP export to the wide CSV contract.
- Selection list: editable notes written to selected.csv. Validate: background model, units, assembly, rank mode and duplicate pairs. Export: placeholder rows per selected individual. docs/keyboard-walkthrough.md.
- `select --per-selected N`; a synthetic BC3F1 fixture generated from the BC2F1 selection proves next_samples.csv loads unchanged as the next generation (docs/adr/0018).
- Advisory QC flag `possible_duplicate` and duplicate pairs in `validate` (docs/adr/0017). IBS is measured over the informative markers called in both individuals, and a pair is reported only when they share calls at half or more of the markers used, so a sample with almost no calls is no longer a duplicate of everyone it overlaps. The fixture pins weighted RPP and staged ranks computed independently.
- `ranking: {mode: staged}` orders survivors lexicographically (docs/adr/0007); `assembly` selects the chromosome-length table (Wm82.a1, a2, a4 or none; docs/adr/0015). Warnings now report a cM request the map cannot honour, chromosomes without an assembly length, and marker positions beyond the assembly length.
- `scripts/read_results.R` reads results.csv with readr and explicit column types; CI verifies it on the fixture and on an empty results file.
- Token profiles (contract 1.4.0 mirrored, docs/adr/0014; backcross docs/adr/0019): the Load screen, the CLI (`--profile`) and backcross read HapMap and wide-CSV cells under a named vocabulary, `tassel`, `soybase-report` (H heterozygous, U missing), `dart` (0/1/2/-), `axiom` (AA/AB/BB, NoCall, or 0/1/2) or `kasp` (X:X/X:Y/?), or a JSON file of the same shape. results.csv and selected.csv gain a trailing `token_profile` column; readers by position keep their columns. A profile with a VCF, or with a wide CSV or HapMap read as A/B/H, is an error, as is an `H` at a marker with an indel allele; a profile file is always recorded as `custom:<id>`; the Load screen's custom profile can be cleared and disables the select while set; error kinds `genotypes.ambiguous_heterozygote` and `genotypes.profile_format`.
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
