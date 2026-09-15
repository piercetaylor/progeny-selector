# Changelog

All notable changes to this project are documented in this file. The format follows Keep a Changelog 1.1.0 (https://keepachangelog.com/en/1.1.0/) and the project uses Semantic Versioning (https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

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
- `scripts/check_contract.py`, a per-commit gate that recomputes `contract/MANIFEST.sha256` and, given `../isoline-browser`, byte-compares the mirror with the canonical copy.

### Changed

- Loaders aligned with contract 1.1.0 (docs/adr/0010): samples.csv, markers.csv and wide CSV may be tab-delimited and quoted; `line_name` is optional; genotype columns not in samples.csv are dropped and samples load in manifest order; a VCF record with ID `.` is named from CHROM as written (`chr13_19000000`, not `Gm13_19000000`); chromosome names accept the `LG` prefix and space or `-` separators, no longer accept `ch6`, and non-soybean names order naturally (`scaffold_2` before `scaffold_10`); `?` is no longer a missing token anywhere, and in HapMap and nucleotide wide CSV a single character outside A, C, G, T and the IUPAC codes (`?`, `B`, `H`, `0`, `+`) is an error naming the line and cell rather than a missing call; the HapMap missing tokens are exactly the contract's eleven (`X`, `XX` added, `?` removed); A/B/H auto-detection reads the whole file instead of its first 200 rows; a leading byte-order mark is accepted on every input.

### Fixed

- Genotype files named `.vcf.bgz` now load, as documented; `.bcf`, which was never documented and was read as text, is rejected.

## [0.1.0] - unreleased

Initial scaffold; no release has been tagged.
