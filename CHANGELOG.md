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

### Fixed

- Genotype files named `.vcf.bgz` now load, as documented; `.bcf`, which was never documented and was read as text, is rejected.

## [0.1.0] - unreleased

Initial scaffold; no release has been tagged.
