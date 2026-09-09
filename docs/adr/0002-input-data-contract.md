# VCF, HapMap and wide CSV inputs with samples.csv and markers.csv as the contract shared with isoline-browser

Status: accepted. Date: 2026-09-04.

## Context and Problem Statement

Progeny genotypes arrive as VCF (array or sequencing pipelines), HapMap (TASSEL GBS) or spreadsheets from KASP software, one generation before the same lines reach the sibling isoline-browser. Which formats are read, and how are parents, generation and family declared so files move between the two tools unchanged?

## Decision Drivers

No reformatting by the user; explicit parent declaration (the recurrent parent is never inferred from REF); a per-generation family structure for navigation and per-family ranks; support for A/B/H matrices without parent columns; one validation boundary.

## Considered Options

1. VCF only.
2. VCF, HapMap and wide CSV, with roles, generation and family in samples.csv and an optional markers.csv map.
3. A bespoke project file (JSON/YAML) embedding genotypes.

## Decision Outcome

Option 2, specified in docs/data-formats.md and mirrored in the sibling repository. `generation` follows the `BCnFk` convention (F1 = no selfing) so expected RPP and heterozygosity can be derived; `family_id` groups progeny for per-family ranking and top-N selection. The VCF reader uses only CHROM, POS, ID, REF, ALT, FORMAT and GT [web] https://samtools.github.io/hts-specs/VCFv4.2.pdf; HapMap uses the 11 fixed columns [web] https://statgen-esalq.github.io/Hapmap-and-VCF-formats-and-its-integration-with-onemap/; the wide CSV auto-detects nucleotide versus A/B/H coding and synthesises parents for coded input.

### Consequences

Good: the program's existing files load; the manifest is a spreadsheet a breeder can write; fixtures and contract tests can be shared with isoline-browser; results.csv and next_samples.csv round-trip into the next generation. Bad: three parsers; content-based detection is not attempted (extension decides), so a mis-named file is an explicit error; coded input hides non-parental alleles that would otherwise flag outcrosses. Deferred: BrAPI loading; `donor_parent_2` for pyramiding.
