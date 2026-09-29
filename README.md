# progeny-selector

progeny-selector ranks candidates in marker-assisted backcross breeding. It checks target and avoid loci, estimates recurrent-parent recovery and donor-segment bounds, applies quality flags and selection criteria, and exports ranked results, selections, and a manifest for the next generation. The analysis runs as a Python package and command-line tool, with a Shiny interface for local use or a browser-hosted site.

Status: version 0.1.0; see [CHANGELOG.md](https://github.com/piercetaylor/progeny-selector/blob/main/CHANGELOG.md) and [Known limitations](#known-limitations).

## Install

Python 3.11 or newer is required.

```sh
python -m pip install progeny-selector
```

or, to keep it out of your other Python projects:

```sh
pipx install progeny-selector
```

On Windows, use `py -3 -m pip install progeny-selector`. If the `progeny-selector` command is not found afterwards, run `py -3 -m progeny_selector ...` instead, or install with pipx (`py -3 -m pip install --user pipx && py -3 -m pipx ensurepath`).

## Open the interface

```sh
progeny-selector app
```

This opens a browser tab at http://127.0.0.1:8000/. Options: `--port`, `--host` and `--no-browser`. Press Ctrl+C in the terminal to stop it.

Or use the [browser site](https://piercetaylor.github.io/progeny-selector/): files are processed in the tab and never uploaded. Both offer **Load example (synthetic BC2F1)** on the Load screen.

## Run the example

```sh
progeny-selector example --out example
progeny-selector validate --genotypes example/synthetic_bc2f1/genotypes.vcf --samples example/synthetic_bc2f1/samples.csv --markers example/synthetic_bc2f1/markers.csv --criteria example/synthetic_bc2f1/criteria.yaml
progeny-selector rank --genotypes example/synthetic_bc2f1/genotypes.vcf --samples example/synthetic_bc2f1/samples.csv --markers example/synthetic_bc2f1/markers.csv --criteria example/synthetic_bc2f1/criteria.yaml --out results.csv
progeny-selector select --results results.csv --top 3 --out selected.csv --next-manifest next_samples.csv --next-generation BC3F1 --samples example/synthetic_bc2f1/samples.csv
```

The example has 500 markers, 2 parents and 40 BC2F1 progeny in 2 families; 11 of the 40 pass the hard filters. It is synthetic, and its rankings are test cases, not breeding recommendations.

For a step-by-step walk through the same data, see the [tutorial](https://github.com/piercetaylor/progeny-selector/blob/main/docs/tutorial.md).

## What the outputs mean

`results.csv` has one row per individual:

- `rank_overall` and `rank_in_family`, `passes_filters` and `exclusion_reason`, and `composite_score`.
- `rpp_total`, `rpp_carrier` and `rpp_noncarrier`: recurrent-parent recovery overall, on chromosomes that carry a target locus, and on those that do not.
- `drag_total_est` and `drag_total_max`, with `drag_unit` naming the unit of both.
- `missing_rate` and `qc_flags`.
- One status column per target and avoid locus, and one `rpp_<chromosome>` column per chromosome.
- `results_schema` names the column layout. `background_model`, `background_unit`, `background_max_marker_coverage`, `assembly` and `crop` record how the numbers were computed.

`selected.csv` holds the chosen individuals, with a `notes` column. `next_samples.csv` is a `samples.csv` for the next genotyping round.

Further reading: the [tutorial](https://github.com/piercetaylor/progeny-selector/blob/main/docs/tutorial.md), the [glossary](https://github.com/piercetaylor/progeny-selector/blob/main/docs/glossary.md), and [selection and output formats](https://github.com/piercetaylor/progeny-selector/blob/main/docs/data-formats.md).

### What `--top N` selects

`select --top N` keeps every passing individual ranked N or better within its family, so up to N per family. With `--overall` it keeps the N best overall instead. Ranks are dense (tied individuals share a rank and the next rank follows without a gap), so ties can return more than N.

On the example, `--top 3` writes 6 rows: 3 in family F1 and 3 in F2. `--top 3 --overall` writes 3. The Selection screen's "add top N per family" does the same as `--top N`.

## Your own data

- Genotypes: VCF, HapMap or wide CSV, each optionally gzip- or bgzip-compressed (`.gz`, `.bgz`); the file extension selects the format. Every call is diploid. On the command line, a BrAPI v2.1 server can be the source instead (`--brapi-url`, `--variant-set`, `--brapi-token-env`).
- `samples.csv`: columns `sample_id` and `role` are required; `line_name`, `generation`, `family_id` and `notes` are optional. `role` is `recurrent_parent`, `donor_parent`, `candidate` or `progeny`. The file needs exactly one recurrent parent, exactly one donor parent, and at least one candidate or progeny.
- `criteria.yaml`: the top-level keys are `name`, `targets`, `avoid`, `flank_window`, `flank_unit`, `assembly`, `background`, `ranking`, `weights` and `filters`. Unknown keys are an error.
- `markers.csv` (optional): `marker_id`, `chrom` and `pos_bp` are required; `cm` is optional and enables cM-weighted recovery and cM segment bounds. Positions here override those in the genotype file.
- Crop: `--crop` selects one of 13 chromosome-naming schemes (soybean, maize, rice, sorghum, wheat, barley, oat, common bean, cotton, cowpea, pea, peanut, sunflower). Soybean is the default. The `assembly` key in `criteria.yaml` (`Wm82.a1`, `Wm82.a2`, `Wm82.a4` or `none`) sets the chromosome-length table used for chromosome ends.
- `--profile` chooses how HapMap and wide CSV cell values are read: one of the built-in token profiles (`tassel`, `soybase-report`, `dart`, `axiom`, `kasp`) or your own JSON file.

Details: the [input contract](https://github.com/piercetaylor/progeny-selector/blob/main/contract/data-contract.md), [selection and output formats](https://github.com/piercetaylor/progeny-selector/blob/main/docs/data-formats.md), and [soybean inputs](https://github.com/piercetaylor/progeny-selector/blob/main/docs/soybean-inputs.md).

## Known limitations

1. `possible_duplicate` is raw pairwise identity by state of at least 0.995 over informative markers, and it does not consider family. It is advisory only ([ADR 0017](https://github.com/piercetaylor/progeny-selector/blob/main/docs/adr/0017-duplicate-flag.md)).
2. The two-generation round trip is verified on synthetic data only. No real linked BC(n) to BC(n+1) dataset has been run.
3. The browser site has no BrAPI source. BrAPI is command-line only until a server sends CORS headers for the site's origin ([ADR 0024](https://github.com/piercetaylor/progeny-selector/blob/main/docs/adr/0024-brapi-allele-matrix-loader.md)).
4. One donor per analysis. An intercross or pyramiding population is not ranked ([contract](https://github.com/piercetaylor/progeny-selector/blob/main/contract/data-contract.md), [ADR 0002](https://github.com/piercetaylor/progeny-selector/blob/main/docs/adr/0002-input-data-contract.md)).
5. Chromosome-length tables exist for soybean only (Wm82.a1, Wm82.a2 and Wm82.a4). Under any other crop a chromosome ends at its last marker ([ADR 0020](https://github.com/piercetaylor/progeny-selector/blob/main/docs/adr/0020-contract-1.5-crop-schemes.md)).
6. The default cap of 2 Mb on base-pair coverage is a soybean translation of 10 cM. Set `background.max_marker_coverage` for other crops ([ADR 0033](https://github.com/piercetaylor/progeny-selector/blob/main/docs/adr/0033-results-schema-1.2-cap-and-end-rule.md)).
7. A VCF 4.4 leading phase indicator in GT is rejected as `genotypes.invalid_gt` (contract 1.10.0).
8. Accessibility: selecting several rows in the Rank grid needs a mouse; three axe rules are recorded exceptions; there is no screen-reader, Firefox, Safari or 320 px reflow testing ([accessibility review](https://github.com/piercetaylor/progeny-selector/blob/main/docs/accessibility.md)).
9. Purdy labels: single crosses only, and no heterozygosity checks, because the label carries no filial generation ([ADR 0034](https://github.com/piercetaylor/progeny-selector/blob/main/docs/adr/0034-purdy-generation-labels.md)).
10. The long-format KASP header names are unconfirmed against a real LGC export ([ADR 0019](https://github.com/piercetaylor/progeny-selector/blob/main/docs/adr/0019-soybean-input-tooling.md)).
11. CI runs the browser tests on Linux only; Windows and macOS run the unit tests.

## Verification

From a source checkout, `ruff check .`, `ruff format --check .`, `mypy` and `pytest -q` check the package. Browser and static-export checks, and the real-data verification with its limits, are described in [the plan](https://github.com/piercetaylor/progeny-selector/blob/main/PLAN.md). Measured run time and memory at several dataset sizes: [limits](https://github.com/piercetaylor/progeny-selector/blob/main/docs/limits.md). The WCAG 2.2 AA review is in the [accessibility review](https://github.com/piercetaylor/progeny-selector/blob/main/docs/accessibility.md). Design decisions are in [docs/adr](https://github.com/piercetaylor/progeny-selector/tree/main/docs/adr).

## Related tool

[backcross](https://github.com/piercetaylor/backcross) is the sibling browser-based tool. backcross characterises finished near-isogenic lines; progeny-selector ranks progeny during the programme. Both read input data contract 1.12.0, so genotype, samples.csv and markers.csv files move between them unchanged. The two tools release v0.1.0 on the same day, and each release note links the other's.

## Citing

Cite the software with [CITATION.cff](https://github.com/piercetaylor/progeny-selector/blob/main/CITATION.cff). GitHub shows a Cite this repository button for it. A Zenodo DOI is added after the first release. In plain text:

Taylor, P. (2026). progeny-selector (version 0.1.0) [Computer software]. https://github.com/piercetaylor/progeny-selector

## Licence

MIT. See [LICENSE](https://github.com/piercetaylor/progeny-selector/blob/main/LICENSE).
