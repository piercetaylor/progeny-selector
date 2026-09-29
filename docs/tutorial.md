# One run on the shipped example

This tutorial takes one backcross population through progeny-selector from start to finish. The population is synthetic data that ships with the package, so nothing here is a breeding recommendation. Every command below is run by the test suite (`tests/test_docs_commands.py`), and every number quoted was produced by that command on version 0.1.0.

Terms are defined in the [glossary](https://github.com/piercetaylor/progeny-selector/blob/main/docs/glossary.md).

## Get the example

Write the example files into a folder called `example` with `progeny-selector example --out example`. The command-line section below starts with this same command, so run it once only. It creates two sub-folders, each holding `genotypes.vcf`, `samples.csv`, `markers.csv` and `criteria.yaml`:

- `example/synthetic_bc2f1`: the BC2F1 population this tutorial ranks.
- `example/synthetic_bc3f1`: a BC3F1 population, the offspring of selected BC2F1 plants, used in the last command section.

In the browser interface the same BC2F1 files load with one button: on screen 1 Load, press "Load example (synthetic BC2F1)". Nothing is uploaded. To open the interface from the command line, run `progeny-selector app`; Ctrl+C stops it.

## Command line

The four basic commands, in order. `validate` checks the inputs and prints QC warnings, `rank` writes `results.csv`, and `select` picks the top 3 in each family and writes the files for the next round.

```sh
progeny-selector example --out example
progeny-selector validate --genotypes example/synthetic_bc2f1/genotypes.vcf --samples example/synthetic_bc2f1/samples.csv --markers example/synthetic_bc2f1/markers.csv --criteria example/synthetic_bc2f1/criteria.yaml
progeny-selector rank --genotypes example/synthetic_bc2f1/genotypes.vcf --samples example/synthetic_bc2f1/samples.csv --markers example/synthetic_bc2f1/markers.csv --criteria example/synthetic_bc2f1/criteria.yaml --out results.csv
progeny-selector select --results results.csv --top 3 --out selected.csv --next-manifest next_samples.csv --next-generation BC3F1 --samples example/synthetic_bc2f1/samples.csv
```

What they print:

- `validate` reports 500 markers by 42 samples, 475 informative markers, 7 QC-flagged individuals and no duplicate pairs.
- `rank` reports 40 individuals, 11 of which pass the hard filters.
- `select` reports 6 selected: `--top 3` means every passing individual ranked 3rd or better within its family, so up to 3 per family, and this population has two families. Ranks are dense, so a tie can return more than 3 in a family. With `--overall`, `--top 3` would instead give the 3 best overall.

It also prints the projection for the next generation: the mean recurrent-parent proportion (RPP) of the selected plants is 0.869 and the expected RPP of their offspring is 0.934.

Run `select` again with ten placeholder rows per selected plant. This replaces `selected.csv` and `next_samples.csv` with a manifest of 60 planted-plant rows (6 selected times 10) plus the two parents:

```sh
progeny-selector select --results results.csv --top 3 --out selected.csv --next-manifest next_samples.csv --next-generation BC3F1 --samples example/synthetic_bc2f1/samples.csv --per-selected 10
```

Finally, rank the next generation. The shipped `synthetic_bc3f1` folder shows the shape of the files for a BC3F1 generation. It is a separate synthetic dataset (4 parents x 10 = 40 progeny), not generated from the manifest you just wrote, so its sample IDs and row count differ from your `next_samples.csv`:

```sh
progeny-selector rank --genotypes example/synthetic_bc3f1/genotypes.vcf --samples example/synthetic_bc3f1/samples.csv --markers example/synthetic_bc3f1/markers.csv --criteria example/synthetic_bc3f1/criteria.yaml --out results_bc3f1.csv
```

This reports 40 individuals, 21 of which pass. Running `validate` on the same files flags 8 individuals as `possible_duplicate`, all in the family descended from `BC2F1-F1-001`. That flag is expected here and is explained in the glossary under `possible_duplicate`.

## What the example contains

The BC2F1 population is a soybean backcross with the following properties.

- 2 parents (`RP_Williams82`, the recurrent parent, and `DONOR_PI_synthetic`, the donor) and 40 BC2F1 progeny in families `F1` and `F2`.
- 500 markers, 25 on each of Gm01 to Gm20, with a genetic map in cM. 475 markers are informative (both parents called, homozygous and different).
- Target locus `T1` at marker `syn_Gm06_13`. The donor allele is required (`required_state: either`, so heterozygous or homozygous donor passes).
- Avoid locus `AV1` at marker `syn_Gm13_10`. The plant must be recurrent-parent homozygous there.
- Flank windows of 6 cM on each side of the target.
- Background model `count`, so RPP is the plain proportion of recurrent-parent alleles.
- Weights 0.5 (non-carrier RPP), 0.2 (carrier RPP), 0.2 (linkage drag) and 0.1 (recombinant flanks). Filters: missing rate at most 0.2; an unknown target status fails; an unknown avoid status passes; QC-excluded plants are excluded.

Five individuals are planted to show specific behaviour:

| individual | what was planted | what results.csv shows |
| --- | --- | --- |
| `BC2F1-F1-001` | Heterozygous at the target with a short donor segment, 28 heterozygous markers elsewhere, no missing calls. | Rank 1, composite score 0.911, RPP 0.971. |
| `BC2F1-F1-002` | Recurrent-parent homozygous at the target. | `exclusion_reason` is `target:T1:fail`; not ranked. |
| `BC2F1-F2-001` | Donor allele at the avoid locus. | `exclusion_reason` is `avoid:AV1:fail`; not ranked. |
| `BC2F1-F2-002` | A selfed plant mislabelled as BC2F1, so it has homozygous-donor calls. | `qc_flags` is `possible_self_or_outcross|family_donor_outlier`; excluded by `qc:possible_self_or_outcross`. |
| `BC2F1-F1-003` | 30 % of calls missing. | `qc_flags` is `high_missing`; `exclusion_reason` is `missing_rate>0.2`. |

Not ranked means `rank_overall` is `NA`. An excluded plant keeps all its other metrics.

## The screens

The interface has seven screens in a navigation bar. Data stays in your browser tab or on your machine; it is not uploaded.

### 1 Load

Choose the genotype file ("Genotypes (VCF, VCF.gz, HapMap, wide CSV)"), `samples.csv`, optionally `markers.csv (optional)`, and `criteria.yaml`, then press "Load and analyse". The "Token profile" select (and "Custom token profile (JSON, optional)", cleared by "Clear custom token profile") tells the reader how to interpret HapMap and wide-CSV cells. The "Crop" select chooses how chromosome names are recognised; soybean is the default. The "Status" card shows what happened, including any error from the data contract, verbatim.

"Load example (synthetic BC2F1)" loads and analyses the shipped example without any upload. The "Criteria (editable)" card shows the criteria in canonical form. Edit it and press "Apply criteria and re-analyse" to rerun; "Download criteria.yaml" saves the applied criteria.

### 2 Validate and QC

The "Summary" card lists markers, samples, progeny, informative markers, uninformative markers by reason, map unit, background model, drag unit, assembly, rank mode, duplicate pairs and all warnings. The "Per-individual QC" table lists each individual's QC flags; excluded and flagged rows are tinted, and the column filters narrow the list.

### 3 Navigate

Pick a family, then a generation, from the accordion. The selection appears as a breadcrumb (cross, family, generation) and filters the Rank table. "(no generation)" appears when some individuals have no generation label. The breadcrumb links step back up.

### 4 Rank

"Ranked individuals (multi-select rows, then open Compare)". Each row is an individual, with `rank_overall`, `rank_in_family`, `sample_id`, `family_id`, `generation`, `passes_filters`, `exclusion_reason`, `composite_score`, `rpp_total`, `rpp_carrier`, `rpp_noncarrier`, `drag_total_est`, `missing_rate` and `qc_flags`. Then, for every target and avoid locus, a status column (`target_T1_status`, `avoid_AV1_status`) and, for targets, the two recombinant columns (`recomb_T1_left`, `recomb_T1_right`).

Status cells are coloured chips with the word kept in the cell: pass, fail or unknown. Colours are from the Okabe-Ito palette, chosen to be distinguishable with common colour blindness, and the text means the colour is never the only signal.

The switch "Show only individuals passing hard filters" is on by default. Turn it off to see excluded individuals and their reasons. The caption above the table says how many individuals are shown and for which family and generation.

Select rows by clicking one, or by Ctrl-clicking (Cmd-clicking on macOS) to add or remove single rows. The selected rows feed screens 5 and 6.

### 5 Compare

"Compare selected individuals" shows up to 6 of the selected rows as cards. Each card gives:

- the individual's rank and score;
- a chip for each target and avoid locus and its status;
- a table of RPP per chromosome;
- a table of donor-segment bounds for each target: `left min`, `left max`, `right min`, `right max`, `total estimate`, `total max`, and `recombinant left` and `recombinant right` (yes or no). The minimum bound is the outermost marker known to carry donor; the maximum bound is the first marker known to be recurrent. The truth lies between them, and the estimate is the midpoint;
- a chromosome strip: one row per chromosome (20 for soybean), coloured by parent-of-origin state in Okabe-Ito colours, with a tick at each target and avoid locus. The legend above the cards names the six states.

If more than six rows are selected, the screen says "showing 6 of N".

### 6 Selection list

Selected individuals appear here with `sample_id`, `line_name`, `family_id`, `rank_overall`, `composite_score`, `rpp_total` and a `notes` column you can type into (click a cell in it). Notes go into `selected.csv`. Avoid the text `NA` as a note, because spreadsheet software and R read it as a missing value.

"Top N per family" with "Add top N per family" adds the passing individuals ranked N or better within their family to the list. With N = 3 on this example the list holds 6 individuals, 3 in each of the two families, the same 6 that `select --top 3` writes. The button adds to the list; it does not clear it.

"Next step" chooses "Backcross to RP" or "Self", and the summary line projects the next generation from the selected group. Backcrossing to the recurrent parent halves the donor content: expected RPP of the offspring is (1 + RPP) / 2 where RPP is the mean of the group. With mean RPP 0.869, the projection is 0.934. The projection assumes Mendelian segregation, unlinked loci and no selection in the next generation.

### 7 Export

Set "Next generation label" (default `BC3F1`) and "Placeholder rows per selected individual" (default 1), then download `results.csv (ranks and statuses)`, `selected.csv`, and `next_samples.csv (manifest for the next genotyping round)`. These are the same files the command line writes.

## Reading results.csv

One row per progeny individual, sorted by rank. Booleans are `TRUE` and `FALSE`. Missing numbers are `NA`; an empty cell means empty text (or an identifier with no source, such as a BrAPI id in a file-loaded dataset). The 41 fixed columns come first, in this order.

| column | meaning |
| --- | --- |
| `rank_overall` | Rank among all ranked individuals; ties share a rank. `NA` if excluded. |
| `rank_in_family` | Rank within the individual's family. `NA` if excluded. |
| `sample_id` | Identifier from `samples.csv`. |
| `line_name` | Line name from `samples.csv`. |
| `family_id` | Family from `samples.csv`. |
| `generation` | Generation label from `samples.csv`. |
| `passes_filters` | Whether the individual passed every hard filter. |
| `exclusion_reason` | Why it failed, joined by `;`, such as `target:T1:fail`, `avoid:AV1:fail`, `missing_rate>0.2`, `qc:possible_self_or_outcross`. |
| `composite_score` | Weighted mean of the score components, from 0 to 1. Computed for excluded individuals too. |
| `foreground_all_pass` | Every target passed. |
| `avoid_all_pass` | Every avoid locus passed. |
| `rpp_total` | Recurrent-parent proportion over all informative markers. |
| `rpp_carrier` | RPP on chromosomes that hold a target. |
| `rpp_noncarrier` | RPP on all other chromosomes. |
| `expected_rpp` | RPP expected for the generation label, with no selection. |
| `drag_total_est` | Estimated total donor segment length around the targets (left plus right, midpoint of min and max). |
| `drag_total_max` | Maximum total donor segment length around the targets. |
| `drag_unit` | Unit of the two drag columns: `cm` or `bp`. |
| `ibs_rp` | Allele sharing with the recurrent parent. |
| `ibs_donor` | Allele sharing with the donor parent. |
| `missing_rate` | Fraction of all markers with no call. |
| `het_rate` | Heterozygous fraction of called informative markers. |
| `expected_het` | Heterozygous fraction expected for the generation label. |
| `qc_flags` | QC flags joined by `|`. |
| `role` | `progeny` or `candidate`. |
| `n_informative_called` | Informative markers with a call. |
| `frac_a`, `frac_h`, `frac_b` | Fractions of called informative markers that are A, H and B. |
| `background_model` | `count` or `weighted`. |
| `background_unit` | Unit of the background model: `cm` or `bp`. |
| `rank_mode` | `weighted` or `staged`. |
| `assembly` | Chromosome-length table used for chromosome ends. |
| `results_schema` | Version of this file layout (`1.2.0`). |
| `crop` | Chromosome scheme used to read chromosome names. |
| `call_set_db_id` | BrAPI call-set id; empty for a file-loaded dataset. |
| `sample_db_id` | BrAPI sample id; empty for a file-loaded dataset. |
| `background_max_marker_coverage` | Cap the weighted model applied, in `background_unit`; `NA` under the count model. |
| `tool` | `progeny-selector`. |
| `tool_version` | Version that wrote the file. |
| `tool_commit` | Source commit: `g` plus seven hex digits, `-dirty` if the working tree differed from it, `NA` if unknown. |

After the fixed columns come the dynamic columns, named by locus id and chromosome. The prefixes are those in [docs/data-formats.md](https://github.com/piercetaylor/progeny-selector/blob/main/docs/data-formats.md):

| column family | meaning |
| --- | --- |
| `target_<id>_status` | `pass`, `fail` or `unknown` for target `<id>`. |
| `drag_<id>_left_max`, `drag_<id>_right_max` | Maximum donor segment on each side of target `<id>`, in `drag_unit`. |
| `recomb_<id>_left`, `recomb_<id>_right` | Whether a recombination lies within the flank window on that side of target `<id>`. |
| `avoid_<id>_status` | `pass`, `fail` or `unknown` for avoid locus `<id>`. |
| `rpp_<chrom>` | RPP on one chromosome, such as `rpp_Gm06`. |

The last column is `token_profile`, the genotype vocabulary the file was read with (`default` for a VCF). In this example the target is `T1`, so the dynamic columns are `target_T1_status`, `drag_T1_left_max`, `drag_T1_right_max`, `recomb_T1_left`, `recomb_T1_right`, `avoid_AV1_status` and `rpp_Gm01` to `rpp_Gm20`. For `BC2F1-F1-001`, `drag_T1_left_max` is 5.094 cM and `drag_T1_right_max` is 10.19 cM, and `recomb_T1_left` is `TRUE` (5.094 is inside the 6 cM window) while `recomb_T1_right` is `FALSE`.

`selected.csv` carries `sample_id, line_name, family_id, generation, rank_overall, rank_in_family, composite_score, rpp_total, notes, results_schema, crop, call_set_db_id, sample_db_id, tool, tool_version, tool_commit, token_profile`. The complete column definitions are in docs/data-formats.md.

## Next generation

`next_samples.csv` is the manifest for the next round. It lists both parents and one placeholder row for each planted plant: `<sample_id>-<next generation>-001`, role `progeny`, `family_id` set to the selected parent's id. In the first `select` run above it holds 6 placeholder rows; in the run with `--per-selected 10` it holds 60.

Once you have planted and tagged your plants, edit the `sample_id` values to your real plant tags and genotype the plants. The file already satisfies the `samples.csv` contract, so it loads unchanged as the next generation's `samples.csv`. `generation` is `BC3F1` because the run passed `--next-generation BC3F1`. If you are selecting a self-generation instead, the Selection list screen's "Next step" projects that case, and `select --project self` does the same on the command line.

`example/synthetic_bc3f1/samples.csv` is a separate, ready-made BC3F1 example (40 progeny; its sample IDs do not match your manifest, which has 60 rows), with rows such as `BC2F1-F1-001-BC3F1-001` with role `progeny`, generation `BC3F1` and `family_id` `BC2F1-F1-001`. Use it with the `criteria.yaml` from the previous round to rank the new generation, as in the last command above.
