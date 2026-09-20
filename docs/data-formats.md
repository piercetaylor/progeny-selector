# Data formats

This document is the contract for every file progeny-selector reads or writes. The input contract is `contract/data-contract.md`; criteria.yaml and the outputs are specific to this tool. Validation happens once, in `progeny_selector.io`; violations raise `DataContractError` or `CriteriaError` naming the file, line and column. A change to this document is a breaking change (CONTRIBUTING.md).

## Input contract

The input contract is `contract/data-contract.md`, version 1.4.0, mirrored byte for byte from backcross, where the canonical copy lives. It defines chromosome names (soybean-only in this version), the genotype file (VCF, HapMap, wide CSV; diploid calls only), the cell vocabulary per format (IUPAC codes expand to the heterozygote; `?`, stray `B`/`H`, `0`, `+`, `A?` are errors; `X`/`XX` are missing in HapMap, and in a wide CSV only under a token profile that lists them, such as `tassel`), samples.csv, markers.csv and the class codes; `contract/README.md` gives the version rules, `tests/test_contract_cases.py` runs every case under `contract/cases/` through `load_dataset`, and `scripts/check_contract.py` recomputes `MANIFEST.sha256` and, given `../backcross`, byte-compares the mirror with the canonical copy. backcross's docs/input-coding.md tabulates the accepted, missing and rejected codes per format and applies here unchanged. What this tool adds on top of the contract: the genotype format is chosen by file extension alone (`.vcf`, `.hmp.txt`, `.hmp`, `.hapmap`, `.csv`, `.tsv`, `.txt`, each optionally `.gz` or `.bgz`), which is the contract's minimum; for an A/B/H-coded file whose parents have no column the loader creates them internally (recurrent = all A, donor = all B, with a warning) and lists them in `Dataset.synthetic_sample_ids`; a pair of one nucleotide and one of N, `-`, `.` (`AN`, `A-`) is read as missing, which the contract has not yet decided; the chromosome lengths of the assembly named by `assembly` (`constants.py`) are the chromosome ends for drag bounds and marker weights when positions are in bp; `generation` values such as `BC2F1`, `BC3F2`, `F2`, `BC1`, `BC2S1` are parsed for expected values and unparseable strings are kept and flagged `generation_unparsed`; `family_id` groups progeny for per-family ranks and top-N selection; and markers absent from markers.csv have no cM, which disables cM mode for the whole dataset (all markers need cM).

## criteria.yaml (selection criteria)

Top-level keys: `name`, `targets`, `avoid`, `flank_window`, `flank_unit`, `assembly`, `background`, `ranking`, `weights`, `filters`. Unknown keys anywhere are an error. The UI's Download criteria.yaml writes this schema with every key explicit (defaults included); an optional key with no value is written as `null`, and regions are written as `chrom`/`start_bp`/`end_bp`; the file reloads unchanged.

### Locus definition (targets and avoid)

Exactly one of three forms per locus:

| form | keys |
|---|---|
| marker | `marker_id` |
| region | `chrom`, `start_bp`, `end_bp` (or shorthand `region: "Gm06:24,000,000-27,000,000"`) |
| flanking | `left_marker`, `right_marker` (both must satisfy the state) |

Common keys: `locus_id` (or `id`; unique across targets and avoid), `rule` (`all` default, or `any`: over called informative markers in the locus; or `run` on region targets, below), `min_markers` (default 1; fewer called markers gives status `unknown`), `notes`.

Targets on a region locus also accept `rule: run` (docs/adr/0011), with `min_run` (integer >= 1, default 3), `anchor_bp` (integer within the region, default `(start_bp + end_bp) // 2`) and `tolerate_isolated` (`true` default); these three keys are an error on avoid loci and under any other rule, and Download criteria.yaml writes all three for every run target. Over the called informative markers in position order, the locus passes when the counted markers nearest the anchor on each side, and every counted marker tied exactly on the anchor, lie in one contiguous run of markers meeting `required_state` and that run holds at least `min_run` of them; with `tolerate_isolated`, a single non-matching call between two matching calls joins the run without adding to its length. A counted marker exactly on the anchor serves as the neighbour on both sides; otherwise a sample with no counted marker on one side of the anchor fails. `min_markers` keeps its meaning under `run`.

### targets

| key | values | meaning |
|---|---|---|
| required_state | `hom_donor`, `het`, `either` (default) | donor state each counted marker must show; `either` = H or B |
| flank_left, flank_right | number in flank_unit | recombinant window per side; default `flank_window` |

### avoid

| key | values | meaning |
|---|---|---|
| allow_het | `false` (default) | when true, H passes; B still fails |

### flank_window and flank_unit

`flank_window` (default 5) in `flank_unit` (`cm` default, or `bp`). When `cm` is requested but the map has no cM, windows are interpreted in bp with a warning.

### assembly

The chromosome-length table used for chromosome ends: `Wm82.a1`, `Wm82.a2`, `Wm82.a4` (default) or `none`. Chromosome ends for RPP weights, drag bounds and strips when positions are in bp; `none` uses the last marker. A chromosome the assembly does not cover also falls back to its last marker, and both that fallback and a marker position beyond the assembly length are reported as warnings (docs/adr/0015).

### ranking

| key | values | default |
|---|---|---|
| mode | `weighted`, `staged` | `weighted` |

Under `staged`, hard filters still apply first, then the survivors are ordered: recombinant flanks (count) descending, `rpp_carrier` descending, `rpp_noncarrier` descending, `drag_total_est` ascending, `missing_rate` ascending, then `sample_id`, with exact ties at each key and no bins. `composite_score` is still written (docs/adr/0007).

### background

| key | values | default |
|---|---|---|
| model | `count`, `weighted` | `weighted` |
| map_unit | `auto`, `bp`, `cm` | `auto` (cM when available) |
| max_marker_coverage | number in map_unit | 10 cM or 4,000,000 bp |

### weights

Non-negative numbers; at least one positive. Components are defined in PLAN.md, algorithm 7.

| key | default |
|---|---|
| rpp_noncarrier | 0.5 |
| rpp_carrier | 0.2 |
| drag | 0.2 |
| recombinant | 0.1 |
| similarity_rp | 0.0 |
| completeness | 0.0 |

### filters

| key | default | meaning |
|---|---|---|
| max_missing_rate | 0.2 | individuals above this are excluded |
| unknown_target_is | `fail` | how an `unknown` target status is treated by the hard filter (`fail` or `pass`) |
| unknown_avoid_is | `pass` | same for avoid loci |
| exclude_qc_flagged | true | exclude individuals flagged possible_self_or_outcross, possible_outcross, possible_rp_sample, possible_donor_sample |
| max_hom_donor_rate_bcf1 | 0.02 | B-call rate above which a BCnF1 individual is flagged possible_self_or_outcross |
| max_nonparental_rate | 0.02 | X-call rate above which possible_outcross is flagged |
| het_rate_tolerance | 0.15 | absolute deviation from the expected het rate that flags het_rate_deviates (advisory) |
| parent_max_het_rate | 0.02 | parent heterozygosity above which a warning is issued |

Example:

```yaml
name: Rag1 introgression, BC2F1
targets:
  - locus_id: Rag1
    region: "Gm07:5,800,000-6,200,000"
    required_state: either
    rule: all
avoid:
  - locus_id: E1_donor
    marker_id: ss715607145
flank_window: 5
flank_unit: cm
assembly: Wm82.a4
ranking:
  mode: weighted
background:
  model: weighted
  max_marker_coverage: 10
weights:
  rpp_noncarrier: 0.5
  rpp_carrier: 0.2
  drag: 0.2
  recombinant: 0.1
filters:
  max_missing_rate: 0.2
```

The marker id in the example is illustrative; the tool does not ship marker lists.

## Outputs

### results.csv

One row per progeny/candidate, sorted by rank. The fixed prefix, in order (34 columns): `rank_overall, rank_in_family, sample_id, line_name, family_id, generation, passes_filters, exclusion_reason, composite_score, foreground_all_pass, avoid_all_pass, rpp_total, rpp_carrier, rpp_noncarrier, expected_rpp, drag_total_est, drag_total_max, drag_unit, ibs_rp, ibs_donor, missing_rate, het_rate, expected_het, qc_flags, role, n_informative_called, frac_a, frac_h, frac_b, background_model, background_unit, rank_mode, assembly, results_schema`. Then the dynamic columns: per target `target_<id>_status, drag_<id>_left_max, drag_<id>_right_max, recomb_<id>_left, recomb_<id>_right`, per avoid locus `avoid_<id>_status`, per chromosome `rpp_<chrom>`. Only the three summaries `rpp_total`, `rpp_carrier` and `rpp_noncarrier` share that prefix, so a reader may select the per-chromosome columns by `^rpp_` minus those three; a chromosome whose name would collide with a fixed column is an error naming the chromosome. Last `token_profile` (the token profile the genotype file was read with: `default`, a built-in id, or `custom:<id>`; contract 1.4.0). A dataset loaded without a profile writes `default` here.

The metadata columns record the run: `background_model` (`count` or `weighted`), `background_unit` (`bp` or `cm`, the unit the background analysis ran in: the unit of the marker weights under `model: weighted`, and under `model: count`, where markers are unweighted, the unit the map was read in and the one the drag bounds and flank windows use), `rank_mode` (`weighted` or `staged`), `assembly` (the criteria.yaml `assembly` value) and `results_schema`, the SemVer version of this format, `1.0.0`. A column added bumps the minor version; a column removed or renamed bumps the major.

Booleans are `TRUE`/`FALSE`; missing values are `NA`, while empty text (`exclusion_reason`, `qc_flags`, `notes`) is an empty cell, not `NA`; ranks are `NA` for excluded individuals; `exclusion_reason` is `;`-joined (`target:<id>:fail`, `avoid:<id>:fail`, `missing_rate>0.2`, `qc:<flag>`). The header row is written when there are no results, and then holds exactly the fixed prefix followed by `token_profile` (35 names, no dynamic columns).

### selected.csv

`sample_id, line_name, family_id, generation, rank_overall, rank_in_family, composite_score, rpp_total, notes, results_schema, token_profile`. Missing values are `NA` (a missing `family_id` or `generation`, an excluded individual's rank or score); `notes` is an empty cell when the breeder typed none. A note whose text is exactly `NA` is written as `NA` and reads back as a missing value in pandas and readr; quoting does not help (pandas reads a quoted `NA` as missing too, and `csv.writer` cannot quote one field without quoting every field), so avoid `NA` as a note. `results_schema` is the results.csv schema version, `1.0.0`. `token_profile` is copied from the results.csv row, so it is `default` for a dataset loaded without a profile and empty when the results.csv predates contract 1.4.0 and has no such column.

### next_samples.csv

A samples.csv skeleton for the next genotyping round: both parents, then one placeholder progeny row per selected individual (`<sample_id>-<next_generation>-001`, role `progeny`, `family_id` = the selected parent's id). Edit ids after planting; the file already satisfies this contract. Columns and empty cells are the samples.csv contract's, so this file carries no `results_schema` column and writes an empty cell, never `NA`. When more than one placeholder per selected individual is requested (the Export count), one row per planted plant is written and the file loads unchanged as the next generation's samples.csv.

## State and status codes

Parent-of-origin states: `A` recurrent homozygous, `H` heterozygous, `B` donor homozygous, `X` non-parental allele, `N` missing, `U` uninformative marker (numeric 0–5 in `constants.py`; identical definitions to the backcross classes rp_hom, het, donor_hom, nonparental, missing, uninformative). Locus statuses: `pass`, `fail`, `unknown`.
