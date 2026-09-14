# Data formats

This document is the contract for every file progeny-selector reads or writes. The input contract (genotype file, samples.csv, markers.csv) is shared verbatim with the sibling isoline-browser project so files move between the tools unchanged; criteria.yaml and the outputs are specific to this tool. Validation happens once, in `progeny_selector.io`; violations raise `DataContractError` or `CriteriaError` naming the file, line and column. A change to this document is a breaking change (CONTRIBUTING.md).

## Chromosome names

`Gm01`..`Gm20`, `Gm1`..`Gm20`, `chr01`/`Chr1`/`chr1`, `chromosome6`, and bare `1`..`20` or `01`..`20` are normalised to `Gm01`..`Gm20`. Other names are kept unchanged and sorted after Gm20. Positions are 1-based bp on whichever Williams 82 assembly the files use; SoyBase lists Wm82.a1.v1.1, Wm82.a2.v1 and Wm82.a4.v1 with Gm01–Gm20 naming [web] https://www.soybase.org/resources/genome_info/. The Wm82.a4.v1 chromosome lengths in `constants.py` are used as chromosome ends for drag bounds and marker weights when positions are in bp.

## Genotype file

Format is chosen by extension: `.vcf` / `.vcf.gz` / `.vcf.bgz` (VCF), `.hmp.txt` / `.hmp` / `.hapmap` (HapMap), `.csv` / `.tsv` / `.txt` (wide CSV). Gzip and bgzip are read through Python's gzip module (bgzip is concatenated gzip members).

### VCF 4.2 or later

Columns `#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT` then one column per sample [web] https://samtools.github.io/hts-specs/VCFv4.2.pdf. Only CHROM, POS, ID, REF, ALT, FORMAT and GT are read; INFO, QUAL and FILTER are ignored (filter upstream with bcftools). GT indices refer to the REF,ALT list; `.` is missing; `/` and `|` are equivalent; haploid GT is read as homozygous; multiallelic ALT is supported. Records with ID `.` get `<chrom>_<pos>`.

### HapMap (TASSEL style)

Eleven fixed columns `rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode`, then one column per taxon [web] https://statgen-esalq.github.io/Hapmap-and-VCF-formats-and-its-integration-with-onemap/. Cells: two nucleotides (`AA`, `AT`), a slash pair (`A/T`) or one IUPAC letter (A, C, G, T; R, Y, S, W, K, M). `N`, `NN`, `-`, `--`, empty are missing. Tab-delimited.

### Wide CSV

| column | type | rule |
|---|---|---|
| marker_id | text | unique |
| chrom | text | any accepted chromosome spelling |
| pos_bp | integer | 1-based position |
| `<sample_id>` ... | text | one column per sample; the header is the sample_id used in samples.csv |

| coding | homozygous | heterozygous | missing |
|---|---|---|---|
| nucleotide | `A`, `AA` | `A/T`, `A\|T`, `AT`, IUPAC `W` | empty, `N`, `NA`, `-`, `.`, `./.` |
| abh | `A` (recurrent-parent allele), `B` (donor allele) | `H` | empty, `N`, `NA` |

`--coding auto` (default) selects `abh` when every non-missing cell in the first 200 rows is A, B or H and at least one B or H occurs; otherwise `nucleotide`. In `abh` coding allele 0 is A and allele 1 is B at every marker; if the parents named in samples.csv are absent from the file they are synthesised (recurrent = all A, donor = all B) with a warning.

Example (nucleotide):

```
marker_id,chrom,pos_bp,RP_Williams82,DONOR_PI,BC2F1-F1-001
syn_Gm06_12,Gm06,23435098,C,T,C
syn_Gm06_13,6,25472933,G,A,G/A
syn_Gm06_14,chr6,27510769,A,T,AT
```

## samples.csv (sample manifest)

| column | required | values |
|---|---|---|
| sample_id | yes | must match a genotype column (except parents of an abh-coded file); unique |
| line_name | no | display name; defaults to sample_id |
| role | yes | `recurrent_parent`, `donor_parent`, `candidate`, `progeny` |
| generation | no | `BC2F1`, `BC3F2`, `F2`, `BC1`, `BC2S1`; unparseable strings are kept and flagged `generation_unparsed` |
| family_id | no | grouping label used for per-family ranks and top-N selection |
| notes | no | free text |

Exactly one `recurrent_parent` and exactly one `donor_parent`; at least one `candidate` or `progeny` (both roles are ranked). Column names are case-insensitive. Genotype columns absent from the manifest are ignored with a warning; manifest samples absent from the genotype file are an error. Pyramiding with two donors is a later extension that adds `donor_parent_2` and per-target donor assignment; it does not change the meaning of existing columns.

Example:

```
sample_id,line_name,role,generation,family_id,notes
RP_Williams82,Williams 82 (synthetic),recurrent_parent,,,
DONOR_PI_synthetic,PI synthetic donor,donor_parent,,,
BC2F1-F1-001,LF1-001,progeny,BC2F1,F1,planted best
```

## markers.csv (marker map, optional)

| column | required | values |
|---|---|---|
| marker_id | yes | matches the genotype file |
| chrom | yes | any accepted spelling |
| pos_bp | yes | integer |
| cm | no | genetic position; enables cM windows, cM-weighted RPP and cM drag |

Map values override genotype-file positions when they differ (reported as warnings); markers absent from the map keep their genotype-file position and have no cM, which disables cM mode for the whole dataset (all markers need cM).

## criteria.yaml (selection criteria)

Top-level keys: `name`, `targets`, `avoid`, `flank_window`, `flank_unit`, `background`, `weights`, `filters`. Unknown keys anywhere are an error. The UI's Download criteria.yaml writes this schema with every key explicit (defaults included) and regions as `chrom`/`start_bp`/`end_bp`; the file reloads unchanged.

### Locus definition (targets and avoid)

Exactly one of three forms per locus:

| form | keys |
|---|---|
| marker | `marker_id` |
| region | `chrom`, `start_bp`, `end_bp` (or shorthand `region: "Gm06:24,000,000-27,000,000"`) |
| flanking | `left_marker`, `right_marker` (both must satisfy the state) |

Common keys: `locus_id` (or `id`; unique across targets and avoid), `rule` (`all` default, or `any`: over called informative markers in the locus), `min_markers` (default 1; fewer called markers gives status `unknown`), `notes`.

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

One row per progeny/candidate, sorted by rank. Leading columns: `rank_overall, rank_in_family, sample_id, line_name, family_id, generation, passes_filters, exclusion_reason, composite_score, foreground_all_pass, avoid_all_pass, rpp_total, rpp_carrier, rpp_noncarrier, expected_rpp, drag_total_est, drag_total_max, drag_unit, ibs_rp, ibs_donor, missing_rate, het_rate, expected_het, qc_flags`, then `role, n_informative_called, frac_a, frac_h, frac_b`, then per target `target_<id>_status, drag_<id>_left_max, drag_<id>_right_max, recomb_<id>_left, recomb_<id>_right`, per avoid locus `avoid_<id>_status`, and per chromosome `rpp_<chrom>`. Booleans are `TRUE`/`FALSE`; missing values are empty; ranks are empty for excluded individuals; `exclusion_reason` is `;`-joined (`target:<id>:fail`, `avoid:<id>:fail`, `missing_rate>0.2`, `qc:<flag>`).

### selected.csv

`sample_id, line_name, family_id, generation, rank_overall, rank_in_family, composite_score, rpp_total, notes`.

### next_samples.csv

A samples.csv skeleton for the next genotyping round: both parents, then one placeholder progeny row per selected individual (`<sample_id>-<next_generation>-001`, role `progeny`, `family_id` = the selected parent's id). Edit ids after planting; the file already satisfies this contract.

## State and status codes

Parent-of-origin states: `A` recurrent homozygous, `H` heterozygous, `B` donor homozygous, `X` non-parental allele, `N` missing, `U` uninformative marker (numeric 0–5 in `constants.py`; identical definitions to the isoline-browser classes rp_hom, het, donor_hom, nonparental, missing, uninformative). Locus statuses: `pass`, `fail`, `unknown`.
