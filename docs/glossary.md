# Glossary

Terms as progeny-selector uses them, in alphabetical order. Where a definition comes from the literature, the source is named; the ADR number is the decision record in [docs/adr](https://github.com/piercetaylor/progeny-selector/tree/main/docs/adr). The [tutorial](https://github.com/piercetaylor/progeny-selector/blob/main/docs/tutorial.md) shows most of these terms on the shipped example.

## Assembly

The reference genome release whose chromosome lengths mark the ends of chromosomes when marker positions are in base pairs: `Wm82.a1`, `Wm82.a2`, `Wm82.a4` or `none`. Under the soybean scheme the default is `Wm82.a4`; under any other crop it is `none`, because the length tables are soybean's. Lengths matter for the last marker's weight in the background score and for where a donor segment can end (docs/adr/0015, 0020, 0033). A warning appears when markers run past the recorded length, which usually means the data are on a different assembly than the one selected.

## Avoid locus

A marker, region or flanking pair where the donor allele is unwanted, such as a known undesirable gene near the target. By default the individual must be homozygous recurrent-parent there; `allow_het: true` lets a heterozygote pass, and homozygous donor always fails. The status is `pass`, `fail` or `unknown`. Flapjack, the tool this analysis follows, likewise lets a breeder declare for each locus whether the desired allele comes from the donor or the recurrent parent (docs/reference-repos.md; PLAN.md, section 5). By default an unknown avoid status passes, because missing data cannot show the donor allele is present (`unknown_avoid_is`).

## Background selection and recurrent-parent proportion (RPP)

Background selection picks, among plants that carry the target, those with the most recurrent-parent genome elsewhere. RPP is the measure: at informative markers, each A call counts 1, each H call 0.5 and each B call 0; X, N and U calls are left out. The count model is the plain average over called informative markers. The weighted model (Flapjack's MABC analysis, https://flapjack.hutton.ac.uk/en/latest/mabc.html) weights each marker by the map distance it covers, so a marker in a sparse region counts for more than one in a dense cluster: the weight is the smaller of half the gap to the neighbouring informative marker and half the maximum marker coverage, on each side, added together. The maximum coverage is `background.max_marker_coverage`, by default 10 cM or 2,000,000 bp, and the value applied is written to `results.csv`.

At chromosome ends the outer side of the first marker is the smaller of its position and half the cap when positions are in bp, and half the cap in cM. The last marker's outer side is the smaller of the distance to the assembly end and half the cap when the length is known, else half the cap. Interior weighting follows Flapjack; totals differ from Flapjack's at chromosome ends by design (docs/adr/0006, 0033). Background theory: Hospital and Charcosset 1997, Genetics 147:1469-1485.

## Carrier and non-carrier chromosome

A carrier chromosome holds at least one target locus; every other chromosome is non-carrier. RPP is reported for each group separately (`rpp_carrier`, `rpp_noncarrier`) because recovery of the recurrent parent on the carrier chromosome is held back by linkage to the target (Hospital and Charcosset 1997, Genetics 147:1469-1485). Frisch, Bohn and Melchinger 1999 (Crop Science 39:1295-1301) select recombinants on the carrier chromosome before the rest of the genome. See PLAN.md, section 3.

## Composite score

The weighted mean of up to six components, each between 0 and 1, computed for every individual: non-carrier RPP, carrier RPP, linkage drag, recombinant flanks, similarity to the recurrent parent and completeness (1 minus the missing rate). The default weights are 0.5, 0.2, 0.2, 0.1, 0 and 0. A component that cannot be computed drops out and the remaining weights are rescaled. The score orders only the individuals that pass the hard filters (docs/adr/0007).

## Crop scheme

The rule set that recognises chromosome names for a crop and puts them in order, so that `chr1`, `Chr01` and `1` are the same chromosome. Thirteen are built in: soybean (the default), maize, rice, sorghum, wheat, barley, oat, common-bean, cotton, cowpea, pea, peanut and sunflower. Choosing a crop changes names and order, never a coordinate; chromosome lengths are recorded for soybean only (docs/adr/0020).

## Duplicate

See `possible_duplicate`.

## Expected RPP by generation

The RPP a backcross plant is expected to have with no selection at all: 1 - (1/2)^(n+1) after n backcrosses, so BC1 is 0.75, BC2 0.875 and BC3 0.9375. Selfing after a backcross leaves expected RPP unchanged and halves heterozygosity each generation. Source: the Iowa State molecular plant breeding chapter on marker-assisted backcrossing cited in PLAN.md, section 3 (https://iastate.pressbooks.pub/molecularplantbreeding/chapter/marker-assisted-backcrossing/). Selection makes real plants deviate from this, which is the point of selecting.

## `family_donor_outlier`

A QC flag on an individual whose donor fraction (the rate of homozygous-donor calls plus half the heterozygous rate, at informative markers) is unusually high compared with the rest of its family. The cut is the family median plus 2.5 times a robust spread, in families of at least 6 individuals. It is advisory: a high donor fraction can be a real off-target introgression, so it never excludes an individual (docs/adr/0012; robust spread after Leys et al. 2013).

## Flanking window and recombinant flags

The flanking window is the distance on each side of a target within which a recombination counts as close (`flank_window`, default 5, in `cm` or `bp`; it can be set per target and side). For each side, the recombinant flag (`recomb_<id>_left`, `recomb_<id>_right`) is true when the maximum donor bound on that side falls inside the window, meaning the donor segment was cut short there. This uses Flapjack's "first recombination on each side" as the boundary of the donor segment (docs/adr/0006; PLAN.md, section 4).

## Foreground selection and `required_state`

Foreground selection keeps plants that carry the target allele. `required_state` says which state each counted marker must show: `hom_donor` (B), `het` (H) or `either` (H or B, the default). The `rule` says how to combine markers in a region: `all` (default), `any`, or `run` (docs/adr/0011). Under `run`, the markers nearest an anchor position must lie in one unbroken stretch of matching calls at least `min_run` long (default 3, following the "at least three markers" of Hospital and Charcosset 1997), with a single stray call tolerated between two matching calls. If fewer than `min_markers` markers are called the status is `unknown`. See PLAN.md, section 2.

## Generation labels

`BCnFm` means n backcrosses to the recurrent parent and m generations of selfing since the last one: BC2F1 is the offspring of the second backcross. The label determines expected RPP and expected heterozygosity. A single-cross Purdy label with one `/` and one `*` is also read, such as `RP*3/DONOR`; it gives expected RPP but no filial generation, so heterozygosity checks are skipped for it (docs/adr/0034). A label that cannot be read is kept and flagged `generation_unparsed`.

## `generation_unparsed`

A QC flag: the `generation` text could not be read as a generation label, so no expected RPP or heterozygosity is available for the individual. Fix the label in `samples.csv` (PLAN.md, section 8).

## Genotype states

Each progeny call, at a marker where the two parents differ, is classified as one of six states:

- **A**: homozygous recurrent parent (backcross label `rp_hom`).
- **H**: heterozygous, one allele from each parent (`het`).
- **B**: homozygous donor (`donor_hom`).
- **X**: an allele found in neither parent (`nonparental`).
- **N**: missing call (`missing`).
- **U**: uninformative marker (`uninformative`).

The letters A, H and B follow the ABHgenotypeR coding convention (docs/reference-repos.md); in a file already coded that way, allele 0 means recurrent and 1 means donor. X is kept apart from N because non-parental alleles point to an outcross. The same definitions are used by the sibling tool backcross (docs/adr/0006; PLAN.md, section 1).

## Hard filters

Conditions an individual must meet to be ranked at all, applied before any scoring: every target passes, every avoid locus passes, the missing rate is at most `max_missing_rate` (default 0.2), and no excluding QC flag (`possible_self_or_outcross`, `possible_outcross`, `possible_rp_sample`, `possible_donor_sample`) is present. An excluded individual keeps its metrics and a reason in `exclusion_reason`, so no one is ranked above another because a high background score hid a failed target (docs/adr/0007).

## `het_rate_deviates`

An advisory QC flag: the observed heterozygous rate differs from the rate expected for the generation label by more than 0.15. It does not exclude, because a plant selected for high recurrent-parent content legitimately departs from the expectation (PLAN.md, section 8).

## `high_missing`

A QC flag on an individual with too many missing calls. It uses the same rule as the missing-rate hard filter (`filters.max_missing_rate`, default 0.2): one rule gives a QC flag, which does not exclude, and a hard-filter reason, which does (PLAN.md, section 8).

## IBS (identity by state)

The fraction of alleles two individuals share, averaged over markers called in both: 1 for two matching genotypes, 0.5 when they share one allele, 0 when none. The tool compares each progeny with each parent (`ibs_rp`, `ibs_donor`) and progeny with each other. For biallelic markers it equals the measure of the same name in the SNPRelate package (PLAN.md, section 6).

## Informative and uninformative marker

A marker is informative when both parents are called, both are homozygous, and their alleles differ. Every call at any other marker is U, uninformative, with a recorded reason: the parents are identical, the recurrent or donor parent is missing, or one parent is heterozygous. The Validate and QC screen counts uninformative markers by reason. PLAN.md, section 1.

## Linkage drag, minimum and maximum bounds, and estimate

Linkage drag is the length of donor chromosome that comes along with the target because it is linked to it. Starting at the target and moving outward over informative markers, skipping N and X calls, the maximum bound on a side is the distance to the first A (recurrent) marker, or to the chromosome end if there is none. The minimum bound is the distance to the outermost H or B marker before that A marker. The truth lies between the two; the estimate is their midpoint, and the total is the left side plus the right side. The maximum bound follows Flapjack's "first recombination on each side" (PLAN.md, section 4).

## `possible_duplicate`

A QC flag, set on both members of a pair of progeny whose IBS is at least 0.995, measured over the informative markers called in both. It suggests one tube may have been sampled twice. It is advisory and never excludes, because the tool cannot tell which of the two is the real plant. Every pair of progeny is compared regardless of family, so the flag does not know the two are supposed to be related. The 0.995 cut is deliberately strict; nonetheless very high recurrent-parent content makes real siblings similar. In the shipped BC3F1 example, 8 individuals descended from `BC2F1-F1-001`, a plant with very high recurrent-parent content, carry this flag, while the other three families have none. The scan is skipped above 2,000 individuals (docs/adr/0017).

## `possible_self_or_outcross`, `possible_outcross`, `possible_rp_sample`, `possible_donor_sample`

QC flags that exclude an individual by default (`exclude_qc_flagged: true`).

- `possible_self_or_outcross`: a BCnF1 plant shows more than 2 % homozygous-donor calls, which a true backcross cannot produce apart from error, so it may be a selfed or outcrossed plant.
- `possible_outcross`: more than 2 % of calls are non-parental (X).
- `possible_rp_sample`, `possible_donor_sample`: IBS to that parent is above 0.995 and the heterozygous rate is below 0.5 %, so the sample may be the parent itself.

PLAN.md, section 8.

## `results_schema`

The version of the `results.csv` column layout, written on every row and currently `1.2.0`. Adding a column raises the minor number; removing or renaming one raises the major number, so a reader that accepts `1.x` keeps working when columns are added (docs/adr/0016, 0033).

## Staged and weighted ranking

Weighted ranking (the default) orders passing individuals by the composite score. Staged ranking (`ranking: {mode: staged}`) orders them by a fixed sequence of criteria: recombinant flanks (more first), carrier RPP, non-carrier RPP, estimated drag (less first), missing rate (less first), then `sample_id`, with no bins. The staged order reflects Frisch, Bohn and Melchinger 1999, whose four-stage scheme selects carrier-chromosome recombinants before the background. Their staging concerns which markers are genotyped in each generation, so it saves genotyping only when that is done stage by stage; here the full file is read either way (docs/adr/0007, amendment of 2026-09-16). `rank_mode` in `results.csv` records which was used.

## Token profile

The vocabulary used to read cells of a HapMap or wide-CSV genotype file: which text means missing, which pairs mean heterozygous. Built-in profiles are `tassel`, `soybase-report`, `dart`, `axiom` and `kasp`, or a JSON file of your own. The profile used is written to the `token_profile` column of `results.csv` (docs/adr/0014).
