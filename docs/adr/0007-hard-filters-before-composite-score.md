# Hard filters first, then a weighted composite of bounded components, with a fixed tie-break order

Status: accepted. Date: 2026-09-04.

## Context and Problem Statement

Breeders select progeny by a sequence of decisions (must carry the target, must not carry the avoid allele, then best background and shortest drag), not by one number. A ranking tool nevertheless needs a total order. How are hard requirements, soft preferences and ties combined?

## Considered Options

1. A single weighted score including penalties for failing targets.
2. Hard filters (foreground, avoid, missing rate, QC exclusions) applied first; a weighted mean of components in [0, 1] ranks the survivors; excluded individuals are listed with reasons.
3. Lexicographic ranking on fixed criteria without weights.

## Decision Outcome

Option 2 (`core/score.py`), which is the four-stage logic of Frisch, Bohn and Melchinger 1999 (recombinants on the carrier chromosome before genome-wide background) expressed as filters and weights [web] https://experts.illinois.edu/en/publications/comparison-of-selection-strategies-for-marker-assisted-backcrossi/. Components: rpp_noncarrier, rpp_carrier, drag (1 − estimated donor segment / carrier chromosome length), recombinant (fraction of flanks with a recombination inside the window), similarity_rp, completeness. NaN components drop out with weight renormalisation. Ties break by RPP total, then shorter drag, then lower missing rate, then sample_id. Ranks are dense, overall and within family. Policies for `unknown` statuses are explicit in criteria.yaml because a missing call at the target is a genotyping decision, not a scoring one.

### Consequences

Good: an individual is never ranked above another because a high background score masked a failed target; exclusion reasons are auditable in results.csv; weights are visible and versioned. Bad: the drag normalisation by chromosome length makes the component small for most individuals (weights must be read with that in mind); users who want lexicographic behaviour approximate it with extreme weights. Neutral: the projected next-generation expectation assumes Mendelian segregation with unlinked loci and no selection, stated in `core/selection.py`.

## Amendments, 2026-09-16

**Staged ranking as an optional second mode; weighted stays the default.** Maintainer decision, recorded before implementation and planned for M2. criteria.yaml gains `ranking: {mode: weighted|staged}`, default `weighted`, so every existing criteria file ranks unchanged. Under `staged`, hard filters still apply first, and survivors are ordered lexicographically: recombinant flanks (count) descending, `rpp_carrier` descending, `rpp_noncarrier` descending, `drag_total_est` ascending, `missing_rate` ascending, then `sample_id`, with exact ties at each key and no bins. `composite_score` is still written in staged mode, and results.csv gains a `rank_mode` column so `rank_overall` can be read downstream.

Why weighted stays the default: Frisch, Bohn and Melchinger 1999 found the staged scheme cut marker data points by about 75 % against a genome-wide index, a saving that applies when progeny are genotyped stage by stage. This tool reads a complete genotype file, so the saving does not arise here, and index selection is also current practice (Herzog, Falke and Frisch 2014, https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0092429). Flapjack's MABC view computes no composite; users sort on several columns and assign ranks by hand (https://flapjack.hutton.ac.uk/en/latest/mabc_tutorial.html). The staged mode serves breeders who want that order without approximating it through extreme weights, which the Consequences above name as a limitation.

The Decision Outcome's reading of option 2 as "the four-stage logic of Frisch, Bohn and Melchinger 1999 expressed as filters and weights" is loose: their staging concerns which markers are genotyped in which generation, not how a fully genotyped population is ordered.
