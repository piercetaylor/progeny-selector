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
