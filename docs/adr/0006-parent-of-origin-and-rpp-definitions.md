# Six-state parent-of-origin classification and Flapjack-style RPP shared with isoline-browser

Status: accepted. Date: 2026-09-04.

## Context and Problem Statement

Every metric depends on classifying each progeny call relative to the two parents and on summarising those classes into recurrent-parent proportion. The definitions must be identical in both sibling tools, handle heterozygous or missing parents, multiallelic markers and coded input, and be robust to uneven marker density.

## Considered Options

1. Three states (A/H/B) with everything else treated as missing.
2. Six states: A, H, B, X (non-parental allele), N (missing), U (uninformative marker), with count-based RPP.
3. Option 2 plus a map-weighted RPP with a coverage cap.

## Decision Outcome

Option 3 (`core/classify.py`, `core/background.py`). A marker is informative only when both parents are called, homozygous and different; all calls at other markers are U. X is kept separate from N because non-parental alleles are the signal for outcrosses. RPP contribution is 1/0.5/0 for A/H/B with X, N and U excluded from both sums. The weighted model follows Flapjack: per side min(half the gap to the neighbouring informative marker, half the maximum coverage), chromosome ends bounded by assembly length when positions are in bp [web] https://flapjack.hutton.ac.uk/en/latest/mabc.html; defaults 10 cM or 4 Mb. Carrier and non-carrier chromosome RPP are reported separately because background recovery on the carrier chromosome is depressed by linkage to the target (Hospital and Charcosset 1997 [web] https://academic.oup.com/genetics/article-abstract/147/3/1469/6054126).

### Consequences

Good: identical class semantics in both tools (numeric codes differ, labels map one-to-one and are documented); weighted RPP is comparable with Flapjack output; expected values by generation give context. Bad: two RPP variants can confuse users, so the criteria file fixes one (`background.model`) and results.csv names it; without a cM map, windows and weights fall back to bp with a warning.
