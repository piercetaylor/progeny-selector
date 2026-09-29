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

## Amendment, 2026-09-27

Decided by research the maintainer delegated to Fable (Q-A and Q-B of 2026-09-27); the full record is docs/adr/0033. Two sentences above no longer hold.

- "chromosome ends bounded by assembly length when positions are in bp" now reads (revised 2026-09-29 by a second delegated Fable research question, recorded in docs/adr/0033):
  - the first marker's outer side: min(p, c/2) when positions are in bp, since the assembly origin is a chromosome end (VCF 4.3 places telomeres at POS 0 and N+1); c/2 when positions are in cM, since a linkage map's 0 cM is its first marker, not the telomere, and the distance to the telomere is unknown at both ends;
  - the last marker's outer side: min(max(L − p, 0), c/2) when positions are in bp and an assembly length L is known, else c/2.

  In bp the first marker used to get `c/2` whenever no length was known, which disagreed with backcross; in cM nothing changes.
- "defaults 10 cM or 4 Mb" now reads: defaults 10 cM or 2 Mb. 2 Mb is a soybean euchromatic translation of 10 cM (about 197 kb/cM, Schmutz et al. 2010, via backcross docs/adr/0006), not a Flapjack value; other crops should set the cap.
- In Consequences, "weighted RPP is comparable with Flapjack output" now reads: interior weighting follows Flapjack; totals differ at chromosome ends by design. Flapjack's ends are asymmetric: the first marker gets the full `min(pos, c)` and the last `min(mapLength - pos, c)`, with `mapLength` defaulting to the last marker's position, so its last marker normally gets zero.
