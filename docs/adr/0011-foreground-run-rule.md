# Foreground rule `run`: an anchored contiguous run of donor calls

Status: accepted. Date: 2026-09-15.

## Context and Problem Statement

On the SoySNP50K Clark isolines, neither existing foreground rule separates an introgression from array noise. `rule: any` over gene ±1 Mb passed lines on 1–6 scattered donor calls, and `rule: all` in a tight window failed carriers whose donor segment covered only part of it. The maintainer (2026-09-15) asked for a rule that rejects scattered donor calls, proposed "at least 5 donor markers", and delegated the design to what standard academic breeding workflows use. What rule should a region target use?

## Considered Options

- A count threshold (`at_least` N donor markers in the region).
- Redefining `min_markers` under `any`.
- A fraction threshold (at least a given share of counted markers are donor).
- An anchored contiguous run: a run of donor calls through an anchor at the gene, of at least `min_run` calls.

Only contiguity plus an anchor at the gene separates array noise from an introgression, because scattered calls can reach any count or fraction.

## Decision Outcome

Chosen option: the anchored contiguous run, as `rule: run` on region targets with `min_run` (default 3), `anchor_bp` (default the region midpoint, `(start_bp + end_bp) // 2`) and `tolerate_isolated` (default true). N, U and X calls are removed first; fewer than `min_markers` counted calls gives unknown; the counted calls nearest the anchor on each side (a call exactly on the anchor serves as both) and every counted call tied on the anchor position must lie in the run, and bridged calls do not add to its length. `all` and `any` are unchanged. Semantics: `core/foreground.py`; keys: docs/data-formats.md.

Why `min_run` defaults to 3 and not the proposed 5: at Lf1 (Gm08) the ±1 Mb window has 4 informative markers, and 3 of the 4 Lf1 carriers carry exactly 4 donor markers even at ±2 Mb (the next informative marker is 1.3 Mb away), so 5 fails real carriers wherever array density is low. 3 is Hospital & Charcosset's (1997, Genetics 147:1469) "at least three markers per QTL". Measured against GRIN gene lists over 80 NIL × target calls in both families, with gene ±1 Mb and tolerance on: min_run 2 → 77/80, 3 → 76/80, 4 → 75/80, 5 → 71/80; `rule: any` was 67/80. At 3, every remaining disagreement is PI547634 (24 % donor genome-wide, an off-type: donor segments at T and R, none at pa1) or PI547592 at R. PI547592 has scattered single donor calls (11/81, longest run 2) and about 1 % donor genome-wide, which suggests the Higan seed lot used for the isoline differs from the genotyped PI548342 haplotype there; no rule recovers it.

Why isolated-mismatch tolerance is on by default, departing from the design consult's "no tolerance in v1": at pa1 on PI547481 and PI547482, a single recurrent-parent call next to the anchor splits a 28-marker donor segment. PLINK `--homozyg-window-het` is the precedent for tolerating isolated discordant calls in a run. Only single calls flanked by predicate calls are bridged, so two scattered donor calls separated by one A still form a run of 2, not a segment.

### Consequences

Good: scattered donor calls in a wide window no longer pass a target; real carriers pass at low array density; existing criteria files load and rank exactly as before. Bad: on strictly alternating calls, tolerance bridges every inner mismatch, so such a sample can reach `min_run` (a hand-built case in `tests/test_run_rule.py` documents it). Not done: a gap break (see isoline-browser docs/adr/0008 if added); a run-length column in results.csv, an output-format change that needs the maintainer. criteria.yaml is specific to this tool (docs/data-formats.md), not part of `contract/`, so isoline-browser is unaffected.

## More Information

- Hospital & Charcosset 1997: https://academic.oup.com/genetics/article-abstract/147/3/1469/6054126
- PLINK 1.9 `--homozyg`: https://www.cog-genomics.org/plink/1.9/ibd
- Flapjack MABC scoring: https://raw.githubusercontent.com/cropgeeks/flapjack/master/src/jhi/flapjack/analysis/MabcAnalysis.java
- Gilbert et al. 2023, Plant Genome: https://acsess.onlinelibrary.wiley.com/doi/10.1002/tpg2.20310
