# Foreground `rule: run` (one phase)

Maintainer decision, 2026-09-15: target regions need a rule that rejects scattered donor calls. The maintainer proposed "at least 5 donor markers" and delegated the design to what standard academic breeding workflows use. The decision, its evidence from the Clark isoline runs, and what was rejected go in `docs/adr/0011-foreground-run-rule.md`.

## Acceptance criteria

1. A target locus of kind `region` accepts `rule: run` with optional keys `min_run` (int >= 1, default 3), `anchor_bp` (int, default the region midpoint) and `tolerate_isolated` (bool, default true).
2. `all` and `any` behave exactly as before. Every existing criteria file loads and ranks with byte-identical results; the synthetic fixture and `tests/fixtures` are unchanged.
3. The run status follows the semantics below and is covered by hand-built unit tests, plus a randomised test that compares `locus_status` with an independent loop implementation written in the test file.
4. Criteria round-trip: `dump_criteria_yaml` then `read_criteria_text` gives an equal `Criteria`. The canonical dump writes `min_run`, `anchor_bp` and `tolerate_isolated` only for run loci, and always writes all three for them, with `anchor_bp` resolved to its default.
5. `docs/data-formats.md` documents the rule additively. `CHANGELOG.md` gets an Added line. ADR 0011 is written. All gates pass.

## Semantics (`core/foreground.py`)

Inputs per sample: the locus markers in position order (`resolved.marker_idx`, already sorted by position; assert or sort). The predicate comes from `marker_predicate(states, required_state)` (`hom_donor`: B; `het`: H; `either`: H or B). Counted calls are the called informative calls (A, H, B); N, U and X are removed before anything else, so they neither count nor break a run.

For one sample, let `c` be the counted calls in position order, with positions `p`, and let `n = len(c)`:

1. If `n < min_markers`, the status is **unknown** (the existing meaning of `min_markers`).
2. Let `sat[k]` = predicate true at counted call k. Let `joined[k] = sat[k]`. If `tolerate_isolated` is true, also set `joined[k] = True` for every k with `0 < k < n-1`, `not sat[k]`, `sat[k-1]` and `sat[k+1]`. Compute this from `sat` alone, never from already-joined neighbours.
3. Anchor neighbours: `a` = the largest k with `p[k] <= anchor_bp`; `b` = the smallest k with `p[k] >= anchor_bp`. They are the same index when one counted marker sits exactly on the anchor, and `a > b` when several counted markers share the anchor position. If either does not exist, the status is **fail**.
4. If `joined[k]` is true for every k from `min(a, b)` to `max(a, b)` inclusive, the run is the maximal contiguous stretch of `joined` containing that span. Run length = number of `sat` calls in that stretch; bridged calls are not counted. Otherwise the status is **fail**.
5. **pass** if the run exists and its length >= `min_run`. Otherwise **fail**.

Marker positions come from `gm.positions("bp")` (the region is in bp). `anchor_bp` must satisfy `start_bp <= anchor_bp <= end_bp`.

Implement it vectorised or with a per-sample loop over only the locus markers. Loci are small; clarity wins. Suggested signature change: `locus_status(states, resolved, predicate, rule="all", min_markers=1, *, positions=None, min_run=3, anchor_bp=None, tolerate_isolated=True)`, with `positions` required when `rule == "run"`. `foreground_status` gains the same keyword pass-through. `pipeline.py:85` passes them from the `TargetSpec` along with `gm.positions("bp")`. Keep the flanking override (`all`, 2) untouched.

## Validation (`model/criteria.py`, `io/criteria.py`)

- `LOCUS_RULES` stays `("all", "any")` for `LocusSpec` and `AvoidSpec`. Add `TARGET_RULES = ("all", "any", "run")`. `TargetSpec.validate` checks the target's rule against `TARGET_RULES`; restructure `LocusSpec.validate` so the base class does not reject `run` for targets. An avoid locus with `rule: run` raises `CriteriaError` naming the locus and stating that `run` is for targets only.
- `TargetSpec` gains `min_run: int = 3`, `anchor_bp: int | None = None` and `tolerate_isolated: bool = True`.
- `rule: run` on a target whose kind is not `region` raises `CriteriaError` ("rule run needs a region locus").
- `min_run < 1` raises. An `anchor_bp` outside `[start_bp, end_bp]` raises.
- In YAML, `min_run`, `anchor_bp` or `tolerate_isolated` present while the rule is not `run` raises `CriteriaError`, so a typo in `rule` never silently ignores them. Type checks use the existing `_check_int` and `_check_bool`. `anchor_bp` accepts thousands separators only if `start_bp` does in the region shorthand; otherwise plain int.
- These keys are target-only: reject them on avoid loci with the existing unknown-key mechanism, whatever it is.
- `_locus_to_dict` writes the three keys only for run targets, with `anchor_bp` resolved to `(start_bp + end_bp) // 2` when None.
- The UI criteria editor: if it offers a rule choice, add `run` only if it can edit the three keys. Otherwise leave the UI alone and make sure a loaded run criteria file is displayed and re-applied unchanged. Report which case applies.

## Tests (new `tests/test_run_rule.py`, plus additions to the existing criteria IO tests)

Hand-built cases (one sample column each; positions 1..N Mb; anchor between markers unless stated):
- `BBBBB` spanning the anchor, min_run 5: pass. Same with min_run 6: fail.
- `AABB|BAA` (run of 3 across the anchor), min_run 3: pass. `ABA|BAB` (alternating, anchor between the 3rd and 4th calls): with tolerance off, min_run 2 gives fail. With tolerance on, both inner A calls bridge, the run spans calls 2–6 with length 3 (B only), so min_run 3 gives pass and min_run 4 gives fail. This documents the known cost of tolerance on strictly alternating calls.
- `BBB|ABBB` with tolerance on: the run bridges the single A at the anchor neighbour b, length 6, pass at min_run 5. With tolerance off: fail.
- `BBAA|BBB` (two mismatches): no bridge; the anchor neighbour a is A, so fail.
- Run entirely right of the anchor (`AAA|ABBBBBB`): fail.
- N, U and X inside a run are ignored: `BBNBB` = length 4 run.
- An anchor exactly on a marker that is B, with B on both sides: pass. An anchor exactly on an A marker flanked by B, tolerance on: bridged, pass if the length is enough.
- No counted marker on one side of the anchor: fail. `n < min_markers`: unknown.
- `required_state: het` with H calls, and `either` mixing H and B.
- Randomised: 500 samples over random A/H/B/N/U/X strings of length 0–40, random anchors, min_run 1–6, both tolerance settings, all three required_states. Compare with a straightforward per-sample reference loop written independently in the test.
- Validation: run on a marker locus, a flanking locus and an avoid locus; min_run 0; anchor outside the region; `min_run` given with `rule: any`; round-trip equality; the dump omits the keys for `all` and `any` loci.
- Pipeline: a small in-memory dataset through `run_analysis` with one run target gives the expected target status column.

Do not change `scripts/make_fixture.py` or `tests/fixtures`. The run rule is a new status rule, not a new metric column, so the hand-built cases and the independent reference loop carry the weight (CLAUDE.md).

## Docs

- `docs/data-formats.md` line 23 area: add `run` to the rule list with its three keys, defaults, the anchor semantics and the isolated-mismatch bridge, in two or three sentences matching the file's style. State that `min_markers` keeps its meaning under `run`. Do not change any other documented key.
- `CHANGELOG.md` [Unreleased] > Added: one line.
- `docs/adr/0011-foreground-run-rule.md`, in the same MADR shape as 0010 (Status accepted, Date 2026-09-15; Context and Problem Statement; Considered Options; Decision Outcome; Consequences; More Information). Content to record:
  - Problem: on SoySNP50K Clark isolines, `rule: any` over gene ±1 Mb passed lines on 1–6 scattered donor calls; `rule: all` in a tight window failed carriers whose segment covered only part of it.
  - Options considered: a count threshold (`at_least` N donor markers); redefining `min_markers` under `any`; a fraction threshold; an anchored contiguous run. Only contiguity plus an anchor at the gene separates array noise from an introgression, because scattered calls can reach any count or fraction.
  - Why `min_run` defaults to 3 and not the proposed 5: at Lf1 (Gm08) the ±1 Mb window has 4 informative markers, and 3 of the 4 Lf1 carriers carry exactly 4 donor markers even at ±2 Mb (the next informative marker is 1.3 Mb away), so 5 fails real carriers wherever array density is low. 3 is Hospital & Charcosset's (1997, Genetics 147:1469) "at least three markers per QTL". Measured against GRIN gene lists over 80 NIL × target calls in both families, with gene ±1 Mb and tolerance on: min_run 2 → 77/80, 3 → 76/80, 4 → 75/80, 5 → 71/80; `rule: any` was 67/80. At 3, every remaining disagreement is PI547634 (24 % donor genome-wide, an off-type: donor segments at T and R, none at pa1) or PI547592 at R. PI547592 has scattered single donor calls (11/81, longest run 2) and about 1 % donor genome-wide, which suggests the Higan seed lot used for the isoline differs from the genotyped PI548342 haplotype there; no rule recovers it.
  - Why isolated-mismatch tolerance is on by default, and why this departs from the design consult's "no tolerance in v1": at pa1 on PI547481 and PI547482, a single recurrent-parent call next to the anchor splits a 28-marker donor segment. PLINK `--homozyg-window-het` is the precedent for tolerating isolated discordant calls in a run. Only single calls flanked by predicate calls are bridged, so two scattered donor calls separated by one A still form a run of 2, not a segment.
  - Not done: a gap break (see isoline-browser ADR 0008 if added); a run-length column in results.csv (an output-format change that needs the maintainer).
  - Sources: Hospital & Charcosset 1997 https://academic.oup.com/genetics/article-abstract/147/3/1469/6054126; PLINK 1.9 `--homozyg` https://www.cog-genomics.org/plink/1.9/ibd; Flapjack MABC scoring https://raw.githubusercontent.com/cropgeeks/flapjack/master/src/jhi/flapjack/analysis/MabcAnalysis.java; Gilbert et al. 2023 Plant Genome https://acsess.onlinelibrary.wiley.com/doi/10.1002/tpg2.20310.
  - criteria.yaml is specific to this tool (data-formats.md), not part of `contract/`, so isoline-browser is unaffected.
