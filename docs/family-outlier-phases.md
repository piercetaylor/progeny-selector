# QC flag `family_donor_outlier` (one phase)

Maintainer decision, 2026-09-16: flag an individual whose donor-genome fraction is far above the rest of its family. Advisory only. The decision and its evidence go in `docs/adr/0012-family-donor-outlier.md`.

## Acceptance criteria

1. `sample_qc` adds the flag `family_donor_outlier` by the rule below. No other flag changes.
2. The flag is advisory: it is not in `QC_EXCLUDING_FLAGS`, so `exclude_qc_flagged` and the hard filters behave exactly as before, and `qc_excluded` stays false for an individual carrying only this flag.
3. No new criteria.yaml key, no new results.csv column, no change to any documented column, parameter or filter. The thresholds are module constants in `core/qc.py`, like the 0.995 IBS and 0.005 het thresholds already there.
4. The fixture is regenerated from `scripts/make_fixture.py` with an independent implementation of the statistic, and the regenerated files are deterministic (running the generator twice gives identical bytes).
5. Hand-built unit tests cover the rule. All gates pass.

## Rule (`core/qc.py`)

Measure, per individual, genome-wide:

    donor_fraction = hom_donor_rate + het_rate / 2

The denominator is the A + H + B calls at informative markers; X (non-parental) calls are excluded, and `possible_outcross` covers those individuals. Both rates are already on `SampleQC`. Genome-wide, from the count model, not the map-weighted RPP: the weighted measure inflates individuals whose donor calls are scattered singletons, which is exactly the noise `rule: run` was added to ignore (docs/adr/0011).

Grouping: by `Sample.family_id` exactly as stored, with the empty string or None its own group, so a dataset with no families is one group.

Per group, over the individuals whose `donor_fraction` is not NaN (an individual with no called informative marker has NaN rates and is never flagged and never contributes; nor does an individual already flagged `high_missing`, whose few calls make its fraction unreliable):

1. If the group has fewer than `FAMILY_OUTLIER_MIN_N = 6` such individuals, flag nobody in it.
2. `med` = median of `donor_fraction`. `mad` = median of `abs(donor_fraction - med)`. `scale = max(1.4826 * mad, FAMILY_OUTLIER_MIN_SCALE)` with `FAMILY_OUTLIER_MIN_SCALE = 0.01`, the floor so a very uniform family does not flag small deviations.
3. Flag `family_donor_outlier` when `donor_fraction > med + FAMILY_OUTLIER_Z * scale`, with `FAMILY_OUTLIER_Z = 2.5` (Leys et al. 2013). One-sided: a low donor fraction is never flagged here; `possible_rp_sample` covers that case.

Comparison is within the family, not against the generation expectation: the generation gives a mean but no usable spread, and selection shrinks the spread below the unselected value (Frisch and Melchinger 2005).

Implementation: a second pass inside `sample_qc`, after the per-sample loop, so the flag order stays deterministic (append it last). Keep it a separate helper, e.g. `_family_donor_outliers(qc, dataset) -> set[str]`, so the tests can drive the statistic directly. Use plain `statistics.median` or numpy; be explicit about the even-length median (the ordinary mean of the two middle values).

## Fixture (`scripts/make_fixture.py`, `tests/fixtures/synthetic_bc2f1/`)

The generator computes expectations independently of the package; keep that. Add the flag to whatever structure holds expected QC flags, computed from the generator's own per-individual donor counts with its own median/MAD code written in the generator's style. Do not import anything from `progeny_selector.core.qc`.

Expected effect: in family `F2` (20 BC2F1 progeny) the planted selfed contaminant `BC2F2-F2-002` (check the exact id; the fixture's selfed contaminant) is flagged, because its B calls put it far above its family. Verify that the ids that gain the flag are exactly the ones the rule predicts from the generator's own numbers, and report them. Then regenerate the fixture, run the full suite, and confirm a second regeneration is byte-identical.

If the fixture's expected file does not carry QC flags at all, say so and leave it unchanged rather than adding a column.

## Tests (`tests/test_qc*.py`, follow the existing file layout)

- A family of 8 with donor fractions 0.03 repeated and one 0.14: only the high one is flagged.
- A family of 5 with the same shape: nobody is flagged (below the minimum size).
- A perfectly uniform family of 10 (all 0.02) plus one at 0.03: MAD is 0, so the floor applies and 0.03 is not flagged (0.03 < 0.02 + 2.5 * 0.01 = 0.045). One at 0.05 is flagged.
- One-sided: a family where one individual is far below the median flags nobody.
- Two families in one dataset are assessed separately; a value that is an outlier in one family is not flagged in the other where it is typical.
- Individuals with no called informative markers (NaN rates) neither count toward the minimum size nor get flagged.
- Heterozygous calls count half: an individual that is all H has donor_fraction 0.5.
- The flag does not exclude: build a result through `run_analysis` where one individual carries only `family_donor_outlier` with `exclude_qc_flagged: true`, and assert it still passes the filters and that `qc_excluded` is false in `qc_table_rows`.
- Flag order is deterministic and the flag is last on an individual that also carries another flag.

## Docs

- `CHANGELOG.md` [Unreleased] > Added: one line.
- `docs/adr/0012-family-donor-outlier.md`, MADR shape as in 0011 (Status accepted, Date 2026-09-16). Record:
  - Problem: the two real SoySNP50K NIL families each contain individuals whose donor fraction is several times their family's, and nothing flagged them: PI547454 (14.3 %), PI547634 (24.1 %), PI547646 (8.7 %), against 1-5 % for the rest. NIL families are genotyped with `generation` empty, so the generation expectation cannot be the baseline.
  - Options: a fixed donor-fraction threshold; expectation from the backcross generation; robust within-family outlier detection; making it a hard filter.
  - Decision: robust one-sided within-family rule as above, advisory only, because a high donor fraction can be a legitimate off-target introgression: Gilbert et al. 2023 found off-target introgressions from 2.6 kb to 54.8 Mb in this very collection, and 368 of 611 NILs at 80 % or more identity to the recurrent parent.
  - Effect on the real families, to state as evidence: Clark x PI86024 (n = 8, median 2.9 %, robust scale 1.2 %) flags PI547454 only, with PI547429 at 5.3 % (z 2.0) not flagged; Clark x Higan (n = 14, median 4.1 %, robust scale 1.5 %) flags PI547634 (z 13.1) and PI547646 (z 3.0), with PI547457 at 7.0 % (z 1.9) not flagged. Record that PI547646 is near the cut: at z 3.0 it would not be flagged, and only 22 individuals have been tested.
  - Why count and not weighted donor fraction (PI547457 is 7.0 % by count and 13.6 % weighted).
  - Consequences: the fixture's expected QC flags change; the flag appears in `qc_flags` in results.csv and in the Validate screen's flags column, an added value in an existing column, not a new column; no criteria key, so the thresholds change only by editing `qc.py`; isoline-browser does not read `qc_flags`.
  - Sources: Leys et al. 2013 https://www.sciencedirect.com/science/article/pii/S0022103113000668; plinkQC check_het_and_miss https://meyer-lab-cshl.github.io/plinkQC/reference/check_het_and_miss.html; Gilbert et al. 2023 https://acsess.onlinelibrary.wiley.com/doi/10.1002/tpg2.20310; Frisch and Melchinger 2005 https://academic.oup.com/genetics/article-abstract/170/2/909/6059341.
- Do not edit `docs/data-formats.md`: no documented column, parameter or filter changes. If you find that the flag vocabulary is documented there after all, stop and report instead of editing it.
