# QC flag `family_donor_outlier`: a robust one-sided within-family rule

Status: accepted. Date: 2026-09-16.

## Context and Problem Statement

The two real SoySNP50K NIL families each contain individuals whose donor fraction is several times their family's, and no QC flag marked them: PI547454 (14.3 %), PI547634 (24.1 %) and PI547646 (8.7 %), against 1-5 % for the rest. NIL families are genotyped with `generation` empty, so the generation expectation cannot be the baseline. The maintainer (2026-09-16) asked for a flag on an individual whose donor-genome fraction is far above the rest of its family. How should it be detected, and should it exclude?

## Considered Options

- A fixed donor-fraction threshold.
- The expectation from the backcross generation.
- Robust within-family outlier detection.
- Making the flag a hard filter.

## Decision Outcome

Chosen option: robust one-sided within-family outlier detection, advisory only. Per individual, genome-wide, `donor_fraction = hom_donor_rate + het_rate / 2`, whose denominator is the A + H + B calls at informative markers; X (non-parental) calls are excluded, and `possible_outcross` covers those individuals. Individuals are grouped by `family_id` as stored, with None and the empty string one group, so a dataset with no families is one group; progeny and candidates in a family are grouped together. Individuals with no called informative marker (NaN rates), or already flagged `high_missing`, neither count nor get flagged. In a group of at least `FAMILY_OUTLIER_MIN_N = 6`, with `med` the median and `mad` the median absolute deviation from it (an even count takes the mean of the two middle values), `scale = max(1.4826 * mad, FAMILY_OUTLIER_MIN_SCALE = 0.01)`, and the flag is set when `donor_fraction > med + FAMILY_OUTLIER_Z * scale` with `FAMILY_OUTLIER_Z = 2.5` (Leys et al. 2013). The floor keeps a very uniform family from flagging small deviations. A low donor fraction is never flagged here; `possible_rp_sample` covers that case. Semantics and constants: `core/qc.py`.

The flag is advisory, not in `QC_EXCLUDING_FLAGS`, because a high donor fraction can be a legitimate off-target introgression: Gilbert et al. 2023 found off-target introgressions from 2.6 kb to 54.8 Mb in this very collection, and 368 of 611 NILs at 80 % or more identity to the recurrent parent.

The comparison is within the family, not against the generation expectation: the generation gives a mean but no usable spread, and selection shrinks the spread below the unselected value (Frisch and Melchinger 2005).

Effect on the real families, as evidence: Clark x PI86024 (n = 8, median 2.9 %, robust scale 1.2 %) flags PI547454 only, with PI547429 at 5.3 % (z 2.0) not flagged; Clark x Higan (n = 14, median 4.1 %, robust scale 1.5 %) flags PI547634 (z 13.1) and PI547646 (z 3.0), with PI547457 at 7.0 % (z 1.9) not flagged. PI547646 is near the cut: at z 3.0 it would not be flagged, and only 22 individuals have been tested.

Why the count model and not the map-weighted donor fraction: the weighted measure inflates individuals whose donor calls are scattered singletons, which is exactly the noise `rule: run` was added to ignore (docs/adr/0011). PI547457 is 7.0 % by count and 13.6 % weighted.

### Consequences

- The fixture's expected QC flags change: `expected_results.csv` gains a `family_donor_outlier` column computed by the generator's own implementation.
- The flag appears in `qc_flags` in results.csv and in the Validate screen's flags column: an added value in an existing column, not a new column.
- No criteria.yaml key, so the thresholds change only by editing `core/qc.py`.
- isoline-browser does not read `qc_flags`, so it is unaffected.

## More Information

- Leys et al. 2013: https://www.sciencedirect.com/science/article/pii/S0022103113000668
- plinkQC check_het_and_miss: https://meyer-lab-cshl.github.io/plinkQC/reference/check_het_and_miss.html
- Gilbert et al. 2023, Plant Genome: https://acsess.onlinelibrary.wiley.com/doi/10.1002/tpg2.20310
- Frisch and Melchinger 2005: https://academic.oup.com/genetics/article-abstract/170/2/909/6059341
