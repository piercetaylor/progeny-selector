# A duplicate pair is an advisory flag on both members, not an exclusion

Status: accepted. Date: 2026-09-21.

## Context and Problem Statement

`core/qc.py` has carried `duplicate_pairs` since the initial scaffold, and both PLAN.md's algorithm 8 and the module docstring say duplicate pairs are reported, but `run_analysis` never called it: a plate re-run, a doubled tube or a mislabelled well reached the ranking unremarked, and the breeder learned of it only when two "different" selections turned out to be one plant. Three questions were open. Where does a detected pair surface — in results.csv, in the UI, or both? Does it exclude? And over which markers is identity by state measured, given that backcross siblings share most of their genome by descent?

## Considered Options

1. Report pairs in the UI only, leaving results.csv silent.
2. A new results.csv column carrying the duplicate partner.
3. An added value in the existing `qc_flags` column, on both members of every pair.

## Decision Outcome

Option 3, the value `possible_duplicate`, exactly as docs/adr/0012 added `family_donor_outlier`: an added value in an existing column is not a schema change (docs/adr/0016 — these columns are open vocabularies), while a new column is a minor bump and a UI-only report does not survive the download the breeder actually archives. `AnalysisResult.duplicates` keeps the pairs themselves, `(a, b, ibs)`, for the `validate` command and the Validate screen, which state the IBS a flag cannot carry.

The flag is advisory and is not in `QC_EXCLUDING_FLAGS`. A duplicate is a bookkeeping error the breeder resolves — one of the two rows is a tube that should not exist, and which one is the real plant is knowledge the tool does not have. Dropping both would discard a genuine selection candidate; dropping one arbitrarily would be worse. The other advisory flags set the precedent: the tool says what it sees and the ranking policy stays the breeder's.

The threshold is IBS >= 0.995, measured over the **informative markers called in both** individuals — the markers where the recurrent and donor parent are called, homozygous and carry different alleles. This is the part of the decision that matters. The uninformative markers are shared by construction, so where the informative fraction is small they dominate the average and backcross siblings reach 0.995 over all markers; restricting the denominator to the informative markers is what keeps the measure discriminating in a backcross. (This project's synthetic fixture is not such a panel — 475 of its 500 markers are informative, and its closest progeny pair sits at 0.931 either way — so the fixture cannot demonstrate the difference; a sparse real panel can.) `core/similarity.py`'s `pairwise_ibs` therefore takes a `marker_idx` pool, defaulting to every marker so `ibs_to_sample` and the parent comparisons are unchanged, and `run_analysis` passes `np.flatnonzero(classification.informative)`. At most 2,000 informative markers are used: above that a seed-0 subsample is drawn, which keeps the pairwise scan bounded and is reproducible across runs.

0.995 is deliberately stricter than PLINK's and KING's duplicate practice, which flags far lower. In a backcross programme BC5 sibs and selfed siblings legitimately approach a very high IBS, and a threshold tuned for an unrelated-sample cohort would flag whole families. A false negative here costs one unnoticed re-run; a false positive costs the breeder's trust in every other flag in the column.

The scan is skipped entirely above 2,000 individuals, with the warning `duplicate detection skipped above 2000 individuals`, because it is quadratic in the sample count.

### Consequences

Good: a doubled sample is visible in results.csv, in the Validate screen and in `progeny-selector validate` without changing the schema version; the measure is meaningful in a backcross rather than saturated; the default behaviour of every other IBS caller is untouched. Bad: `qc_flags` content changes for any dataset that does contain duplicates, so a reader matching that column exactly must be updated; with few informative markers the restricted denominator is small and the estimate noisy. Neutral: a progeny identical to a parent is already `possible_rp_sample` or `possible_donor_sample` and is not in the progeny-only pairwise set, so the two rules do not overlap; a pair where one member is `high_missing` is still flagged, since IBS uses only the markers called in both.

## Amendments, 2026-09-21

**A minimum overlap between the two members.** Maintainer decision, taken after the phase 6 review, which found the false positive: IBS accepted any denominator at or above one, so a sample called at almost none of its markers reached IBS 1.0 against everyone it happened to agree with on its handful of calls. With 30 markers, a sample called at 2 of them matched two plants that differ from each other at 16 markers, and all three were flagged.

A pair is now reported only when the two members are called in common at `min_overlap_frac` of the markers actually used for that comparison — the marker set after the `marker_idx` restriction and after subsampling, not the full panel — with the floor computed as `ceil(frac * n_markers_used)` and never below one. `min_overlap_frac = 0.5`, as a constant (`core/qc.py` `MIN_DUPLICATE_OVERLAP_FRAC`) and a keyword default on `duplicate_pairs`. `core/similarity.py` `pairwise_ibs` gains `return_counts`, returning the per-pair overlap matrix the loop already computes; its default single-array return is unchanged, so `ibs_to_sample` and the parent-similarity callers are untouched. `scripts/make_fixture.py` applies the same floor, so the generator's independent implementation and the pipeline agree on the rule; the fixture's expectation is still an empty duplicate set, since its missingness is far too low for the floor to bite.

A fraction rather than an absolute count, because the panels this tool reads run from a few hundred markers to 50K, and a count that is generous on a 50K array would reject every honest pair on a 500-marker KASP panel.

0.5 is a judgement about how much agreement is worth believing, not a figure derived from the hard filters, and the two are not comparable: `Filters.max_missing_rate` is a fraction of *all* markers, while the floor is a fraction of the *informative* pool the comparison actually uses. Where missingness is spread evenly across the panel the two do line up — an individual 20 % missing overall is about 20 % missing among the informative markers, two such individuals overlap at roughly two thirds of that pool, and the floor never bites. Where it is not spread evenly they come apart, and the trade-off is real: on 500 markers of which 100 are informative, a genuine duplicate whose partner is missing at 60 of those 100 informative markers is only 12 % missing overall, so it passes the hard filters and is ranked, and the floor still drops the pair from duplicate detection. Half the pool was judged the right place to stand: below it a match rests on so few markers that a moderately similar sibling can reach 0.995 by chance, and an advisory flag that fires on chance agreements is one breeders stop reading.

The consequence is stated rather than hidden: an individual whose missing calls concentrate in the informative markers can be ranked and still be invisible to duplicate detection. `validate` and the Validate screen report the pairs that were found, not the pairs that were skipped, so this is a false negative the tool does not announce.

Not a criteria key. Exposing it would change a documented parameter in `docs/data-formats.md` and pull the shared input contract into the question, which is not worth doing for a threshold nobody should be tuning; unlike the assembly or the ranking mode, there is no run in which a breeder wants a different value.

This supersedes nothing above, but it narrows the original Consequences: the Bad clause "with few informative markers the restricted denominator is small and the estimate noisy" described a sparse *panel*, where every sample has few informative markers and the floor scales down with them. It did not describe a sparse *sample* on a dense panel, where the denominator collapses for one pair while the panel stays large — the case the floor closes.
