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

The threshold is IBS >= 0.995, measured over the **informative markers called in both** individuals — the markers where the recurrent and donor parent are called, homozygous and carry different alleles. This is the part of the decision that matters. The uninformative markers are shared by construction, so where the informative fraction is small they dominate the average and backcross siblings reach 0.995 over all markers; restricting the denominator to the informative markers is what keeps the measure discriminating in a backcross. (This project's synthetic fixtures are not such a panel — 475 of the 500 markers are informative — so they cannot demonstrate the difference between the two denominators; a sparse real panel can. The BC2F1 fixture's closest progeny pair sits at 0.931 either way. What the BC3F1 fixture added later does demonstrate is the threshold itself: see "What the BC3F1 fixture shows" below.) `core/similarity.py`'s `pairwise_ibs` therefore takes a `marker_idx` pool, defaulting to every marker so `ibs_to_sample` and the parent comparisons are unchanged, and `run_analysis` passes `np.flatnonzero(classification.informative)`. At most 2,000 informative markers are used: above that a seed-0 subsample is drawn, which keeps the pairwise scan bounded and is reproducible across runs.

0.995 is deliberately stricter than PLINK's and KING's duplicate practice, which flags far lower. In a backcross programme BC5 sibs and selfed siblings legitimately approach a very high IBS, and a threshold tuned for an unrelated-sample cohort would flag whole families. A false negative here costs one unnoticed re-run; a false positive costs the breeder's trust in every other flag in the column.

The scan is skipped entirely above 2,000 individuals, with the warning `duplicate detection skipped above 2000 individuals`, because it is quadratic in the sample count.

### Consequences

Good: a doubled sample is visible in results.csv, in the Validate screen and in `progeny-selector validate` without changing the schema version; the measure is meaningful in a backcross rather than saturated; the default behaviour of every other IBS caller is untouched. Bad: `qc_flags` content changes for any dataset that does contain duplicates, so a reader matching that column exactly must be updated; with few informative markers the restricted denominator is small and the estimate noisy. Neutral: a progeny identical to a parent is already `possible_rp_sample` or `possible_donor_sample` and is not in the progeny-only pairwise set, so the two rules do not overlap; a pair where one member is `high_missing` is still flagged, since IBS uses only the markers called in both.

## Amendments, 2026-09-21

**A minimum overlap between the two members.** Maintainer decision, taken after the phase 6 review, which found the false positive: IBS accepted any denominator at or above one, so a sample called at almost none of its markers reached IBS 1.0 against everyone it happened to agree with on its handful of calls. With 30 markers, a sample called at 2 of them matched two plants that differ from each other at 16 markers, and all three were flagged.

A pair is now reported only when the two members are called in common at `min_overlap_frac` of the markers actually used for that comparison — the marker set after the `marker_idx` restriction and after subsampling, not the full panel — with the floor computed as `ceil(frac * n_markers_used)` and never below one. `min_overlap_frac = 0.5`, as a constant (`core/qc.py` `MIN_DUPLICATE_OVERLAP_FRAC`) and a keyword default on `duplicate_pairs`. `core/similarity.py` `pairwise_ibs` gains `return_counts`, returning the per-pair overlap matrix the loop already computes; its default single-array return is unchanged, so `ibs_to_sample` and the parent-similarity callers are untouched. `scripts/make_fixture.py` applies the same floor, so the generator's independent implementation and the pipeline agree on the rule; the BC2F1 fixture's expectation is still an empty duplicate set, since its missingness is far too low for the floor to bite. (Written before the BC3F1 fixture existed; that fixture does contain flagged pairs, and its overlaps are far above the floor too. See below.)

A fraction rather than an absolute count, because the panels this tool reads run from a few hundred markers to 50K, and a count that is generous on a 50K array would reject every honest pair on a 500-marker KASP panel.

0.5 is a judgement about how much agreement is worth believing, not a figure derived from the hard filters, and the two are not comparable: `Filters.max_missing_rate` is a fraction of *all* markers, while the floor is a fraction of the *informative* pool the comparison actually uses. Where missingness is spread evenly across the panel the two do line up — an individual 20 % missing overall is about 20 % missing among the informative markers, two such individuals overlap at roughly two thirds of that pool, and the floor never bites. Where it is not spread evenly they come apart, and the trade-off is real: on 500 markers of which 100 are informative, a genuine duplicate whose partner is missing at 60 of those 100 informative markers is only 12 % missing overall, so it passes the hard filters and is ranked, and the floor still drops the pair from duplicate detection. Half the pool was judged the right place to stand: below it a match rests on so few markers that a moderately similar sibling can reach 0.995 by chance, and an advisory flag that fires on chance agreements is one breeders stop reading.

The consequence is stated rather than hidden: an individual whose missing calls concentrate in the informative markers can be ranked and still be invisible to duplicate detection. `validate` and the Validate screen report the pairs that were found, not the pairs that were skipped, so this is a false negative the tool does not announce.

Not a criteria key. Exposing it would change a documented parameter in `docs/data-formats.md` and pull the shared input contract into the question, which is not worth doing for a threshold nobody should be tuning; unlike the assembly or the ranking mode, there is no run in which a breeder wants a different value.

This supersedes nothing above, but it narrows the original Consequences: the Bad clause "with few informative markers the restricted denominator is small and the estimate noisy" described a sparse *panel*, where every sample has few informative markers and the floor scales down with them. It did not describe a sparse *sample* on a dense panel, where the denominator collapses for one pair while the panel stays large — the case the floor closes.

## What the BC3F1 fixture shows, 2026-09-21

Recorded when the two-generation round trip (docs/adr/0018) added `tests/fixtures/synthetic_bc3f1/`. It corrects two statements above that were true when this ADR shipped and are no longer true of every fixture: the duplicate set is empty only in the BC2F1 fixture.

Eight of the forty BC3F1 individuals carry `possible_duplicate`, and **all eight are in the family of `BC2F1-F1-001`**. The pairs at or above the threshold are 003/010 at 0.99892 over 464 shared markers, 006/009 at 0.99786 over 467, 001/005 at 0.99784 over 463, 002/008 at 0.99679 over 467 and 002/009 at 0.99572 over 467. The overlap floor for this panel is 238 of the 475 informative markers and every overlap is above 460, so the floor plays no part in this result; the threshold alone produces it.

These flags are **not a normal BC3F1 result, and nothing here should be read as a claim that backcross sibs generally reach 0.995**. `BC2F1-F1-001` is the deliberately extreme individual planted in the BC2F1 fixture (`scripts/make_fixture.py`, the "best possible passer"): 28 heterozygous markers in four blocks, RPP 0.9705 and heterozygosity 0.059 — BC4/BC5-equivalent genome recovery wearing a BC2F1 label. Its progeny inherit half of what little donor genome it has, so they are nearly identical stretches of recurrent parent and differ at a handful of markers. The other three families descend from realistic parents (RPP 0.78, 0.90 and 0.85) and their maximum within-family IBS is 0.936, 0.982 and 0.963 respectively, with not one flag between them. Eighty per cent of one deliberately extreme family and zero per cent of the three normal ones is the honest summary.

That is the trade-off this ADR already recorded, now visible in a fixture rather than argued from first principles: at BC4/BC5-equivalent recovery, full sibs legitimately reach the duplicate threshold, and the flag is advisory precisely so that a breeder — not the tool — decides what to do about it. The fixture's `README.md` says the same in the place a reader of the data will look.
