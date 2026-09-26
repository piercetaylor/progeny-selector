# progeny-selector: plan

## Problem statement

In a marker-assisted backcross (MABC) program a breeder must pick, from hundreds to thousands of segregating progeny per generation across several families, the few plants that carry the donor allele at the target locus, do not carry donor alleles at loci or regions the program wants to avoid, have recovered the most recurrent-parent genome (especially on the non-carrier chromosomes), and carry the shortest donor segment around the target with a recombination close to it on each flank. Today the genotype file is recoded in a spreadsheet or an ad-hoc R script, percent recurrent parent is computed by counting, drag is judged from a heatmap, and the selection list is typed by hand; each generation the script is rewritten. progeny-selector reads the genotype file, a sample manifest and a criteria file, computes every status and metric with documented formulas, applies hard filters before a weighted composite score, ranks overall and within family, and exports the ranked table, the selected IDs and the sample manifest for the next genotyping round. The user navigates cross → family → generation → individual, compares candidates side by side, and never uploads unreleased genotype data to a third party.

## Users and workflow today vs. target

Users: the breeder who decides which plants advance, and the graduate student or technician who runs genotyping and prepares lists; both fluent in R, working on a laptop in a field office during the selection window, with a Snakemake/HPC environment available for batch work.

Today: KASP or array calls arrive as a spreadsheet or VCF; someone recodes calls as A/B/H against the parents, filters carriers by eye, sorts by percent recurrent parent, checks the avoid marker by hand, and writes tag numbers into a planting list. Recombinant selection on the carrier chromosome is rarely done because it needs positions and a window. The BC generation and family structure live in file names.

Target: load genotypes, samples.csv and criteria.yaml; the Validate screen reports contract errors and QC flags; the Navigate screen shows cross → family → generation; the Rank screen lists individuals with pass/fail chips for each target and avoid locus, RPP (total, carrier, non-carrier), drag bounds, recombinant flags, composite score and rank, with column filters and multi-select; the Compare screen shows selected individuals side by side with compact chromosome strips; the Selection screen holds the list with notes and a top-N-per-family helper and shows the expected next-generation composition; the Export screen writes results.csv, selected.csv and next_samples.csv. The CLI does the same in one command on an HPC node; the CSV outputs are read in R with `scripts/read_results.R`.

## Scope and non-goals

In scope: MABC selection among segregating progeny per generation; foreground selection at target loci defined by marker, region or flanking markers with a required donor state; background selection with count- and map-weighted RPP overall, per chromosome, carrier versus non-carrier; recombinant selection on carrier-chromosome flanks; negative selection at avoid loci and regions; identity-by-state similarity to each parent; QC (missing, heterozygosity, expected values by generation, selfs, outcrosses, swaps, duplicates); hard filters, composite score, ranking, selection lists, top N per family, projected next-generation expectation; navigation cross → family → generation → individual; CSV exports; CLI; Shiny UI runnable locally or as a static Shinylive site.

Non-goals: rich graphical-genotype browsing and QC of finished isolines (the sibling backcross; this tool shows only a compact per-individual chromosome strip); phenotype or trial analysis; marker design; imputation or error correction; multi-donor pyramiding (documented extension); mobile layouts; server-side data storage.

## Reference projects and what is reused

Details, licences and fetched URLs are in docs/reference-repos.md. The weighted RPP model, the linkage-drag definition and the output-table layout (per-chromosome RPP, total, QTL status, drag, selected state, rank, comment) follow Flapjack's MABC analysis [web] https://flapjack.hutton.ac.uk/en/latest/mabc.html and its tutorial [web] https://flapjack.hutton.ac.uk/en/latest/mabc_tutorial.html; Flapjack's QTL "Source" column (DP or RP) is generalised into `required_state` for targets and RP-required avoid loci. The compact strip in the Compare screen borrows the one-row-per-line colour-by-state idea from flapjack-bytes [web] https://github.com/cropgeeks/flapjack-bytes. The A/B/H coding vocabulary follows ABHgenotypeR [web] https://github.com/StefanReuscher/ABHgenotypeR/ and GenoSee [web] https://github.com/hashimotoshumpei/GenoSee. The IBS definition is SNPRelate's [web] https://rdrr.io/bioc/SNPRelate/man/snpgdsIBS.html. Package layout (src package, docs, tests, pytest) follows scikit-allel [web] https://github.com/cggh/scikit-allel and py-shiny [web] https://github.com/posit-dev/py-shiny; the logic/view split follows rhino [web] https://github.com/Appsilon/rhino. BrAPI's Genotyping module [web] https://github.com/plantbreeding/API is the intended later data source. No code was copied from any reference.

## Architecture

```mermaid
flowchart LR
  subgraph Inputs
    G[VCF / HapMap / wide CSV]
    S[samples.csv]
    M[markers.csv optional]
    C[criteria.yaml]
  end
  G --> IO[progeny_selector.io<br/>parsers + boundary validation]
  S --> IO
  M --> IO
  C --> IO
  IO --> DS[(Dataset<br/>numpy int8 calls)]
  IO --> CR[(Criteria dataclasses)]
  DS --> CORE[progeny_selector.core<br/>classify → foreground / background / drag / avoid / similarity / qc → score → selection]
  CR --> CORE
  CORE --> ROWS[(AnalysisResult.rows<br/>one dict per individual)]
  ROWS --> APP[progeny_selector.app<br/>Shiny screens 1–7]
  ROWS --> CLI[progeny_selector.cli<br/>validate / rank / select]
  APP --> EXP[progeny_selector.io.export<br/>results.csv, selected.csv, next_samples.csv]
  CLI --> EXP
  EXP --> R[R reader<br/>scripts/read_results.R]
  APP -. shiny run or shinylive export .-> HOST[Laptop / university server / GitHub Pages]
```

`core` is pure: arrays and dataclasses in, arrays and dataclasses out, no I/O, no pandas, no Shiny. `io` is the only place that validates. `app` renders rows and calls `run_analysis`; it computes nothing. The same package runs under CPython locally and under Pyodide in the browser.

## Data model and file contracts

`GenotypeMatrix` stores calls as an int8 array of shape (markers, samples, 2) of allele indices (−1 missing) into a per-marker allele table, so VCF GT, HapMap nucleotides and A/B/H codes share one structure. `Sample` carries sample_id, line_name, role, generation, family_id, notes; `Dataset` joins the two and asserts exactly one recurrent and one donor parent. `Criteria` holds `TargetSpec` and `AvoidSpec` loci (marker, region or flanking pair), flank windows, `BackgroundOptions`, `Weights` and `Filters`. `Classification` is an int8 (markers, samples) state matrix with six states (A, H, B, X, N, U) plus an informative mask. `AnalysisResult.rows` is a list of flat dicts whose column names are the results.csv contract. Every file (genotypes, samples.csv, markers.csv, criteria.yaml, results.csv, selected.csv, next_samples.csv) is specified with column tables and examples in docs/data-formats.md; the input contract is shared with backcross.

## Core algorithms and metrics

1. Parent-of-origin classification (`core/classify.py`). A marker is informative when both parents are called, both homozygous and their alleles differ; otherwise every call at that marker is U (uninformative) with a recorded reason (parents identical; recurrent or donor missing; recurrent or donor heterozygous). At an informative marker a progeny call is A (both alleles recurrent), B (both donor), H (one of each), X (any allele in neither parent), or N (missing). Coded A/B/H input defines allele 0 = recurrent, 1 = donor. This definition is identical to backcross's (labels rp_hom, donor_hom, het, nonparental, missing, uninformative).

2. Foreground status (`core/foreground.py`). A target is a marker, a region (all markers with chrom and start_bp ≤ pos ≤ end_bp) or a flanking pair. Per marker the predicate for `required_state` is hom_donor → B; het → H; either → H or B. Over the locus, counted markers are those with A/H/B calls; with fewer than `min_markers` counted the status is unknown; rule all requires every counted marker to satisfy the predicate, rule any at least one; flanking loci use rule all with both markers required. Result per individual: pass, fail or unknown.

3. Background (`core/background.py`). Contribution s = 1 for A, 0.5 for H, 0 for B; X, N and U leave both sums. Count model: RPP = Σs / n_called. Weighted model: RPP = Σ w s / Σ w with w = min(d_left/2, c/2) + min(d_right/2, c/2), d the distance to the neighbouring informative marker (bp or cM), c the maximum marker coverage (default 10 cM or 4 Mb), chromosome ends bounded by the Wm82.a4.v1 length when positions are in bp [web] https://www.soybase.org/resources/genome_info/; this is Flapjack's model [web] https://flapjack.hutton.ac.uk/en/latest/mabc.html. Reported overall, per chromosome, and over carrier chromosomes (those holding any target) versus non-carrier chromosomes. Expected RPP after n backcrosses without selection is 1 − (1/2)^(n+1): BC1 0.75, BC2 0.875, BC3 0.9375 [web] https://iastate.pressbooks.pub/molecularplantbreeding/chapter/marker-assisted-backcrossing/; selfing generations leave expected RPP unchanged and halve heterozygosity each generation (`core/generation.py`). Theory and context: Hospital and Charcosset 1997, Genetics 147:1469–1485, DOI 10.1093/genetics/147.3.1469 [web] https://academic.oup.com/genetics/article-abstract/147/3/1469/6054126 (foreground marker positioning, background selection, carrier-chromosome recovery); Frisch, Bohn and Melchinger 1999, Crop Science 39:1295–1301, DOI 10.2135/cropsci1999.3951295x [web] https://experts.illinois.edu/en/publications/comparison-of-selection-strategies-for-marker-assisted-backcrossi/ (four-stage selection prioritising carrier-chromosome recombinants cut marker data points by about 75 %); Frisch and Melchinger 2005, Genetics 170:909–917, DOI 10.1534/genetics.104.035451 [web] https://academic.oup.com/genetics/article-abstract/170/2/909/6059341 (selection theory accounting for variance reduction after selection).

4. Linkage drag and recombinant selection (`core/drag.py`). For each target and individual, walk outward from the locus boundary over informative markers on the carrier chromosome, skipping N and X. The max bound is the distance to the first A marker (or to the chromosome end: 0 on the left; the assembly length in bp, or the last mapped marker in cM, on the right), matching Flapjack's "first recombination on each side" [web] https://flapjack.hutton.ac.uk/en/latest/mabc.html; the min bound is the distance to the outermost H/B marker before that A marker (0 if none). Estimate = (min + max) / 2 per side; total = left + right. recombinant_left is true when left_max ≤ the left window (default `flank_window`, 5 cM, per-target override); likewise right. Distances are in cM when every marker has cM and `flank_unit` is cm, else bp with a warning; bp locus boundaries are converted to cM by linear interpolation on the chromosome's map.

5. Avoid criteria (`core/avoid.py`). Same locus resolution and rule machinery with the predicate A (or A/H when `allow_het`); pass, fail or unknown per individual. Policy for unknown is `filters.unknown_avoid_is` (default pass, because no data cannot prove the donor allele is present; the row keeps the unknown chip).

6. Similarity (`core/similarity.py`). IBS to each parent over all markers called in both: per marker, shared alleles / 2 (multiset intersection: 1, 0.5 or 0), averaged; equals SNPRelate's 1 − |g1 − g2| / 2 for biallelic markers [web] https://rdrr.io/bioc/SNPRelate/man/snpgdsIBS.html. Pairwise IBS among progeny on a random subset of at most 2,000 markers finds duplicates (≥ 0.995), skipped above 2,000 individuals.

7. Composite score, ranking, selection (`core/score.py`, `core/selection.py`). Hard filters first: every target must pass (unknown → `unknown_target_is`, default fail), every avoid locus must pass, missing rate ≤ `max_missing_rate` (0.2), and no excluding QC flag (possible_self_or_outcross, possible_outcross, possible_rp_sample, possible_donor_sample) when `exclude_qc_flagged`. Excluded individuals keep their metrics and a `;`-joined `exclusion_reason`. Components in [0, 1]: rpp_noncarrier (falls back to rpp_total when no non-carrier markers), rpp_carrier, drag = 1 − min(1, total_est / carrier chromosome span) averaged over targets, recombinant = flagged flanks / (2 × targets), similarity_rp, completeness = 1 − missing rate. Score = Σ w_k c_k / Σ w_k over non-NaN components; default weights 0.5, 0.2, 0.2, 0.1, 0, 0. Ties break by rpp_total desc, drag_est asc, missing_rate asc, sample_id asc; ranks are dense, overall and within family. `select_top_n(rows, n, per_family)` returns rank_in_family ≤ n (or rank_overall ≤ n). Projection from the selected group's mean A/H/B fractions: backcross → A' = A + H/2, H' = H/2 + B, B' = 0, so expected RPP' = (1 + RPP)/2 and half the offspring carry a heterozygous target; self → A' = A + H/4, H' = H/2, B' = B + H/4, RPP unchanged, 3/4 carry the target (1/4 homozygous). Assumptions: Mendelian segregation, unlinked loci, no selection within the next generation, no distortion, no genotyping error.

8. QC (`core/qc.py`). Per individual: missing rate over all markers; het, hom-donor and non-parental rates over called informative markers; expected het and RPP from the parsed generation. Flags: high_missing; generation_unparsed; possible_self_or_outcross when a BCnF1 shows B calls above 2 % (homozygous donor is impossible in a true backcross apart from error); het_rate_deviates when |observed − expected| > 0.15 (advisory; selected high-RPP plants legitimately deviate); possible_outcross when non-parental calls exceed 2 %; possible_rp_sample or possible_donor_sample when IBS to that parent > 0.995 with het < 0.5 %; family_donor_outlier (advisory, never excludes; docs/adr/0012) when the count-model donor fraction (hom-donor rate + het rate / 2) exceeds the family median + 2.5 × max(1.4826 × MAD, 0.01), in families of at least 6 individuals that are neither high_missing nor without called informative markers; possible_duplicate (advisory, never excludes; docs/adr/0017) when pairwise IBS over informative markers called in both individuals is at least 0.995, on both members of the pair. Parents: heterozygosity above 2 %, missing above 10 %, fewer than 20 % informative markers, or none, are warnings. Duplicates from pairwise IBS are reported as pairs.

Edge cases: a marker with a heterozygous or missing parent is uninformative rather than guessed; multiallelic markers work by allele index; an individual with no called informative markers has NaN RPP and score and is excluded by the missing-rate filter; a target inside a region with zero called markers is unknown; targets on different chromosomes make several carrier chromosomes; when the map lacks cM for any marker, all cM options fall back to bp with a warning.

## UI walkthrough

Screen 1 Load: file inputs for genotypes, samples.csv, markers.csv (optional) and criteria.yaml, a Load and analyse button, and a status panel quoting contract errors verbatim or the counts and warnings on success (implemented). Screen 2 Validate and QC: informative-marker count, parent warnings, per-individual QC table with flags and expected values, QC-excluded and flagged rows tinted (implemented on the real QC objects). Screen 3 Navigate: an accordion of families, each holding a radio group of its generations, with breadcrumbs whose ancestor crumbs step back up; the selection filters every later screen (implemented; resolution 4, docs/m1-phases.md). Screen 4 Rank: DataGrid of individuals with rank, family, pass/fail, exclusion reason, composite score, RPP total/carrier/non-carrier, drag, missing rate, QC flags and one status column per target and avoid locus; header filters, sortable columns, multi-row selection, a switch to hide excluded individuals (implemented on the real rows). Screen 5 Compare: one column per selected individual with status chips, per-chromosome RPP, drag bounds and recombinant flags, and a compact 20-row chromosome strip in Okabe-Ito colours, capped at the first six selected individuals (implemented). Screen 6 Selection list: the selected individuals with notes, "add top N per family", and the projected next-generation composition for backcross or self (implemented on the real rows). Screen 7 Export: results.csv, selected.csv and next_samples.csv downloads (implemented against the writers). Keyboard: Shiny's navbar and inputs are focusable in order; the Rank grid needs a mouse for multi-row selection, a limit of the third-party DataGrid component recorded as an exception in `docs/accessibility.md` and docs/adr/0023; status chips carry text labels so colour is never the only cue. Palette: Okabe-Ito, A blue #0072B2, H bluish green #009E73, B vermilion #D55E00, X reddish purple #CC79A7, N grey #999999, U yellow #F0E442; pass/fail/unknown bluish green/vermilion/grey [web] https://raw.githubusercontent.com/wch/r-source/trunk/src/library/grDevices/R/colorstuff.R.

## Technology decisions

Python 3.11+ package with numpy and PyYAML as the only core dependencies; Shiny for Python UI (pandas only for DataGrid); two deployments from one codebase, `shiny run` and a Shinylive static export (docs/adr/0001, 0003). Shared input contract with backcross (0002). MIT licence (0004). Criteria as strict YAML (0005). Six-state classification and Flapjack-style RPP (0006). Hard filters before a bounded weighted composite with fixed tie-breaks (0007). Tooling: ruff (lint and format), mypy, pytest with a reproducible fixture, GitHub Actions on Python 3.11 and 3.12, Conventional Commits, SemVer, Keep a Changelog, MADR.

## Repository layout

```
progeny-selector/
├── PLAN.md                         this document
├── README.md                       what it does, quickstart, status
├── CHANGELOG.md                    Keep a Changelog, 0.1.0 unreleased
├── CONTRIBUTING.md                 conventions
├── LICENSE                         MIT
├── pyproject.toml                  package metadata, ruff, mypy, pytest config
├── .env.example                    twelve-factor variables read only by config.py
├── .gitignore                      caches, site/, user genotype files
├── .github/workflows/ci.yml        ruff, mypy, pytest, fixture reproducibility, shinylive export, Pages deploy
├── docs/
│   ├── data-formats.md             input, criteria and output contracts
│   ├── reference-repos.md          fetched repositories and sources
│   └── adr/0001..0027-*.md         MADR records
├── scripts/make_fixture.py         deterministic synthetic BC2F1 fixture generator
├── src/progeny_selector/
│   ├── __init__.py, __main__.py, cli.py, config.py, constants.py
│   ├── model/    dataset.py (GenotypeMatrix, Sample, Dataset), criteria.py (specs, weights, filters)
│   ├── io/       calls.py, vcf.py, hapmap.py, wide_csv.py, manifest.py, criteria.py, export.py, __init__.py (loader)
│   ├── core/     chrom, generation, classify, foreground, avoid, background, drag, similarity, qc, score, selection, pipeline
│   └── app/      app.py (navbar shell, AppState), requirements.txt (shinylive), screens/ load, qc, navigate, rank, compare, select, export
└── tests/
    ├── conftest.py, test_classify.py, test_status_metrics.py, test_io.py, test_smoke_pipeline.py
    └── fixtures/synthetic_bc2f1/   genotypes.vcf, samples.csv, markers.csv, criteria.yaml, expected_results.csv, README.md
```

## Testing and CI

Unit tests cover classification on all six states and every uninformative reason, generation parsing and expectations, chromosome normalisation, foreground rules (marker, region all/any, flanking, unknown), avoid with and without allow_het, count and weighted RPP with hand-computed weights, drag bounds and recombinant flags in bp and cM including chromosome ends, IBS, composite renormalisation and tie-breaks, every parser (VCF plain and gzip, HapMap, wide nucleotide and coded CSV, manifest rules, marker map, criteria validation). The smoke test loads the fixture (2 parents, 40 BC2F1 progeny in 2 families, 500 markers, 475 informative, planted target, avoid locus, self contaminant, high-missing individual) and compares statuses, exclusion reasons, RPP (total, carrier, non-carrier), drag bounds, recombinant flags, composite scores, advisory flags and ranks against expected_results.csv produced by the generator's independent implementation, then exercises top-N selection, projection, the results writer and the CLI end to end. CI runs ruff check, ruff format --check, mypy, pytest with coverage, and regenerates the fixture and fails on any diff, in the `check` job on Python 3.11 and 3.12; the `e2e` job installs Chromium and runs the Playwright screen tests under `pytest -m e2e`; the `export` job runs `scripts/build_shinylive.py` and the export smoke test, then uploads the Pages artifact; `deploy-pages` deploys it to GitHub Pages on pushes to main when `ENABLE_PAGES` is set; the `r-reader` job installs R and readr and reads the fixture's results.csv, plus an empty results.csv, with `scripts/read_results.R`, not a per-commit gate. Playwright tests of the screens arrived in M1 (docs/adr/0008); the shared contract test suite with backcross (`tests/test_contract_cases.py`, `scripts/check_contract.py`) shipped in M2 and is a gate in both repositories.

## Deployment and cost

Local: `pip install progeny-selector` (or the repository in editable mode) and `shiny run`, or the CLI inside a Snakemake rule on the HPC cluster; no server, no cost. Shared: the CI-built Shinylive site on GitHub Pages, static files only, zero cost, genotype data processed in the browser tab; the exported site measures 44.5 MB, built 2026-09-26 from `829f247` (docs/limits.md). A university VM running `shiny run` behind the campus network is the alternative for large datasets shared by several users; a Dockerfile is a later addition. Posit Connect Cloud and shinyapps.io are not used because unreleased data would leave the university.

## Milestones

M0 scaffold (complete): package layout, all parsers, criteria reader, six-state classification, foreground/avoid/background/drag/similarity/QC, hard filters, composite, ranking, selection, projection, writers, CLI, Shiny shell with working Load, Rank, Selection and Export screens, fixture and tests, CI, documents. Acceptance: ruff, mypy and pytest pass; the CLI reproduces expected_results.csv on the fixture.

M1 vertical slice: Navigate tree with breadcrumbs driving the Rank filter; Validate screen table; status chips and per-locus columns in the Rank grid; Compare screen with side-by-side statuses and compact chromosome strips; criteria editable in the UI and downloadable; Shinylive export verified in CI. Acceptance: a breeder loads a real BC2F1 KASP or SoySNP6K dataset and produces selected.csv without leaving the browser; the export runs on GitHub Pages.

M2 usable (complete): weighted background model default with map warnings surfaced in the UI; per-target windows; notes persisted in selected.csv; next_samples.csv used as input for the next generation without editing; duplicate detection surfaced; keyboard walkthrough documented; CSV columns frozen under `results_schema`. Acceptance: two consecutive generations processed with the round-tripped manifest; `scripts/read_results.R` reads results.csv unchanged in CI.

M3 hardened (complete): streaming VCF parse to bound memory; measured performance and a documented limits table for CPython and Shinylive; accessibility review; the screen-test gaps M1 deferred; BrAPI loader. Acceptance, measured in docs/limits.md: the 50K x 2,000 parse peak is 241.0 MiB, below the matrix plus 64 MiB (190.9 + 64 = 254.9 MiB); the 50K x 2,000 analysis peak is 481.0 MiB, 2.52 times the 190.9 MiB matrix, within the 4x threshold.

### M3 scope, decided 2026-09-22

The milestone was stated two ways. This line named Playwright screen tests and shared contract tests as M3 work, and both had already shipped in M2: `tests/e2e/` holds six screen tests and `tests/test_contract_cases.py` with `scripts/check_contract.py` are gates. Both are struck above. The `docs/m2-phases.md` decisions of 2026-09-16 pushed the crop work past M2, which read as a second and conflicting M3.

The maintainer delegated the scope to research under the standing criterion, what is best for academic and open-source plant-breeding users, meaning reproducible with standard tools, tolerant of spreadsheet exports, and never silently wrong. The outcome follows.

The crop work is not this milestone. What M2 deferred is the crop chromosome schemes, and not the token profiles, which landed in M2 as contract 1.4.0. The schemes are contract 1.5.0, specified in `../backcross/docs/m4-phases.md` section 8 with this repository's half as S6. A contract version originates in backcross and is copied here, so S6 cannot precede backcross phase 6 and is not a design task for this repository. S6 lands before the M3 phases, because it threads a crop scheme through `io/vcf.py`, which the streaming parse then rewrites.

Multi-donor pyramiding is deferred past M3. `contract/data-contract.md` line 65 already reserves `donor_parent_2` and a per-marker donor assignment, and states that they will be added and will not change the meaning of existing columns. Adding a second donor is a contract version starting in backcross, and it changes the six-state classifier, per-target `required_state`, the meaning of `frac_b` and `ibs_donor`, and the results columns. That is a feature carrying a contract cost and not hardening. The forward commitment at line 65 means nothing shipped now is invalidated by waiting. Published marker-assisted backcrossing runs the two-gene case as parallel single-donor series intercrossed at the end, as in the Xa21 and Xa33 pyramiding into DRR17B (Balachiranjeevi et al., 2018, https://pmc.ncbi.nlm.nih.gov/articles/PMC6201878/), where every backcross generation is a one-donor problem this tool already ranks. The maintainer answered the one fact this turned on, on 2026-09-22: the two-gene step is an intercross of finished single-donor lines, so no backcross generation genotyped this season has two donors segregating in the same family. Every generation the program genotypes is a one-donor problem the tool already ranks, and the deferral is settled and not provisional. What the answer leaves open is the intercross itself. Selecting double homozygotes among its progeny is foreground selection at two target loci with no recurrent-parent-proportion question, and whether this tool should rank that population at all is a separate question from `donor_parent_2`, which is not raised here.

The performance figure this line previously carried, 2,000 individuals under 60 s in CPython, was marked speculation in this document and has been replaced by a measurement phase. At 50K markers and 2,000 individuals the working set is 100 million calls, and in a browser the binding constraint is memory and not throughput. Duplicate detection is skipped above 2,000 individuals under `docs/adr/0017`, so the quadratic scan still runs at exactly the target and the measurement must cover that path.

## Risks and open questions

Marker density on the carrier chromosome limits recombinant selection: with KASP panels of a few markers per chromosome the flank windows must be wider than the marker spacing, and the min/max drag bounds will be far apart; the UI must show both bounds. Missing data at the target marker is the single most consequential input problem; policy is explicit but users may not read it. Weighted versus count RPP disagree on uneven maps; results.csv names the model. Pyodide performance and memory for large VCFs are measured for CPython and Shinylive at three dataset sizes (docs/limits.md). Open questions for the user: which platform the program genotypes BC progeny on (KASP panel size, SoySNP6K, or 50K) and on which assembly positions are reported (BARCSoySNP6K positions are on Wm82.a2.v1 [web] https://digitalcommons.unl.edu/cgi/viewcontent.cgi?article=2402&context=agronomyfacpub; SoySNP50K on Glyma1.01 [web] https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0054985); whether a genetic map with cM is available for the panel or bp windows suffice; whether the composite weights should default to Frisch et al.'s four-stage lexicographic order instead of a weighted mean; how the intercross step of a two-gene pyramid should be handled, given that the backcross generations themselves carry one donor each (answered in part on 2026-09-22; see the M3 scope note); whether the GitHub Pages site may be public.

## Assumptions

Parents are inbred and nearly homozygous; both parents are genotyped on the same panel as the progeny, or the input is A/B/H-coded; one donor per analysis set; generation labels follow BCnFk; positions in the genotype file and markers.csv are on one assembly; progeny of one generation are analysed together per run; users have Python 3.11+ locally or a current browser for the static site; the sibling projects keep the samples.csv and wide-CSV contracts unchanged.

## Verification status

Run on 2026-09-04 in the scaffold container (Python 3.11.15, ruff 0.15.11, pytest 9.1.1, mypy 2.3.1, numpy 2.4.4, shiny 1.7.0). Commands executed from the repository root with `PYTHONPATH=src` (the package was not pip-installed).

```
$ PYTHONPATH=src python3 -m pytest
........................                                                 [100%]
24 passed in 0.49s

$ ruff check . && ruff format --check .
All checks passed!
46 files already formatted

$ PYTHONPATH=src python3 -m mypy
Success: no issues found in 39 source files

$ python3 scripts/make_fixture.py && diff -r <previous copy> tests/fixtures/synthetic_bc2f1
wrote fixture to .../tests/fixtures/synthetic_bc2f1: 500 markers, 40 progeny, 11 pass; informative 475
fixture unchanged

$ PYTHONPATH=src python3 -m progeny_selector validate --genotypes tests/fixtures/synthetic_bc2f1/genotypes.vcf \
    --samples tests/fixtures/synthetic_bc2f1/samples.csv --markers tests/fixtures/synthetic_bc2f1/markers.csv \
    --criteria tests/fixtures/synthetic_bc2f1/criteria.yaml
genotypes: 500 markers x 42 samples; coded=False; cM map=yes
samples.csv: 42 rows; progeny/candidates=40
informative markers: 475
QC-flagged individuals: 7

$ PYTHONPATH=src python3 -m shiny run src/progeny_selector/app/app.py --port 8123   (20 s) ; curl -s -o /dev/null -w "%{http_code}" http://127.0.0.1:8123/
200   (page title "progeny-selector")
```

Not verified: `shinylive export` (pip install shinylive failed in this container while building the lzstring wheel; the CI step is written but has not run); the Shiny screens beyond serving the shell (no browser automation here); the GitHub Actions workflow itself (no repository yet); behaviour on real SoySNP50K-scale data (fixture is 500 markers). File count and size are recorded in the final report.

### M1 verification, 2026-09-14, from the maintainer's laptop

Windows, Python 3.12.4, ruff 0.16.7, pytest 9.1.1, mypy 2.3.1, shiny 1.7.0, shinylive 0.8.11 (shinylive web assets 0.10.14), pytest-playwright 0.9.0, playwright 1.62.0, Chromium 151.0.7922.34. Commands run from the repository root with `.venv\Scripts\python.exe` (editable install, `pip install -e ".[dev,export,e2e]"`).

```
> ruff check . ; ruff format --check .
All checks passed!
84 files already formatted

> mypy
Success: no issues found in 42 source files

> pytest
95 passed, 17 deselected in 5.50s

> pytest -m e2e
16 passed, 1 skipped, 95 deselected in 46.99s
(the 1 skip is tests/e2e/test_shinylive_export.py without PS_SITE_DIR; Chromium 151.0.7922.34)

> python scripts/build_shinylive.py
exported .../site (44.4 MB) in 0.9s
(warm shinylive asset cache; the first run on this laptop, which also downloaded
the ~400 MB shinylive web-assets archive into that cache, took 47.8 s)

> $env:PS_SITE_DIR="site"; pytest -m e2e tests/e2e/test_shinylive_export.py -s
time to #load-run: 8.7s
time to #load-status: 0.3s
1 passed in 15.10s

> python scripts/make_fixture.py ; git diff --exit-code -- tests/fixtures
wrote fixture to .../tests/fixtures/synthetic_bc2f1: 500 markers, 40 progeny, 11 pass; informative 475
(no diff)
```

Packaging path in use: resolution 5's primary path (docs/m1-phases.md; docs/adr/0009). `scripts/build_shinylive.py` stages a copy of `src/progeny_selector` next to a stub `app.py` and exports that; the browser test recorded every request the exported page made while loading Pyodide, uploading the fixture and running the pipeline, and all of them started with the page's own origin (`http://127.0.0.1:8008/`) — no request to PyPI, npm or any CDN. The pipeline output (`500 markers, 40 progeny; 11 pass hard filters`) requires numpy, which only runs if the staged `progeny_selector` package imported inside Pyodide. The wheel-URL fallback was not needed and is not in use.

`python -m shinylive` does not work on this install (shinylive 0.8.11 ships no `__main__.py`; `python -m shinylive` raises "No module named shinylive.__main__"). `scripts/build_shinylive.py` detects that and falls back to the `shinylive` console script installed next to the interpreter, which is the documented form.

The exported app renders inside a same-origin `iframe` (path `/app_<random id>/`), discovered while writing the smoke test; it was not anticipated by resolution 6's description of the Playwright controllers, which assume the app is the top-level document. `tests/e2e/test_shinylive_export.py` scopes every locator through `page.frame_locator("iframe")` rather than reusing `shiny.playwright.controller`, and waits for each file upload's progress bar to reach `width: 100%` in its `style` attribute (not just to attach) before clicking Run, because Pyodide's upload round trip is slow enough to race a plain click otherwise.

GitHub Pages: enabled 2026-09-14 with the GitHub Actions source and `ENABLE_PAGES=true`. First deploy on 278536a (run 34888105283): check on 3.11 and 3.12, e2e, export and deploy-pages all passed. URL: https://piercetaylor.github.io/progeny-selector/ (returned HTTP 200).

Real-data protocol (Q7): run 2026-09-14 to 2026-09-16 on SoySNP50K Clark isolines through the CLI; results are in the Handoff block. The browser run on real segregating data moved to M2.

GitHub Actions on d1adc20 (run 34884283693, 2026-09-14): check on 3.11 and 3.12, e2e, and export including the smoke test all passed; deploy-pages skipped because ENABLE_PAGES is unset.

### M2 verification, 2026-09-21

Real segregating data, not the fixture: the Clark x PI86024 SoySNP50K NILs in the gitignored `data/nils/`, recurrent parent
PI548533, 41,358 markers x 10 samples, 8 NILs. Target dt1 (Glyma.19g194300, Gm19 Wm82.a2.v1) +/- 1 Mb, `hom_donor` under
`rule: run`, weighted background model, `exclude_qc_flagged: false` — the program's own `criteria_dt1_run.yaml`, with
`assembly: Wm82.a2` added in a scratch copy because that file predates phase 3. Nothing under `data/` was written or committed.

The CLI, the local `shiny run` and the Pages site agree exactly: 41,358 markers, 8 progeny, 4 pass hard filters.

Browser run, local `shiny run` on port 8009, Chromium 151 driven by Playwright: uploads complete 4.5 s, status 8.5 s. Validate
reports background model weighted, RPP unit bp, drag unit bp, assembly Wm82.a2, rank mode weighted, duplicate pairs none, and
14,538 informative markers against 26,820 uninformative because the parents are identical there. Two individuals were ranked,
selected and exported with a note typed on the Selection screen; selected.csv carried `results_schema` 1.0.0 and the note on the
right `sample_id` (PI547472, rank 1, rpp_total 0.946929). All 48 requests went to `localhost:8009`; none left the origin.

Browser run, https://piercetaylor.github.io/progeny-selector/ , same files: 20.7 s to the Run button, 21.8 s with the three
uploads complete, 29.1 s to the status line. All 81 requests went to `piercetaylor.github.io`; none left the origin, so the
Pyodide runtime and the package both came from the page. Shinylive export 44.4 MB, unchanged from M1.

Phase 3 on real data, which is the point of the assembly selector: with the criteria as the program wrote them in September,
and so with the `Wm82.a4` default, the run warns that marker positions on Gm04, Gm06, Gm08, Gm13 and Gm20 reach past the
recorded a4 lengths — Gm04 to 52,367,723 against 51,203,390. With `assembly: Wm82.a2`, which is what the data is, there are no
warnings. The dataset's own README had recorded that mismatch as an accepted approximation; it is now detected and fixable.

Phase 6 on real data: no pair of the 8 NILs reaches the 0.995 duplicate threshold. The closest is PI547429 x PI547472 at 0.99100
over 2,000 informative markers subsampled from 14,538, then 0.98300 and 0.97000. These are sibling isolines from one backcross
family, which is the population the threshold was most likely to fire on spuriously, and it does not. This run is also the only
exercise the marker-subsampling path has had on real data: both fixtures sit under the 2,000-marker bound, so neither reaches it.

The two-generation real-data round trip (Question B) was not run: no linked consecutive-generation dataset has been found, and
the maintainer is still sourcing one. The synthetic round trip covers the writer and the reader (phase 8, `tests/test_round_trip.py`).


### M3 verification, 2026-09-26

Versions: Python 3.12.4, numpy 2.5.3, shiny 1.7.0, shinylive 0.8.11, Playwright 1.62.0, Chromium
151.0.7922.34, axe-core 4.12.1 via axe-playwright-python 0.1.8 (docs/limits.md, docs/accessibility.md).

Limits (docs/limits.md): the 50K x 2,000 parse peak is 241.0 MiB, below the matrix plus 64 MiB
(190.9 + 64 = 254.9 MiB, docs/adr/0022). The 50K x 2,000 analysis peak is 481.0 MiB, 2.52 times the
190.9 MiB matrix, within the 4x threshold (docs/adr/0025). Wall-clock and memory are recorded for
6K x 2,000, 50K x 200 and 50K x 2,000 in CPython, and for all three in Shinylive; the browser
completed the 50K x 2,000 case as well as the two the acceptance line named.

Accessibility: axe-core checked nine page states against WCAG 2.2 AA tags. `html-has-lang` is fixed
(`lang="en"`); text contrast on parent-of-origin state A is fixed (black on `#0072B2` measured
4.05:1, below AA's 4.5:1; white measures 5.19:1). Three exceptions are documented in
docs/accessibility.md and docs/adr/0023: `aria-allowed-attr` and `label` on the DataGrid, and
`scrollable-region-focusable` on bslib card bodies, each a limit of a third-party component this
app does not author. A Tab walk from the Load screen reached the navbar tabs and `#load-run`, and
never focused an element with no layout box.

BrAPI: `tests/test_brapi.py` has 28 tests, all passing, with no network. `test_brapi_equals_vcf_on_the_fixture`
loads the fixture through BrAPI and through the VCF and gets the same `Dataset`.
`test_positions_are_one_based` confirms `pos_bp = start + 1`; `test_region_resolves_the_shifted_marker`
proves a foreground region resolves the shifted marker and not its neighbour. `test_cli_brapi` runs the
CLI end to end against a stub server. The bearer token appears only in the `Authorization` header and
is dropped on a cross-host redirect (`test_token_only_in_the_header_and_never_in_errors`).
`brapi-callsets` writes the call-set table.

No request left the origin: `PS_SITE_DIR=site pytest -m e2e tests/e2e/test_shinylive_export.py`
passed on 2026-09-26 against a site built from `829f247`; the test asserts that every request the
exported page makes stays on the page's own origin.

Gates at `829f247`: `pytest` 617 passed, 44 deselected; `ruff check .`, `ruff format --check .` and
`mypy` clean; `scripts/check_contract.py ../backcross` reports "contract 1.7.0: manifest ok, 273
files" and "contract mirror: 274 files identical"; the fixture regenerates byte-identical
(`git diff --exit-code -- tests/fixtures` clean); `pytest -m e2e` 43 passed, 1 skipped (the export
test above, run separately). CI run 36265695866 green on all six jobs against `1f3247c`: `check (3.11)`, `check (3.12)`, `e2e`, `export`, `r-reader` and `deploy-pages`.

M3 acceptance checklist (docs/m3-phases.md, "M3 acceptance checklist"):

- [x] An empty or blank next-generation label is refused in `io.export` for the CLI and the Export
  screen; a padded label is trimmed.
- [x] `docs/limits.md` records measured wall-clock and peak memory for 6K x 2,000, 50K x 200 and
  50K x 2,000 in CPython, and for 6K x 2,000 and 50K x 200 in Shinylive, with the 2,000-progeny
  rows marked "ran" and the 50K x 2,000 browser attempt recorded. The browser attempt completed.
- [x] `read_vcf` parses within 1.5x the matrix on the bound test and below matrix + 64 MiB at
  50K x 2,000; every VCF contract case, `test_io.py`, `test_crops.py` and the smoke test pass
  unchanged; the fixture regenerates without a diff.
- [x] `classify`, `rpp` and `ibs_to_sample` are bit-identical to their pre-phase implementations on
  random inputs and stay under their memory thresholds; the 50K x 2,000 analysis peak is at or
  below 4x the matrix.
- [x] Navigate offers "(no generation)"; the Rank caption names "(no family)" and "(no
  generation)"; the reload guard has a browser test. That test pins only the final crumb state
  (`tests/e2e/test_navigate_gaps.py:97-102`); `all([])` is vacuously true if the mutation observer
  never fires, so the test would also pass with the guard removed. Commit `59b1170` recorded this
  as an accepted limit, carried in the Handoff below.
- [x] axe-core finds no WCAG 2.2 AA violation outside `KNOWN_A11Y_EXCEPTIONS` on every screen
  state; every exception is justified in docs/accessibility.md and docs/adr/0023; a Tab walk from the Load
  screen never focuses an element with no layout box; statuses carry text.
- [x] A BrAPI variant set loads through the CLI into the same `Dataset` as the VCF of the same
  data; positions are `start + 1`, proven by `test_positions_are_one_based` and
  `test_region_resolves_the_shifted_marker`; the token never appears outside the `Authorization`
  header; `brapi-callsets.csv` is written; no test touches the network.
- [ ] `contract/` untouched: false as written. No M3 phase touched `contract/`; it moved to 1.6.0
  and 1.7.0 through sibling-led changes in backcross (`9f0cbd6`, `cbd8fe8`), and
  `scripts/check_contract.py ../backcross` now reports 274 files identical, not the count at
  `fc21019`. (The spec's `../isoline-browser` names the same directory under its old name.)
  `constants.RESULTS_SCHEMA` is confirmed `"1.1.0"`, unchanged, so that half of the item holds.
- [x] ADRs 0021-0025 are present, from this milestone's phases; 0026 and 0027 also exist, the
  contract 1.6.0 and 1.7.0 sibling halves.
  CHANGELOG entries under `[Unreleased]` were written by the main session; this PLAN.md
  verification block and the Handoff are updated; CLAUDE.md State is updated; the end-of-milestone
  reviewer pass is done.

## Handoff (updated after every phase)

- Milestone / phase: **M3 is complete**, 2026-09-26, at `829f247`. All nine phases of `docs/m3-phases.md` are done: phase 1 `next_generation` validation (`08e5021`), phase 2 the measurement harness and `docs/limits.md` baseline (`c3e54ab`), phase 3 streaming `read_vcf` (`25852a2`, docs/adr/0022), phase 4 the screen-test gaps (`59b1170`), phase 5 the accessibility review (`3d529cd`, docs/adr/0023), phase 6 the analysis-memory bound (`b2bd711`, docs/adr/0025), phase 7 the BrAPI loader (`4d93907`, docs/adr/0024), and phase 9's re-measurement and documents. Phase 8 (BrAPI on the Load screen under `shiny run`) is struck from M3 (Q6). The end-of-milestone review landed two fix commits, `937ba0b` (BrAPI) and `567f110` (CHANGELOG shapes, accessibility test). A Load-screen upload fix followed at `74c9f99`, and benchmark fixes at `829f247`. See PLAN.md, "M3 verification, 2026-09-26", for the measured acceptance, the gate counts and the checklist. M4 is not yet planned.

  Open items carried into M4, none blocking this milestone's acceptance:

  - Q14: whether `possible_duplicate` should mean something narrower than raw IBS, and the silent skip of pairs with insufficient marker overlap, are both still open (raised by the M2 phase 8 reviewer, not resolved in M3's scope).
  - Question B: no linked two-generation real dataset has been found; the round trip that needs one has not run on real data (the synthetic round trip shipped in M2).
  - Shinylive BrAPI waits on a server that sends `Access-Control-Allow-Origin` for the Pages origin; the browser build has no BrAPI source until then (Q6). The Load-screen source card (the struck phase 8) is carried forward for that milestone.
  - results.csv and selected.csv carry no BrAPI provenance columns in M3; `Dataset.call_sets` and `brapi-callsets.csv` keep the server ids so a `results_schema` bump can add `call_set_db_id` and `sample_db_id` later without re-fetching (Q7).
  - Three pre-existing defects in `read_vcf`, unchanged by phase 3 and found by its reviewer, deliberately left alone: a GT allele index of 128 or more raises a raw `OverflowError` from numpy before the ALT-count check can report it; `int()` leniency accepts `-5/0` as -5, `1_0/0` as 10 and `" 1/0"` as 1; a Latin-1 VCF escapes as a raw `UnicodeDecodeError`. Each is a contract-error candidate for a later milestone.
  - `rpp_contribution` in `core/classify.py` is unused in `src/` and holds a second copy of the 1/0.5/0 table that `constants.RPP_LOOKUP` now carries; docs/adr/0025 records its removal as a follow-up.
  - A review low not fixed: `tests/e2e/test_navigate_gaps.py:97-102` pins only the final crumb state after a reload. `all([])` is vacuously true when the mutation observer never fires, so the reload-guard test would still pass with the guard removed. Commit `59b1170` recorded this as an accepted limit rather than adding a wait or a retry.
  - Benchmark caveats from `docs/limits.md`: the `parse peak MiB` column brackets `load_dataset`, so it includes markers.csv as well as the VCF, and is not the parser's own peak; the duplicate scan dominates analysis memory and time at 2,000 individuals (175.6 of 201.7 MiB and 204.50 of 305.95 s at 6K x 2,000); the 2026-09-25 and 2026-09-26 CPython re-run timings were taken while a WSL virtual machine that grew from 1.9 to 4.5 GB was left running on the same machine, so their wall-clock direction is more trustworthy than their size.

- Milestone / phase (previous): **M2 is complete** at 6b1d3a4, 2026-09-21; CI run 35666253372 green on all six jobs; the 15-box acceptance checklist in `docs/m2-phases.md` is ticked with its evidence; the end-of-milestone reviewer pass over `65121bf..HEAD` found no high or medium finding and its six lows are closed bar one recorded below. M3 is next: non-soybean crops and the crop profiles, which is a contract version after 1.4.0. Phases as shipped: 1 (0bb7712), 2 (b8382f2), 3 (19c2770, adr 0015), 4 (419a4c5, adr 0016), 5 (e97d566 with 2b966bc), 6 (66f968c, adr 0017), 7 (ed4c4ea), 8 (e73d7cc, adr 0018), 9 (19a9cfc), 10 (b5136b2, adr 0019), 11 (the verification block above, f991b1c, dd9840e, 6b1d3a4). Wm82.a5 and a6 lengths remain deliberately absent from ASSEMBLIES; `docs/m2-phases.md` carries a dated note because the decisions record says otherwise.
- Next:
  1. Done 2026-09-14: end-of-milestone reviewer pass over `git diff 1ef748e..65121bf`. 11 findings, none high, listed under "Deferred reviewer findings"; they form the first M2 phase.
  2. Done 2026-09-14 to 2026-09-16: real-data runs on SoySNP50K Clark isolines (recurrent PI548533), converted by `scripts/soysnp50k_nils.py` into the gitignored `data/nils/`; never commit anything from `data/`. Targets are Wm82.a2.v1 gene models ±1 Mb, `hom_donor`, `exclude_qc_flagged: false`, generation empty. Against GRIN gene lists under `rule: run`: Clark x PI86024 24/24, Clark x Higan 52/56, 76/80 in total against 67/80 under `rule: any`. The four remaining disagreements are PI547634 at T, R and pa1 (an off-type, 24 % donor genome-wide) and PI547592 at R (scattered single donor calls, probably a donor seed-lot haplotype difference). `family_donor_outlier` flags exactly PI547454, PI547634 and PI547646 across both families; `possible_rp_sample` fires on PI547435 and PI547592 as expected. Open: whether Purdy dose n = 5 means BC4 or BC5; whether chromosome normalisation should accept SoyBase's `glyma.Wm82.gnmN.` prefix (a contract question that starts in backcross; the conversion script strips it for now).
  3. M2 planning: the planner writes `docs/m2-phases.md`; the main session puts its flagged questions to the maintainer before dispatch. Decided inputs, 2026-09-16:
     - Scope: PLAN.md "Milestones" M2, plus the M1 gap (a breeder takes real segregating data to selected.csv in the browser).
     - results.csv freeze: missing cells written `NA` (as backcross does); add `background_model`, `rpp_unit`, `rank_mode` and `results_schema` (SemVer); route `het_rate`, `expected_het` and `expected_rpp` through `_num` (a literal `nan` can be written today); write the header row when there are no results.
     - Consumer (question A, answered): no external R dashboard exists; the "R Shiny dashboards" lines were a scaffold assumption and now name `scripts/read_results.R` (readr, explicit column types), read against the fixture's results.csv by a CI job with R installed, plus a Python header/schema test.
     - Ranking: `ranking: {mode: weighted|staged}`, weighted default, staged order as docs/adr/0007's amendment states.
     - Soybean inputs: a SoySNP50K/6K marker position table across assemblies joined on marker id from the SoyBase Data Store GFF3s (`Wm82.gnm{1,2,4,5,6}.mrk.*`, license Open); chromosome-length tables per assembly with an assembly selector (`core/chrom.py` hardcodes Wm82.a4); a converter from Song et al. 2016 BMC Genomics 17:33 Table S1 (CC BY 4.0) to markers.csv with monotone clamped Marey-map cM interpolation for unmapped markers; a KASP export to wide CSV converter. No cM map exists for the program otherwise.
     - Two-generation data (question B, open): no public linked consecutive-generation dataset found; Kim et al. 2021 Plants 10:804 Table S3 (in `data/kim2021/`) has parent calls only. The maintainer is finding a fallback; the round-trip phase takes whatever dataset is supplied.
     - Out of M2: non-soybean crops and the crop profiles (M3; now a contract version after 1.3.0, since 1.3.0 is the lines-and-rows change); the `glyma.Wm82.gnmN.` prefix.
- Model for next phase: M3 has no plan yet. It starts the way M2 did — a planner (fable) writes `docs/m3-phases.md`, the main session puts its flagged questions to the maintainer before any dispatch, then one implementer per phase with a reviewer on every core, io, model or contract diff. What M2 taught, and what the next planner should be told: demand a discriminating input for every new test. Five assertions this milestone would have passed with the behaviour they covered removed, and no gate caught any of them: the informative-marker restriction had no failing test at all; the dedup test compared two samples identical at every indexed marker; `test_no_duplicates_in_fixture` asserted an empty list that is also the default; the KASP `Uncallable` case reached the same `N` by both paths until a stderr count was added (the doer caught that one itself); and the KASP grid test checked one cell, so a transposed sample assignment would have passed. A sixth, the Rank keyboard test, is a different class — it pinned third-party behaviour that differs by runner, and it went red in CI; and push before calling a phase finished, because `e2e`, `export` and `r-reader` are CI-only and one of them went red on a runner difference the local run could not reproduce.
- Review required: on any diff touching core, io, model or contract paths; one pass at the end of M2.
- Open maintainer questions: B (two-generation dataset; it blocked only the real-data round trip, which the M2 verification block records as not run — the synthetic round trip shipped in phase 8); `io/export.py` does not validate `next_generation`, so an empty label writes ids like `BC2F1-F1-001--001` and a blank `generation` column on both the CLI and the screen, the same shape as the placeholder-row count closed in 6b1d3a4 and the obvious first M3 fix; whether `possible_duplicate` should mean something narrower than raw IBS (raised 2026-09-21 by the phase 8 reviewer, probably M3): the measure saturates, since sib IBS over informative markers has mean about 1 - 0.25*(parent het), so segment-identical sibs of a nearly homozygous parent sit above any threshold below 1.0 and a real duplicate with ordinary call error sits below 1.0. What would discriminate for a breeder is whether the pair shares `family_id` — a same-family late-generation pair is expected, a cross-family or cross-generation pair at 0.995 is the alarming one — and whether discordance is exactly zero. ADR 0017's BC3F1 section has the numbers; whether `validate` should report pairs skipped for insufficient overlap, since ADR 0017's amendment records that this false negative is currently silent (raised 2026-09-21, not urgent: the flag is advisory); and whatever `docs/m2-phases.md` flags. Known limits deferred from M1 (`docs/m1-phases.md`, phase 4 section): no "(no generation)" radio choice; Rank caption omits an unassigned family; the Navigate stale-input guard has no browser test.
- Deferred reviewer findings, for the first M2 phase (earlier two first, then the M1 milestone review of 2026-09-14; locate by function name, some line numbers were diff positions):
  - `io/criteria.py` `yaml.safe_load` resolves anchors and aliases, so a hostile pasted document could exhaust memory (local effect only); `app/screens/load.py` Apply catches only CriteriaError and DataContractError.
  - medium, do first: `core/qc.py` `qc_table_rows` marks `qc_excluded` from flags alone, while `core/score.py` excludes only when `filters.exclude_qc_flagged` is true, so Validate and Rank disagree under `exclude_qc_flagged: false`. Needs a design choice: pass `filters` in, or derive from `exclusion_reason`. `tests/test_qc_table.py` only runs the fixture criteria.
  - medium (test gap): Compare drag values are asserted only for BC2F1-F1-001, which is manifest index 0 and rank 1, so indexing by row position would pass. Assert BC2F1-F2-015 (manifest index 34, rank 11, `drag_total_max_cm = 109.534`).
  - low: `io/criteria.py` `_normalise_locus` does not coerce `marker_id`/`left_marker`/`right_marker` to str (YAML int IDs report "not in genotype file"), and silently truncates float `start_bp`; `_check_number` accepts `.nan`/`.inf` (NaN composite scores).
  - low: Load `accept=` omits `.bgz`, `.tsv`, `.hapmap`; `criteria_to_dict` omits None keys although `docs/data-formats.md` says every key is explicit; Export results CSV with no result has no header row, while the manifest CSV has one.
  - low: CI shinylive asset cache key ignores the unpinned shinylive version; the `export` extra still lists `build`; the export smoke test records `page.on("request")`, which may miss service-worker fetches (confirm, or record on the context); the Navigate "(no family)" path has no UI test.
- Coordination: the sibling is `backcross` (GitHub piercetaylor/backcross, package name backcross), at `../backcross` since the folder was renamed on 2026-09-23, so `python3 scripts/check_contract.py ../backcross` is the working form. On 2026-09-16 it reported 159 files identical, both trees at uncommitted 1.3.0. Any contract change starts in backcross and is mirrored here by copying `contract/`. Contract 1.6.0 (backcross docs/adr/0021, here docs/adr/0026), committed 2026-09-23: a pair of one nucleotide and one of `N`, `-`, `.` (`AN`, `A-`) is read as missing, which both tools already did. A pair of two missing characters (`N/N`, `..`) remains undefined by the contract and is read as missing by both. The contract text says so explicitly, because a draft that omitted it would have let the stricter-reading rule break such files in a minor bump. Contract 1.7.0, chromosome schemes for cowpea, pea and peanut, is committed in backcross and mirrored here at `cbd8fe8` (docs/adr/0027); `scripts/check_contract.py ../backcross` reports 274 files identical. Two spec inaccuracies found by doing the work, both resolved against the fixture rather than the text: phase 7's "Add top N per family" arithmetic assumed two selected individuals where `core/selection.py` `select_top_n` gives three, because it unions across every family and not only the ones already selected; and phase 8's "selected parent with a missing target call" edge case is vacuous, since the generator never punches a missing call at the target index and such a parent could not pass the foreground filter to be selected. Check phase 11's acceptance numbers against the same possibility. Unexercised after phase 10: `soysnp_positions.py --download` verifies each directory's `CHECKSUM.*.md5`, but only the `.gff3.gz` files were ever fetched, so the checksum path has never run against a live checksum file. Decision 9's long-format KASP headers (`SubjectID`, `SNPID`, `Call`) are still unconfirmed against a primary source; no real LGC export has been seen, and `docs/adr/0019`'s addendum says so. Other Claude sessions may write in this working tree at the same time, as the backcross session did while contract 1.6.0 and 1.7.0 were mirrored here during M3. Commit only by explicit pathspec, never `git add -A`; treat a gate failure in a file outside the phase's list as another writer's edit in progress and re-run before reporting it; and edit a shared file such as `CHANGELOG.md` only after the other writer's change has been committed. `ListAgents` shows which sessions are live.
- Gate baseline, 2026-09-26 at `829f247`: 617 unit tests pass, 44 deselected; `ruff check .`, `ruff format --check .` and `mypy` clean; `scripts/check_contract.py ../backcross` reports "contract 1.7.0: manifest ok, 273 files" and "contract mirror: 274 files identical"; the fixture regenerates byte-identical; `pytest -m e2e` 43 passed, 1 skipped (the export smoke test, run separately with `PS_SITE_DIR`). CI run 36265695866 green on all six jobs against `1f3247c`: `check (3.11)`, `check (3.12)`, `e2e`, `export`, `r-reader` and `deploy-pages`. Run `pytest` with no extra `-q`: the project `addopts` already carries one, and a second reaches quiet level 2 and suppresses the summary line. The Shinylive site measures 44.5 MB, built 2026-09-26 from `829f247` (docs/limits.md), superseding the M1 site figure. Browser waits default to 30 s (PS_E2E_TIMEOUT_MS). Note that `git diff --exit-code -- tests/fixtures` cannot see an untracked fixture directory, so a new fixture's determinism rests on `git hash-object` until its first commit.
