# The assembly is a criteria.yaml key, and every fallback from it is a warning

Status: accepted. Date: 2026-09-16.

## Context and Problem Statement

Chromosome ends are a run parameter: they normalise the drag component, cap the terminal marker weights of the weighted RPP model and set the width of every chromosome strip. Until now `core/chrom.py` hardcoded the Wm82.a4.v1 table, so a2-positioned data (the SoySNP50K NIL families, whose positions are Wm82.a2.v1) was measured against the wrong ends with only a README caveat, and a marker beyond the assumed end was silently absorbed. Where does the assembly belong, which assemblies are offered, and what happens when a length is missing or contradicted by the data?

## Considered Options

1. A column in markers.csv.
2. A control on the Load screen only.
3. A top-level `assembly` key in criteria.yaml.

## Decision Outcome

Option 3, `assembly` (default `Wm82.a4`). markers.csv is contract-owned, so a column there starts in isoline-browser and binds both tools (docs/adr/0010); a UI-only control would not reach the CLI and would not round-trip through the criteria download. criteria.yaml is this repo's file, is versioned with the run, and already carries the other run parameters. Values: `Wm82.a1`, `Wm82.a2`, `Wm82.a4` and `none`; `constants.py` keys one table per assembly and `core/chrom.py` reads `chrom_length_bp(name, fallback, assembly)`. The three tables are SoyBase's published lengths (https://www.soybase.org/resources/genome_info/), 1-based chromosome ends as listed. Checked 2026-09-19 against NCBI's sequence reports for GCF_000004515.3, .4 and .6: a1 and a4 are uniformly one greater than NCBI's lengths, SoyBase's end convention rather than an error here, and ten a2 chromosomes are 300 bp to 3.6 kb longer because SoyBase's Wm82.a2.v1 differs from NCBI's copy; SoyBase governs, since SoySNP50K positions are published against it. `Wm82.a5` and `Wm82.a6` are left out and wait on verified lengths: SoyBase publishes no length table for them, and the figures circulated with the 2026-09-16 decision were never checked against the Data Store `genome_main.fna.gz` files. Adding an assembly is a table plus a name in `ASSEMBLIES`.

`none` and any chromosome the chosen assembly does not cover (a scaffold, a non-soybean name) fall back to the last marker on that chromosome, as before, but the fallback is now a warning rather than silence, and so is a chromosome whose markers run past the assembly length — the signature of data positioned on a different assembly than the one selected. At most five missing-length chromosomes are listed, then `... and N more`, so a scaffold-rich file does not bury the other warnings. Two more warnings state when `cm` was requested but the map has no cM (separately for `background.map_unit` and `flank_unit`) and when the weighted model therefore runs in bp with a bp cap. The UI already renders `AnalysisResult.warnings`.

### Consequences

Good: a2 data can be analysed correctly; a mismatch between the data's assembly and the setting is visible in the UI, the CLI and the warning list rather than hidden in a slightly wrong drag normalisation; the CLI and the browser behave identically because the setting travels in the criteria file. Bad: one more key to set per run, and a file written against the default is silently assumed to be a4 — the beyond-length warning only catches mismatches large enough to push a marker past the end. Neutral: under `none`, and on any chromosome the assembly does not place, the drag right bound stops at the last marker, while the weighted RPP leaves that chromosome out of its length table so the terminal markers keep the coverage cap's half as their outer weight rather than nothing; a terminal marker must not drop out of the background score because a length is unknown.

## Amendments, 2026-09-21

**A chromosome whose markers run past the recorded length is treated as having no usable length.** Maintainer decision, taken after the end-of-milestone M2 review, which found the inconsistency: `core/pipeline.py` `_chrom_lengths` set the chromosome end to the last marker when markers ran beyond the assembly length, but did not add that chromosome to the set of chromosomes without a length. Under `model: weighted` with `map_unit: bp` the terminal marker's outer weight was therefore `min(length − p[-1], cap)` against the *last-marker* end, which is zero whenever the last marker on that chromosome is informative, while a chromosome the assembly gives no length for at all kept the coverage cap's half. The same data — markers whose positions the chosen assembly cannot account for — was scored two opposite ways depending only on whether the assembly happened to list the chromosome name. Gm01 markers at 10 Mb to 60 Mb under `assembly: Wm82.a4`, weighted bp, gave weights `[4e6 ×5, 2e6]`; the same markers under `assembly: none` gave `[4e6 ×6]`.

A beyond-length chromosome is now added to the same set, so the weighted model leaves it out of its length table and its terminal markers weigh the cap's half, exactly as the existing fallback path does. This also moves the *first* marker's outer side from `min(p[0], cap)` to the cap, which matters only when that marker sits inside half the cap of the chromosome start; that follows from the decision as written, "the same as an unknown-length one", rather than being special-cased. Read it that way and not more narrowly: **both** outer sides of the chromosome take the cap, not only the terminal one, because the chromosome has no usable recorded length at either end. The chromosome end recorded for drag bounds and the chromosome strips is unchanged: still the last marker. The beyond-length warning still fires, and remains the way the breeder learns the assembly setting looks wrong; nothing here makes a mismatch quieter.

The real-data case that motivated it: SoySNP50K NIL material carrying Wm82.a2.v1 positions, run on this tool's `Wm82.a4` default, warns on Gm04, Gm06, Gm08, Gm13 and Gm20 — the chromosomes whose a2 positions run past the a4 end. That is precisely the run in which the old behaviour bit.

This narrows the Consequences' Neutral clause above, which reads "under `none`, and on any chromosome the assembly does not place". Read it now as: under `none`, on any chromosome the assembly does not place, **and on any chromosome whose markers run past the assembly length**.
