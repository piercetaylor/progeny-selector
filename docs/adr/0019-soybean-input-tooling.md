# Soybean input tooling: genetic maps from Song et al. 2016

Status: accepted. Date: 2026-09-17.

## Context and Problem Statement

Real soybean data reaches this tool with physical positions on a Williams 82 assembly and, usually, no genetic map, so cM-weighted RPP and cM segment bounds stay off. Song et al. 2016 (BMC Genomics 17:33) publish two dense SoySNP50K linkage maps, Williams 82 x PI 479752 (WP) and Essex x Williams 82 (EW), with each marker's position on Glyma1.01 (a1) and Wm82.a2.v1 (a2), in Table S1 (https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12864-015-2344-0/MediaObjects/12864_2015_2344_MOESM1_ESM.xls, CC BY 4.0). How should that table become markers.csv, including cM for markers the maps do not carry, without writing a wrong genetic position silently?

The maintainer delegated this question (docs/m2-phases.md, question 8 and decision 8) to research, with the criterion: whatever serves academic and open-source plant-breeding users best, reproducible with standard tools, tolerant of spreadsheet exports, never silently wrong.

## Scope

This record covers `scripts/song2016_map.py` and `src/progeny_selector/core/genetic_map.py`. The position-table and genotype-conversion scripts are recorded in this file when they are added.

## Decision Outcome

- **Hand export, no download.** Table S1 is `.xls`; the maintainer exports it to CSV by hand into gitignored `data/song2016/`. No spreadsheet dependency is added and the script never fetches it. The export keeps the title row, so the header is row 2, and the script reads the column names published in the table: `ss ID` (marker id `"ss" + int`), `Glyma1.01 Chromosome`/`Coordinate`, `Wm82.a2.v1 Chromosome`/`Coordinate` (`Chr01`; scaffold cells may carry a leading line break), `WP Linkage Group`, `WP linkage position`, `EW Linkage Group ` (trailing space), `EW linkage position`. `--assembly a2|a1` (default a2) chooses the coordinates.
- **One map at a time.** `--map WP|EW`, default `WP` (WP: Williams 82 x PI 479752; EW: Essex x Williams 82; per the article, not verified here). The two maps come from different crosses and are never averaged.
- **Clean by invalidation, not smoothing.** Per chromosome, in order: drop rows without a position on the chosen assembly (scaffolds included) or without a position on the chosen map; drop rows whose linkage group is not the chromosome number (Table S1 writes linkage groups as numbers, `1.0` in the WP column and `1` in the EW column, equal to the chromosome number; the cell is read as a whole number and anything else is a mismatch); then, with candidates sorted by (bp, cM), keep only the markers in the longest subsequence of cM that is non-decreasing along bp (`longest_nondecreasing_mask`; ties allowed; among equal-length subsequences the one keeping earlier markers). This follows MareyMap practice, where markers inconsistent with the Marey curve are invalidated before interpolation (https://cran.r-project.org/web/packages/MareyMap/vignettes/vignette.pdf). It replaces the running maximum first proposed, which flattened local inversions into plateaus and so kept misplaced markers with made-up positions.
- **Interpolate, clamp, never extrapolate.** `interpolate_cm` sorts the mapped markers by (bp, cM), applies the mask, collapses kept markers at equal bp to one anchor at their mean cM so the abscissae are strictly increasing, and uses `np.interp`: linear between kept markers, clamped to the end cM beyond them. One mapped marker gives its cM to every query on its chromosome; a chromosome with none leaves `cm` empty and the script warns that cM mode will be disabled, since the contract requires cM on every marker for cM mode.
- **One assembly.** `--markers-in` must be on the assembly chosen. Its markers are matched to the table by `ss ID`, else by `SNP ID`; when none matches, the assembly cannot be checked and the script stops; a matched marker at a different position stops it too, because interpolating a4 positions against a2 anchors would be silently wrong. A `--markers-in` marker that already has a cM keeps it; one matched to a mapped row takes that row's cM. Kept table markers at a shared bp are written at their mean cM, the anchor value, so cM is non-decreasing along the file; a marker id is written once, under the `--markers-in` id when matched by `SNP ID`.
- **Nothing kept stops.** When no marker is kept on any chromosome the script exits with the largest drop count and writes no markers.csv; a single chromosome with none kept is a warning.
- **Counted.** A `.log` beside `--out` records rows read, no position on assembly, no linkage position on map, LG mismatch, non-monotone dropped, mapped, cM kept from --markers-in, interpolated, clamped (a subset of interpolated) and unplaced. Column names are fixed by decision 8; there is no `--columns` option.
- **Scripts, not io.** These are one-off conversions of third-party tables outside the input contract; the loader keeps reading markers.csv only. The pure computation lives in `core/genetic_map.py` so it is unit-tested on hand-computed cases.

### Consequences

Good: markers.csv built this way is reproducible from a public CC BY 4.0 table and one command; every dropped marker is counted; interpolation is monotone by construction.

Bad: Table S1 carries a1 and a2 positions only, so data on Wm82.a4 or later must first be placed on a2 (or a1); markers dropped as non-monotone lose their published cM and receive an interpolated one when they are in `--markers-in`; the hand export is a manual step, and a change in the exported header names stops the script with the missing column named.

## Addendum: SoySNP position table and KASP converter (Phase 10)

This section records `scripts/soysnp_positions.py` and `scripts/kasp_to_wide.py`, added under decisions 7 and 9 (docs/m2-phases.md).

- **Download, verified.** `--download` (never run in CI, never exercised by the test suite) fetches the SoyBase Data Store GFF3s named in `SOYBASE_URLS` (`https://data.soybase.org/Glycine/max/markers/Wm82.gnm{1,2,4,5,6}.mrk.SoySNP{50K,6K}/glyma.<dir>.gff3.gz`, licence Open) with a browser User-Agent, and verifies each directory's `CHECKSUM.*.md5` before use. There are no `##sequence-region` pragmas.
- **`marker_id` from `Name=`, not `ID=`.** The GFF3 `ID=` attribute is `glyma.Wm82.gnmN.ss...`, an assembly-specific id; `Name=` is the marker id shared across assemblies, which is what the joined table keys on.
- **Chromosome by stripping, then normalising.** The `seqid`'s `glyma.Wm82.gnmN.` prefix is stripped the way `scripts/soysnp50k_nils.py::strip_chrom` does (text after the last `.`), then passed through `core.chrom.normalize_chrom`; this needs no gnm5-specific code, because gnm5's raw seqids are `Chr01`..`Chr20` and `normalize_chrom` already accepts a `Chr` prefix. Scaffolds keep their stripped name.
- **Panel from the file name.** The phase names a `--panel` option "per file via name" without defining its syntax; since SoyBase's own directory names carry SoySNP50K/SoySNP6K, the panel is read directly from each `--gff3` path instead of a separate flag.
- **No allele column.** Decision 9's correction that allele attributes are `alleles=` (gnm1/2 50K and all 6K) or `ref_allele=` (gnm4/5/6 50K) is noted, but Phase 10's output table has no allele column to apply it to, so nothing reads those attributes yet.
- **KASP: two shapes.** `kasp_to_wide.py` accepts long (`SubjectID`/`SNPID`/`Call` headers, case-insensitive; unverified against a primary source) and grid (SNPviewer's export, orientation detected by matching `--markers` ids against the header row or the first column). A call matches `^[ACGT][:\-.]?[ACGT]$`; the phase table's narrower `^[ACGT]:[ACGT]$` for the long shape is superseded by this decision-9 pattern for both shapes, so a call with no separator (`AG`) or a `-`/`.` separator is also accepted. `Uncallable`, `Missing`, `?`, `Bad`, `Dupe`, `NTC`, empty and any other non-matching text become `N`, counted per reason.
- **Conflicts and controls.** Two rows for the same sample x marker with different resulting calls raise, naming both source lines; identical repeats do not. `--controls` (default `NTC`) drops sample rows/columns before conversion, in both shapes.

### Consequences

Good: one converter handles both a QC lab's long export and SNPviewer's grid export; every missing call is attributed to a reason instead of being silently coalesced.

Bad: the long-format header names are not confirmed against a KASP/LGC primary source (decision 9's own caveat); grid-orientation detection can be ambiguous when `--markers` matches equally well (or not at all) on both axes, which the phase table does not define a tie-break for beyond "in either orientation" (this implementation prefers markers-as-columns on a tie and errors when neither axis matches).
