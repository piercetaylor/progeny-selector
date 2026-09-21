# Soybean inputs

Scripts that turn public soybean resources into files this tool reads. Everything they read or write under `data/` is gitignored and never committed. Decisions: docs/adr/0019.

## Song et al. 2016 genetic maps: `scripts/song2016_map.py`

Source: Song et al. 2016, BMC Genomics 17:33, Table S1 (https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12864-015-2344-0/MediaObjects/12864_2015_2344_MOESM1_ESM.xls). Licence: CC BY 4.0; cite the article when you publish results that use the map.

Preparing the input: open the `.xls` in a spreadsheet program and save the sheet as CSV, unchanged, to `data/song2016/table_s1.csv`. Keep the title row; the header is row 2. The script reads the columns `ss ID`, `Wm82.a2.v1 Chromosome`, `Wm82.a2.v1 Coordinate` (or `Glyma1.01 Chromosome`, `Glyma1.01 Coordinate`), and `WP Linkage Group`, `WP linkage position` (or `EW Linkage Group `, `EW linkage position`).

```
python3 scripts/song2016_map.py --table data/song2016/table_s1.csv --map WP --assembly a2 --out data/song2016/markers_wp_a2.csv
python3 scripts/song2016_map.py --table data/song2016/table_s1.csv --map WP --assembly a2 --out data/song2016/markers.csv --markers-in data/my_markers.csv
```

- `--map WP|EW` (default `WP`): the two maps of Table S1 (Williams 82 x PI 479752 and Essex x Williams 82, per the article; not verified here). The two maps are never averaged.
- `--assembly a2|a1` (default `a2`): Wm82.a2.v1 or Glyma1.01 coordinates.
- `--markers-in`: a markers.csv (`marker_id,chrom,pos_bp`, optional `cm`) for your panel, on the same assembly. Markers are matched to the table by `ss` id or by `SNP ID`; at least one must match, and a matched marker at a different position stops the script. A marker with a `cm` keeps it; one matched to a kept table marker takes its cM; other markers get cM interpolated between the kept markers of their chromosome (markers at the same bp count once, at their mean cM), clamped to the end cM beyond them; markers on a chromosome with no kept marker get an empty `cm`, which disables cM mode for a dataset that includes them.

Output: `marker_id,chrom,pos_bp,cm`, sorted by chromosome, position and cM, and a `.log` with the same name beside it counting rows read, no position on assembly, no linkage position on map, LG mismatch, non-monotone dropped, mapped, cM kept from --markers-in, interpolated, clamped (a subset of interpolated) and unplaced. When no marker is kept on any chromosome, the script stops with the largest drop count and writes no markers.csv. Kept table markers that share a position are written at their mean cM, so cM never decreases along the file. Each marker_id appears once: a `--markers-in` marker with the same id as a table marker, or matched to it by `SNP ID`, is written once under the `--markers-in` id, with its own cM when it has one. Column names are fixed; there is no `--columns` option.

Cleaning: rows on scaffolds or without a coordinate are dropped; rows whose linkage group (a number such as `1` or `1.0`, equal to the chromosome number in Table S1) is not the chromosome number are dropped; within a chromosome, with markers ordered by position and then cM, only the longest run of markers whose cM never decreases along bp is kept, so a marker placed out of order on the map is dropped instead of bending the curve.

Assembly caveat: Table S1 has a1 and a2 positions only. Genotype data on Wm82.a4 or later must be placed on a2 (or a1) first, and `--markers-in` positions must be on the assembly passed with `--assembly`; a marker found in both files at different positions stops the script. Set `assembly` in criteria.yaml to the same assembly when you run the analysis.

## SoySNP marker positions: `scripts/soysnp_positions.py`

Source: SoyBase Data Store, `https://data.soybase.org/Glycine/max/markers/<dir>/glyma.<dir>.gff3.gz` with `<dir>` in `Wm82.gnm{1,2,4,5,6}.mrk.SoySNP50K` and `Wm82.gnm{1,2,4,5,6}.mrk.SoySNP6K` (licence Open, per each directory's `README.*.yml`). Each GFF3 row is `seqid source type start end score strand phase attributes`; `marker_id` is the `Name=` attribute, not `ID=` (`ID=` is `glyma.Wm82.gnmN.ss...`). The chromosome is `seqid` with the `glyma.Wm82.gnmN.` prefix stripped (the text after the last `.`, as `scripts/soysnp50k_nils.py::strip_chrom` does) and then normalised; scaffolds keep their stripped name. There are no `##sequence-region` pragmas.

```
python3 scripts/soysnp_positions.py --gff3 Wm82.a2=data/soybase/glyma.Wm82.gnm2.mrk.SoySNP50K.gff3.gz \
    --gff3 Wm82.a2=data/soybase/glyma.Wm82.gnm2.mrk.SoySNP6K.gff3.gz \
    --out data/soybase/soysnp_positions.csv
python3 scripts/soysnp_positions.py --gff3 Wm82.a4=data/soybase/glyma.Wm82.gnm4.mrk.SoySNP50K.gff3.gz \
    --out data/soybase/soysnp_positions.csv --emit-markers-csv Wm82.a4 --markers-out data/soybase/markers.csv
```

- `--gff3 ASSEMBLY=PATH` (repeatable): `ASSEMBLY` is one of `Wm82.a1`, `Wm82.a2`, `Wm82.a4`, `Wm82.a5`, `Wm82.a6`. The panel (SoySNP50K or SoySNP6K) is read from the file name, since SoyBase's directory names carry it and there is no separate `--panel` flag.
- `--out`: the joined wide table `marker_id, in_SoySNP50K, in_SoySNP6K, chrom_Wm82.a1, pos_bp_Wm82.a1, chrom_Wm82.a2, pos_bp_Wm82.a2, chrom_Wm82.a4, pos_bp_Wm82.a4, chrom_Wm82.a5, pos_bp_Wm82.a5, chrom_Wm82.a6, pos_bp_Wm82.a6`, one row per `marker_id`, blank where an assembly's file was not given or does not carry the marker; booleans are `TRUE`/`FALSE`.
- `--emit-markers-csv ASSEMBLY --markers-out PATH`: a `marker_id,chrom,pos_bp` file for one assembly, Gm01..Gm20 only (scaffolds excluded), sorted by chromosome then position; loads directly through `io.manifest.read_markers`.
- `--download`: fetches every file in `SOYBASE_URLS` into `data/soybase/` with a browser User-Agent and verifies each directory's `CHECKSUM.*.md5`; never run in CI or by the test suite.

`--emit-markers-csv Wm82.a5` and `Wm82.a6` produce a markers.csv on positions no `assembly` value in criteria.yaml can name: `constants.py` carries no length table for a5 or a6 (docs/adr/0015), because no published SoyBase table exists and the figures once circulated for them were never checked against the Data Store FASTA files. A run against a markers.csv emitted this way falls back to the last marker on each chromosome for the drag right bound and the weighted RPP's terminal weight, and reports beyond-length warnings against whichever assembly `assembly` is set to, since every a5 or a6 position beyond the a1/a2/a4 tables' shorter ends is read as past the chosen assembly's length.

The script also prints each given assembly's maximum position per chromosome, to sanity-check against `constants.py`'s chromosome-length tables.

## KASP export conversion: `scripts/kasp_to_wide.py`

Converts an LGC KASP export to the wide-CSV genotype contract. Two input shapes are accepted (docs/adr/0019, decision 9):

- **long**: one row per (sample, marker) call, detected case-insensitively by the headers `SubjectID`, `SNPID`, `Call`. This shape is unverified against a primary source.
- **grid**: SNPviewer's export, one row per sample or one row per marker; detected by matching the header row or the first column against the marker ids in `--markers`, in either orientation.

```
python3 scripts/kasp_to_wide.py --kasp export.csv --markers markers.csv --out genotypes.csv
python3 scripts/kasp_to_wide.py --kasp export.csv --markers markers.csv --out genotypes.csv --drop-unplaced --controls NTC,H2O
```

- A call matches `^[ACGT][:\-.]?[ACGT]$` (case-insensitive; separator `:`, `-`, `.` or none) and becomes the nucleotide pair of its first and last character (`A:G` -> `AG`). `Uncallable`, `Missing`, `?`, `Bad`, `Dupe`, `NTC`, empty, and any call that does not match the pattern become `N`, counted per reason and printed to stderr.
- `--controls` (default `NTC`): sample ids dropped before conversion.
- Two rows for the same sample x marker with different resulting calls are an error naming both source lines; identical repeats are not a conflict.
- A marker id absent from `--markers` is an error listing the ids, unless `--drop-unplaced` (counted and dropped).
- Output: `marker_id,chrom,pos_bp,<samples...>`, chrom/pos from `--markers`; sample columns in first-seen order, marker rows sorted by chromosome then position. The output loads through `io.load_genotypes`/`io.load_dataset` with nucleotide coding auto-detected.
