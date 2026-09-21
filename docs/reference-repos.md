# Reference repositories

Every repository and document below was fetched on 2026-09-04 at the URL given; facts are as reported by the fetched page. No code was copied from any of them; what was borrowed is definitions, conventions and layout.

## cropgeeks/flapjack

URL: https://github.com/cropgeeks/flapjack (fetched). Licence: BSD-2-Clause. Language: Java; Ant build; GitHub Actions in .github/workflows; tests/ directory; docs on Read the Docs. Desktop genotype visualisation with a marker-assisted backcrossing analysis documented at https://flapjack.hutton.ac.uk/en/latest/mabc.html (fetched; also the raw source https://raw.githubusercontent.com/cropgeeks/flapjack/master/docs/mabc.rst) and a tutorial at https://flapjack.hutton.ac.uk/en/latest/mabc_tutorial.html (fetched). Borrowed: the weighted RPP model (marker weight = per side min(half the gap to the neighbour, half the maximum coverage); heterozygous calls weigh 0.5), the linkage-drag definition (distance to the first recombination on each side, or to the chromosome end), the QTL file idea of declaring whether the desired allele is from the donor or the recurrent parent (our `required_state` and avoid loci), and the output-table layout of per-chromosome RPP, total, QTL status, drag, selected state, rank and comment. Not borrowed: Java desktop architecture and native formats.

## cropgeeks/flapjack-bytes

URL: https://github.com/cropgeeks/flapjack-bytes (fetched). Licence: BSD-2-Clause. Language: JavaScript; Rollup, Babel, ESLint; test/ and sample-data/; no CI file visible. Canvas graphical-genotype library with BrAPI/file/URL loading. Borrowed: the compact-strip idea for the Compare screen (colour by state per marker on one row per line). Not borrowed: rendering code; detailed browsing is the sibling project's scope.

## plantbreeding/API (BrAPI)

URL: https://github.com/plantbreeding/API (fetched). Licence: MIT. Specification for the Breeding API; V2.1 with Core, Phenotyping, Genotyping (samples, markers, variant sets, variants, call sets, calls) and Germplasm modules. Borrowed: vocabulary for a later loader and the sample/germplasm distinction behind `sample_id` versus `line_name`. Not borrowed: nothing implemented.

## Breeding-Insight (organisation)

URL: https://github.com/Breeding-Insight (fetched). bi-web (TypeScript, Apache-2.0), bi-api (Java, Apache-2.0), brapi (Java, Apache-2.0), deltabreedquery (R, BrAPI to R), familia and AlloMate (Shiny). Borrowed: the pattern of Shiny apps and R helpers around a breeding database, which is how this tool's CSV exports are meant to be consumed. Not borrowed: server architecture.

## zhengxwen/SNPRelate

URL: https://github.com/zhengxwen/SNPRelate (fetched). Licence: GPL-3. Language: R with C/C++; Bioconductor; GitHub Actions; tests/ and vignettes/. IBS definition from https://rdrr.io/bioc/SNPRelate/man/snpgdsIBS.html (fetched): the average over SNPs of 1 − |g1 − g2| / 2. Borrowed: this definition, implemented as shared alleles / 2 per marker (equal for biallelic markers, defined for multiallelic ones). Not borrowed: code (GPL-3).

## cggh/scikit-allel

URL: https://github.com/cggh/scikit-allel (fetched). Licence: MIT. Language: Python; allel/ package, docs/, notebooks/; pytest; Travis/AppVeyor; PyPI; maintenance-only, successor sgkit. Distance documentation at https://scikit-allel.readthedocs.io/en/stable/stats/distance.html (fetched). Borrowed: the src-package/docs/tests layout, pytest, and computing on allele-count encodings; the decision not to depend on it (maintenance mode, compiled extensions unavailable in Pyodide). Not borrowed: code.

## knausb/vcfR

URL: https://github.com/knausb/vcfR (fetched). Language: R with Rcpp; CRAN; GitHub Actions R-CMD-check, AppVeyor, Coveralls. Licence not shown on the fetched page. Borrowed: the expectation that R users will consume results.csv rather than VCF, and the fixed-columns/genotype-matrix separation. Not borrowed: code.

## solgenomics/sgn (Breedbase)

URL: https://github.com/solgenomics/sgn (fetched). Licence: MIT. Language: Perl; lib/, mason/, db/, js/, t/, selenium/. Evaluated; nothing borrowed (server-side database outside scope).

## hashimotoshumpei/GenoSee

URL: https://github.com/hashimotoshumpei/GenoSee (fetched). Licence: MIT. Language: Python 3.6+ (matplotlib, numpy, pandas). Graphical genotype figures from A/B/H/N or VCF-style calls; multiallelic sites unsupported; no CI or tests visible. Borrowed: the A/B/H/N vocabulary for the coded wide CSV. Not borrowed: code.

## StefanReuscher/ABHgenotypeR

URL: https://github.com/StefanReuscher/ABHgenotypeR/ (fetched). Language: R; CRAN; R/, man/, vignettes/; no CI visible; licence not shown on the fetched page. ABH coding with imputation and error correction for parent-based populations. Borrowed: the A/B/H coding convention. Not borrowed: correction functions (out of scope; would be an explicit step).

## posit-dev/py-shiny

URL: https://github.com/posit-dev/py-shiny (fetched). Licence: MIT. Language: Python; shiny/ package, examples/, docs/, tests/, js/; GitHub Actions; pytest; PyPI and conda-forge. Shinylive documentation https://shiny.posit.co/py/docs/shinylive.html (fetched) and the shinylive package https://pypi.org/project/shinylive/ (fetched, version 0.8.11): apps run in the browser on Pyodide 0.27.3 with pure-Python wheels installable from requirements.txt; about 13 MB base download plus NumPy 7.5 MB. DataGrid documentation https://shiny.posit.co/py/components/outputs/data-grid/ (fetched): `selection_mode="rows"`, `filters=True`, `.cell_selection()["rows"]`. Borrowed: the UI framework and deployment model (docs/adr/0001, 0003); the module pattern for screens. Pyodide package list https://pyodide.org/en/stable/usage/packages-in-pyodide.html (fetched): numpy, pyyaml and pandas are included.

## ThinkR-open/golem and Appsilon/rhino

URLs: https://github.com/ThinkR-open/golem (fetched; MIT; R package framework for Shiny with R/, inst/app/www, dev/ scripts, tests/), https://github.com/Appsilon/rhino (fetched; LGPL-3; app/logic and app/view separation, Cypress e2e, renv). Evaluated as the R Shiny alternatives in docs/adr/0001; borrowed the logic/view separation as the core/app split. Not used.

## SoyBase Data Store (data.soybase.org)

URL: `https://data.soybase.org/Glycine/max/markers/` (referenced 2026-09-21). `Wm82.gnm{1,2,4,5,6}.mrk.SoySNP{50K,6K}` directory paths per docs/m2-phases.md decision 7. Fetched and verified 2026-09-21: `glyma.Wm82.gnm1.mrk.SoySNP50K.gff3.gz` (60,800 rows), `glyma.Wm82.gnm2.mrk.SoySNP50K.gff3.gz` (60,556 rows), `glyma.Wm82.gnm4.mrk.SoySNP50K.gff3.gz` (58,394 rows) all returned HTTP 200 at `https://data.soybase.org/Glycine/max/markers/<dir>/<file>`; attributes confirmed `alleles=` on gnm1 and gnm2, `ref_allele=` on gnm4, matching decision 9's correction. `CHECKSUM.*.md5` verification and gnm5/gnm6 were not fetched this session. Used by `scripts/soysnp_positions.py --download` to fetch GFF3 marker files and their `CHECKSUM.*.md5`.

## Searches that found no reusable MABC tool

Searches for open marker-assisted backcrossing or background-selection packages ("marker assisted backcrossing R package github", "background selection software recurrent parent genome") returned PLABSIM (simulation software described in a 2000 article, no maintained repository found), Flapjack, and journal articles describing analyses done in Flapjack or spreadsheets. No maintained open-source repository implementing foreground/background/recombinant selection lists other than Flapjack was found.

## Literature and platform sources cited in PLAN.md

Hospital and Charcosset 1997, Genetics 147:1469–1485, DOI 10.1093/genetics/147.3.1469, https://academic.oup.com/genetics/article-abstract/147/3/1469/6054126 (fetched). Frisch, Bohn and Melchinger 1999, Crop Science 39:1295–1301, DOI 10.2135/cropsci1999.3951295x, https://experts.illinois.edu/en/publications/comparison-of-selection-strategies-for-marker-assisted-backcrossi/ (fetched). Frisch and Melchinger 2005, Genetics 170:909–917, DOI 10.1534/genetics.104.035451, https://academic.oup.com/genetics/article-abstract/170/2/909/6059341 (fetched). Lübberstedt, Beavis and Suza, Molecular Plant Breeding chapter 6, https://iastate.pressbooks.pub/molecularplantbreeding/chapter/marker-assisted-backcrossing/ (fetched). Song et al. 2013 SoySNP50K, PLOS ONE 8(1):e54985, https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0054985 (fetched). Song et al. 2020 BARCSoySNP6K, Plant Journal 104:800–811, https://digitalcommons.unl.edu/cgi/viewcontent.cgi?article=2402&context=agronomyfacpub (fetched). SoyBase genome information https://www.soybase.org/resources/genome_info/ (fetched). VCF 4.2 specification https://samtools.github.io/hts-specs/VCFv4.2.pdf (fetched). HapMap format https://statgen-esalq.github.io/Hapmap-and-VCF-formats-and-its-integration-with-onemap/ (fetched). Okabe-Ito hex values from R grDevices source https://raw.githubusercontent.com/wch/r-source/trunk/src/library/grDevices/R/colorstuff.R (fetched). MADR https://adr.github.io/madr/, Keep a Changelog https://keepachangelog.com/en/1.1.0/, Conventional Commits https://www.conventionalcommits.org/en/v1.0.0/ (all fetched). PubMed Central pages for the two Genetics papers returned a browser-check page and were not used.
