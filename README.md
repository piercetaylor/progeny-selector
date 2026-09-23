# progeny-selector

Status: scaffold / pre-alpha (0.1.0 unreleased).

progeny-selector ranks and selects progeny in marker-assisted backcross (MABC) programs: foreground status at target loci, negative selection at avoid loci, recurrent-parent genome recovery overall and on carrier versus non-carrier chromosomes, donor-segment (linkage drag) bounds with recombinant flags on each flank, similarity to each parent, QC flags for selfs, outcrosses and swaps, and a weighted composite score applied after hard filters. It reads VCF, HapMap or wide CSV genotypes with a samples.csv manifest and a criteria.yaml file, and writes results.csv, selected.csv and a sample manifest for the next genotyping round. The compute core is a pure Python package; the UI is Shiny for Python, run locally or exported with Shinylive as a static site so genotype data never leaves the browser tab.

## What exists now

Parsers for all input formats and the criteria file, with an `assembly` selector (Wm82.a1, a2, a4 or none) for chromosome lengths and beyond-length warnings; classification, foreground (including per-target windows), avoid, background (count and weighted RPP, weighted the default), drag, similarity, QC (including the advisory `possible_duplicate` flag and duplicate pairs), scoring, ranking (weighted or staged), selection and projection; results/selection/manifest writers with results.csv and selected.csv columns frozen under `results_schema` and `NA` for missing cells; a CLI (`validate`, `rank`, `select --per-selected N`, each taking `--crop ID` for one of nine crop chromosome schemes); a Shiny app with all seven screens running the real pipeline — Load, Validate/QC (background model, units, assembly, rank mode and duplicate pairs shown), Navigate (tree and breadcrumbs), Rank (status chips, filters, multi-select), Compare (side-by-side cards and chromosome strips), Selection list (editable notes) and Export; a documented keyboard walkthrough; a Shinylive static export that stages the package into the bundle and is smoke-tested in a real browser; a synthetic BC2F1 fixture, a BC3F1 fixture proving next_samples.csv round-trips unchanged as the next generation's samples.csv, and a passing unit and browser test suite. See PLAN.md, "Milestones" and "Verification status".

## Quickstart

Python 3.11 or later.

```
pip install -e ".[dev]"
pytest
ruff check . && ruff format --check . && mypy
progeny-selector validate --genotypes tests/fixtures/synthetic_bc2f1/genotypes.vcf \
  --samples tests/fixtures/synthetic_bc2f1/samples.csv \
  --markers tests/fixtures/synthetic_bc2f1/markers.csv \
  --criteria tests/fixtures/synthetic_bc2f1/criteria.yaml
progeny-selector rank --genotypes tests/fixtures/synthetic_bc2f1/genotypes.vcf \
  --samples tests/fixtures/synthetic_bc2f1/samples.csv \
  --markers tests/fixtures/synthetic_bc2f1/markers.csv \
  --criteria tests/fixtures/synthetic_bc2f1/criteria.yaml --out results.csv
# --crop maize, rice, sorghum, wheat, barley, oat, common-bean or cotton reads chromosome
# names under that crop's convention; soybean is the default.
progeny-selector select --results results.csv --top 3 --out selected.csv \
  --next-manifest next_samples.csv --next-generation BC3F1 \
  --samples tests/fixtures/synthetic_bc2f1/samples.csv
shiny run src/progeny_selector/app/app.py
Rscript scripts/read_results.R results.csv
```

`python3 scripts/build_shinylive.py` produces the static site in `site/` (see docs/adr/0003, docs/adr/0009).

## Input files

Genotypes (VCF 4.2+, plain or bgzip; HapMap; wide CSV with nucleotide or A/B/H calls), samples.csv with exactly one `recurrent_parent` and one `donor_parent`, optional markers.csv with cM, and criteria.yaml. The contract is in docs/data-formats.md and is shared with the sibling backcross project.

## Documents

PLAN.md (problem, algorithms, UI walkthrough, milestones, verification status), docs/data-formats.md, docs/reference-repos.md, docs/adr/ (MADR records), CHANGELOG.md, CONTRIBUTING.md.

## Limits

Measured wall-clock and peak memory for CPython and Shinylive at several dataset sizes: docs/limits.md.

## Licence

MIT (LICENSE).
