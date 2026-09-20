# progeny-selector

Status: scaffold / pre-alpha (0.1.0 unreleased).

progeny-selector ranks and selects progeny in marker-assisted backcross (MABC) programs: foreground status at target loci, negative selection at avoid loci, recurrent-parent genome recovery overall and on carrier versus non-carrier chromosomes, donor-segment (linkage drag) bounds with recombinant flags on each flank, similarity to each parent, QC flags for selfs, outcrosses and swaps, and a weighted composite score applied after hard filters. It reads VCF, HapMap or wide CSV genotypes with a samples.csv manifest and a criteria.yaml file, and writes results.csv, selected.csv and a sample manifest for the next genotyping round. The compute core is a pure Python package; the UI is Shiny for Python, run locally or exported with Shinylive as a static site so genotype data never leaves the browser tab.

## What exists now

Parsers for all input formats and the criteria file; classification, foreground, avoid, background, drag, similarity, QC, scoring, ranking, selection and projection; results/selection/manifest writers; a CLI (`validate`, `rank`, `select`); a Shiny app with all seven screens running the real pipeline — Load, Validate/QC, Navigate (tree and breadcrumbs), Rank (status chips, filters, multi-select), Compare (side-by-side cards and chromosome strips), Selection list and Export; a Shinylive static export that stages the package into the bundle and is smoke-tested in a real browser; a synthetic BC2F1 fixture with independently computed expected results and a passing unit and browser test suite. See PLAN.md, "Milestones" and "Verification status".

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

## Licence

MIT (LICENSE).
