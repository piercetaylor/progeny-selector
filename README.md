# progeny-selector

progeny-selector ranks candidates in marker-assisted backcross breeding. It checks target and avoid loci, estimates recurrent-parent recovery and donor-segment bounds, applies quality flags and selection criteria, and exports ranked results, selections, and a manifest for the next generation. The analysis runs as a Python package and CLI, with a Shiny interface for local use or a [browser-hosted site](https://piercetaylor.github.io/progeny-selector/).

**Status:** Pre-release software (version 0.1.0 is untagged). The synthetic BC2F1 fixture exercises the full selection pipeline; its rankings are test cases, not breeding recommendations.

## Data and interpretation

Inputs are VCF, HapMap, or wide CSV genotypes; `samples.csv` identifies one recurrent and one donor parent; `criteria.yaml` defines the selection rules. `markers.csv` can supply genetic-map positions. The shared [input contract](contract/data-contract.md) and [selection and output formats](docs/data-formats.md) document the fields. Crop-specific chromosome conventions cover soybean, maize, rice, sorghum, wheat, barley, oat, common bean, cotton, cowpea, pea, and peanut.

Marker spacing limits how tightly the software can bound donor segments. Criteria, assembly choice, missing calls, and quality flags affect rankings and should be reviewed before selecting plants. The browser-hosted Shiny export processes uploaded files in the browser tab. No real genotype dataset is distributed here.

## Run locally

Python 3.11 or newer is required.

```sh
python -m pip install -e ".[dev]"
progeny-selector rank --genotypes tests/fixtures/synthetic_bc2f1/genotypes.vcf --samples tests/fixtures/synthetic_bc2f1/samples.csv --markers tests/fixtures/synthetic_bc2f1/markers.csv --criteria tests/fixtures/synthetic_bc2f1/criteria.yaml --out results.csv
progeny-selector rank --brapi-url https://host/brapi/v2 --variant-set VS1 --samples samples.csv --criteria criteria.yaml --out results.csv
progeny-selector select --results results.csv --top 3 --out selected.csv --next-manifest next_samples.csv --next-generation BC3F1 --samples tests/fixtures/synthetic_bc2f1/samples.csv
```

`progeny-selector validate` checks inputs and reports quality warnings. `shiny run src/progeny_selector/app/app.py` opens the local interface. The [keyboard walkthrough](docs/keyboard-walkthrough.md) covers its screens, and the [accessibility review](docs/accessibility.md) records the WCAG 2.2 AA check against every screen.

## Verification and documentation

`pytest -q`, `ruff check .`, `ruff format --check .`, and `mypy` check the package. Browser and static-export checks are described in [the plan](PLAN.md), along with real-data verification and its limits. [Design decisions](docs/adr/) and the [archived README](docs/legacy-readme.md) retain more detail.

## Limits

Measured wall-clock and peak memory for CPython and Shinylive at several dataset sizes: [docs/limits.md](docs/limits.md). The Shinylive build loaded and analyzed 50,000 markers by 2,000 individuals in the browser in 226 s.

## Licence

The software is available under the [MIT license](LICENSE). There is no associated paper; cite this repository with the commit or version used.
