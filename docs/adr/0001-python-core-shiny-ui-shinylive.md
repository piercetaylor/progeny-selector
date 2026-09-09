# Python compute core with a Shiny for Python UI, deployable locally or as a Shinylive static site

Status: accepted. Date: 2026-09-04. Format: MADR 4.0.0 [web] https://adr.github.io/madr/.

## Context and Problem Statement

progeny-selector ranks hundreds to thousands of segregating backcross progeny per generation against target, avoid and background criteria and produces selection lists. It must run on a laptop in a field office, cost nothing to host at a public university, keep unreleased genotype data off third-party servers, be maintainable by an R/Python-fluent breeder-bioinformatician who builds Shiny dashboards downstream, and be verifiable in CI. Which stack?

## Decision Drivers

Data sensitivity (client-side or self-hosted processing); zero hosting cost; maintainer fluency (R, Python, Snakemake, HPC); table-centric UI (ranked tables, filters, multi-select, chips) rather than pixel rendering; batch use from the command line on an HPC node; interoperability with the sibling isoline-browser through files; the scaffold container has Python 3.11 and Node 22 but no R, so an R choice could not be verified here.

## Considered Options

1. R Shiny structured with golem or rhino.
2. Python FastAPI back end with a TypeScript front end.
3. Client-only TypeScript (Vite) with in-browser parsing.
4. Python package (numpy) with a Shiny for Python UI, run locally with `shiny run` or exported with Shinylive to a static site.

## Decision Outcome

Option 4. The compute core is a pure Python package (`progeny_selector.core`, numpy only) with boundary parsers, a CLI and pytest coverage; the UI is Shiny for Python, whose module and reactive model matches the maintainer's R Shiny experience. Two deployment modes use the same code: `shiny run` on the laptop or a university server, and `shinylive export` to static files served from GitHub Pages, where the app runs on Pyodide in the browser and data never leaves the tab [web] https://shiny.posit.co/py/docs/shinylive.html. numpy, pyyaml and pandas are available in Pyodide [web] https://pyodide.org/en/stable/usage/packages-in-pyodide.html; the package itself is pure Python, so it installs from a wheel URL listed in requirements.txt.

### Consequences

Good: one language the maintainer already uses; the core is importable from Snakemake pipelines and callable from Python notebooks; CI verifies everything in this container; the CSV outputs feed R Shiny dashboards. Bad: Shinylive downloads about 13 MB plus NumPy before the first screen and computes slower than native Python, so very large VCFs (50K markers × thousands of individuals) are better run locally or via the CLI; R users cannot call the core directly; the DataGrid is less keyboard-rich than a purpose-built table component. Neutral: pandas is used only in the UI layer for DataGrid; the core never imports it.

## Pros and Cons of the Options

R Shiny (golem [web] https://github.com/ThinkR-open/golem, MIT; rhino [web] https://github.com/Appsilon/rhino, LGPL-3): the maintainer's home stack and direct reuse in downstream dashboards. Against: every user needs R or a Shiny Server; shinyapps.io free tier would receive unreleased data; no R in the scaffold environment means no verified tests; vcfR-based parsing of 50K × 2000 genotypes is memory-heavy [inference].

Python FastAPI + TypeScript: suited to a multi-user database (the Breeding Insight pattern [web] https://github.com/Breeding-Insight) but requires a server that receives data, two languages, and deployment effort out of proportion for a single-program tool.

Client-only TypeScript: best on data sensitivity and cost and the right choice for the sibling isoline-browser's canvas-heavy UI; for this table-centric tool it would move the maintainer entirely out of R/Python and split the core from the Snakemake/HPC workflows where batch ranking will also run.

Python core + Shiny for Python + Shinylive (chosen): see Decision Outcome. Verified in this container: pytest, ruff, mypy pass and `shiny run` serves the shell; `shinylive export` could not be installed here (lzstring wheel build failed) and is marked unverified in PLAN.md.
