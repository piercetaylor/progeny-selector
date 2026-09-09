# Run locally with `shiny run`, publish as a Shinylive static export on GitHub Pages

Status: accepted. Date: 2026-09-04.

## Context and Problem Statement

The tool must be available to program staff at no recurring cost and must not transmit unreleased genotype data. Where does it run?

## Considered Options

1. `shiny run` on each user's laptop (pip install) or on a university VM.
2. Shinylive static export on GitHub Pages (Pyodide in the browser).
3. Posit Connect Cloud or shinyapps.io.
4. A Docker image on a university server.

## Decision Outcome

Options 1 and 2 from the same code. The CI workflow builds the wheel, runs `shinylive export src/progeny_selector/app site`, and deploys `site/` to GitHub Pages on pushes to main when Pages is enabled; the app's requirements.txt names the package wheel so Pyodide installs it in the browser [web] https://pypi.org/project/shinylive/. For large datasets or batch runs the CLI and `shiny run` on a laptop or HPC login node are the documented path. Option 3 is rejected because unreleased data would be uploaded to a third party; option 4 remains available (a Dockerfile is a later addition) for a shared university instance.

### Consequences

Good: zero cost; a URL the breeder can open anywhere; data stays in the tab; no server to patch. Bad: the first load downloads about 13 MB plus NumPy and pandas [web] https://shiny.posit.co/py/docs/shinylive.html; Pyodide is slower than native Python; GitHub Pages is public, so the code is public (the data never is). Neutral: pure-Python constraint on dependencies (no cyvcf2, no scikit-allel) is already satisfied by the core.
