# MIT licence

Status: accepted. Date: 2026-09-04.

## Context and Problem Statement

The repository will be public on the maintainer's GitHub. Which licence, given the licences of the projects whose ideas were reused?

## Considered Options

MIT; BSD-2-Clause; GPL-3; Apache-2.0.

## Decision Outcome

MIT. No code was copied from any reference repository (docs/reference-repos.md). The RPP weighting and linkage-drag definitions come from Flapjack's documentation (BSD-2-Clause project), the IBS definition from SNPRelate's manual (GPL-3 project); both are formulas, not code. Runtime dependencies: numpy (BSD-3), pyyaml (MIT), shiny (MIT), pandas (BSD-3). If GPL code were ever incorporated the derivative would have to be GPL and this record would be superseded.

### Consequences

Good: compatible with the sibling isoline-browser (MIT), BrAPI (MIT) and Breeding Insight (Apache-2.0) ecosystems; no obstacle to inclusion in a university pipeline. Bad: no patent grant; not a concern here.
