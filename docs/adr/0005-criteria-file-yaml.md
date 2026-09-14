# Selection criteria as one YAML file with strict keys

Status: accepted. Date: 2026-09-04.

## Context and Problem Statement

A selection run needs target loci with required donor states, avoid loci, flank windows, background options, weights and filter policies. Loci can be a marker, a region or a flanking pair. How is this expressed so that it is versionable, reviewable by a breeder, and validated before any computation?

## Considered Options

1. Two CSVs (loci table and a key-value settings table).
2. One YAML document with nested sections.
3. UI-only configuration saved as JSON.

## Decision Outcome

Option 2 (docs/data-formats.md, "criteria.yaml"; `progeny_selector.io.criteria`). Unknown keys anywhere are errors so typos fail at the boundary; the `region: "Gm06:24,000,000-27,000,000"` shorthand keeps region loci to one line; defaults are defined once in `progeny_selector.model.criteria` dataclasses. PyYAML is pure Python and ships with Pyodide, so the same reader runs in the browser. The UI edits the same structure and offers it for download, so a run is reproducible from files alone.

### Consequences

Good: criteria live next to the data in the project folder and in git; the CLI and UI share one reader; per-target windows and per-locus rules are expressible. Bad: breeders who prefer spreadsheets must learn a small YAML schema (a CSV loci import is a possible later extension); YAML's implicit typing means quoted region strings are required.

## Amendment, 2026-09-12

The in-memory `Criteria` object is the single source of truth for a run; any editor, whether the M1 YAML text area or a later field-by-field form, is a view that reads from it and writes back to it. Serialisation goes through one path, `criteria_to_dict` and `dump_criteria_yaml` in `progeny_selector.io.criteria`, the inverse of `criteria_from_dict`, and `criteria_from_dict(criteria_to_dict(c)) == c` is tested. The emitted form is canonical: every scalar key explicit with defaults included, loci written with only the keys of their kind, and the `region:` shorthand accepted on input and normalised to `chrom`/`start_bp`/`end_bp` on output. Download always serialises the applied `Criteria`, never the raw editor text. Saved presets are named criteria.yaml files in the project folder, not an in-app store.
