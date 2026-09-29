# Release decisions and provenance columns: version, cap column name, tool/tool_version/tool_commit, citation, joint release

Status: accepted. Date: 2026-09-29. Decided by research the maintainer delegated to Fable, run jointly with backcross, whose record is its docs/adr/0030.

## Context and Problem Statement

Both tools were about to be released for the first time, and five questions were open, taken for the two tools together. The criterion was the academic and open-source breeder: files that readr, pandas and a spreadsheet read without special handling, output that is never silently wrong, and a release that can be cited.

- **Q1, version.** What number the first release carries, what bumps it, and how the app version relates to the input contract version.
- **Q3, the cap column.** results_schema 1.2.0 (docs/adr/0033) drafted the weighted model's coverage cap as `background_max_coverage`, while backcross names its caps `max_gap_bp` and `max_gap_cm`. The two tools should use one stem.
- **Q4, provenance.** A results.csv or selected.csv said which schema it followed but not which build of the tool wrote it, so two files differing in weighted RPP could not be traced to a code change.
- **Q5, citation.** How a user cites the software, and which metadata file Zenodo reads.
- **Q6, release.** Whether the two tools release together.

## Considered Options

- Q1: 0.1.0 with a minor bump per release; 1.0.0 now.
- Q3: keep `background_max_coverage`; the shared stem `max_marker_coverage`; the stem `max_gap`.
- Q4: a `#` comment line before the header; a sidecar file; constant-per-row columns; nothing.
- Q5: CITATION.cff only; .zenodo.json only; both.
- Q6: release together; release independently.

## Decision Outcome

**Q1.** Both tools tag v0.1.0 and bump the minor version per release. 1.0.0 waits until the columns and flags are frozen and a second cross or crop has been validated on real data. The app version and the input contract version are independent; each release note says "implements input data contract 1.12.0". The classifier at 0.1.0 is "Development Status :: 3 - Alpha".

**Q3.** The shared stem is `max_marker_coverage`, after Flapjack's `mabcMaxMrkrCoverage` and this tool's criteria key `background.max_marker_coverage`. The results.csv column is `background_max_marker_coverage`, with `background_unit` beside it giving its unit. backcross renames `max_gap_bp` and `max_gap_cm` to `max_marker_coverage_bp` and `max_marker_coverage_cm`. "Gap" is avoided because it collides with segments.csv's `gap_criterion`, the run-breaking gap of backcross docs/adr/0008.

**Q4.** Three constant-per-row columns, `tool`, `tool_version` and `tool_commit`. No `#` comment line and no sidecar file.

- `tool` is `progeny-selector`.
- `tool_version` is `progeny_selector.__version__`, which equals the release tag.
- `tool_commit` is `g` plus the first seven hex digits of the commit (for example `g8473292`), suffixed `-dirty` when a tracked file differs from that commit (`git status --porcelain --untracked-files=no` is non-empty; an untracked file does not count), and `NA` when the commit cannot be known: an installed wheel, Pyodide without a recorded build commit, or no git.
- results.csv: directly after `background_max_marker_coverage` and before the dynamic columns, so the fixed prefix is 41 names and the header of an empty file 42. selected.csv: after `sample_db_id`. Not in next_samples.csv, which is the input contract's samples.csv, nor in the `brapi-callsets` output. The contract version is not a column.
- The values name the build that writes the file, so the export writers fill them at write time. selected.csv written by `select` names the build that ran `select`, not the one that wrote the results.csv it read, and a results.csv of schema 1.1.0 without these columns still feeds `select`.
- `scripts/read_results.R` reads all three with `col_character()`.

How `tool_commit` is found (`src/progeny_selector/provenance.py`, cached per process), first hit wins:

1. A file `_build_commit.txt` next to `provenance.py`, stripped and used verbatim. Only `scripts/build_shinylive.py` writes it, into its staged copy of the package, by the rule in 3 run against the source tree; the browser never calls git. The file is gitignored in the source tree.
2. Under Pyodide (`sys.platform == "emscripten"`), `NA`.
3. git, run in the package directory with a 5 s timeout, only when `git rev-parse --show-toplevel` names a work tree whose `src/progeny_selector` is the package directory itself. A wheel installed in a venv that sits inside some other git repository must give `NA`, never that repository's commit. Any OSError, timeout or non-zero exit gives `NA`.

backcross agreed the `-dirty` rule and the build-time staging on 2026-09-29: its CI uses `GITHUB_SHA` on a clean checkout, its Pages site and bundled CLI bake the value in at build time, and a source run asks git at run time.

**Q5.** CITATION.cff only, no .zenodo.json: Zenodo uses .zenodo.json and ignores the CFF when both exist (help.zenodo.org, checked by the backcross session). The CFF is committed before tagging, with `version` equal to the tag. The maintainer enables Zenodo's GitHub integration before the tag, and the DOI goes into `identifiers` afterwards. The release is a normal release, not a pre-release. No ORCID unless the maintainer has one. This tool's CFF: `cff-version: 1.2.0`; `message: "If you use this software, please cite it as below."`; `type: software`; `title: "progeny-selector: marker-assisted backcross progeny ranking and selection"`; `authors: [{family-names: Taylor, given-names: Pierce}]`; `license: MIT`; `repository-code: https://github.com/piercetaylor/progeny-selector`; `url: https://piercetaylor.github.io/progeny-selector/`; `version: 0.1.0`; `date-released:` the tag day; `keywords: [plant breeding, marker-assisted backcrossing, progeny selection, recurrent parent genome, soybean]`; `references:` backcross as `type: software`, same author, `repository-code: https://github.com/piercetaylor/backcross`, `url: https://piercetaylor.github.io/backcross/`, `notes: "Sibling tool; both read input data contract 1.12.0"`. The README gains a "Sibling tool" paragraph: backcross characterises finished NILs; progeny-selector ranks progeny during the programme; both read contract 1.12.0, so genotype, samples.csv and markers.csv files move between them unchanged.

**Q6.** Both tools release v0.1.0 on the same day, and each release note names contract 1.12.0 and links the sibling's release. Later releases are independent. Tagging, the GitHub Release and Zenodo happen only on the maintainer's explicit yes.

### Consequences

Good: a results.csv or selected.csv says which build wrote it, down to an uncommitted working tree; the two tools share one name for the coverage cap; the release can be cited with a DOI. Bad: the header grows by three columns; `tool_commit` is `NA` for a pip-installed wheel, so a wheel user is traced by `tool_version` alone; a committed file must never hold a `tool_version` or `tool_commit` that changes per commit, so tests set `tool_commit` explicitly rather than reading it from the checkout.
