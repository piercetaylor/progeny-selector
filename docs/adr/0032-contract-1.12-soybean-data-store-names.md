# Contract 1.12.0: the soybean scheme reads the SoyBase / LIS Data Store names

Status: accepted. Date: 2026-09-27. Mirrors backcross docs/adr/0027, the canonical record for this contract version.

## Context and Problem Statement

A VCF downloaded from the SoyBase / LIS Data Store names its chromosomes `glyma.Wm82.gnmN.GmNN`, or `glyma.Wm82.gnm5.ChrNN` for the gnm5 assembly, and Williams 82 V1.1 files use `GLYMAchr_01`. Up to contract 1.11.0 the soybean scheme kept these names as written, so they were ordered after the canonical names and could not match a markers.csv or a target region written as `Gm01`. Contract 1.12.0 adds them to the soybean pattern. This record is the sibling half: what changed in this repository and what it checks. The decision, the verified facts (the NCBI CM000834-CM000853 lineage that carries chromosomes 1-20 in every Wm82 version, and the Data Store's per-assembly `chromosome_prefix`), the options weighed and the reasons each spelling stays kept as written are in backcross docs/adr/0027 and are not restated here.

## Decision Outcome

The canonical record's decision applies here unchanged: under `soybean` the pattern is `^(?:glyma\.wm82\.gnm[0-9]+\.(?:gm|chr)|glymachr|gm|chr|chromosome|lg)?[_\s-]?0*([1-9]|1[0-9]|20)$`, read case-insensitively, with the same keys and canonical names.

- **The literal changes.** `SOYBEAN_SCHEME` in `src/progeny_selector/core/chrom.py` carries the new pattern as a raw string equal to the decoded pattern of `contract/crops/soybean.json`, and the LIS Data Store README for Wm82.gnm5.NRKG is appended to its sources (the NCBI GCF_000004515.6 assembly report was already listed). `normalize_chrom`, `chrom_sort_key` and `chrom_length_bp` are unchanged; the new spellings reach them through the pattern alone.
- **The docstrings no longer say the soybean scheme "reproduces the 1.2.0 rule exactly".** `io/crops.py`, `io/__init__.py` and `model/dataset.py` now say it reads every spelling of the 1.2.0 rule and, since 1.12.0, the Data Store names of the Wm82 reference assemblies and the V1.1 spelling. docs/adr/0020 keeps its dated wording; this record supersedes that sentence.
- **The rule is pinned by hand-built tests, not only by the shared cases.** `tests/test_crops.py` asserts, through the Python literal, that the ten spellings the old pattern read map as before; that the six newly read spellings (`glyma.Wm82.gnm4.Gm01`, `glyma.Wm82.gnm5.Chr13`, `glyma.Wm82.gnm2.Gm20`, `GLYMAchr_01`, `GLYMAchr01`, `glyma.wm82.gnm12.gm05`) map to their Gm names; that the twelve spellings the canonical record keeps as written (`21`, `Gm00`, `scaffold_22`, `glyma.Lee.gnm1.Gm01`, `glyma.Wm82.gnm4.scaffold_22`, `glyma.Wm82.gnm2.Gm21`, `glyma.Wm82.gnm4.LG7`, `glyma.Wm82.gnm6.01`, `NC_016088.4`, `MT`, `Pltd`, `ChrUn`) stay unchanged; and that the names of the shared case order as `Gm07`, `Gm13`, then `glyma.Wm82.gnm4.scaffold_22`.
- **The contract case comes from the sibling unchanged.** `contract/cases/crop-soybean-data-store-spellings` is generated in backcross and runs here through `tests/test_contract_cases.py`; `test_literal_equals_file` compares the literal with the copied `soybean.json`, and `scripts/check_contract.py ../backcross` byte-compares the mirror.
- **`docs/data-formats.md`** moves its version reference to 1.12.0 and describes the soybean scheme in the new words.

## Consequences

A minor bump: no input read before changes meaning, since every newly read name was kept as written and unordered before. A dataset carrying Data Store names now merges them onto `Gm01`-`Gm20`, orders them with the canonical names, matches markers.csv rows and target regions written as `Gm01`, and takes the Williams 82 chromosome lengths for drag bounds and marker weights. Scripts that strip the Data Store prefix before normalising (`scripts/soysnp_positions.py`, `scripts/soysnp50k_nils.py`) give the same result as before. `results_schema` does not move and no error kind is added.

## Revisit when

- A RefSeq accession table (`NC_016088.4` and its siblings) is wanted, which backcross docs/adr/0027 defers.
- A Wm82 assembly renumbers a chromosome, breaking the lineage the prefix relies on.
- backcross docs/adr/0027 is amended or superseded.
