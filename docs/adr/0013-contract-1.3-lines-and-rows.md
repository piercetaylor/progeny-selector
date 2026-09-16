# Contract 1.3.0 mirrored: a blank line before the header, the sniff on the first non-blank line, and `genotypes.column_count`

Status: accepted. Date: 2026-09-16.

## Context and Problem Statement

Contract 1.3.0 is decided in the canonical repository, backcross docs/adr/0018 ("Contract 1.3.0: blank lines, `#` as data, quoted line breaks, two-fault rows and `genotypes.column_count`"), and mirrored here byte for byte under `contract/`. That record carries the decisions and their reasons: blank and whitespace-only lines and all-empty delimited rows are skipped everywhere, including before the header; `#` is not a comment marker; a quoted field may contain a line break; a data line with the wrong field count is the error kind `genotypes.column_count`; a row with two faults reports one of them, unspecified. What does this repository change to read files that way?

## Decision Outcome

This repository already behaved as 1.3.0 states on every point but two: a blank line before the header, and the delimiter sniff on line 1.

- `io/delimited.py`, `sniff_delimiter`: reads the first line that is not blank (`line.strip() == ""` is blank); text that is all blank sniffs a comma.
- `io/hapmap.py`: reads lines until one is not blank; that line is the header, and the `rs#` check and its message are unchanged, so a `#` line before the header is an error. The body already skipped blank lines; body line numbers now count from the header's physical line.
- `io/wide_csv.py`: the header is the first row with any non-blank cell, instead of the first row.
- `io/manifest.py`, `_read_rows` (samples.csv and markers.csv): a `csv.reader` replaces `DictReader`, whose field-name read stops at the first non-empty row even when that row is whitespace-only; the first row with any non-blank cell supplies the field names, and each following row with any non-blank cell is mapped to a dict by those names. Line numbers still come from `reader.line_num`, the row's last physical line.
- `tests/test_contract_cases.py` maps `genotypes.column_count` to `expected \d+ columns, found`, the message the VCF, HapMap and wide-CSV readers already raised.

Unchanged by decision: `#` lines were already data here, VCF and HapMap bodies already skipped blank lines, and quoted line breaks were already read by the `csv` module.

### Consequences

Good: a HapMap, wide CSV, samples.csv or markers.csv with a blank or whitespace-only first line now loads, and a tab-delimited file whose first line is blank is no longer sniffed as comma-delimited. Neutral: which fault a two-fault row reports stays as this repository checks it.

## Amendment, 2026-09-16: review findings (still contract 1.3.0)

Following backcross docs/adr/0018's amendment. `io/delimited.py` gains `is_blank` (empty or only spaces and tabs, so a no-break space is not blank; used by the HapMap, VCF, wide-CSV and manifest readers in place of `str.strip()`) and `csv_rows`, which the wide-CSV and manifest readers now share: a linear pre-scan with the contract's quote automaton raises `DataContractError` naming the physical line where a quoted field left open at the end of the file opened (error kind `delimited.unterminated_quote`), line breaks inside quoted fields are read as a single LF, and rows whose every cell is blank are skipped. The `csv` module's non-strict reading is kept because it already treats a `"` that is not at the start of a field as a literal. `io/vcf.py` raises on a line beginning with a single `#` after the `#CHROM` line, including a second `#CHROM` line, which it previously accepted by replacing the sample list (error kind `genotypes.repeated_header`). `io/manifest.py` raises `DataContractError` naming the line for an unparseable `cm`.
