# Contract 1.10.0: the VCF GT grammar

Status: accepted. Date: 2026-09-27. Mirrors backcross docs/adr/0025, the canonical record for this contract version.

## Context and Problem Statement

The contract's VCF section never stated the grammar of a GT allele index, and `_parse_gt` (`src/progeny_selector/io/vcf.py`) read each side with Python `int()` into an int8 matrix. Contract 1.10.0 states the grammar and adds the error kind `genotypes.invalid_gt`. This record is the sibling half: what this repository did before, what changed, and what did not. The decision, the options weighed, the VCF 4.2 and htslib evidence and the criterion ("never silently wrong, tolerant of standard-tool exports") are in backcross docs/adr/0025 and are not restated here.

Measured here before the change, one sample and one ALT unless stated:

| GT                     | progeny-selector before 1.10.0                                                 |
| ---------------------- | ------------------------------------------------------------------------------ |
| `-5/0`                 | int8 -5 stored; the ALT-count check tested only the maximum, so it passed      |
| `1e1/0`                | raw `ValueError`                                                               |
| `/0`, `0/`             | raw `ValueError`                                                               |
| empty GT               | raw `ValueError`                                                               |
| `+1/0`, ` 0/1`, `01/0` | read as 1/0, 0/1, 1/0                                                          |
| `1_0/0`                | 10/0 (`int()` accepts `_`), given 11 or more alleles                             |
| `0/0/1`, `\|0\|1`      | `non-diploid GT` error, without a line number                                  |
| `2/0`                  | `line N: GT allele index exceeds ALT count`                                    |
| `128/0`                | raw `OverflowError` from numpy, before the ALT-count check could run           |

## Decision Outcome

The canonical record's decision applies here unchanged: the GT sub-field is `.` or an allele index, or two of them separated by `/` or `|`; an index is `0` or `[1-9][0-9]*`, at most 127, and less than the number of alleles in REF and ALT together; anything else is `genotypes.invalid_gt`.

- **The same regular expression.** `_GT_RE` is backcross's `GT_PATTERN` character for character, `^(\.|0|[1-9][0-9]*)(?:[/|](\.|0|[1-9][0-9]*))?$`, applied with `re.fullmatch` so that a trailing newline cannot pass `$`, and written with `[0-9]` because Python's `\d` matches Unicode digits.
- **The same messages.** A grammar failure is `line N: invalid GT "<value>"`; a range failure adds `, allele index above 127` or `, allele index I but the record has K alleles`, judged first side then second, as backcross does. The prefix is this reader's `line N:`, as in its `invalid position` message, where backcross writes `VCF line N:`. `tests/test_contract_cases.py` maps the kind to `invalid GT`.
- **The range is judged per record, before the int8 store.** `_parse_gt` judges the grammar only and its result is cached by GT text; the reader compares each cached pair's larger index with the smaller of the record's allele count and 128 before writing the row. A cached `2/2` from a record with ALT `G,T` is therefore still refused on a record with one ALT, and no index can reach the int8 store out of range, so the `OverflowError` at 128 and the stored -5 are both gone. The old post-store check (`GT allele index exceeds ALT count`) is removed; its case is now the `allele index I but the record has K alleles` message.
- **An empty GT is read as missing.** It raised a raw `ValueError` here; backcross already read it as missing. It names no allele, so reading it as missing cannot be silently wrong, and it is pinned in `vcf-gt-haploid-and-missing`.
- **A sample field that ends before its GT sub-field is read as missing.** VCF lets a sample drop trailing sub-fields, so with FORMAT `DP:GT` a bare `12` has no GT. The reader indexed the split field and raised a raw `IndexError`; it now reads the missing sub-field as `.`, as backcross does (`?? '.'`). Pinned by record `h7` of `vcf-gt-haploid-and-missing` and by `tests/test_io.py`.
- **Kept as before:** `./.`, `.|.` and `.` are missing; a half-missing `./1` is stored as (-1, 1) and read as missing downstream (contract 1.6.0); haploid `1` is read as 1/1; phased `1|0` is the pair 1/0, unordered downstream; a multiallelic index such as `2/2` with ALT `G,T` is read. The GT is still the text before the first `:` when GT leads FORMAT, and the `GT` sub-field otherwise.
- **The rule is pinned by hand-built tests.** `tests/test_io.py` covers every row above and more spellings (a trailing space, `/`, `0//1`, a backslash, an Arabic-Indic digit, two out-of-range sides in both orders), the 127 and 128 boundary with 128 ALT alleles, the cache case, and every reading kept. The shared cases `err-vcf-gt-negative`, `-exponent`, `-underscore`, `-empty-side`, `-leading-zero`, `-triploid`, `-exceeds-alt`, `-leading-phase`, `vcf-gt-haploid-and-missing` and `vcf-gt-multiallelic-index` are generated in backcross and run here through `tests/test_contract_cases.py`; `scripts/check_contract.py ../backcross` byte-compares the mirror.
- **`docs/data-formats.md`** moves its version reference to 1.10.0.

## Consequences

No stated input changes meaning, no output changes, and `results_schema` does not move. One error kind is added. A file carrying a sign, an exponent, an underscore, whitespace, a leading zero, an empty side, three alleles or the VCF 4.4 leading phase indicator in GT now fails with a contract error naming the line and the value, where before it was read, or failed with a raw Python exception. A file with an empty GT now loads, with that call missing.

## Revisit when

- The deferred item of backcross docs/adr/0025 is taken up: reading the VCF 4.4 leading phase indicator when the file declares `##fileformat=VCFv4.4`. This reader does not gate on the version either.
- backcross docs/adr/0025 is amended or superseded.
