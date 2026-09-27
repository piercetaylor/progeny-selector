# Contract 1.11.0: text inputs are UTF-8

Status: accepted. Date: 2026-09-27. Mirrors backcross docs/adr/0026, the canonical record for this contract version.

## Context and Problem Statement

Contract 1.0.0 never stated an encoding. Contract 1.11.0 states that every text input is UTF-8, as VCF 4.3 section 1.2 states for VCF, and adds the error kind `text.invalid_utf8`: a byte sequence that is not well-formed UTF-8 (Unicode Table 3-7) is an error naming the file and the physical line, in the genotype file (plain, gzip or bgzip), samples.csv and markers.csv. This record is the sibling half: what this repository did before, what changed, and what did not. The decision, the options weighed and the evidence are in backcross docs/adr/0026 and are not restated here.

Measured here before the change, each file carrying the Latin-1 byte `0xE9`:

| Input                         | progeny-selector before 1.11.0                                                 |
| ----------------------------- | ------------------------------------------------------------------------------ |
| VCF (plain, `.gz`, `.bgz`)    | raw `UnicodeDecodeError` from `open_text`, with a stream offset and no line   |
| HapMap                        | raw `UnicodeDecodeError`, no line                                              |
| wide CSV                      | raw `UnicodeDecodeError` from `read_text`, no line                             |
| samples.csv                   | raw `UnicodeDecodeError` from `read_text`, no line                             |
| markers.csv                   | raw `UnicodeDecodeError` from `read_text`, no line                             |
| token profile JSON            | already `DataContractError` (`except (OSError, ValueError)`; `UnicodeDecodeError` is a `ValueError`) |
| criteria.yaml                 | raw `UnicodeDecodeError` from `read_criteria`                                  |
| results.csv (`select`)        | raw `UnicodeDecodeError` from `cli._read_results`                              |

## Decision Outcome

The canonical record's decision applies here unchanged.

What changed, all in `src/progeny_selector/io`:

- **`delimited.py`** gains `open_binary(path)` (`gzip.open(path, "rb")` for `.gz` and `.bgz`, `open(path, "rb")` otherwise), `invalid_utf8_message(path)` and `invalid_utf8_error(path)`. `invalid_utf8_message` re-reads the file as bytes, line by line, and returns `{path}: line N: not valid UTF-8 (byte 0xHH at position P)` for the first line whose bytes fail a strict `utf-8` decode, where `P` is `UnicodeDecodeError.start + 1`, the 1-based byte offset of the ill-formed sequence's lead byte within that line. It returns `{path}: not valid UTF-8` when no line fails or the re-read itself raises `EOFError` or `OSError`. `read_text` turns a `UnicodeDecodeError` into that `DataContractError`, which covers the wide CSV, samples.csv and markers.csv readers.
- **`vcf.py`** `read_vcf` turns a `UnicodeDecodeError` from either pass into the same error; **`hapmap.py`** does the same around its `open_text` block.
- **`criteria.py`** `read_criteria` turns a `UnicodeDecodeError` into `CriteriaError` with the same message. criteria.yaml is outside the contract, so no case covers it.
- **`cli.py`** `_read_results`, the results.csv reader of `select`, turns a `UnicodeDecodeError` into `DataContractError` with the same message, which `main` reports as an error with exit status 1. results.csv is outside the contract, so no case covers it; `tests/test_utf8.py` pins it through `select`.
- **The custom token profile JSON** (`cli._profile_ref`, `app/screens/load.read_profile_json`) was already a `DataContractError`, but it carried the codec's stream offset; a `UnicodeDecodeError` is now caught before the generic `ValueError` and reported as `token profile file <path>: line N: not valid UTF-8 (byte 0xHH at position P)` on the command line (`token profile file: line N: ...` on the Load screen, whose upload path is a temporary file), the line and position backcross gives its `token profile`. `invalid_utf8_message` takes an optional `label` for this. Outside the contract (D8), so no case; `tests/test_utf8.py` pins the command-line path.
- **`tests/test_contract_cases.py`** maps `text.invalid_utf8` to `not valid UTF-8`, and `tests/test_utf8.py` pins every reader's line and position, the three VCF compressions, the overlong, surrogate and truncated sequences, and the positives.

What did not change:

- The decode is still a strict `TextIOWrapper` with `utf-8-sig`, so the happy path does no extra work; the byte re-read happens only after a `UnicodeDecodeError`.
- The byte-order mark stays ignored (contract 1.1.0). On line 1 the re-read decodes a BOM as U+FEFF, so its three bytes count in the position, as the contract message states.
- The strictness table is CPython's strict `utf-8` codec, which agrees with Unicode Table 3-7 and with backcross's `TextDecoder {fatal: true}`: an overlong encoding, an encoded surrogate, a code point above U+10FFFF, a lone continuation byte and a sequence cut short by the end of the file are all errors. A well-formed character split across a gzip member boundary is read whole, because the incremental decoder sees one stream.

**Precedence (backcross docs/adr/0026 D6).** Decoding is judged per `TextIOWrapper` buffer, and `read_vcf`'s pass 1 reads the whole file (docs/adr/0022), so a UTF-8 error can outrank a row error on an earlier line. Contract 1.3.0 D4.4 allows this, and no case carries two faults. docs/adr/0022 carries a dated amendment saying so.

**The message prefix is `{path}: line N:`**, as the `cannot read` message names the path; backcross labels the file (`genotype file`, `samples.csv`, `markers.csv`). The position and label are tool wording, not contract vocabulary; the stable substring is `not valid UTF-8`.

## Consequences

No well-formed input changes meaning, no output changes, and `results_schema` does not move. One error kind is added. A file carrying a byte sequence that is not UTF-8, such as a Windows-1252 export of a `notes` column, now fails with a contract error naming the file, line and byte, where before it failed with a raw Python exception; the fix is to save it as "CSV UTF-8".

## Revisit when

- backcross docs/adr/0026 is amended or superseded.
