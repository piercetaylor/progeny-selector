"""VCF 4.2+ reader (plain or gzip/bgzip) producing a GenotypeMatrix.

Responsibility: parse the ``#CHROM`` header for sample ids and the GT field of
every record into allele indices (REF = 0, ALT_k = k, '.' = -1); keep REF/ALT
strings as the per-marker allele table; normalise chromosome names. Only GT is
read; phasing is ignored; haploid GT is duplicated.
The GT grammar is contract 1.10.0: `.` or an allele index, or two of them joined by
`/` or `|`; an index is `0` or `[1-9][0-9]*`, at most 127 and less than the REF,ALT
allele count. An empty GT is missing. Anything else is an `invalid GT` error naming
the line and the value (docs/adr/0030).
Records with ID "." or empty get "<CHROM>_<POS>" from CHROM as written in the file, before normalisation (contract 1.1.0).
POS goes through position.py with the "digits" grammar, and an ID-less record is named from the parsed POS (contract 1.2.0).
The file is read twice (docs/adr/0022): pass 1 (`_scan`) takes the sample ids and
counts data lines, judging nothing, so every message and physical line number is
pass 2's as it always was; pass 2 (`_fill`) writes each record into a matrix
preallocated from that count; the emptiness of a VCF is judged after pass 2, so a file
whose records precede its header is told that rather than called empty. The peak is therefore the genotype matrix plus one
line's tokens plus the marker and allele tables, instead of about twice the
matrix, and the cost is one extra decompression of a `.gz`. Chromosome names are
normalised with `scheme` in pass 2 only; pass 1 tokenises nothing.
A byte sequence that is not valid UTF-8 in either pass raises DataContractError naming the
file and the physical line (contract 1.11.0, docs/adr/0031); pass 1 reads the whole file, so it
can outrank a row error on an earlier line.
Format facts: VCF 4.2 specification [web]
https://samtools.github.io/hts-specs/VCFv4.2.pdf.

Interface:
    read_vcf(path: str | Path, scheme: CompiledScheme = SOYBEAN) -> GenotypeMatrix
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from progeny_selector.core.chrom import SOYBEAN, CompiledScheme, normalize_chrom
from progeny_selector.io.delimited import invalid_utf8_error, is_blank, open_text
from progeny_selector.io.position import parse_position
from progeny_selector.model.dataset import DataContractError, GenotypeMatrix, Marker

# Contract 1.10.0 GT grammar, the same pattern as backcross's GT_PATTERN. `[0-9]`, never `\d`,
# which matches Unicode digits in Python; applied with fullmatch, because `$` accepts a trailing newline.
_GT_RE = re.compile(r"^(\.|0|[1-9][0-9]*)(?:[/|](\.|0|[1-9][0-9]*))?$")

# Highest allele index accepted (contract 1.10.0), so an index always fits the int8 store.
MAX_ALLELE_INDEX = 127


def _index(text: str) -> int:
    """One allele index already matched by `_GT_RE`; '.' is -1.

    More than three digits without a leading zero is at least 1000, so it is returned as
    MAX_ALLELE_INDEX + 1 rather than converted: `int()` raises ValueError past 4300 digits
    (Python 3.11+), which must reach the caller as the range error, not as a raw exception.
    """
    if text == ".":
        return -1
    return MAX_ALLELE_INDEX + 1 if len(text) > 3 else int(text)


def _parse_gt(gt: str, line_no: int) -> tuple[int, int, int]:
    """One GT value as (a, b, max(a, b)); '.' is -1, a haploid call is homozygous, an empty GT is missing.

    Judges the grammar only. The range against the record's alleles is judged per record by the
    caller, because a cached pair is shared by records with different ALT counts.
    """
    if gt == "":
        return (-1, -1, -1)
    m = _GT_RE.fullmatch(gt)
    if m is None:
        raise DataContractError(f'line {line_no}: invalid GT "{gt}"')
    a = _index(m[1])
    b = a if m[2] is None else _index(m[2])
    return (a, b, max(a, b))


def _range_error(gt: str, pair: tuple[int, int, int], n_alleles: int, line_no: int) -> DataContractError:
    """The message for the first side of `pair` out of range, in backcross's order (first side, then second)."""
    for i in pair[:2]:
        if i > MAX_ALLELE_INDEX:
            return DataContractError(f'line {line_no}: invalid GT "{gt}", allele index above {MAX_ALLELE_INDEX}')
        if i >= n_alleles:
            return DataContractError(f'line {line_no}: invalid GT "{gt}", allele index {i} but the record has {n_alleles} alleles')
    raise AssertionError("no allele index out of range")


def _scan(path: Path) -> tuple[list[str], int]:
    """Pass 1: the sample ids and the number of data lines, judging nothing at all.

    The header's shape, a line beginning with '#' after the header, a data line before
    it and a second #CHROM line are all left to pass 2, which raises the existing
    message at the same physical line, so precedence between errors inside the file is
    unchanged. `seen_header` is an explicit flag, never the truthiness of `sample_ids`:
    a malformed #CHROM line yields no sample ids, and counting records only while a
    non-empty list was in hand would silently shorten the matrix of a file that must
    instead be refused in pass 2.
    """
    sample_ids: list[str] = []
    seen_header = False
    n_records = 0
    with open_text(path) as fh:
        for _line_no, line in enumerate(fh, start=1):
            if is_blank(line):
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                if not seen_header and line.startswith("#CHROM"):
                    seen_header = True
                    sample_ids = line.rstrip("\r\n").split("\t")[9:]
                continue
            if not seen_header:
                continue
            n_records += 1
    return sample_ids, n_records


def _fill(path: Path, sample_ids: list[str], calls: np.ndarray, scheme: CompiledScheme) -> tuple[list[Marker], list[list[str]]]:
    """Pass 2: every validation of the single-pass reader, writing record i into `calls[i]`."""
    markers: list[Marker] = []
    alleles: list[list[str]] = []
    gt_cache: dict[str, tuple[int, int, int]] = {}
    seen_header = False
    filled = 0
    with open_text(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            if is_blank(line):
                continue
            if line.startswith("##"):
                continue
            if seen_header and line.startswith("#"):
                # A second #CHROM line or any other single-# line after the header (contract 1.3.0).
                raise DataContractError(f"line {line_no}: line beginning with '#' after the #CHROM header")
            if line.startswith("#CHROM"):
                # The shape is judged here, where the single-pass reader judged it, so an error
                # on an earlier line still outranks a malformed header.
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 10 or fields[8] != "FORMAT":
                    raise DataContractError("VCF header must have FORMAT and at least one sample column")
                seen_header = True
                continue
            if not seen_header:
                raise DataContractError("VCF data line before #CHROM header")
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9 + len(sample_ids):
                raise DataContractError(f"line {line_no}: expected {9 + len(sample_ids)} columns, found {len(fields)}")
            chrom, pos, mid, ref, alt, _qual, _filt, _info, fmt = fields[:9]
            fmt_keys = fmt.split(":")
            if fmt_keys[0] != "GT":
                if "GT" not in fmt_keys:
                    raise DataContractError(f"line {line_no}: FORMAT has no GT")
                gt_index = fmt_keys.index("GT")
                # A sample field that ends before its GT sub-field (VCF lets trailing sub-fields be
                # dropped) has a missing GT, read as '.' as backcross does (contract 1.10.0).
                tokens = []
                for f in fields[9:]:
                    subs = f.split(":")
                    tokens.append(subs[gt_index] if gt_index < len(subs) else ".")
            else:
                tokens = fields[9:]
            try:
                pos_bp = parse_position(pos, grammar="digits")
            except ValueError as exc:
                raise DataContractError(f"line {line_no}: {exc}") from exc
            marker_id = mid if mid not in (".", "") else f"{chrom}_{pos_bp}"
            alt_alleles = [] if alt in (".", "") else alt.split(",")
            allele_list = [ref, *alt_alleles]
            n_alleles = len(allele_list)
            # Every index must be below both the allele count and 128 (contract 1.10.0); checked
            # before the int8 store, so no index can overflow it.
            limit = min(n_alleles, MAX_ALLELE_INDEX + 1)
            gts: list[tuple[int, int]] = []
            for token in tokens:
                # Keyed on the GT sub-field alone, so a GT:DP file with varying depths cannot grow the cache.
                key = token.split(":", 1)[0] if ":" in token else token
                pair = gt_cache.get(key)
                if pair is None:
                    pair = _parse_gt(key, line_no)
                    gt_cache[key] = pair
                if pair[2] >= limit:
                    raise _range_error(key, pair, n_alleles, line_no)
                gts.append((pair[0], pair[1]))
            if filled >= calls.shape[0]:
                raise DataContractError("VCF changed while it was being read")
            calls[filled] = gts
            markers.append(Marker(marker_id=marker_id, chrom=normalize_chrom(chrom, scheme), pos_bp=pos_bp))
            alleles.append(allele_list)
            filled += 1
    if filled != calls.shape[0]:
        raise DataContractError("VCF changed while it was being read")
    return markers, alleles


def read_vcf(path: str | Path, scheme: CompiledScheme = SOYBEAN) -> GenotypeMatrix:
    path = Path(path)
    try:
        sample_ids, n_records = _scan(path)
        # A zero-row allocation, so pass 2 still runs on a file with no countable records and
        # raises the structural message its first line deserves ("VCF data line before #CHROM
        # header", or the '#'-after-header message) instead of the false "no variant records".
        calls = np.empty((n_records, len(sample_ids), 2), dtype=np.int8)
        markers, alleles = _fill(path, sample_ids, calls, scheme)
    except UnicodeDecodeError:
        # Contract 1.11.0: strict UTF-8, located by a byte re-read (delimited.invalid_utf8_message).
        raise invalid_utf8_error(path) from None
    except FileNotFoundError:
        # A missing file fails at open, not while reading, and the CLI catches it by name.
        raise
    except (EOFError, OSError) as exc:
        raise DataContractError(f"{path}: cannot read: {exc}") from exc
    if n_records == 0:
        raise DataContractError("VCF contains no variant records")
    return GenotypeMatrix(markers=markers, sample_ids=sample_ids, alleles=alleles, calls=calls, coded=False)
