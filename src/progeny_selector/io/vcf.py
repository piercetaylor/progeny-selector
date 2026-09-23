"""VCF 4.2+ reader (plain or gzip/bgzip) producing a GenotypeMatrix.

Responsibility: parse the ``#CHROM`` header for sample ids and the GT field of
every record into allele indices (REF = 0, ALT_k = k, '.' = -1); keep REF/ALT
strings as the per-marker allele table; normalise chromosome names. Only GT is
read; phasing is ignored; haploid GT is duplicated.
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
Format facts: VCF 4.2 specification [web]
https://samtools.github.io/hts-specs/VCFv4.2.pdf.

Interface:
    read_vcf(path: str | Path, scheme: CompiledScheme = SOYBEAN) -> GenotypeMatrix
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from progeny_selector.core.chrom import SOYBEAN, CompiledScheme, normalize_chrom
from progeny_selector.io.delimited import is_blank, open_text
from progeny_selector.io.position import parse_position
from progeny_selector.model.dataset import DataContractError, GenotypeMatrix, Marker


def _parse_gt(token: str) -> tuple[int, int]:
    gt = token.split(":", 1)[0]
    parts = gt.replace("|", "/").split("/")
    if len(parts) == 1:
        parts = [parts[0], parts[0]]
    if len(parts) != 2:
        raise DataContractError(f"non-diploid GT {gt!r}")
    a = -1 if parts[0] == "." else int(parts[0])
    b = -1 if parts[1] == "." else int(parts[1])
    return (a, b)


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
    gt_cache: dict[str, tuple[int, int]] = {}
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
                tokens = [f.split(":")[gt_index] for f in fields[9:]]
            else:
                tokens = fields[9:]
            try:
                pos_bp = parse_position(pos, grammar="digits")
            except ValueError as exc:
                raise DataContractError(f"line {line_no}: {exc}") from exc
            marker_id = mid if mid not in (".", "") else f"{chrom}_{pos_bp}"
            alt_alleles = [] if alt in (".", "") else alt.split(",")
            allele_list = [ref, *alt_alleles]
            gts: list[tuple[int, int]] = []
            for token in tokens:
                # Keyed on the GT sub-field alone, so a GT:DP file with varying depths cannot grow the cache.
                key = token.split(":", 1)[0] if ":" in token else token
                pair = gt_cache.get(key)
                if pair is None:
                    pair = _parse_gt(key)
                    gt_cache[key] = pair
                gts.append(pair)
            if filled >= calls.shape[0]:
                raise DataContractError("VCF changed while it was being read")
            calls[filled] = gts
            if calls[filled].max() >= len(allele_list):
                raise DataContractError(f"line {line_no}: GT allele index exceeds ALT count")
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
    except FileNotFoundError:
        # A missing file fails at open, not while reading, and the CLI catches it by name.
        raise
    except (EOFError, OSError) as exc:
        raise DataContractError(f"{path}: cannot read: {exc}") from exc
    if n_records == 0:
        raise DataContractError("VCF contains no variant records")
    return GenotypeMatrix(markers=markers, sample_ids=sample_ids, alleles=alleles, calls=calls, coded=False)
