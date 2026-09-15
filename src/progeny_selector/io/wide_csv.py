"""Wide CSV genotype reader (marker_id, chrom, pos_bp, then one column per sample).

Responsibility: parse nucleotide calls ("A", "AT", "A/T", IUPAC codes expanded; any
other cell such as "?", a stray "B" or "H", "X", "0", "+" or "A?" is an error) or A/B/H
coding with auto-detection over every row; the first three columns are marker_id, chrom,
pos_bp in that order (contract 1.1.0); comma or tab delimited (sniffed from the header),
RFC 4180 quoting, CRLF and a leading BOM accepted; when coded, the parents may be absent
from the file and are synthesised as all-A (recurrent) and all-B (donor) by build_dataset,
which records them in Dataset.synthetic_sample_ids.

Interface:
    read_wide_csv(path, coding='auto') -> GenotypeMatrix
    add_synthetic_parents(gm, rp_id, donor_id) -> GenotypeMatrix   (coded matrices only)
"""

from __future__ import annotations

import csv
import io
from pathlib import Path

import numpy as np

from progeny_selector.core.chrom import normalize_chrom
from progeny_selector.io.calls import detect_coding, encode_marker, parse_coded_call, parse_nucleotide_call
from progeny_selector.io.delimited import read_text, sniff_delimiter
from progeny_selector.model.dataset import DataContractError, GenotypeMatrix, Marker

FIXED = ("marker_id", "chrom", "pos_bp")


def read_wide_csv(path: str | Path, coding: str = "auto") -> GenotypeMatrix:
    path = Path(path)
    text = read_text(path)
    reader = csv.reader(io.StringIO(text, newline=""), delimiter=sniff_delimiter(text))
    header = [h.strip() for h in next(reader, [])]
    if tuple(h.lower() for h in header[:3]) != FIXED:
        raise DataContractError(f"wide CSV must start with columns {FIXED}, found {header[:3]}")
    sample_ids = header[3:]
    if not sample_ids:
        raise DataContractError("wide CSV has no sample columns")
    records = [row for row in reader if any(cell.strip() for cell in row)]
    if not records:
        raise DataContractError("wide CSV has no marker rows")
    if coding == "auto":
        coding = detect_coding(cell for row in records for cell in row[3:])
    markers: list[Marker] = []
    alleles: list[list[str]] = []
    rows: list[np.ndarray] = []
    for line_no, row in enumerate(records, start=2):
        if len(row) != len(header):
            raise DataContractError(f"line {line_no}: expected {len(header)} columns, found {len(row)}")
        try:
            markers.append(Marker(marker_id=row[0].strip(), chrom=normalize_chrom(row[1]), pos_bp=int(float(row[2]))))
            if coding == "abh":
                pairs = np.array([parse_coded_call(c) for c in row[3:]], dtype=np.int8)
                alleles.append(["A", "B"])
            else:
                allele_list, pairs = encode_marker([parse_nucleotide_call(c) for c in row[3:]])
                alleles.append(allele_list)
        except ValueError as exc:
            raise DataContractError(f"line {line_no}: {exc}") from exc
        rows.append(pairs)
    return GenotypeMatrix(markers=markers, sample_ids=sample_ids, alleles=alleles, calls=np.stack(rows), coded=(coding == "abh"))


def add_synthetic_parents(gm: GenotypeMatrix, rp_id: str, donor_id: str) -> GenotypeMatrix:
    """For A/B/H-coded input without parent columns, append RP (all A) and donor (all B)."""
    if not gm.coded:
        raise DataContractError("synthetic parents are only valid for A/B/H-coded input")
    new_ids = [s for s in (rp_id, donor_id) if s not in gm.sample_ids]
    if not new_ids:
        return gm
    cols = []
    for s in new_ids:
        value = 0 if s == rp_id else 1
        cols.append(np.full((gm.n_markers, 1, 2), value, dtype=np.int8))
    return gm.with_samples_added(new_ids, np.concatenate(cols, axis=1))
