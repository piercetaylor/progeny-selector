"""Wide CSV genotype reader (marker_id, chrom, pos_bp, then one column per sample).

Responsibility: parse nucleotide calls ("A", "AT", "A/T", IUPAC codes expanded; any
other cell such as "?", a stray "B" or "H", "X", "0", "+" or "A?" is an error) or A/B/H
coding with auto-detection over every row; a row whose every cell is empty or whitespace
is skipped, an empty marker_id or an invalid pos_bp (position.py, contract 1.2.0) is an
error naming the physical line; the first three columns are marker_id, chrom, pos_bp in
that order (contract 1.1.0); the header is the first row with a non-blank cell, so
blank lines before it are skipped (contract 1.3.0); comma or tab delimited (sniffed from the header),
RFC 4180 quoting, CRLF and a leading BOM accepted; when coded, the parents may be absent
from the file and are synthesised as all-A (recurrent) and all-B (donor) by build_dataset,
which records them in Dataset.synthetic_sample_ids.
Under a token profile (contract 1.4.0) nucleotide cells go through the profile first and a
heterozygote token that names no alleles resolves to the row's two alleles; detection is
skipped (nucleotide) under a base "none" profile, the profile's tokens do not vote under a
base "nucleotide" profile, and a profile on a file whose coding is "abh", requested or detected,
is an error naming the profile and how the coding was reached (contract 1.4.0).

Interface:
    read_wide_csv(path, coding='auto', profile=None, scheme=SOYBEAN) -> GenotypeMatrix
    add_synthetic_parents(gm, rp_id, donor_id) -> GenotypeMatrix   (coded matrices only)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from progeny_selector.core.chrom import SOYBEAN, CompiledScheme, normalize_chrom
from progeny_selector.io.calls import detect_coding, encode_marker, parse_coded_call, parse_nucleotide_call
from progeny_selector.io.delimited import csv_rows, read_text
from progeny_selector.io.position import parse_position
from progeny_selector.io.profiles import TokenProfile, compile_profile
from progeny_selector.model.dataset import DataContractError, GenotypeMatrix, Marker

FIXED = ("marker_id", "chrom", "pos_bp")


def read_wide_csv(
    path: str | Path, coding: str = "auto", profile: TokenProfile | None = None, scheme: CompiledScheme = SOYBEAN
) -> GenotypeMatrix:
    path = Path(path)
    requested = coding
    compiled = compile_profile(profile) if profile is not None else None
    text = read_text(path)
    all_rows = csv_rows(text, "wide CSV")
    header = [h.strip() for h in all_rows[0][1]] if all_rows else []
    if tuple(h.lower() for h in header[:3]) != FIXED:
        raise DataContractError(f"wide CSV must start with columns {FIXED}, found {header[:3]}")
    sample_ids = header[3:]
    if not sample_ids:
        raise DataContractError("wide CSV has no sample columns")
    records = all_rows[1:]
    if not records:
        raise DataContractError("wide CSV has no marker rows")
    if coding == "auto" and compiled is not None and compiled.base == "none":
        coding = "nucleotide"
    elif coding == "auto":
        claimed = set(compiled.missing) | set(compiled.homozygous) | set(compiled.heterozygous) if compiled is not None else set()
        coding = detect_coding(cell for _, row in records for cell in row[3:] if cell.strip().upper() not in claimed)
    if profile is not None and coding == "abh":
        how = "requested" if requested == "abh" else "detected"
        raise DataContractError(
            f'token profile "{profile.id}" applies to HapMap and wide CSV nucleotide calls; the genotype file is coded A/B/H ({how})'
        )
    markers: list[Marker] = []
    alleles: list[list[str]] = []
    rows: list[np.ndarray] = []
    for line_no, row in records:
        if len(row) != len(header):
            raise DataContractError(f"line {line_no}: expected {len(header)} columns, found {len(row)}")
        marker_id = row[0].strip()
        if not marker_id:
            raise DataContractError(f"line {line_no}: empty marker_id")
        try:
            markers.append(Marker(marker_id=marker_id, chrom=normalize_chrom(row[1], scheme), pos_bp=parse_position(row[2])))
            if coding == "abh":
                pairs = np.array([parse_coded_call(c) for c in row[3:]], dtype=np.int8)
                alleles.append(["A", "B"])
            else:
                cells = row[3:]
                allele_list, pairs = encode_marker([parse_nucleotide_call(c, profile=compiled) for c in cells], (), cells)
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
