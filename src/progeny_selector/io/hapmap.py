"""HapMap (TASSEL-style) reader producing a GenotypeMatrix.

Responsibility: parse the 11 fixed columns (rs#, alleles, chrom, pos, strand,
assembly#, center, protLSID, assayLSID, panelLSID, QCcode) followed by one
column per sample; calls go through calls.parse_nucleotide_call with the HapMap
missing list "", N, NN, NA, -, --, ., ./., .|., X, XX (contract/data-contract.md
1.1.0, "HapMap"); IUPAC codes expand; any other cell is an error naming the line.
`pos` goes through position.py (contract 1.2.0). Blank lines are skipped
wherever they occur; the header is the first line that is not blank and must
start with rs#, so a `#` line before it is an error (contract 1.3.0).
Under a token profile (contract 1.4.0) cells go through the profile first; a
heterozygote token that names no alleles resolves to the row's two alleles (the
`alleles` column upper-cased, minus N, plus the symbols the cells show), and any other number of
alleles is an error naming the line.
Column facts [web]
https://statgen-esalq.github.io/Hapmap-and-VCF-formats-and-its-integration-with-onemap/.

Interface:
    read_hapmap(path: str | Path, profile: TokenProfile | None = None) -> GenotypeMatrix
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from progeny_selector.constants import HAPMAP_MISSING
from progeny_selector.core.chrom import normalize_chrom
from progeny_selector.io.calls import encode_marker, parse_nucleotide_call
from progeny_selector.io.delimited import is_blank, open_text
from progeny_selector.io.position import parse_position
from progeny_selector.io.profiles import TokenProfile, compile_profile
from progeny_selector.model.dataset import DataContractError, GenotypeMatrix, Marker

N_FIXED = 11


def read_hapmap(path: str | Path, profile: TokenProfile | None = None) -> GenotypeMatrix:
    path = Path(path)
    compiled = compile_profile(profile) if profile is not None else None
    markers: list[Marker] = []
    alleles: list[list[str]] = []
    rows: list[np.ndarray] = []
    with open_text(path) as fh:
        header_line = fh.readline()
        header_no = 1
        while header_line and is_blank(header_line):
            header_line = fh.readline()
            header_no += 1
        header = header_line.rstrip("\r\n").split("\t")
        if len(header) <= N_FIXED or not header[0].lower().startswith("rs#"):
            raise DataContractError("HapMap header must start with rs# and have 11 fixed columns plus samples")
        sample_ids = header[N_FIXED:]
        for line_no, line in enumerate(fh, start=header_no + 1):
            if is_blank(line):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(header):
                raise DataContractError(f"line {line_no}: expected {len(header)} columns, found {len(fields)}")
            try:
                cells = fields[N_FIXED:]
                calls = [parse_nucleotide_call(t, HAPMAP_MISSING, compiled) for t in cells]
                seed = [a.strip().upper() for a in fields[1].split("/") if a.strip() and a.strip().upper() != "N"]
                allele_list, pairs = encode_marker(calls, seed, cells)
            except ValueError as exc:
                raise DataContractError(f"line {line_no}: {exc}") from exc
            try:
                pos_bp = parse_position(fields[3])
            except ValueError as exc:
                raise DataContractError(f"line {line_no}: {exc}") from exc
            markers.append(Marker(marker_id=fields[0], chrom=normalize_chrom(fields[2]), pos_bp=pos_bp))
            alleles.append(allele_list)
            rows.append(pairs)
    if not markers:
        raise DataContractError("HapMap contains no marker rows")
    return GenotypeMatrix(markers=markers, sample_ids=sample_ids, alleles=alleles, calls=np.stack(rows), coded=False)
