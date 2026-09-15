"""HapMap (TASSEL-style) reader producing a GenotypeMatrix.

Responsibility: parse the 11 fixed columns (rs#, alleles, chrom, pos, strand,
assembly#, center, protLSID, assayLSID, panelLSID, QCcode) followed by one
column per sample; calls go through calls.parse_nucleotide_call with the HapMap
missing list "", N, NN, NA, -, --, ., ./., .|., X, XX (contract/data-contract.md
1.1.0, "HapMap"); IUPAC codes expand; any other cell is an error naming the line.
Column facts [web]
https://statgen-esalq.github.io/Hapmap-and-VCF-formats-and-its-integration-with-onemap/.

Interface:
    read_hapmap(path: str | Path) -> GenotypeMatrix
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from progeny_selector.constants import HAPMAP_MISSING
from progeny_selector.core.chrom import normalize_chrom
from progeny_selector.io.calls import encode_marker, parse_nucleotide_call
from progeny_selector.io.delimited import open_text
from progeny_selector.model.dataset import DataContractError, GenotypeMatrix, Marker

N_FIXED = 11


def read_hapmap(path: str | Path) -> GenotypeMatrix:
    path = Path(path)
    markers: list[Marker] = []
    alleles: list[list[str]] = []
    rows: list[np.ndarray] = []
    with open_text(path) as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        if len(header) <= N_FIXED or not header[0].lower().startswith("rs#"):
            raise DataContractError("HapMap header must start with rs# and have 11 fixed columns plus samples")
        sample_ids = header[N_FIXED:]
        for line_no, line in enumerate(fh, start=2):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(header):
                raise DataContractError(f"line {line_no}: expected {len(header)} columns, found {len(fields)}")
            try:
                calls = [parse_nucleotide_call(t, HAPMAP_MISSING) for t in fields[N_FIXED:]]
            except ValueError as exc:
                raise DataContractError(f"line {line_no}: {exc}") from exc
            allele_list, pairs = encode_marker(calls)
            try:
                pos_bp = int(fields[3])
            except ValueError as exc:
                raise DataContractError(f"line {line_no}: invalid position {fields[3]!r}") from exc
            markers.append(Marker(marker_id=fields[0], chrom=normalize_chrom(fields[2]), pos_bp=pos_bp))
            alleles.append(allele_list)
            rows.append(pairs)
    if not markers:
        raise DataContractError("HapMap contains no marker rows")
    return GenotypeMatrix(markers=markers, sample_ids=sample_ids, alleles=alleles, calls=np.stack(rows), coded=False)
