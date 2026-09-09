"""VCF 4.2+ reader (plain or gzip/bgzip) producing a GenotypeMatrix.

Responsibility: parse the ``#CHROM`` header for sample ids and the GT field of
every record into allele indices (REF = 0, ALT_k = k, '.' = -1); keep REF/ALT
strings as the per-marker allele table; normalise chromosome names. Only GT is
read; phasing is ignored; haploid GT is duplicated. Records without an ID get
``<chrom>_<pos>``. Format facts: VCF 4.2 specification [web]
https://samtools.github.io/hts-specs/VCFv4.2.pdf.

Interface:
    read_vcf(path: str | Path) -> GenotypeMatrix
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np

from progeny_selector.core.chrom import normalize_chrom
from progeny_selector.model.dataset import DataContractError, GenotypeMatrix, Marker


def _open_text(path: Path):
    if str(path).endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


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


def read_vcf(path: str | Path) -> GenotypeMatrix:
    path = Path(path)
    sample_ids: list[str] = []
    markers: list[Marker] = []
    alleles: list[list[str]] = []
    rows: list[np.ndarray] = []
    with _open_text(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            if not line.strip():
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10 or fields[8] != "FORMAT":
                    raise DataContractError("VCF header must have FORMAT and at least one sample column")
                sample_ids = fields[9:]
                continue
            if not sample_ids:
                raise DataContractError("VCF data line before #CHROM header")
            fields = line.rstrip("\n").split("\t")
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
            marker_id = mid if mid not in (".", "") else f"{normalize_chrom(chrom)}_{pos}"
            alt_alleles = [] if alt in (".", "") else alt.split(",")
            allele_list = [ref, *alt_alleles]
            pairs = np.array([_parse_gt(t) for t in tokens], dtype=np.int8)
            if pairs.max() >= len(allele_list):
                raise DataContractError(f"line {line_no}: GT allele index exceeds ALT count")
            markers.append(Marker(marker_id=marker_id, chrom=normalize_chrom(chrom), pos_bp=int(pos)))
            alleles.append(allele_list)
            rows.append(pairs)
    if not markers:
        raise DataContractError("VCF contains no variant records")
    calls = np.stack(rows, axis=0)
    return GenotypeMatrix(markers=markers, sample_ids=sample_ids, alleles=alleles, calls=calls, coded=False)
