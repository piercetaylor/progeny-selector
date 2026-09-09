"""Chromosome name normalisation and ordering.

Responsibility: map the chromosome spellings accepted by the data contract
("Gm06", "gm6", "chr6", "Chr06", "6", "06") onto the canonical soybean form "Gm06",
and provide a stable sort key. Names that do not match the soybean pattern are
kept unchanged so the tool still runs on non-soybean data.

Interface:
    normalize_chrom(name: str) -> str
    chrom_sort_key(name: str) -> tuple[int, str]
    chrom_length_bp(name: str, fallback: int | None) -> int | None
"""

from __future__ import annotations

import re

from progeny_selector.constants import SOYBEAN_CHROM_LENGTHS_BP_WM82A4, SOYBEAN_CHROMOSOMES

_SOY_PATTERN = re.compile(r"^(?:gm|chr|chromosome|ch)?_?0*([1-9]|1[0-9]|20)$", re.IGNORECASE)


def normalize_chrom(name: str) -> str:
    """Return the canonical chromosome name ("Gm01".."Gm20") when recognisable, else the input stripped."""
    raw = str(name).strip()
    match = _SOY_PATTERN.match(raw)
    if match:
        return f"Gm{int(match.group(1)):02d}"
    return raw


def chrom_sort_key(name: str) -> tuple[int, str]:
    """Sort soybean chromosomes numerically first, then any other names alphabetically."""
    canonical = normalize_chrom(name)
    if canonical in SOYBEAN_CHROMOSOMES:
        return (SOYBEAN_CHROMOSOMES.index(canonical), "")
    return (len(SOYBEAN_CHROMOSOMES), canonical)


def chrom_length_bp(name: str, fallback: int | None = None) -> int | None:
    """Assembly length for a soybean chromosome (Wm82.a4.v1) or ``fallback`` when unknown."""
    return SOYBEAN_CHROM_LENGTHS_BP_WM82A4.get(normalize_chrom(name), fallback)
