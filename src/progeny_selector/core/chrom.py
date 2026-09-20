"""Chromosome name normalisation and ordering.

Responsibility: map the chromosome spellings accepted by the contract (prefix Gm, Chr,
Chromosome or LG, optional _, space or - separator, 1..20 with leading zeros;
contract/data-contract.md 1.4.0) onto the canonical "Gm06", keep other names unchanged,
and provide a sort key that puts Gm01..Gm20 first and everything else after them in
natural (numeric-aware) order, the same order as backcross's compareChromosomes.

Interface:
    normalize_chrom(name: str) -> str
    chrom_sort_key(name: str) -> tuple[int, tuple[tuple[int, str], ...]]
    chrom_length_bp(name: str, fallback: int | None, assembly: str) -> int | None
"""

from __future__ import annotations

import re

from progeny_selector.constants import DEFAULT_ASSEMBLY, SOYBEAN_CHROM_LENGTHS_BP, SOYBEAN_CHROMOSOMES

_SOY_PATTERN = re.compile(r"^(?:gm|chr|chromosome|lg)?[_\s-]?0*([1-9]|1[0-9]|20)$", re.IGNORECASE)
_DIGIT_RUN = re.compile(r"([0-9]+)")


def normalize_chrom(name: str) -> str:
    """Return the canonical chromosome name ("Gm01".."Gm20") when recognisable, else the input stripped."""
    raw = str(name).strip()
    match = _SOY_PATTERN.match(raw)
    if match:
        return f"Gm{int(match.group(1)):02d}"
    return raw


def _natural_key(name: str) -> tuple[tuple[int, str], ...]:
    """Split on digit runs: even parts are text, odd parts are numbers, so aligned parts share a type."""
    parts = _DIGIT_RUN.split(name)
    return tuple((int(p), "") if i % 2 else (0, p) for i, p in enumerate(parts))


def chrom_sort_key(name: str) -> tuple[int, tuple[tuple[int, str], ...]]:
    """Gm01..Gm20 numerically first, then any other name in natural order (scaffold_2 before scaffold_10)."""
    canonical = normalize_chrom(name)
    if canonical in SOYBEAN_CHROMOSOMES:
        return (SOYBEAN_CHROMOSOMES.index(canonical), ())
    return (len(SOYBEAN_CHROMOSOMES), _natural_key(canonical))


def chrom_length_bp(name: str, fallback: int | None = None, assembly: str = DEFAULT_ASSEMBLY) -> int | None:
    """Chromosome length in the named assembly, or ``fallback`` when the chromosome or the assembly is unknown.

    ``assembly`` is the criteria.yaml key (contract 1.4.0 leaves chromosome lengths to the tool); "none"
    and any name without a table return ``fallback``.
    """
    return SOYBEAN_CHROM_LENGTHS_BP.get(assembly, {}).get(normalize_chrom(name), fallback)
