"""Shared call-string parsing for HapMap and wide-CSV inputs.

Responsibility: turn per-sample call strings into (allele_a, allele_b) pairs,
build a per-marker allele table, and detect A/B/H coding. VCF has its own
allele indices and does not use this module.

Interface:
    parse_nucleotide_call(text) -> tuple[str, str] | None   (None = missing)
    parse_coded_call(text) -> tuple[int, int]               (-1, -1 = missing)
    detect_coding(values: Iterable[str]) -> 'abh' | 'nucleotide'
    encode_marker(calls: list[tuple[str, str] | None]) -> (alleles, int8 array (n, 2))
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from progeny_selector.constants import MISSING_TOKENS

IUPAC_HET: dict[str, tuple[str, str]] = {
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("C", "G"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
}
NUCLEOTIDES = frozenset("ACGT")
CODED_TOKENS = frozenset({"A", "B", "H"})


def parse_nucleotide_call(text: str) -> tuple[str, str] | None:
    """Accept "A", "AA", "AT", "A/T", "A|T", IUPAC het codes; return sorted allele pair or None."""
    t = str(text).strip().upper()
    if t in MISSING_TOKENS:
        return None
    t = t.replace("/", "").replace("|", "")
    if len(t) == 1:
        if t in NUCLEOTIDES:
            return (t, t)
        if t in IUPAC_HET:
            return IUPAC_HET[t]
        return None
    if len(t) == 2 and all(c in NUCLEOTIDES or c in ("N", "-", ".") for c in t):
        if any(c not in NUCLEOTIDES for c in t):
            return None
        return tuple(sorted(t))  # type: ignore[return-value]
    raise ValueError(f"unrecognised nucleotide call {text!r}")


def parse_coded_call(text: str) -> tuple[int, int]:
    """A -> (0,0), B -> (1,1), H -> (0,1), missing -> (-1,-1)."""
    t = str(text).strip().upper()
    if t in MISSING_TOKENS:
        return (-1, -1)
    if t == "A":
        return (0, 0)
    if t == "B":
        return (1, 1)
    if t == "H":
        return (0, 1)
    raise ValueError(f"unrecognised coded call {text!r} (expected A, B, H, N or NA)")


def detect_coding(values: Iterable[str]) -> str:
    """'abh' when every non-missing call is A, B or H and at least one B or H occurs; else 'nucleotide'."""
    seen: set[str] = set()
    for v in values:
        t = str(v).strip().upper()
        if t in MISSING_TOKENS:
            continue
        seen.add(t)
        if len(seen) > 8:
            break
    if seen and seen <= CODED_TOKENS and seen & {"B", "H"}:
        return "abh"
    return "nucleotide"


def encode_marker(calls: list[tuple[str, str] | None]) -> tuple[list[str], np.ndarray]:
    """Allele table (sorted alphabetically) and int8 (n_samples, 2) indices for one marker."""
    alleles = sorted({a for c in calls if c is not None for a in c})
    lookup = {a: i for i, a in enumerate(alleles)}
    out = np.full((len(calls), 2), -1, dtype=np.int8)
    for i, c in enumerate(calls):
        if c is not None:
            out[i, 0] = lookup[c[0]]
            out[i, 1] = lookup[c[1]]
    return alleles, out
