"""Shared call-string parsing for HapMap and wide-CSV inputs (contract/data-contract.md 1.1.0; mirrors isoline-browser src/io/calls.ts).

Responsibility: turn per-sample call strings into (allele_a, allele_b) pairs,
build a per-marker allele table, and detect A/B/H coding. VCF has its own
allele indices and does not use this module.

Interface:
    parse_nucleotide_call(text, missing=WIDE_NUCLEOTIDE_MISSING) -> tuple[str, str] | None
        (None = missing; ValueError for any other cell, including "?")
    parse_coded_call(text) -> tuple[int, int]               (-1, -1 = missing)
    detect_coding(values: Iterable[str]) -> 'abh' | 'nucleotide'
    encode_marker(calls: list[tuple[str, str] | None]) -> (alleles, int8 array (n, 2))
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from progeny_selector.constants import WIDE_CODED_MISSING, WIDE_NUCLEOTIDE_MISSING

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


def parse_nucleotide_call(text: str, missing: frozenset[str] = WIDE_NUCLEOTIDE_MISSING) -> tuple[str, str] | None:
    """Accept "A", "AA", "AT", "A/T", "A|T" (a separator only between two characters) and the IUPAC codes R Y S W K M (expanded);
    return the sorted allele pair,
    None for a token in ``missing`` or a pair with N, - or . (half-missing, undecided in contract 1.1.0);
    ValueError for any other cell ("?", "B", "H", "X", "0", "+", "A?")."""
    t = str(text).strip().upper()
    if t in missing:
        return None
    if len(t) == 1:
        if t in NUCLEOTIDES:
            return (t, t)
        if t in IUPAC_HET:
            return IUPAC_HET[t]
        raise ValueError(f"unrecognised nucleotide call {text!r}")
    if len(t) == 2:
        pair = (t[0], t[1])
    elif len(t) == 3 and t[1] in ("/", "|"):
        pair = (t[0], t[2])
    else:
        raise ValueError(f"unrecognised nucleotide call {text!r}")
    if not all(c in NUCLEOTIDES or c in ("N", "-", ".") for c in pair):
        raise ValueError(f"unrecognised nucleotide call {text!r}")
    if any(c not in NUCLEOTIDES for c in pair):
        return None
    return (min(pair), max(pair))


def parse_coded_call(text: str) -> tuple[int, int]:
    """A -> (0,0), B -> (1,1), H -> (0,1), missing -> (-1,-1)."""
    t = str(text).strip().upper()
    if t in WIDE_CODED_MISSING:
        return (-1, -1)
    if t == "A":
        return (0, 0)
    if t == "B":
        return (1, 1)
    if t == "H":
        return (0, 1)
    raise ValueError(f"unrecognised coded call {text!r} (expected A, B, H, N, NA or empty)")


def detect_coding(values: Iterable[str]) -> str:
    """'abh' when every cell outside the nucleotide missing set is A, B or H and at least one B or H occurs;
    else 'nucleotide'. Scans every value given."""
    seen: set[str] = set()
    for v in values:
        t = str(v).strip().upper()
        if t in WIDE_NUCLEOTIDE_MISSING:
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
