"""Shared call-string parsing for HapMap and wide-CSV inputs (contract/data-contract.md 1.1.0; mirrors backcross src/io/calls.ts).

Responsibility: turn per-sample call strings into (allele_a, allele_b) pairs,
build a per-marker allele table, and detect A/B/H coding. VCF has its own
allele indices and does not use this module.

Under a token profile (profiles.py, contract 1.4.0) a cell is resolved in this order:
the profile's missing tokens; ``missing`` when the profile's base is nucleotide; the
profile's homozygous tokens; its heterozygous tokens (a ``*`` token returns
HET_OF_MARKER, resolved by encode_marker once the row's alleles are known); the grammar
above when the base is nucleotide; otherwise ValueError.

Interface:
    parse_nucleotide_call(text, missing=WIDE_NUCLEOTIDE_MISSING, profile=None) -> tuple[str, str] | None | HET_OF_MARKER
        (None = missing; ValueError for any other cell, including "?")
    resolve_het_of_marker(alleles, cell) -> tuple[str, str]   (ValueError unless exactly two alleles)
    parse_coded_call(text) -> tuple[int, int]               (-1, -1 = missing)
    detect_coding(values: Iterable[str]) -> 'abh' | 'nucleotide'
    encode_marker(calls, seed_alleles=(), cells=None) -> (alleles, int8 array (n, 2))
        (a HET_OF_MARKER call is resolved against the seed alleles plus the symbols of the row's other calls;
        ``cells`` are the raw cell texts, named in that error)
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from enum import Enum

import numpy as np

from progeny_selector.constants import WIDE_CODED_MISSING, WIDE_NUCLEOTIDE_MISSING
from progeny_selector.io.profiles import CompiledProfile


class _HetOfMarker(Enum):
    """Sentinel for a profile heterozygote token that names no alleles (``*``)."""

    HET_OF_MARKER = "het-of-marker"


HET_OF_MARKER = _HetOfMarker.HET_OF_MARKER
Call = tuple[str, str] | None | _HetOfMarker

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


def parse_nucleotide_call(text: str, missing: frozenset[str] = WIDE_NUCLEOTIDE_MISSING, profile: CompiledProfile | None = None) -> Call:
    """Accept "A", "AA", "AT", "A/T", "A|T" (a separator only between two characters) and the IUPAC codes R Y S W K M (expanded);
    return the sorted allele pair,
    None for a token in ``missing`` or a pair with N, - or . (half-missing, undecided in contract 1.1.0);
    ValueError for any other cell ("?", "B", "H", "X", "0", "+", "A?")."""
    t = str(text).strip().upper()
    if profile is not None:
        if t in profile.missing or (profile.base == "nucleotide" and t in missing):
            return None
        if t in profile.homozygous:
            symbol = profile.homozygous[t]
            return (symbol, symbol)
        if t in profile.heterozygous:
            pair = profile.heterozygous[t]
            if pair == "*":
                return HET_OF_MARKER
            return (min(pair[0], pair[1]), max(pair[0], pair[1]))
        if profile.base == "none":
            raise ValueError(
                f'unrecognised nucleotide call {text!r} (not a token of the "{profile.label}" token profile; '
                "an A/B/H-coded file is read without a profile)"
            )
    elif t in missing:
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


def resolve_het_of_marker(alleles: Sequence[str], cell: str) -> tuple[str, str]:
    """The pair a ``*`` heterozygote stands for: the marker's two alleles; ValueError with any other number or an indel ``-``."""
    if len(alleles) != 2:
        listed = ", ".join(alleles)
        raise ValueError(f'"{cell}" is a heterozygote token but the marker shows {len(alleles)} allele(s) ({listed}); it needs exactly two')
    if "-" in alleles:
        listed = ", ".join(alleles)
        raise ValueError(f'"{cell}" is a heterozygote token but the marker shows an indel allele ({listed}); it needs two nucleotides')
    a, b = sorted(alleles)
    return (a, b)


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


def encode_marker(
    calls: Sequence[Call], seed_alleles: Sequence[str] = (), cells: Sequence[str] | None = None
) -> tuple[list[str], np.ndarray]:
    """Allele table (sorted alphabetically) and int8 (n_samples, 2) indices for one marker.

    A HET_OF_MARKER call resolves to the row's two alleles: ``seed_alleles`` (HapMap's alleles column)
    plus the symbols of the other calls; the allele table holds only symbols the resolved calls use."""
    if any(c is HET_OF_MARKER for c in calls):
        row_alleles = sorted(set(seed_alleles) | {a for c in calls if isinstance(c, tuple) for a in c})
        resolved: list[tuple[str, str] | None] = []
        for i, c in enumerate(calls):
            if c is HET_OF_MARKER:
                cell = str(cells[i]).strip().upper() if cells is not None else "*"
                resolved.append(resolve_het_of_marker(row_alleles, cell))
            else:
                resolved.append(c if isinstance(c, tuple) else None)
    else:
        resolved = [c if isinstance(c, tuple) else None for c in calls]
    alleles = sorted({a for c in resolved if c is not None for a in c})
    lookup = {a: i for i, a in enumerate(alleles)}
    out = np.full((len(calls), 2), -1, dtype=np.int8)
    for i, c in enumerate(resolved):
        if c is not None:
            out[i, 0] = lookup[c[0]]
            out[i, 1] = lookup[c[1]]
    return alleles, out
