"""Background selection: recurrent-parent proportion (RPP) overall and by chromosome.

Responsibility: compute count-based and map-weighted RPP from the state matrix
(PLAN.md, algorithm 3). The weighted model follows Flapjack's MABC definition:
each marker represents at most ``max_marker_coverage`` of map, split half to
each side, and heterozygous calls contribute half [web]
https://flapjack.hutton.ac.uk/en/latest/mabc.html.

Interface:
    marker_weights(gm, informative, unit, max_coverage, chrom_lengths) -> float (n_markers,)
    rpp(states, weights=None, marker_mask=None) -> float (n_samples,)
    rpp_per_chromosome(states, gm, weights=None) -> dict[chrom, float (n_samples,)]
    carrier_mask(gm, carrier_chroms) -> bool (n_markers,)
"""

from __future__ import annotations

import numpy as np

from progeny_selector.core.chrom import chrom_sort_key
from progeny_selector.core.classify import is_called_informative, rpp_contribution
from progeny_selector.model.dataset import GenotypeMatrix

DEFAULT_MAX_COVERAGE = {"cm": 10.0, "bp": 4_000_000.0}


def marker_weights(
    gm: GenotypeMatrix,
    informative: np.ndarray,
    unit: str = "bp",
    max_coverage: float | None = None,
    chrom_lengths: dict[str, float] | None = None,
) -> np.ndarray:
    """Map-interval weight per informative marker; zero for uninformative markers.

    For each informative marker on a chromosome, weight = min(d_left/2, c/2) + min(d_right/2, c/2)
    where d_left/d_right are distances to the neighbouring informative markers and c is
    ``max_coverage``. At chromosome ends the outer side uses min(distance to the chromosome
    end, c/2) when the chromosome length is known, else c/2. Positions must be sorted.
    """
    cap = (max_coverage if max_coverage is not None else DEFAULT_MAX_COVERAGE[unit]) / 2.0
    pos = gm.positions(unit)
    chroms = gm.chroms()
    weights = np.zeros(gm.n_markers, dtype=float)
    for chrom in sorted(set(chroms.tolist()), key=chrom_sort_key):
        idx = np.where((chroms == chrom) & informative & ~np.isnan(pos))[0]
        if idx.size == 0:
            continue
        p = pos[idx]
        length = None if chrom_lengths is None else chrom_lengths.get(chrom)
        left = np.empty(idx.size)
        right = np.empty(idx.size)
        if idx.size == 1:
            left[0] = cap if length is None else min(p[0], cap)
            right[0] = cap if length is None else min(length - p[0], cap)
        else:
            gaps = np.diff(p) / 2.0
            left[1:] = np.minimum(gaps, cap)
            right[:-1] = np.minimum(gaps, cap)
            left[0] = cap if length is None else min(p[0], cap)
            right[-1] = cap if length is None else min(length - p[-1], cap)
        weights[idx] = np.maximum(left, 0) + np.maximum(right, 0)
    return weights


def rpp(states: np.ndarray, weights: np.ndarray | None = None, marker_mask: np.ndarray | None = None) -> np.ndarray:
    """RPP per sample: sum(w * contribution) / sum(w) over called informative markers in the mask.

    contribution: A=1, H=0.5, B=0. X, N and U calls are excluded from both sums, which is
    equivalent to imputing each individual's own mean. NaN when no marker qualifies.
    """
    contrib = rpp_contribution(states)
    counted = is_called_informative(states)
    if marker_mask is not None:
        counted = counted & marker_mask[:, None]
    w = np.ones(states.shape[0], dtype=float) if weights is None else weights
    w2 = np.where(counted, w[:, None], 0.0)
    denom = w2.sum(axis=0)
    numer = np.nansum(w2 * np.nan_to_num(contrib), axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = numer / denom
    out[denom == 0] = np.nan
    return out


def rpp_per_chromosome(states: np.ndarray, gm: GenotypeMatrix, weights: np.ndarray | None = None) -> dict[str, np.ndarray]:
    chroms = gm.chroms()
    return {c: rpp(states, weights, chroms == c) for c in sorted(set(chroms.tolist()), key=chrom_sort_key)}


def carrier_mask(gm: GenotypeMatrix, carrier_chroms: set[str]) -> np.ndarray:
    """Boolean marker mask for chromosomes that carry at least one target locus."""
    chroms = gm.chroms()
    return np.isin(chroms, list(carrier_chroms))


def n_called_informative(states: np.ndarray, marker_mask: np.ndarray | None = None) -> np.ndarray:
    counted = is_called_informative(states)
    if marker_mask is not None:
        counted = counted & marker_mask[:, None]
    return counted.sum(axis=0)
