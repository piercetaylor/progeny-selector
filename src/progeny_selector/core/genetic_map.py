"""Marey-map cleaning and interpolation of genetic positions along physical positions.

Responsibility: given markers with both a physical position (bp) and a genetic
position (cM) on one chromosome, keep the markers whose cM is consistent with
their bp order and interpolate cM for other markers on the same assembly.
Outlier markers are invalidated, not smoothed (MareyMap practice; docs/adr/0019):
the kept set is the longest subsequence of cM that is non-decreasing along bp.
Pure numpy and the standard library; no I/O.

Interface:
    longest_nondecreasing_mask(cm_sorted_by_bp) -> np.ndarray[bool]
    interpolate_cm(mapped_bp, mapped_cm, query_bp) -> np.ndarray[float]
"""

from __future__ import annotations

from bisect import bisect_right

import numpy as np


def longest_nondecreasing_mask(cm_sorted_by_bp: np.ndarray) -> np.ndarray:
    """Boolean mask of a longest non-decreasing subsequence (equal values allowed) of finite cM values in bp order.

    Tie-break: among equal-length subsequences, the one keeping earlier markers (lexicographically smallest indices).
    """
    cm = np.asarray(cm_sorted_by_bp, dtype=float)
    n = cm.shape[0]
    mask = np.zeros(n, dtype=bool)
    if n == 0:
        return mask
    if not np.all(np.isfinite(cm)):
        raise ValueError("cm values must be finite")
    # length of the longest non-decreasing subsequence starting at each index, by a right-to-left pass
    # over -cm (a subsequence non-decreasing in cm going right is non-decreasing in -cm going left)
    start_len = np.zeros(n, dtype=int)
    tails: list[float] = []
    for i in range(n - 1, -1, -1):
        v = -float(cm[i])
        pos = bisect_right(tails, v)
        if pos == len(tails):
            tails.append(v)
        else:
            tails[pos] = v
        start_len[i] = pos + 1
    need = len(tails)
    last = -np.inf
    for i in range(n):
        if need == 0:
            break
        if start_len[i] == need and cm[i] >= last:
            mask[i] = True
            last = cm[i]
            need -= 1
    return mask


def interpolate_cm(mapped_bp: np.ndarray, mapped_cm: np.ndarray, query_bp: np.ndarray) -> np.ndarray:
    """cM at each query position, in query order, by linear interpolation over the kept mapped markers.

    Mapped markers are sorted by (bp, cM) and reduced to ``longest_nondecreasing_mask``; kept markers
    at equal bp collapse to one anchor at their mean cM, so the interpolation abscissae are strictly
    increasing. Queries outside the mapped range are clamped to the end cM (no extrapolation). One
    mapped marker gives its cM for every query; none gives NaN.
    """
    bp = np.asarray(mapped_bp, dtype=float)
    cm = np.asarray(mapped_cm, dtype=float)
    query = np.asarray(query_bp, dtype=float)
    if bp.shape != cm.shape:
        raise ValueError("mapped_bp and mapped_cm must have the same length")
    if bp.shape[0] == 0:
        return np.full(query.shape, np.nan)
    order = np.lexsort((cm, bp))
    bp_sorted = bp[order]
    cm_sorted = cm[order]
    keep = longest_nondecreasing_mask(cm_sorted)
    anchors_bp, inverse = np.unique(bp_sorted[keep], return_inverse=True)
    anchors_cm = np.bincount(inverse, weights=cm_sorted[keep]) / np.bincount(inverse)
    return np.interp(query, anchors_bp, anchors_cm)
