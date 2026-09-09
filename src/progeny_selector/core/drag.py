"""Linkage drag and recombinant selection around a target locus.

Responsibility: for each individual, bound the donor segment that carries the
target on its chromosome and flag recombination on each flank inside a user
window (PLAN.md, algorithm 4). The "max" bound matches Flapjack's linkage-drag
definition (distance from the locus to the first recombination on each side,
or to the chromosome end) [web] https://flapjack.hutton.ac.uk/en/latest/mabc.html.
The "min" bound is the distance to the outermost marker still carrying donor.

Interface:
    donor_segment(states, gm, resolved, unit, chrom_length, window_left, window_right) -> DragResult
    DragResult fields are float arrays of length n_samples (distances in ``unit``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from progeny_selector.constants import STATE_A, STATE_B, STATE_H
from progeny_selector.core.foreground import ResolvedLocus
from progeny_selector.model.dataset import GenotypeMatrix


@dataclass
class DragResult:
    left_min: np.ndarray
    left_max: np.ndarray
    right_min: np.ndarray
    right_max: np.ndarray
    recombinant_left: np.ndarray
    recombinant_right: np.ndarray
    unit: str

    @property
    def total_max(self) -> np.ndarray:
        return self.left_max + self.right_max

    @property
    def total_min(self) -> np.ndarray:
        return self.left_min + self.right_min

    @property
    def total_est(self) -> np.ndarray:
        """Midpoint estimate: the recombination is equally likely anywhere in the bounding interval."""
        return (self.total_min + self.total_max) / 2.0


def _walk(sub_states: np.ndarray, distances: np.ndarray, end_distance: float) -> tuple[np.ndarray, np.ndarray]:
    """Walk outward from the locus over markers ordered by increasing distance.

    Returns (min_bound, max_bound) per sample. max_bound = distance to the first A marker
    (or ``end_distance`` when none). min_bound = distance to the outermost H/B marker before
    that A marker (0 when none). N, X and U calls are skipped without terminating the walk.
    """
    n_samples = sub_states.shape[1]
    if sub_states.shape[0] == 0:
        return np.zeros(n_samples), np.full(n_samples, end_distance, dtype=float)
    is_a = sub_states == STATE_A
    has_a = is_a.any(axis=0)
    first_a = np.where(has_a, is_a.argmax(axis=0), sub_states.shape[0])
    max_bound = np.where(has_a, distances[np.minimum(first_a, len(distances) - 1)], end_distance)
    order = np.arange(sub_states.shape[0])[:, None]
    donor = ((sub_states == STATE_H) | (sub_states == STATE_B)) & (order < first_a[None, :])
    has_donor = donor.any(axis=0)
    last_donor = np.where(has_donor, sub_states.shape[0] - 1 - donor[::-1].argmax(axis=0), 0)
    min_bound = np.where(has_donor, distances[last_donor], 0.0)
    return min_bound.astype(float), max_bound.astype(float)


def donor_segment(
    states: np.ndarray,
    gm: GenotypeMatrix,
    resolved: ResolvedLocus,
    unit: str = "bp",
    chrom_length: float | None = None,
    window_left: float | None = None,
    window_right: float | None = None,
) -> DragResult:
    """Bound the donor segment around ``resolved`` on its chromosome for every sample.

    Distances are measured from the locus boundaries (start_bp / end_bp, converted to ``unit``)
    to marker positions outside the locus. Markers inside the locus define foreground status
    and are excluded from the walk. Without a chromosome length the left end is 0 and the
    right end is the last marker position on the chromosome.
    """
    pos = gm.positions(unit)
    chroms = gm.chroms()
    on_chrom = (chroms == resolved.chrom) & ~np.isnan(pos)
    bp = gm.positions("bp")
    locus_members = np.zeros(gm.n_markers, dtype=bool)
    locus_members[resolved.marker_idx] = True
    left_start = _to_unit(resolved.start_bp, bp, pos, on_chrom)
    right_start = _to_unit(resolved.end_bp, bp, pos, on_chrom)

    left_idx = np.where(on_chrom & ~locus_members & (bp < resolved.start_bp))[0]
    right_idx = np.where(on_chrom & ~locus_members & (bp > resolved.end_bp))[0]
    left_idx = left_idx[np.argsort(-pos[left_idx])]
    right_idx = right_idx[np.argsort(pos[right_idx])]

    chrom_pos = pos[on_chrom]
    right_end = (chrom_length if chrom_length is not None else float(chrom_pos.max())) - right_start
    left_end = left_start - 0.0
    left_min, left_max = _walk(states[left_idx], left_start - pos[left_idx], max(left_end, 0.0))
    right_min, right_max = _walk(states[right_idx], pos[right_idx] - right_start, max(right_end, 0.0))

    rec_left = np.zeros(states.shape[1], dtype=bool) if window_left is None else left_max <= window_left
    rec_right = np.zeros(states.shape[1], dtype=bool) if window_right is None else right_max <= window_right
    return DragResult(left_min, left_max, right_min, right_max, rec_left, rec_right, unit)


def _to_unit(value_bp: float, bp: np.ndarray, pos: np.ndarray, on_chrom: np.ndarray) -> float:
    """Convert a bp coordinate to the working unit by linear interpolation on this chromosome's map."""
    if pos is bp or np.array_equal(pos[on_chrom], bp[on_chrom]):
        return float(value_bp)
    xs = bp[on_chrom]
    ys = pos[on_chrom]
    order = np.argsort(xs)
    return float(np.interp(value_bp, xs[order], ys[order]))
