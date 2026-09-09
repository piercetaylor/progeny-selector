"""Locus resolution and foreground status at target loci.

Responsibility: map a LocusSpec (marker, region, or flanking pair) onto marker
indices, then decide pass/fail/unknown per individual for the required donor
state. Pure functions over the classification matrix (PLAN.md, algorithm 2).

Interface:
    resolve_locus(spec, gm) -> ResolvedLocus
    marker_predicate(states, required_state) -> bool array (True where the call meets the state)
    locus_status(states, resolved, predicate, rule, min_markers) -> int8 status per sample
    foreground_status(states, resolved, required_state, rule, min_markers) -> int8 per sample
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from progeny_selector.constants import STATE_A, STATE_B, STATE_H, STATUS_FAIL, STATUS_PASS, STATUS_UNKNOWN
from progeny_selector.core.chrom import normalize_chrom
from progeny_selector.core.classify import is_called_informative
from progeny_selector.model.criteria import CriteriaError, LocusSpec
from progeny_selector.model.dataset import GenotypeMatrix


@dataclass(frozen=True)
class ResolvedLocus:
    locus_id: str
    kind: str
    chrom: str
    start_bp: int
    end_bp: int
    marker_idx: np.ndarray  # indices into gm.markers, sorted by position

    @property
    def mid_bp(self) -> float:
        return (self.start_bp + self.end_bp) / 2.0


def resolve_locus(spec: LocusSpec, gm: GenotypeMatrix) -> ResolvedLocus:
    """Resolve a locus definition to marker indices; raise CriteriaError when nothing matches."""
    kind = spec.kind()
    if kind == "marker":
        i = gm.marker_index(spec.marker_id)  # type: ignore[arg-type]
        m = gm.markers[i]
        return ResolvedLocus(spec.locus_id, kind, m.chrom, m.pos_bp, m.pos_bp, np.array([i]))
    if kind == "flanking":
        i = gm.marker_index(spec.left_marker)  # type: ignore[arg-type]
        j = gm.marker_index(spec.right_marker)  # type: ignore[arg-type]
        mi, mj = gm.markers[i], gm.markers[j]
        if mi.chrom != mj.chrom:
            raise CriteriaError(f"locus {spec.locus_id!r}: flanking markers are on different chromosomes")
        if mi.pos_bp > mj.pos_bp:
            i, j, mi, mj = j, i, mj, mi
        return ResolvedLocus(spec.locus_id, kind, mi.chrom, mi.pos_bp, mj.pos_bp, np.array([i, j]))
    chrom = normalize_chrom(spec.chrom)  # type: ignore[arg-type]
    idx = [
        k
        for k, m in enumerate(gm.markers)
        if m.chrom == chrom and spec.start_bp <= m.pos_bp <= spec.end_bp  # type: ignore[operator]
    ]
    if not idx:
        raise CriteriaError(f"locus {spec.locus_id!r}: no markers in {chrom}:{spec.start_bp}-{spec.end_bp}")
    idx.sort(key=lambda k: gm.markers[k].pos_bp)
    return ResolvedLocus(spec.locus_id, kind, chrom, int(spec.start_bp), int(spec.end_bp), np.array(idx))  # type: ignore[arg-type]


def marker_predicate(states: np.ndarray, required_state: str) -> np.ndarray:
    """Boolean array: does each call satisfy the required donor state?"""
    if required_state == "hom_donor":
        return states == STATE_B
    if required_state == "het":
        return states == STATE_H
    if required_state == "either":
        return (states == STATE_H) | (states == STATE_B)
    if required_state == "hom_rp":
        return states == STATE_A
    if required_state == "rp_or_het":
        return (states == STATE_A) | (states == STATE_H)
    raise ValueError(f"unknown required_state {required_state!r}")


def locus_status(
    states: np.ndarray,
    resolved: ResolvedLocus,
    predicate: np.ndarray,
    rule: str = "all",
    min_markers: int = 1,
) -> np.ndarray:
    """Combine per-marker predicates over the locus markers into pass/fail/unknown per sample.

    Only called, informative markers (A/H/B) count. With fewer than ``min_markers`` such
    markers the status is unknown. rule='all': every counted marker must satisfy the
    predicate; rule='any': at least one must.
    """
    sub = states[resolved.marker_idx, :]
    pred = predicate[resolved.marker_idx, :]
    counted = is_called_informative(sub)
    n_counted = counted.sum(axis=0)
    n_pass = (pred & counted).sum(axis=0)
    status = np.full(states.shape[1], STATUS_UNKNOWN, dtype=np.int8)
    enough = n_counted >= min_markers
    if rule == "all":
        ok = n_pass == n_counted
    elif rule == "any":
        ok = n_pass >= 1
    else:
        raise ValueError(f"unknown rule {rule!r}")
    status[enough & ok] = STATUS_PASS
    status[enough & ~ok] = STATUS_FAIL
    return status


def foreground_status(
    states: np.ndarray,
    resolved: ResolvedLocus,
    required_state: str,
    rule: str = "all",
    min_markers: int = 1,
) -> np.ndarray:
    """Foreground status per individual at one target locus."""
    if resolved.kind == "flanking":
        rule, min_markers = "all", 2
    return locus_status(states, resolved, marker_predicate(states, required_state), rule, min_markers)
