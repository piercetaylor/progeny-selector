"""Hard filters, composite score and ranking.

Responsibility: apply hard filters (foreground must pass, avoid must pass,
missing-rate and QC policies), combine bounded [0, 1] components with user
weights, and assign ranks with a fixed tie-break order (PLAN.md, algorithm 7).
Pure functions on per-sample arrays; no I/O.

Interface:
    qc_excluding_hits(flags, filters) -> list[str]   excluding QC flags under the applied filters
    hard_filters(target_status, avoid_status, missing_rate, qc_flags, filters) -> (passes, reasons)
    composite_score(components, weights) -> (score, weight_sum_used)
    rank_rows(score, rpp_total, drag_est, missing_rate, sample_ids, passes, family_ids) -> (rank_overall, rank_in_family)
    rank_rows_staged(rec_count, rpp_carrier, rpp_noncarrier, drag_est, missing_rate, sample_ids, passes, family_ids)
        -> (rank_overall, rank_in_family)   lexicographic order (docs/adr/0007, amendment 2026-09-16)
"""

from __future__ import annotations

from typing import SupportsFloat

import numpy as np

from progeny_selector.constants import STATUS_FAIL, STATUS_PASS, STATUS_UNKNOWN
from progeny_selector.model.criteria import Filters, Weights

QC_EXCLUDING_FLAGS = ("possible_self_or_outcross", "possible_outcross", "possible_rp_sample", "possible_donor_sample")


def qc_excluding_hits(flags: list[str], filters: Filters) -> list[str]:
    """The flags in ``QC_EXCLUDING_FLAGS`` when ``filters.exclude_qc_flagged`` is true, else an empty list."""
    if not filters.exclude_qc_flagged:
        return []
    return [f for f in flags if f in QC_EXCLUDING_FLAGS]


def hard_filters(
    target_status: dict[str, np.ndarray],
    avoid_status: dict[str, np.ndarray],
    missing_rate: np.ndarray,
    qc_flags: list[list[str]],
    filters: Filters,
) -> tuple[np.ndarray, list[str]]:
    """Return (passes, reason) per sample. Reasons are ';'-joined, empty when passing."""
    n = len(missing_rate)
    reasons: list[list[str]] = [[] for _ in range(n)]
    for locus_id, status in target_status.items():
        bad = status == STATUS_FAIL
        if filters.unknown_target_is == "fail":
            bad |= status == STATUS_UNKNOWN
        for i in np.where(bad)[0]:
            reasons[i].append(f"target:{locus_id}:{'unknown' if status[i] == STATUS_UNKNOWN else 'fail'}")
    for locus_id, status in avoid_status.items():
        bad = status == STATUS_FAIL
        if filters.unknown_avoid_is == "fail":
            bad |= status == STATUS_UNKNOWN
        for i in np.where(bad)[0]:
            reasons[i].append(f"avoid:{locus_id}:{'unknown' if status[i] == STATUS_UNKNOWN else 'fail'}")
    for i in np.where(missing_rate > filters.max_missing_rate)[0]:
        reasons[i].append(f"missing_rate>{filters.max_missing_rate}")
    for i, flags in enumerate(qc_flags):
        hits = qc_excluding_hits(flags, filters)
        if hits:
            reasons[i].append("qc:" + "|".join(hits))
    passes = np.array([len(r) == 0 for r in reasons], dtype=bool)
    return passes, [";".join(r) for r in reasons]


def composite_score(components: dict[str, np.ndarray], weights: Weights) -> tuple[np.ndarray, np.ndarray]:
    """Weighted mean of components in [0, 1]; NaN components drop out and the weights renormalise.

    Returns (score, weight_sum_used). Score is NaN when every weighted component is NaN.
    """
    w = weights.as_dict()
    n = len(next(iter(components.values())))
    numer = np.zeros(n)
    denom = np.zeros(n)
    for name, values in components.items():
        wk = w.get(name, 0.0)
        if wk <= 0:
            continue
        v = np.asarray(values, dtype=float)
        ok = ~np.isnan(v)
        numer[ok] += wk * np.clip(v[ok], 0.0, 1.0)
        denom[ok] += wk
    with np.errstate(invalid="ignore", divide="ignore"):
        score = numer / denom
    score[denom == 0] = np.nan
    return score, denom


def rank_rows(
    score: np.ndarray,
    rpp_total: np.ndarray,
    drag_est: np.ndarray,
    missing_rate: np.ndarray,
    sample_ids: list[str],
    passes: np.ndarray,
    family_ids: list[str | None],
) -> tuple[np.ndarray, np.ndarray]:
    """Dense ranks (1 = best) among passing samples; NaN for excluded samples.

    Tie-break order: composite desc, RPP total desc, drag estimate asc, missing rate asc, sample_id asc.
    """
    n = len(sample_ids)
    order = sorted(
        range(n),
        key=lambda i: (
            -(score[i] if not np.isnan(score[i]) else -np.inf),
            -(rpp_total[i] if not np.isnan(rpp_total[i]) else -np.inf),
            drag_est[i] if not np.isnan(drag_est[i]) else np.inf,
            missing_rate[i] if not np.isnan(missing_rate[i]) else np.inf,
            sample_ids[i],
        ),
    )
    return _dense_ranks(order, passes, family_ids)


def rank_rows_staged(
    rec_count: np.ndarray,
    rpp_carrier: np.ndarray,
    rpp_noncarrier: np.ndarray,
    drag_est: np.ndarray,
    missing_rate: np.ndarray,
    sample_ids: list[str],
    passes: np.ndarray,
    family_ids: list[str | None],
) -> tuple[np.ndarray, np.ndarray]:
    """Dense ranks (1 = best) among passing samples in staged order; NaN for excluded samples.

    Order: recombinant flanks (count) desc, rpp_carrier desc, rpp_noncarrier desc, drag estimate asc,
    missing rate asc, sample_id asc. Exact ties at each key, no bins. NaN sorts last on every key.
    """
    n = len(sample_ids)
    order = sorted(
        range(n),
        key=lambda i: (
            -_desc(rec_count[i]),
            -_desc(rpp_carrier[i]),
            -_desc(rpp_noncarrier[i]),
            _asc(drag_est[i]),
            _asc(missing_rate[i]),
            sample_ids[i],
        ),
    )
    return _dense_ranks(order, passes, family_ids)


def _desc(value: SupportsFloat) -> float:
    """NaN -> -inf, so a NaN sorts last on a descending key."""
    v = float(value)
    return -np.inf if np.isnan(v) else v


def _asc(value: SupportsFloat) -> float:
    """NaN -> +inf, so a NaN sorts last on an ascending key."""
    v = float(value)
    return np.inf if np.isnan(v) else v


def _dense_ranks(order: list[int], passes: np.ndarray, family_ids: list[str | None]) -> tuple[np.ndarray, np.ndarray]:
    """Number the passing samples 1..k in ``order``, overall and within family; excluded samples stay NaN."""
    n = len(passes)
    rank_overall = np.full(n, np.nan)
    rank_in_family = np.full(n, np.nan)
    counter = 0
    fam_counter: dict[str, int] = {}
    for i in order:
        if not passes[i]:
            continue
        counter += 1
        rank_overall[i] = counter
        fam = family_ids[i] or ""
        fam_counter[fam] = fam_counter.get(fam, 0) + 1
        rank_in_family[i] = fam_counter[fam]
    return rank_overall, rank_in_family


def status_label(status: int) -> str:
    return {STATUS_PASS: "pass", STATUS_FAIL: "fail", STATUS_UNKNOWN: "unknown"}[int(status)]
