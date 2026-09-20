"""Staged ranking: lexicographic order on hand-built arrays (docs/adr/0007, amendment 2026-09-16)."""

from __future__ import annotations

import numpy as np

from progeny_selector.core.score import rank_rows_staged

SAMPLES = ["s1", "s2", "s3", "s4", "s5"]
REC_COUNT = np.array([2.0, 1.0, 1.0, 1.0, 1.0])
RPP_CARRIER = np.array([0.5, 0.9, 0.9, 0.9, np.nan])
RPP_NONCARRIER = np.array([0.5, 0.99, 0.95, 0.95, 0.99])
DRAG = np.array([0.0, 0.0, 5.0, 3.0, 0.0])
MISSING = np.zeros(5)
FAMILIES: list[str | None] = ["a", "a", "b", "b", "b"]


def staged(passes: np.ndarray, families: list[str | None] | None = None):
    return rank_rows_staged(
        REC_COUNT, RPP_CARRIER, RPP_NONCARRIER, DRAG, MISSING, SAMPLES, passes, FAMILIES if families is None else families
    )


def test_staged_order():
    rank_overall, _ = staged(np.ones(5, dtype=bool))
    # s1 leads on recombinant flanks; s2 on rpp_noncarrier; s4 before s3 on drag; s5's NaN carrier sorts last.
    assert list(rank_overall) == [1.0, 2.0, 4.0, 3.0, 5.0]


def test_excluded_row_is_nan_and_others_close_up():
    passes = np.ones(5, dtype=bool)
    passes[1] = False
    rank_overall, rank_in_family = staged(passes)
    assert np.isnan(rank_overall[1])
    assert list(rank_overall[[0, 2, 3, 4]]) == [1.0, 3.0, 2.0, 4.0]
    assert rank_overall[3] == 2.0
    assert np.isnan(rank_in_family[1])
    assert list(rank_in_family[[0, 2, 3, 4]]) == [1.0, 2.0, 1.0, 3.0]


def test_sample_id_breaks_a_full_tie():
    ids = ["b_second", "a_first"]
    ones = np.ones(2)
    rank_overall, _ = rank_rows_staged(ones, ones, ones, ones, np.zeros(2), ids, np.ones(2, dtype=bool), ["f", "f"])
    assert list(rank_overall) == [2.0, 1.0]


def test_all_nan_carrier_ranks_by_the_next_key():
    nan = np.full(3, np.nan)
    ids = ["x", "y", "z"]
    rank_overall, _ = rank_rows_staged(
        np.ones(3), nan, np.array([0.1, 0.9, 0.5]), np.zeros(3), np.zeros(3), ids, np.ones(3, dtype=bool), [None, None, None]
    )
    assert list(rank_overall) == [3.0, 1.0, 2.0]


def test_rpp_carrier_direction():
    """Equal recombinant counts: the higher rpp_carrier wins, even with a worse rpp_noncarrier."""
    ids = ["low_carrier", "high_carrier"]
    rank_overall, _ = rank_rows_staged(
        np.ones(2),
        np.array([0.5, 0.9]),
        np.array([1.0, 0.1]),
        np.zeros(2),
        np.zeros(2),
        ids,
        np.ones(2, dtype=bool),
        ["f", "f"],
    )
    assert list(rank_overall) == [2.0, 1.0]


def test_rpp_noncarrier_direction():
    ids = ["a", "b"]
    rank_overall, _ = rank_rows_staged(
        np.ones(2), np.ones(2), np.array([0.5, 0.9]), np.zeros(2), np.zeros(2), ids, np.ones(2, dtype=bool), ["f", "f"]
    )
    assert list(rank_overall) == [2.0, 1.0]


def test_drag_and_missing_rate_directions():
    """Everything above drag equal: the shorter drag wins; everything above missing_rate equal: the lower wins."""
    ids = ["a", "b"]
    ones = np.ones(2)
    rank_overall, _ = rank_rows_staged(ones, ones, ones, np.array([5.0, 3.0]), np.zeros(2), ids, np.ones(2, dtype=bool), ["f", "f"])
    assert list(rank_overall) == [2.0, 1.0]
    rank_overall, _ = rank_rows_staged(ones, ones, ones, np.zeros(2), np.array([0.2, 0.05]), ids, np.ones(2, dtype=bool), ["f", "f"])
    assert list(rank_overall) == [2.0, 1.0]


def test_recombinant_count_outranks_everything_below_it():
    ids = ["few_flanks", "many_flanks"]
    rank_overall, _ = rank_rows_staged(
        np.array([1.0, 2.0]),
        np.array([0.99, 0.10]),
        np.array([0.99, 0.10]),
        np.array([0.0, 9.0]),
        np.zeros(2),
        ids,
        np.ones(2, dtype=bool),
        ["f", "f"],
    )
    assert list(rank_overall) == [2.0, 1.0]
