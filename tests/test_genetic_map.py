"""Hand-computed cases for core/genetic_map.py (docs/adr/0019)."""

from __future__ import annotations

import numpy as np
import pytest

from progeny_selector.core.genetic_map import interpolate_cm, longest_nondecreasing_mask


def test_mask_drops_the_inversion():
    assert longest_nondecreasing_mask(np.array([0, 5, 4, 10])).tolist() == [True, True, False, True]


def test_mask_allows_ties():
    assert longest_nondecreasing_mask(np.array([0, 5, 5, 10])).tolist() == [True, True, True, True]


def test_mask_prefers_earlier_markers():
    # length-2 subsequences of [1, 0, 3, 2] use indices (0,2), (0,3), (1,2), (1,3); the earliest is (0,2)
    assert longest_nondecreasing_mask(np.array([1, 0, 3, 2])).tolist() == [True, False, True, False]
    assert longest_nondecreasing_mask(np.array([3, 2, 1])).tolist() == [True, False, False]


def test_mask_empty_and_non_finite():
    assert longest_nondecreasing_mask(np.array([])).tolist() == []
    with pytest.raises(ValueError):
        longest_nondecreasing_mask(np.array([0.0, np.nan]))


def test_interpolate_over_kept_markers_with_clamping():
    bp = np.array([1e6, 2e6, 3e6, 4e6])
    cm = np.array([0, 5, 4, 10])
    # marker at 3e6 (cM 4) is dropped: kept (1e6, 0), (2e6, 5), (4e6, 10)
    out = interpolate_cm(bp, cm, np.array([0.5e6, 2.5e6, 3.5e6, 9e6]))
    assert out == pytest.approx([0.0, 6.25, 8.75, 10.0])


def test_interpolate_input_order_preserved():
    bp = np.array([4e6, 2e6, 1e6, 3e6])
    cm = np.array([10, 5, 0, 4])
    out = interpolate_cm(bp, cm, np.array([9e6, 3.5e6, 0.5e6, 2.5e6]))
    assert out == pytest.approx([10.0, 8.75, 0.0, 6.25])


def test_interpolate_equal_bp_uses_mean_anchor():
    out = interpolate_cm(np.array([1e6, 1e6, 2e6]), np.array([5.0, 0.0, 6.0]), np.array([1e6, 0.0, 1.5e6]))
    assert out == pytest.approx([2.5, 2.5, 4.25])


def test_interpolate_single_and_none():
    assert interpolate_cm(np.array([2e6]), np.array([7.5]), np.array([1.0, 5e6])) == pytest.approx([7.5, 7.5])
    out = interpolate_cm(np.array([]), np.array([]), np.array([1.0, 2.0]))
    assert out.shape == (2,)
    assert np.isnan(out).all()
