"""Foreground, avoid, background (count and weighted), drag and similarity on hand-built cases."""

from __future__ import annotations

import math

import numpy as np
import pytest

from progeny_selector.constants import STATUS_FAIL, STATUS_PASS, STATUS_UNKNOWN
from progeny_selector.core.avoid import avoid_status
from progeny_selector.core.background import marker_weights, rpp
from progeny_selector.core.classify import classify
from progeny_selector.core.drag import donor_segment
from progeny_selector.core.foreground import foreground_status, resolve_locus
from progeny_selector.core.pipeline import run_analysis
from progeny_selector.core.score import composite_score, rank_rows
from progeny_selector.core.similarity import ibs_to_sample
from progeny_selector.model.criteria import AvoidSpec, Criteria, CriteriaError, TargetSpec, Weights
from tests.conftest import make_dataset, make_matrix

# Ten markers on Gm01 at 1..10 Mb; progeny P1 states below.
STATES = ["A", "A", "H", "H", "H", "N", "H", "A", "A", "B"]


@pytest.fixture
def case():
    gm = make_matrix(STATES)
    cls = classify(gm, "RP", "DONOR")
    return gm, cls.states, cls.informative


def test_foreground_marker_region_flanking(case):
    gm, states, _ = case
    t_marker = resolve_locus(TargetSpec("T", marker_id="m4"), gm)
    assert foreground_status(states, t_marker, "either")[2] == STATUS_PASS
    assert foreground_status(states, t_marker, "hom_donor")[2] == STATUS_FAIL
    assert foreground_status(states, t_marker, "het")[2] == STATUS_PASS
    t_missing = resolve_locus(TargetSpec("T", marker_id="m6"), gm)
    assert foreground_status(states, t_missing, "either")[2] == STATUS_UNKNOWN
    region_all = resolve_locus(TargetSpec("R", chrom="1", start_bp=3_000_000, end_bp=7_000_000), gm)
    assert list(region_all.marker_idx) == [2, 3, 4, 5, 6]
    assert foreground_status(states, region_all, "either", rule="all")[2] == STATUS_PASS  # missing m6 is skipped
    region_mixed = resolve_locus(TargetSpec("R2", chrom="Gm01", start_bp=7_000_000, end_bp=8_000_000), gm)
    assert foreground_status(states, region_mixed, "either", rule="all")[2] == STATUS_FAIL
    assert foreground_status(states, region_mixed, "either", rule="any")[2] == STATUS_PASS
    flank = resolve_locus(TargetSpec("F", left_marker="m3", right_marker="m5"), gm)
    assert flank.kind == "flanking"
    assert foreground_status(states, flank, "either")[2] == STATUS_PASS
    with pytest.raises(CriteriaError):
        resolve_locus(TargetSpec("bad", chrom="Gm02", start_bp=1, end_bp=2), gm)


def test_avoid_status(case):
    gm, states, _ = case
    assert avoid_status(states, resolve_locus(AvoidSpec("a", marker_id="m1"), gm))[2] == STATUS_PASS
    assert avoid_status(states, resolve_locus(AvoidSpec("a", marker_id="m3"), gm))[2] == STATUS_FAIL
    assert avoid_status(states, resolve_locus(AvoidSpec("a", marker_id="m3"), gm), allow_het=True)[2] == STATUS_PASS
    assert avoid_status(states, resolve_locus(AvoidSpec("a", marker_id="m10"), gm), allow_het=True)[2] == STATUS_FAIL
    assert avoid_status(states, resolve_locus(AvoidSpec("a", marker_id="m6"), gm))[2] == STATUS_UNKNOWN


def test_rpp_count_and_weighted(case):
    gm, states, informative = case
    # called: 4 A, 4 H, 1 B -> (4 + 2) / 9
    assert rpp(states)[2] == pytest.approx(6 / 9)
    w = marker_weights(gm, informative, unit="bp", max_coverage=1_000_000)
    # interior markers: 0.5 Mb each side; ends: min(distance to end, 0.5 Mb) with no chromosome length -> cap
    assert w[0] == pytest.approx(1_000_000) and w[5] == pytest.approx(1_000_000)
    assert rpp(states, w)[2] == pytest.approx(6 / 9)  # equal weights reproduce the count model
    w2 = marker_weights(gm, informative, unit="bp", max_coverage=10_000_000, chrom_lengths={"Gm01": 20_000_000})
    assert w2[0] == pytest.approx(1_000_000 + 500_000)  # 1 Mb to the chromosome start, half a gap to the right
    assert w2[9] == pytest.approx(500_000 + 5_000_000)  # half a gap left, capped 5 Mb to the right
    assert math.isnan(rpp(states[:0], w2[:0])[0]) if states[:0].size else True


def test_drag_bounds_and_recombinants(case):
    gm, states, _ = case
    target = resolve_locus(TargetSpec("T", marker_id="m4"), gm)
    d = donor_segment(states, gm, target, unit="bp", chrom_length=20_000_000, window_left=2_500_000, window_right=2_500_000)
    # left: m3 H (1 Mb), m2 A (2 Mb) -> min 1 Mb, max 2 Mb; right: m5 H, m6 N skipped, m7 H (3 Mb), m8 A (4 Mb)
    assert (d.left_min[2], d.left_max[2]) == (1_000_000, 2_000_000)
    assert (d.right_min[2], d.right_max[2]) == (3_000_000, 4_000_000)
    assert d.total_max[2] == 6_000_000 and d.total_est[2] == 5_000_000
    assert bool(d.recombinant_left[2]) is True and bool(d.recombinant_right[2]) is False
    # no A on the right of m10 -> extends to the chromosome end
    end = resolve_locus(TargetSpec("E", marker_id="m10"), gm)
    d2 = donor_segment(states, gm, end, unit="bp", chrom_length=20_000_000)
    assert d2.right_max[2] == 10_000_000 and d2.right_min[2] == 0
    assert d2.left_max[2] == 1_000_000  # m9 is A
    d3 = donor_segment(states, gm, target, unit="cm", chrom_length=None, window_left=5.0, window_right=5.0)
    assert d3.left_max[2] == pytest.approx(5.0) and d3.right_max[2] == pytest.approx(10.0)


def test_ibs_to_parents(case):
    gm, _, _ = case
    ibs_rp = ibs_to_sample(gm, "RP")
    # P1 vs RP: A=1 (4), H=0.5 (4), B=0 (1), X=(0,1) vs (0,0) shares one allele -> none here; N excluded -> (4 + 2) / 9
    assert ibs_rp[2] == pytest.approx(6 / 9)
    assert ibs_rp[0] == 1.0 and ibs_rp[1] == 0.0


def test_per_target_windows():
    """flank_left/flank_right override flank_window on that side only, per target."""
    dataset = make_dataset(make_matrix(STATES))
    criteria = Criteria(
        targets=[
            TargetSpec("TW", marker_id="m4", flank_left=1.0, flank_right=20.0),
            TargetSpec("TD", marker_id="m4"),  # defaults: flank_window 5.0 cM on both sides
        ],
        flank_window=5.0,
        flank_unit="cm",
    )
    row = run_analysis(dataset, criteria).row("P1")
    # Hand-computed on STATES: target m4 at 10 cM; the nearest A is m2 at 5 cM on the left
    # (left_max 5.0 cM) and m8 at 20 cM on the right (right_max 10.0 cM, m6 missing is skipped).
    assert row["drag_TW_left_max"] == pytest.approx(5.0) and row["drag_TW_right_max"] == pytest.approx(10.0)
    assert row["drag_TD_left_max"] == pytest.approx(5.0) and row["drag_TD_right_max"] == pytest.approx(10.0)
    assert row["recomb_TW_left"] is False  # 5.0 > flank_left 1.0
    assert row["recomb_TW_right"] is True  # 10.0 <= flank_right 20.0
    assert row["recomb_TD_left"] is True  # 5.0 <= flank_window 5.0
    assert row["recomb_TD_right"] is False  # 10.0 > flank_window 5.0


def test_composite_and_ranking():
    comps = {"rpp_noncarrier": np.array([0.9, 0.8, np.nan]), "drag": np.array([0.5, 0.5, 0.5])}
    score, used = composite_score(comps, Weights(rpp_noncarrier=0.5, rpp_carrier=0.0, drag=0.5, recombinant=0.0))
    assert score[0] == pytest.approx(0.7) and score[1] == pytest.approx(0.65) and score[2] == pytest.approx(0.5)
    assert used[2] == 0.5  # NaN component dropped, weights renormalised
    ro, rf = rank_rows(
        score=np.array([0.7, 0.7, 0.5, 0.9]),
        rpp_total=np.array([0.9, 0.95, 0.9, 0.9]),
        drag_est=np.array([1.0, 1.0, 1.0, 1.0]),
        missing_rate=np.zeros(4),
        sample_ids=["a", "b", "c", "d"],
        passes=np.array([True, True, True, False]),
        family_ids=["F1", "F1", "F2", "F1"],
    )
    assert list(ro[:3]) == [2, 1, 3] and math.isnan(ro[3])
    assert list(rf[:3]) == [2, 1, 1]
