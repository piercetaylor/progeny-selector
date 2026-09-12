"""Genome-strip segments and locus ticks on hand-built matrices."""

from __future__ import annotations

import numpy as np
import pytest

from progeny_selector.app.present import strip_rects
from progeny_selector.constants import STATE_A, STATE_B, STATE_H, STATE_N
from progeny_selector.core.classify import classify
from progeny_selector.core.foreground import resolve_locus
from progeny_selector.core.strip import ChromStrip, chromosome_strips, locus_ticks
from progeny_selector.model.criteria import AvoidSpec, TargetSpec
from progeny_selector.model.dataset import GenotypeMatrix, Marker
from tests.conftest import make_matrix

GM01_LEN = 57_932_356


def progeny_states(gm: GenotypeMatrix) -> np.ndarray:
    return classify(gm, "RP", "DONOR").states[:, 2]


def strips_for(states: list[str], **kwargs) -> tuple[GenotypeMatrix, list[ChromStrip]]:
    gm = make_matrix(states, **kwargs)
    return gm, chromosome_strips(progeny_states(gm), gm)


def two_chrom_matrix() -> GenotypeMatrix:
    """Markers on Gm11 and Gm2: alphabetical order would put Gm11 first, chromosome order puts Gm2 first."""
    a, b = make_matrix(["A", "H"], chrom="Gm11"), make_matrix(["B"], chrom="Gm2")
    markers = [Marker(f"{m.chrom}_{m.marker_id}", m.chrom, m.pos_bp, m.cm) for m in a.markers + b.markers]
    return GenotypeMatrix(
        markers=markers,
        sample_ids=a.sample_ids,
        alleles=a.alleles + b.alleles,
        calls=np.concatenate([a.calls, b.calls]),
    ).sorted_by_position()


def positioned_matrix(states: list[str], positions_bp: list[int]) -> GenotypeMatrix:
    """make_matrix with explicit marker positions, kept in the given (possibly unsorted) order."""
    gm = make_matrix(states)
    markers = [Marker(m.marker_id, m.chrom, p) for m, p in zip(gm.markers, positions_bp, strict=True)]
    return GenotypeMatrix(markers=markers, sample_ids=gm.sample_ids, alleles=gm.alleles, calls=gm.calls)


def test_segments_split_at_midpoints():
    _, strips = strips_for(["A", "H", "A"])
    assert len(strips) == 1
    s = strips[0]
    assert (s.chrom, s.length_bp, s.n_markers) == ("Gm01", GM01_LEN, 3)
    expected = [(0.0, 1.5e6 / GM01_LEN, STATE_A), (1.5e6 / GM01_LEN, 2.5e6 / GM01_LEN, STATE_H), (2.5e6 / GM01_LEN, 1.0, STATE_A)]
    assert len(s.segments) == 3
    for seg, (start, end, state) in zip(s.segments, expected, strict=True):
        assert seg.start == pytest.approx(start, abs=1e-9)
        assert seg.end == pytest.approx(end, abs=1e-9)
        assert seg.state == state


def test_equal_states_merge():
    _, strips = strips_for(["A", "A", "A"])
    assert [(seg.start, seg.end, seg.state) for seg in strips[0].segments] == [(0.0, 1.0, STATE_A)]


def test_missing_keeps_its_own_segment():
    _, strips = strips_for(["N", "A"])
    assert [seg.state for seg in strips[0].segments] == [STATE_N, STATE_A]
    assert strips[0].segments[0].end == pytest.approx(1.5e6 / GM01_LEN, abs=1e-9)


def test_input_order_does_not_matter():
    unsorted = positioned_matrix(["A", "H", "B"], [3_000_000, 1_000_000, 2_000_000])
    ordered = positioned_matrix(["H", "B", "A"], [1_000_000, 2_000_000, 3_000_000])
    got = chromosome_strips(progeny_states(unsorted), unsorted)
    want = chromosome_strips(progeny_states(ordered), ordered)
    assert got == want
    assert [seg.state for seg in got[0].segments] == [STATE_H, STATE_B, STATE_A]


def co_located_strips() -> list[ChromStrip]:
    """Three calls at 2 Mb (H, B, H) and two at 70 Mb (A, B), past the Gm01 assembly end."""
    gm = positioned_matrix(
        ["A", "H", "B", "H", "A", "A", "B"],
        [1_000_000, 2_000_000, 2_000_000, 2_000_000, 3_000_000, 70_000_000, 70_000_000],
    )
    return chromosome_strips(progeny_states(gm), gm)


def test_co_located_markers_keep_zero_width_segments():
    strips = co_located_strips()
    s = strips[0]
    assert s.length_bp == 70_000_000
    assert [seg.state for seg in s.segments] == [STATE_A, STATE_H, STATE_B, STATE_H, STATE_A, STATE_B]
    inner, last = s.segments[2], s.segments[-1]
    assert inner.start == inner.end == pytest.approx(2e6 / 70e6)
    assert last.start == last.end == 1.0
    assert all(0.0 <= seg.start <= seg.end <= 1.0 for seg in s.segments)


def test_marker_beyond_assembly_extends_length():
    _, strips = strips_for(["A", "A"], spacing_bp=40_000_000)
    assert strips[0].length_bp == 80_000_000


def test_chromosomes_in_sort_key_order():
    gm = two_chrom_matrix()
    strips = chromosome_strips(progeny_states(gm), gm)
    assert [s.chrom for s in strips] == ["Gm2", "Gm11"]
    assert [s.n_markers for s in strips] == [1, 2]


def test_locus_ticks_region_and_marker():
    gm, strips = strips_for(["A", "H", "A"])
    target = resolve_locus(TargetSpec("T1", chrom="Gm01", start_bp=1_500_000, end_bp=3_000_000), gm)
    avoid = resolve_locus(AvoidSpec("AV1", marker_id="m2"), gm)
    ticks = locus_ticks({"T1": target}, {"AV1": avoid}, strips)
    assert [(t.locus_id, t.kind, t.chrom) for t in ticks] == [("T1", "target", "Gm01"), ("AV1", "avoid", "Gm01")]
    assert ticks[0].start == pytest.approx(1.5e6 / GM01_LEN, abs=1e-12)
    assert ticks[0].end == pytest.approx(3e6 / GM01_LEN, abs=1e-12)
    assert ticks[1].start == ticks[1].end == pytest.approx(2e6 / GM01_LEN, abs=1e-12)

    geom = strip_rects(strips, ticks)
    bar = 320 - 34
    marker_tick = geom.rows[0].ticks[1]
    assert marker_tick.w == pytest.approx(0.004 * bar)
    assert marker_tick.x + marker_tick.w / 2 == pytest.approx(34 + 2e6 / GM01_LEN * bar)
