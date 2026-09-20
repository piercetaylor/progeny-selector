"""Genome-strip segments and locus ticks on hand-built matrices."""

from __future__ import annotations

import numpy as np
import pytest

from progeny_selector.app.present import strip_rects
from progeny_selector.constants import STATE_A, STATE_B, STATE_H, STATE_N
from progeny_selector.core.classify import classify
from progeny_selector.core.foreground import resolve_locus
from progeny_selector.core.pipeline import run_analysis
from progeny_selector.core.strip import ChromStrip, chromosome_strips, locus_ticks
from progeny_selector.model.criteria import AvoidSpec, Criteria, TargetSpec
from progeny_selector.model.dataset import GenotypeMatrix, Marker
from tests.conftest import make_dataset, make_matrix

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


def test_view_marker_order_matches_pipeline_states():
    """Compare strips index result.classification.states with the view's own sorted_by_position() copy.

    run_analysis sorts a copy of the matrix before classifying; the view sorts the loaded matrix the same
    way. With markers out of order across chromosomes and three co-located calls, both orders must agree,
    so the strip equals the one built from a matrix constructed already in position order.
    """
    states = ["B", "A", "H", "B", "H", "A"]
    chroms = ["Gm11", "Gm2", "Gm11", "Gm2", "Gm11", "Gm11"]
    positions = [5_000_000, 3_000_000, 1_000_000, 1_000_000, 5_000_000, 5_000_000]
    base = make_matrix(states)
    markers = [Marker(f"m{i}", c, p) for i, (c, p) in enumerate(zip(chroms, positions, strict=True))]
    loaded = GenotypeMatrix(markers=markers, sample_ids=base.sample_ids, alleles=base.alleles, calls=base.calls)
    dataset = make_dataset(loaded)
    result = run_analysis(dataset, Criteria(targets=[TargetSpec("T1", marker_id="m3")], avoid=[AvoidSpec("AV1", marker_id="m0")]))

    gm = dataset.genotypes.sorted_by_position()
    j = gm.sample_index("P1")
    got = chromosome_strips(result.classification.states[:, j], gm)

    order = [3, 1, 2, 0, 4, 5]  # Gm2 by position, then Gm11 by position with ties in input order
    ordered = GenotypeMatrix(
        markers=[markers[i] for i in order],
        sample_ids=base.sample_ids,
        alleles=[base.alleles[i] for i in order],
        calls=base.calls[order],
    )
    want = chromosome_strips(progeny_states(ordered), ordered)
    assert got == want
    assert [[seg.state for seg in s.segments] for s in got] == [[STATE_B, STATE_A], [STATE_H, STATE_B, STATE_H, STATE_A]]


def test_strip_length_follows_the_assembly():
    """Gm11 markers at 35 Mb and 39 Mb: a4 places both inside the chromosome, a2 (34.8 Mb) does not."""
    gm = make_matrix(["A", "B"], chrom="Gm11", spacing_bp=1)
    gm = GenotypeMatrix(
        markers=[Marker(m.marker_id, m.chrom, bp, m.cm) for m, bp in zip(gm.markers, (35_000_000, 39_000_000), strict=True)],
        sample_ids=list(gm.sample_ids),
        alleles=gm.alleles,
        calls=gm.calls,
    )
    states = progeny_states(gm)
    assert chromosome_strips(states, gm, "Wm82.a4")[0].length_bp == 39_643_746
    assert chromosome_strips(states, gm, "Wm82.a2")[0].length_bp == 39_000_000  # stretched to the last marker
    assert chromosome_strips(states, gm, "none")[0].length_bp == 39_000_000
    assert chromosome_strips(states, gm)[0].length_bp == 39_643_746  # the default is Wm82.a4
