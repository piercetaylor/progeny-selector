"""Compare and QC view helpers that need no browser: strip clamping and stacking, None cells in the QC grid."""

from __future__ import annotations

import pytest

from progeny_selector.app.present import Rect, StripRow, TickShape
from progeny_selector.app.screens.compare import LABEL_PX, _bar_end, _placed_rects, _tick_vertices
from progeny_selector.app.screens.qc import qc_frame

BLUE, GREEN, VERMILION = "#0072B2", "#009E73", "#D55E00"


def _row(rects: list[Rect], ticks: list[TickShape] | None = None) -> StripRow:
    return StripRow(chrom="Gm01", y=0.0, label_x=0.0, rects=rects, ticks=ticks or [])


@pytest.mark.parametrize("fraction", [0.0, 1.0])
def test_tick_at_chromosome_end_stays_on_bar(fraction: float) -> None:
    bar = 100.0
    row = _row([Rect(x=LABEL_PX, w=bar, fill=BLUE, state_label="A")])
    end = _bar_end(row, [False])
    assert end == LABEL_PX + bar
    x = LABEL_PX + fraction * bar
    for kind in ("target", "avoid"):
        for px, _ in _tick_vertices(kind, x - 0.2, x + 0.2, row.y, LABEL_PX, end):
            assert LABEL_PX <= px <= end


def test_only_zero_width_segments_stay_in_bounds() -> None:
    # A bar shorter than the rect floor: both positions are zero-width, one at the start and one at the end.
    rects = [
        Rect(x=LABEL_PX, w=1.0, fill=BLUE, state_label="A"),
        Rect(x=LABEL_PX, w=1.0, fill=GREEN, state_label="H"),
        Rect(x=LABEL_PX + 0.5, w=1.0, fill=VERMILION, state_label="B"),
    ]
    row = _row(rects)
    zero = [True, True, True]
    end = _bar_end(row, zero)
    placed = _placed_rects(row, zero, end)
    assert len(placed) == 3
    for rect, x in placed:
        assert x >= LABEL_PX
        assert x + rect.w <= end + 1e-9


def test_three_states_at_one_position_are_all_visible() -> None:
    bar = 200.0
    pos = LABEL_PX + 50.0
    rects = [
        Rect(x=LABEL_PX, w=50.0, fill=BLUE, state_label="A"),
        Rect(x=pos, w=1.0, fill=GREEN, state_label="H"),
        Rect(x=pos, w=1.0, fill=VERMILION, state_label="B"),
        Rect(x=pos, w=1.0, fill=GREEN, state_label="H"),
        Rect(x=pos, w=bar - 50.0, fill=BLUE, state_label="A"),
    ]
    zero = [False, True, True, True, False]
    row = _row(rects)
    end = _bar_end(row, zero)
    placed = _placed_rects(row, zero, end)
    # Draw order: the positive-width rects first, then the co-located stack.
    assert [r.state_label for r, _ in placed] == ["A", "A", "H", "B", "H"]
    stack = [x for _, x in placed[2:]]
    assert len(set(stack)) == 3
    assert all(x >= LABEL_PX and x + 1.0 <= end for x in stack)


def test_stack_at_bar_end_is_clamped_inside() -> None:
    rects = [
        Rect(x=LABEL_PX, w=100.0, fill=BLUE, state_label="A"),
        Rect(x=LABEL_PX + 100.0, w=1.0, fill=GREEN, state_label="H"),
        Rect(x=LABEL_PX + 100.0, w=1.0, fill=VERMILION, state_label="B"),
    ]
    zero = [False, True, True]
    row = _row(rects)
    end = _bar_end(row, zero)
    xs = [x for _, x in _placed_rects(row, zero, end)[1:]]
    assert xs == [end - 2.0, end - 1.0]


def test_qc_frame_keeps_none_blank() -> None:
    from shiny.render._data_frame_utils._tbl_data import serialize_frame

    rows = [
        {"sample_id": "P1", "generation": None, "expected_het": None, "het_rate": 0.3, "qc_excluded": False},
        {"sample_id": "P2", "generation": "BC2F1", "expected_het": 0.25, "het_rate": None, "qc_excluded": True},
    ]
    frame = qc_frame(rows)
    assert frame.loc[0, "expected_het"] is None
    assert frame.loc[1, "het_rate"] is None
    assert frame["qc_excluded"].dtype == bool
    data = serialize_frame(frame)["data"]
    assert data[0] == ["P1", None, None, 0.3, False]
    assert data[1] == ["P2", "BC2F1", 0.25, None, True]
