"""Shiny-free presentation helpers: status styles, QC row styles, strip geometry."""

from __future__ import annotations

import importlib
import sys

import pytest

from progeny_selector.app import present
from progeny_selector.app.present import chip_style, qc_row_styles, status_cell_styles, status_columns, strip_rects
from progeny_selector.constants import STATE_COLORS, STATUS_COLORS
from progeny_selector.core.strip import chromosome_strips
from tests.test_strip import co_located_strips, progeny_states, strips_for, two_chrom_matrix


def test_status_columns():
    cols = ["sample_id", "target_T1_status", "avoid_AV1_status", "rpp_total", "target_status"]
    assert status_columns(cols) == ["target_T1_status", "avoid_AV1_status"]


def test_status_cell_styles():
    rows = [{"target_T1_status": "pass"}, {"target_T1_status": "fail"}, {"target_T1_status": "pass"}]
    styles = status_cell_styles(["sample_id", "target_T1_status"], rows)
    assert len(styles) == 2
    assert all(s["cols"] == [1] and s["location"] == "body" for s in styles)
    assert [s["rows"] for s in styles] == [[0, 2], [1]]
    assert all(s["style"]["background-color"] in STATUS_COLORS.values() for s in styles)
    assert [s["style"]["background-color"] for s in styles] == [STATUS_COLORS["pass"], STATUS_COLORS["fail"]]


def test_status_styles_without_avoid_columns():
    columns = ["sample_id", "target_T1_status", "recomb_T1_left", "target_T2_status"]
    rows = [
        {"sample_id": "a", "target_T1_status": "pass", "recomb_T1_left": True, "target_T2_status": "unknown"},
        {"sample_id": "b", "target_T1_status": "fail", "recomb_T1_left": False, "target_T2_status": "unknown"},
    ]
    assert status_columns(columns) == ["target_T1_status", "target_T2_status"]
    styles = status_cell_styles(columns, rows)
    assert [(s["cols"], s["rows"]) for s in styles] == [([1], [0]), ([1], [1]), ([3], [0, 1])]
    assert all(columns[s["cols"][0]].startswith("target_") for s in styles)


def test_qc_row_styles():
    rows = [
        {"flags": "possible_outcross", "qc_excluded": True},
        {"flags": "high_missing", "qc_excluded": False},
        {"flags": "", "qc_excluded": False},
    ]
    styles = qc_row_styles(rows)
    assert [s["rows"] for s in styles] == [[0], [1]]
    assert styles[0]["style"]["background-color"] == "rgba(213, 94, 0, 0.25)"  # STATUS_COLORS["fail"] #D55E00
    assert styles[1]["style"]["background-color"] == "rgba(240, 228, 66, 0.35)"  # Okabe-Ito yellow #F0E442
    assert qc_row_styles([{"flags": "", "qc_excluded": False}]) == []


def test_chip_style():
    css = chip_style("unknown")
    assert f"background-color: {STATUS_COLORS['unknown']};" in css
    assert "color: #000000;" in css


@pytest.mark.parametrize("states", [["A", "H", "A"], ["N", "A"], ["A", "A", "A"]])
def test_strip_rects_widths_sum_to_bar(states):
    _, strips = strips_for(states)
    geom = strip_rects(strips, [])
    assert len(geom.rows) == len(strips)
    assert_widths_fit_bar(strips[0], geom.rows[0], 320 - 34)
    assert all(r.fill == STATE_COLORS[r.state_label] for r in geom.rows[0].rects)


def assert_widths_fit_bar(strip, row, bar: float, min_rect_px: float = 1.0) -> None:
    """Rect widths sum to the bar length, plus at most one floor per zero-width segment."""
    n_zero_width = sum(seg.start == seg.end for seg in strip.segments)
    total = sum(r.w for r in row.rects)
    assert bar - 0.5 <= total <= bar + n_zero_width * min_rect_px + 0.5


def test_strip_rects_draws_co_located_states():
    strips = co_located_strips()
    geom = strip_rects(strips, [])
    rects = geom.rows[0].rects
    assert [r.state_label for r in rects] == ["A", "H", "B", "H", "A", "B"]
    assert all(r.w >= 1.0 for r in rects)
    assert rects[2].w == 1.0 and rects[-1].w == 1.0  # inner and chromosome-end zero-width runs
    assert_widths_fit_bar(strips[0], geom.rows[0], 320 - 34)


def test_strip_rects_scales_to_longest_chromosome():
    gm = two_chrom_matrix()
    strips = chromosome_strips(progeny_states(gm), gm)
    geom = strip_rects(strips, [])
    assert [row.chrom for row in geom.rows] == ["Gm2", "Gm11"]
    bar_max = 320 - 34
    longest = max(s.length_bp for s in strips)
    for strip, row in zip(strips, geom.rows, strict=True):
        assert_widths_fit_bar(strip, row, bar_max * strip.length_bp / longest)
    assert geom.rows[1].y > geom.rows[0].y
    assert geom.height == 2 * 8 + 3


def test_import_without_shiny(monkeypatch):
    for name in ("shiny", "pandas", "htmltools"):
        monkeypatch.setitem(sys.modules, name, None)
    reloaded = importlib.reload(present)
    assert reloaded.strip_rects is not None
