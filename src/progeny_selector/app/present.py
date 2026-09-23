"""Shiny-free presentation helpers: status and QC styles, badge CSS, genome-strip geometry.

Responsibility: turn result rows and strip segments into style dicts and pixel
rectangles that the screens hand to Shiny. Colours come only from
``constants.STATE_COLORS`` and ``constants.STATUS_COLORS``; no metric is computed
here. Imports only ``constants`` and ``core`` so the helpers are unit-tested
without Shiny, pandas or a browser (docs/m1-phases.md, invariant 1).

Interface:
    STATUS_COLUMN_RE
    status_columns(columns) -> list[str]
    status_cell_styles(columns, rows) -> list[dict]   styles for the Rank grid's per-locus status cells
    qc_row_styles(rows) -> list[dict]                 row backgrounds for the QC table (rows from core.qc.qc_table_rows)
    chip_style(status) -> str                         inline CSS for a status badge
    strip_rects(strips, ticks, width_px=320, row_px=8, gap_px=3, label_px=34, min_rect_px=1.0) -> StripGeometry
    StripGeometry, StripRow, Rect, TickShape (frozen dataclasses)
    ALL_GENERATIONS, NO_GENERATION, GENERATION_PREFIX
    encode_generation(generation) -> str              Navigate radio value for a generation node
    decode_generation(value) -> str | None            radio value -> breadcrumb generation selector
    generation_choices(fam) -> dict[str, str]          Navigate radio choices for a FamilyNode
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass

from progeny_selector.constants import STATE_COLORS, STATE_LABELS, STATUS_COLORS
from progeny_selector.core.navigation import UNASSIGNED, FamilyNode
from progeny_selector.core.strip import ChromStrip, LocusTick

STATUS_COLUMN_RE = re.compile(r"^(target|avoid)_.+_status$")

# The Navigate radio's "all generations" value is "" today and core.navigation.UNASSIGNED is
# also "", so the two cannot share the value space: real generations are prefixed, sentinels
# are bare words. A stripped manifest cell can never collide with a prefixed value, and a
# generation literally named "all" encodes as "g:all".
ALL_GENERATIONS = "all"
NO_GENERATION = "none"
GENERATION_PREFIX = "g:"


def encode_generation(generation: str | None) -> str:
    """Navigate radio value for a generation node: NO_GENERATION for None, else prefixed."""
    if generation is None:
        return NO_GENERATION
    return GENERATION_PREFIX + generation


def decode_generation(value: str | None) -> str | None:
    """Radio value -> breadcrumb generation selector (None means no filter)."""
    if value is None:
        return None
    if value == "" or value == ALL_GENERATIONS:
        return None
    if value == NO_GENERATION:
        return UNASSIGNED
    if value.startswith(GENERATION_PREFIX):
        return value[len(GENERATION_PREFIX) :]
    return None


def generation_choices(fam: FamilyNode) -> dict[str, str]:
    """Navigate radio choices for a family: all generations, then each named generation, then (no generation)."""
    choices = {ALL_GENERATIONS: f"(all generations) ({fam.n})"}
    for g in fam.generations:
        if g.generation is not None:
            choices[encode_generation(g.generation)] = f"{g.generation} ({g.n})"
    for g in fam.generations:
        if g.generation is None:
            choices[NO_GENERATION] = f"(no generation) ({g.n})"
    return choices


# A zero-width (single-marker) locus tick is drawn this wide, as a fraction of its bar.
MIN_TICK_FRACTION = 0.004


@dataclass(frozen=True)
class Rect:
    x: float
    w: float
    fill: str
    state_label: str


@dataclass(frozen=True)
class TickShape:
    locus_id: str
    kind: str  # 'target' | 'avoid'
    x: float
    w: float


@dataclass(frozen=True)
class StripRow:
    chrom: str
    y: float
    label_x: float
    rects: list[Rect]
    ticks: list[TickShape]


@dataclass(frozen=True)
class StripGeometry:
    width: float
    height: float
    rows: list[StripRow]


def status_columns(columns: Sequence[str]) -> list[str]:
    """Per-locus status columns (target_*_status, avoid_*_status) in their incoming order."""
    return [c for c in columns if STATUS_COLUMN_RE.match(c)]


def status_cell_styles(columns: Sequence[str], rows: list[dict]) -> list[dict]:
    """One style entry per (status column, status value) present, values in STATUS_COLORS order."""
    styles: list[dict] = []
    for j, col in enumerate(columns):
        if not STATUS_COLUMN_RE.match(col):
            continue
        for value, colour in STATUS_COLORS.items():
            hits = [i for i, r in enumerate(rows) if r.get(col) == value]
            if hits:
                styles.append(
                    {
                        "location": "body",
                        "rows": hits,
                        "cols": [j],
                        "style": {"background-color": colour, "color": "#000000", "font-weight": "600"},
                    }
                )
    return styles


def _rgba(hex_colour: str, alpha: float) -> str:
    h = hex_colour.lstrip("#")
    r, g, b = (int(h[k : k + 2], 16) for k in (0, 2, 4))
    return f"rgba({r}, {g}, {b}, {alpha})"


def qc_row_styles(rows: list[dict]) -> list[dict]:
    """QC-excluded rows tinted vermilion at 25 %; flagged rows that are not QC-excluded tinted yellow at 35 %."""
    excluded = [i for i, r in enumerate(rows) if r.get("qc_excluded")]
    flagged = [i for i, r in enumerate(rows) if r.get("flags") and not r.get("qc_excluded")]
    styles: list[dict] = []
    if excluded:
        styles.append({"location": "body", "rows": excluded, "style": {"background-color": _rgba(STATUS_COLORS["fail"], 0.25)}})
    if flagged:
        # Okabe-Ito yellow; STATE_COLORS["U"] is the only palette entry carrying it.
        styles.append({"location": "body", "rows": flagged, "style": {"background-color": _rgba(STATE_COLORS["U"], 0.35)}})
    return styles


def chip_style(status: str) -> str:
    """Inline CSS for a status badge: black text on the status colour."""
    return f"background-color: {STATUS_COLORS[status]}; color: #000000; font-weight: 600; padding: 0 0.4em; border-radius: 0.25em;"


def strip_rects(
    strips: list[ChromStrip],
    ticks: list[LocusTick],
    width_px: int = 320,
    row_px: int = 8,
    gap_px: int = 3,
    label_px: int = 34,
    min_rect_px: float = 1.0,
) -> StripGeometry:
    """Pixel geometry: one row per strip, a label column of ``label_px`` then a bar scaled to the longest chromosome.

    Every rect is at least ``min_rect_px`` wide so a zero-width run (co-located markers with
    different states) stays visible; drawn widths may exceed the bar length only by those floors.
    """
    max_len = max((s.length_bp for s in strips), default=0.0)
    bar_max = float(width_px - label_px)
    rows: list[StripRow] = []
    for k, strip in enumerate(strips):
        bar = bar_max * strip.length_bp / max_len if max_len > 0 else 0.0
        rects = [
            Rect(
                x=label_px + seg.start * bar,
                w=max((seg.end - seg.start) * bar, min_rect_px),
                fill=STATE_COLORS[STATE_LABELS[seg.state]],
                state_label=STATE_LABELS[seg.state],
            )
            for seg in strip.segments
        ]
        shapes: list[TickShape] = []
        for t in ticks:
            if t.chrom != strip.chrom:
                continue
            frac = max(t.end - t.start, MIN_TICK_FRACTION)
            centre = (t.start + t.end) / 2.0
            shapes.append(TickShape(locus_id=t.locus_id, kind=t.kind, x=label_px + (centre - frac / 2.0) * bar, w=frac * bar))
        rows.append(StripRow(chrom=strip.chrom, y=float(k * (row_px + gap_px)), label_x=0.0, rects=rects, ticks=shapes))
    height = len(strips) * row_px + max(len(strips) - 1, 0) * gap_px
    return StripGeometry(width=float(width_px), height=float(height), rows=rows)
