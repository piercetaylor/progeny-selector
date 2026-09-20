"""Screen 5, Compare: side-by-side statuses and compact chromosome strips for selected individuals.

Responsibility: for the first six of ``state.selected_ids`` show one card per
individual with target and avoid status chips, RPP per chromosome (carrier
chromosomes marked), drag bounds and recombinant flags per target, and a
compact chromosome strip drawn as inline SVG from ``result.classification``
states in Okabe-Ito colours (``core.strip`` segments, ``app.present.strip_rects``
geometry). The strip is a summary, not a browser: the sibling backcross
project owns detailed graphical-genotype viewing. Nothing here computes a
metric; every number is read from ``result.rows``.

Interface:
    view(id) -> Tag
    server(id, state) -> None   (reads state.dataset, state.result, state.selected_ids)
"""

from __future__ import annotations

from htmltools import Tag
from shiny import module, render, ui

from progeny_selector.app.present import Rect, StripGeometry, StripRow, chip_style, strip_rects
from progeny_selector.constants import PALETTE_OKABE_ITO, STATE_COLORS
from progeny_selector.core.drag import DragResult
from progeny_selector.core.strip import chromosome_strips, locus_ticks

MAX_CARDS = 6
ROW_PX = 8
LABEL_PX = 34
TICK_COLOUR = PALETTE_OKABE_ITO["black"]
STATE_MEANINGS = (
    ("A", "homozygous recurrent parent"),
    ("H", "heterozygous"),
    ("B", "homozygous donor"),
    ("X", "non-parental allele"),
    ("N", "missing call"),
    ("U", "uninformative marker"),
)
RPP_SUMMARY_KEYS = frozenset({"rpp_total", "rpp_carrier", "rpp_noncarrier"})


@module.ui
def ui_(id: str = "") -> ui.Tag:
    return ui.card(ui.card_header("Compare selected individuals"), ui.output_ui("legend"), ui.output_ui("panels"))


def view(id: str) -> ui.Tag:
    return ui_(id)


def _fmt(value: float | None, spec: str = ".3f") -> str:
    return "-" if value is None else format(value, spec)


def _num(value: float) -> str:
    return f"{value:.2f}"


def _clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(x, hi))


def _tick_vertices(kind: str, x0: float, x1: float, y: float, lo: float, hi: float) -> list[tuple[float, float]]:
    """Target: downward triangle in the top half of the row; avoid: upward triangle in the bottom half.

    A single-marker locus is about 1 px wide, so the glyph is at least one row height wide; every
    vertex is clamped to [lo, hi] so a locus at either chromosome end stays on its bar.
    """
    cx = (x0 + x1) / 2.0
    left, right = min(x0, cx - ROW_PX / 2.0), max(x1, cx + ROW_PX / 2.0)
    mid = y + ROW_PX / 2.0
    edge = y if kind == "target" else y + ROW_PX
    return [(_clamp(left, lo, hi), edge), (_clamp(right, lo, hi), edge), (_clamp(cx, lo, hi), mid)]


def _tick_points(kind: str, x0: float, x1: float, y: float, lo: float, hi: float) -> str:
    return " ".join(f"{_num(px)},{_num(py)}" for px, py in _tick_vertices(kind, x0, x1, y, lo, hi))


def _bar_end(row: StripRow, zero_width: list[bool]) -> float:
    """Right edge of the bar: the end of the last positive-width rect (the last segment ends at 1.0).

    A row with only zero-width segments has no such rect; its bar is taken as wide enough to hold
    the widest stack of co-located rects, starting at its right-most position.
    """
    positive = [r.x + r.w for r, z in zip(row.rects, zero_width, strict=True) if not z]
    if positive:
        return max(positive)
    stacks: dict[float, float] = {}
    for r in row.rects:
        stacks[round(r.x, 6)] = stacks.get(round(r.x, 6), 0.0) + r.w
    right_most = max((x + w for x, w in stacks.items()), default=float(LABEL_PX))
    return max(right_most, LABEL_PX + max(stacks.values(), default=0.0))


def _placed_rects(row: StripRow, zero_width: list[bool], bar_end: float) -> list[tuple[Rect, float]]:
    """Rects with their drawn x, in draw order: positive-width first, then zero-width ones.

    Zero-width rects at one position form a stack: the stack is centred on the position, clamped
    inside [LABEL_PX, bar_end], and each rect sits one width to the right of the previous one, so
    every co-located state stays visible.
    """
    placed = [(r, r.x) for r, z in zip(row.rects, zero_width, strict=True) if not z]
    stacks: dict[float, list[Rect]] = {}
    for r, z in zip(row.rects, zero_width, strict=True):
        if z:
            stacks.setdefault(round(r.x, 6), []).append(r)
    for rects in stacks.values():
        total = sum(r.w for r in rects)
        x = max(LABEL_PX, min(rects[0].x - total / 2.0, bar_end - total))
        for r in rects:
            placed.append((r, x))
            x += r.w
    return placed


def _tick_glyph(kind: str) -> Tag:
    """A legend glyph drawn with the same polygon as the strip ticks."""
    return ui.tags.svg(
        {"viewBox": f"0 0 {ROW_PX + 2} {ROW_PX}", "role": "img", "aria-label": f"{kind} locus tick", "width": "14", "height": "11"},
        ui.tags.title(f"{kind} locus tick"),
        Tag("polygon", {"points": _tick_points(kind, 1.0, 1.0 + ROW_PX, 0.0, 0.0, ROW_PX + 2.0), "fill": TICK_COLOUR}),
    )


def _legend() -> Tag:
    swatches = [
        ui.span(
            f"{label} {meaning}",
            style=f"background-color: {STATE_COLORS[label]}; color: #000000; padding: 0 0.4em; border-radius: 0.25em; margin-right: 0.3em;",
        )
        for label, meaning in STATE_MEANINGS
    ]
    ticks = [ui.span(_tick_glyph(kind), f" {kind} locus", style="margin-right: 0.6em;") for kind in ("target", "avoid")]
    return ui.div(ui.div(*swatches, class_="mb-1"), ui.div(*ticks), class_="mb-2")


def _strip_row(row: StripRow, zero_width: list[bool]) -> Tag:
    """One chromosome: positive-width rects, then zero-width ones centred and clamped, then ticks on top."""
    bar_end = _bar_end(row, zero_width)
    shapes: list[Tag] = []
    for rect, x in _placed_rects(row, zero_width, bar_end):
        shapes.append(
            Tag(
                "rect",
                {"x": _num(x), "y": _num(row.y), "width": _num(rect.w), "height": str(ROW_PX), "fill": rect.fill},
                ui.tags.title(f"{row.chrom} {rect.state_label}"),
            )
        )
    for tick in row.ticks:
        shapes.append(
            Tag(
                "polygon",
                {"points": _tick_points(tick.kind, tick.x, tick.x + tick.w, row.y, LABEL_PX, bar_end), "fill": TICK_COLOUR},
                ui.tags.title(f"{tick.kind} locus {tick.locus_id}"),
            )
        )
    label = Tag("text", {"x": _num(row.label_x), "y": _num(row.y + ROW_PX), "font-size": str(ROW_PX)}, row.chrom)
    return Tag("g", {"data-chrom": row.chrom}, label, *shapes)


def _strip_svg(sid: str, geometry: StripGeometry, zero_width: list[list[bool]]) -> Tag:
    return ui.tags.svg(
        {
            "viewBox": f"0 0 {_num(geometry.width)} {_num(geometry.height)}",
            "role": "img",
            "aria-label": f"chromosome strip for {sid}",
            "style": "width: 100%; height: auto; display: block;",
        },
        ui.tags.title(f"Chromosome strip for {sid}: parent-of-origin state along each chromosome"),
        *[_strip_row(row, z) for row, z in zip(geometry.rows, zero_width, strict=True)],
    )


def _chips(row: dict, result) -> Tag:
    chips = []
    for kind, loci in (("target", result.resolved_targets), ("avoid", result.resolved_avoid)):
        for locus_id in loci:
            status = row[f"{kind}_{locus_id}_status"]
            chips.append(ui.div(f"{kind} {locus_id} ", ui.span(status, class_="status-chip", style=chip_style(status))))
    return ui.div(*chips, class_="mb-2")


def _rpp_table(row: dict, carrier: set[str]) -> Tag:
    body = []
    for key, value in row.items():
        if not key.startswith("rpp_") or key in RPP_SUMMARY_KEYS:
            continue
        chrom = key.removeprefix("rpp_")
        name: ui.TagChild = ui.tags.strong(f"{chrom} (carrier)") if chrom in carrier else chrom
        body.append(ui.tags.tr(ui.tags.td(name), ui.tags.td(_fmt(value))))
    return ui.tags.table(
        {"class": "table table-sm mb-2"},
        ui.tags.thead(ui.tags.tr(ui.tags.th("chromosome"), ui.tags.th("RPP"))),
        ui.tags.tbody(*body),
    )


def _drag_table(drag: DragResult, i: int, locus_id: str) -> Tag:
    """Donor-segment bounds for one target and one individual, read from ``result.drag`` (index ``i`` in progeny order)."""

    def length(values) -> str:
        v = float(values[i])
        return "-" if v != v else f"{v:.2f} {drag.unit}"  # NaN when the bound is undefined

    lines = [
        ("left min", length(drag.left_min)),
        ("left max", length(drag.left_max)),
        ("right min", length(drag.right_min)),
        ("right max", length(drag.right_max)),
        ("total estimate", length(drag.total_est)),
        ("total max", length(drag.total_max)),
        ("unit", drag.unit),
        ("recombinant left", "yes" if drag.recombinant_left[i] else "no"),
        ("recombinant right", "yes" if drag.recombinant_right[i] else "no"),
    ]
    return ui.tags.table(
        {"class": "table table-sm mb-2", "data-locus": locus_id},
        ui.tags.caption(f"drag at {locus_id}", style="caption-side: top;"),
        ui.tags.tbody(*[ui.tags.tr(ui.tags.td(k), ui.tags.td(v)) for k, v in lines]),
    )


@module.server
def server_(input, output, session, state) -> None:
    @render.ui
    def legend() -> ui.TagChild:
        return _legend()

    @render.ui
    def panels() -> ui.TagChild:
        dataset = state.dataset()
        result = state.result()
        criteria = state.criteria()
        ids = state.selected_ids()
        if dataset is None or result is None or criteria is None or not ids:
            return ui.p("Select rows on the Rank screen.")
        by_id = {r["sample_id"]: r for r in result.rows}
        # An id left over from an earlier result is skipped rather than raised on, and counted.
        progeny_ids = [s.sample_id for s in dataset.progeny]  # result.drag arrays follow this order
        known = set(progeny_ids)
        present = [sid for sid in ids if sid in by_id and sid in known]
        missing = len(ids) - len(present)
        notice = ui.p(f"{missing} selected individuals are not in the current results") if missing else None
        if not present:
            return ui.TagList(notice, ui.p("Select rows on the Rank screen."))
        gm = dataset.genotypes.sorted_by_position()  # run_analysis sorted its own copy the same way
        states = result.classification.states
        drag_ok = all(len(d.left_max) == len(progeny_ids) for d in result.drag.values())
        if states.shape != (gm.n_markers, gm.n_samples) or not drag_ok:
            return ui.p("The loaded data changed; run the analysis again.")
        shown = present[:MAX_CARDS]
        cards = []
        for sid in shown:
            row = by_id[sid]
            j = gm.sample_index(sid)
            strips = chromosome_strips(states[:, j], gm, criteria.assembly)
            ticks = locus_ticks(result.resolved_targets, result.resolved_avoid, strips)
            geometry = strip_rects(strips, ticks, row_px=ROW_PX, label_px=LABEL_PX)
            zero_width = [[seg.start == seg.end for seg in strip.segments] for strip in strips]
            cards.append(
                ui.card(
                    ui.card_header(f"{sid} ({row['line_name']})"),
                    ui.p(
                        f"rank {_fmt(row['rank_overall'], '.0f')} overall, {_fmt(row['rank_in_family'], '.0f')} in family; "
                        f"score {_fmt(row['composite_score'])}"
                    ),
                    _chips(row, result),
                    _rpp_table(row, result.carrier_chroms),
                    *[_drag_table(result.drag[locus_id], progeny_ids.index(sid), locus_id) for locus_id in result.resolved_targets],
                    _strip_svg(sid, geometry, zero_width),
                )
            )
        children: list[ui.TagChild] = [notice]
        if len(present) > MAX_CARDS:
            children.append(ui.p(f"showing {MAX_CARDS} of {len(present)}"))
        children.append(ui.layout_columns(*cards, col_widths=[12 // len(cards)]))
        return ui.TagList(*children)


def server(id: str, state) -> None:
    server_(id, state)
