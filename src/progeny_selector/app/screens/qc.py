"""Screen 2, Validate and QC: contract warnings, marker informativeness, per-individual QC.

Responsibility: summarise the loaded dataset and run (marker, sample and progeny
counts, informative markers, uninformative markers by reason from
``core.qc.uninformative_summary``, the map unit, background model, RPP unit,
drag unit, assembly, rank mode, duplicate pairs, warnings verbatim) and show the
``SampleQC`` table from ``core.qc.qc_table_rows`` as a sortable, filterable
DataGrid whose QC-excluded and flagged rows are tinted by
``app.present.qc_row_styles``.

Duplicate pairs (``AnalysisResult.duplicates``, docs/adr/0017 and its
2026-09-21 amendment) are listed as ``sample_a, sample_b: IBS 0.xxx``: the IBS
is measured over informative markers called in both individuals, and a pair
is reported only when the two overlap on half or more of the markers used for
that comparison, not over the whole panel.

Interface:
    view(id) -> Tag
    server(id, state) -> None   (reads state.dataset, state.result, state.criteria)
"""

from __future__ import annotations

from typing import Any

from shiny import module, render, ui

from progeny_selector.app.present import qc_row_styles
from progeny_selector.core.qc import qc_table_rows, uninformative_summary


@module.ui
def ui_(id: str = "") -> ui.Tag:
    return ui.layout_columns(
        ui.card(ui.card_header("Summary"), ui.output_ui("summary")),
        ui.card(ui.card_header("Per-individual QC"), ui.output_data_frame("table")),
        col_widths=(4, 8),
    )


def view(id: str) -> ui.Tag:
    return ui_(id)


def qc_frame(rows: list[dict]) -> Any:
    """The QC rows as a DataFrame in which a missing value stays None, so the grid shows a blank cell.

    Only columns holding a missing value become object dtype; the rest keep their numeric or bool type.
    """
    import pandas as pd  # UI-only dependency; the core never imports pandas

    frame = pd.DataFrame(rows)
    for col in frame.columns[frame.isna().any()]:
        frame[col] = frame[col].astype(object).where(frame[col].notna(), None)
    return frame


@module.server
def server_(input, output, session, state) -> None:
    @render.ui
    def summary() -> ui.TagChild:
        dataset = state.dataset()
        result = state.result()
        criteria = state.criteria()
        if dataset is None or result is None or criteria is None:
            return ui.p("Load data first.")
        reasons = uninformative_summary(result.classification)
        drag_unit = result.rows[0]["drag_unit"] if result.rows else "-"
        dup_items = (
            [ui.tags.li(f"{a}, {b}: IBS {ibs:.3f}") for a, b, ibs in result.duplicates] if result.duplicates else [ui.tags.li("none")]
        )
        items: list[ui.TagChild] = [
            ui.tags.dt("markers"),
            ui.tags.dd(str(dataset.genotypes.n_markers)),
            ui.tags.dt("samples"),
            ui.tags.dd(str(len(dataset.samples))),
            ui.tags.dt("progeny"),
            ui.tags.dd(str(len(dataset.progeny))),
            ui.tags.dt("informative markers"),
            ui.tags.dd(str(result.classification.n_informative)),
            ui.tags.dt("uninformative markers by reason"),
            *([ui.tags.dd(f"{reason}: {count}") for reason, count in reasons] or [ui.tags.dd("none")]),
            ui.tags.dt("map unit"),
            ui.tags.dd(result.unit),
            ui.tags.dt("background model"),
            ui.tags.dd(criteria.background.model),
            ui.tags.dt("RPP unit"),
            ui.tags.dd(result.unit),
            ui.tags.dt("drag unit"),
            ui.tags.dd(drag_unit),
            ui.tags.dt("assembly"),
            ui.tags.dd(criteria.assembly),
            ui.tags.dt("rank mode"),
            ui.tags.dd(criteria.ranking.mode),
            ui.tags.dt("duplicate pairs"),
            ui.tags.dd(ui.tags.ul(*dup_items)),
        ]
        # result.warnings already begins with the dataset's own warnings (core.pipeline.run_analysis).
        warnings = ui.tags.ul(*[ui.tags.li(w) for w in result.warnings]) if result.warnings else ui.p("No warnings.")
        return ui.TagList(ui.tags.dl(*items), ui.tags.h6("Warnings"), warnings)

    @render.data_frame
    def table():
        dataset = state.dataset()
        result = state.result()
        criteria = state.criteria()
        if dataset is None or result is None or criteria is None:
            return None
        rows = qc_table_rows(result.qc, dataset, criteria.filters)
        # present.py returns plain dicts shaped as shiny's StyleInfoBody, which it cannot import (Shiny-free).
        styles: Any = qc_row_styles(rows)
        return render.DataGrid(
            qc_frame(rows),
            filters=True,
            selection_mode="none",
            height="70vh",
            styles=styles,
        )


def server(id: str, state) -> None:
    server_(id, state)
