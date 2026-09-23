"""Screen 4, Rank: ranked individuals table with status chips and filters.

Responsibility: render ``result.rows`` filtered by the breadcrumb (family /
generation, through ``core.navigation.filter_rows``) as a sortable, filterable
DataGrid with multi-row selection, per-locus status and recombinant columns,
and status cells coloured by ``app.present.status_cell_styles`` (Okabe-Ito,
with the status text kept in the cell). The rows selected in the grid's
current view feed the Compare and Selection-list screens.

Interface:
    view(id) -> Tag
    server(id, state) -> None   (reads state.result, state.breadcrumb; writes state.selected_ids)
"""

from __future__ import annotations

from typing import Any

from shiny import module, reactive, render, ui

from progeny_selector.app.present import status_cell_styles, status_columns
from progeny_selector.core.navigation import NavTree, build_tree, crumb_labels, filter_rows

DISPLAY_COLUMNS = (
    "rank_overall",
    "rank_in_family",
    "sample_id",
    "family_id",
    "generation",
    "passes_filters",
    "exclusion_reason",
    "composite_score",
    "rpp_total",
    "rpp_carrier",
    "rpp_noncarrier",
    "drag_total_est",
    "missing_rate",
    "qc_flags",
)


@module.ui
def ui_(id: str = "") -> ui.Tag:
    return ui.card(
        ui.card_header("Ranked individuals (multi-select rows, then open Compare)"),
        ui.output_text("caption"),
        ui.input_switch("only_pass", "Show only individuals passing hard filters", value=True),
        ui.output_data_frame("table"),
    )


def view(id: str) -> ui.Tag:
    return ui_(id)


@module.server
def server_(input, output, session, state) -> None:
    @reactive.calc
    def nav_tree() -> NavTree | None:
        dataset = state.dataset()
        return build_tree(dataset) if dataset else None

    @reactive.calc
    def rows() -> list[dict]:
        result = state.result()
        if result is None:
            return []
        crumb = state.breadcrumb()
        out = filter_rows(result.rows, crumb.get("family"), crumb.get("generation"))
        if input.only_pass():
            out = [r for r in out if r["passes_filters"]]
        return out

    @render.text
    def caption() -> str:
        nav = nav_tree()
        if nav is None:
            where = "no data loaded"
        else:
            crumb = state.breadcrumb()
            where = " > ".join(crumb_labels(nav, crumb.get("family"), crumb.get("generation")))
        return f"{len(rows())} individuals shown ({where}); select rows, then open Compare"

    @render.data_frame
    def table():
        import pandas as pd  # UI-only dependency; the core never imports pandas

        data = rows()
        result = state.result()
        # Columns come from the full result, so an empty view keeps the same columns as a non-empty one.
        keys = list(result.rows[0]) if result is not None and result.rows else []
        columns = list(DISPLAY_COLUMNS) + status_columns(keys) + [k for k in keys if k.startswith("recomb_")]
        frame = pd.DataFrame(data, columns=columns)
        # present.py returns plain dicts shaped as shiny's StyleInfoBody, which it cannot import (Shiny-free).
        styles: Any = status_cell_styles(list(frame.columns), data)
        return render.DataGrid(
            frame,
            selection_mode="rows",
            filters=True,
            height="70vh",
            styles=styles,
        )

    @reactive.effect
    def _publish_selection() -> None:
        view = table.data_view(selected=True)
        state.selected_ids.set(view["sample_id"].tolist() if "sample_id" in view.columns else [])


def server(id: str, state) -> None:
    server_(id, state)
