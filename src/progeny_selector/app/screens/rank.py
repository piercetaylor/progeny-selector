"""Screen 4, Rank: ranked individuals table with status chips and filters.

Responsibility: render ``result.rows`` filtered by the breadcrumb (family /
generation) as a sortable, filterable DataGrid with multi-row selection;
selected rows feed the Compare and Selection-list screens. Colours come from
``constants.STATUS_COLORS`` / ``STATE_COLORS`` (Okabe-Ito).

Interface:
    view(id) -> Tag
    server(id, state) -> None   (reads state.result, state.breadcrumb; writes state.selected_ids)
"""

from __future__ import annotations

from shiny import module, reactive, render, ui

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
        ui.input_switch("only_pass", "Show only individuals passing hard filters", value=True),
        ui.output_data_frame("table"),
    )


def view(id: str) -> ui.Tag:
    return ui_(id)


@module.server
def server_(input, output, session, state) -> None:
    @reactive.calc
    def rows() -> list[dict]:
        result = state.result()
        if result is None:
            return []
        crumb = state.breadcrumb()
        out = []
        for r in result.rows:
            if crumb.get("family") and r.get("family_id") != crumb["family"]:
                continue
            if crumb.get("generation") and r.get("generation") != crumb["generation"]:
                continue
            if input.only_pass() and not r["passes_filters"]:
                continue
            out.append(r)
        return out

    @render.data_frame
    def table():
        import pandas as pd  # UI-only dependency; the core never imports pandas

        data = rows()
        cols = [c for c in DISPLAY_COLUMNS if data and c in data[0]]
        extra = [c for c in (data[0] if data else {}) if c.startswith(("target_", "avoid_"))]
        frame = pd.DataFrame(data, columns=cols + extra) if data else pd.DataFrame(columns=list(DISPLAY_COLUMNS))
        return render.DataGrid(frame, selection_mode="rows", filters=True, height="70vh")

    @reactive.effect
    def _publish_selection() -> None:
        sel = table.cell_selection()
        idx = sel.get("rows", ()) if sel else ()
        data = rows()
        state.selected_ids.set([data[i]["sample_id"] for i in idx if i < len(data)])


def server(id: str, state) -> None:
    server_(id, state)
