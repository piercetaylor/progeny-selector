"""Screen 6, Selection list: chosen individuals, editable notes and top-N helper.

Responsibility: maintain the selection list (from Rank multi-select or from
``select_top_n`` per family), free-text notes per individual edited in place
in the grid, and the projected next-generation expectation for the list.
Notes are keyed by ``sample_id`` in ``state.notes``, never by row, so a note
survives sorting the grid or changing the selection.

Interface:
    view(id) -> Tag
    server(id, state) -> None   (reads state.result, state.selected_ids; writes state.notes)
"""

from __future__ import annotations

from shiny import module, reactive, render, ui

from progeny_selector.core.selection import project_next_generation, select_top_n

NOTES_COLUMNS = ("sample_id", "line_name", "family_id", "rank_overall", "composite_score", "rpp_total", "notes")


@module.ui
def ui_(id: str = "") -> ui.Tag:
    return ui.card(
        ui.card_header("Selection list"),
        ui.layout_columns(
            ui.input_numeric("top_n", "Top N per family", value=5, min=1),
            ui.input_action_button("apply_top_n", "Add top N per family"),
            ui.input_radio_buttons("step", "Next step", choices={"backcross": "Backcross to RP", "self": "Self"}),
            col_widths=(3, 3, 6),
        ),
        ui.output_data_frame("table"),
        ui.output_text_verbatim("summary"),
    )


def view(id: str) -> ui.Tag:
    return ui_(id)


@module.server
def server_(input, output, session, state) -> None:
    @reactive.effect
    @reactive.event(input.apply_top_n)
    def _top_n() -> None:
        result = state.result()
        if result is None:
            return
        chosen = select_top_n(result.rows, int(input.top_n()), per_family=True)
        state.selected_ids.set(sorted(set(state.selected_ids()) | {r["sample_id"] for r in chosen}))

    @render.data_frame
    def table():
        import pandas as pd  # UI-only dependency; the core never imports pandas

        result = state.result()
        ids = state.selected_ids()
        notes = state.notes()
        if result is None or not ids:
            frame = pd.DataFrame(columns=list(NOTES_COLUMNS))
        else:
            records = []
            for sid in ids:
                row = result.row(sid)
                records.append({**{c: row.get(c) for c in NOTES_COLUMNS if c != "notes"}, "notes": notes.get(sid, "")})
            frame = pd.DataFrame(records, columns=list(NOTES_COLUMNS))
        return render.DataGrid(frame, editable=True, selection_mode="none", height="50vh")

    @table.set_patch_fn
    async def _(*, patch: render.CellPatch) -> render.CellValue:
        if patch["column_index"] != NOTES_COLUMNS.index("notes"):
            return table.data().iat[patch["row_index"], patch["column_index"]]
        sid = table.data().iat[patch["row_index"], 0]
        state.notes.set({**state.notes(), sid: str(patch["value"])})
        return patch["value"]

    @render.text
    def summary() -> str:
        result = state.result()
        ids = state.selected_ids()
        if result is None or not ids:
            return "No individuals selected."
        rows = [result.row(s) for s in ids]
        proj = project_next_generation(rows, input.step())
        return (
            f"{len(rows)} selected; mean RPP {proj.mean_rpp_selected:.3f}; expected RPP next {proj.expected_rpp_next:.3f}; "
            f"expected het next {proj.expected_het_next:.3f}; target carriers next {proj.expected_target_carrier_fraction:.2f}\n{proj.note}"
        )


def server(id: str, state) -> None:
    server_(id, state)
