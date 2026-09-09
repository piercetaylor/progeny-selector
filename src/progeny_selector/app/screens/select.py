"""Screen 6, Selection list: chosen individuals with notes and top-N helper.

Responsibility: maintain the selection list (from Rank multi-select or from
``select_top_n`` per family), free-text notes per individual, and the
projected next-generation expectation for the list. Placeholder in M0.

Interface:
    view(id) -> Tag
    server(id, state) -> None   (reads state.result, state.selected_ids; writes state.notes)
"""

from __future__ import annotations

from shiny import module, reactive, render, ui

from progeny_selector.core.selection import project_next_generation, select_top_n


@module.ui
def ui_(id: str = "") -> ui.Tag:
    return ui.card(
        ui.card_header("Selection list"),
        ui.input_numeric("top_n", "Top N per family", value=5, min=1),
        ui.input_action_button("apply_top_n", "Add top N per family"),
        ui.input_radio_buttons("step", "Next step", choices={"backcross": "Backcross to RP", "self": "Self"}),
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
