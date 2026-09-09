"""Screen 5, Compare: side-by-side statuses and compact chromosome strips for selected individuals.

Responsibility: for ``state.selected_ids`` show one column per individual with
foreground/avoid chips, RPP per chromosome, drag bounds, recombinant flags,
and a compact 20-row chromosome strip drawn from ``result.classification``
states in Okabe-Ito colours. The strip is a summary, not a browser: the
sibling isoline-browser project owns detailed graphical-genotype viewing.
Placeholder in M0; interface fixed.

Interface:
    view(id) -> Tag
    server(id, state) -> None   (reads state.result, state.selected_ids)
"""

from __future__ import annotations

from shiny import module, render, ui


@module.ui
def ui_(id: str = "") -> ui.Tag:
    return ui.card(ui.card_header("Compare selected individuals"), ui.output_text_verbatim("panel"))


def view(id: str) -> ui.Tag:
    return ui_(id)


@module.server
def server_(input, output, session, state) -> None:
    @render.text
    def panel() -> str:
        result = state.result()
        ids = state.selected_ids()
        if result is None or not ids:
            return "Select rows on the Rank screen."
        lines = []
        for sid in ids:
            r = result.row(sid)
            drag = f"{r['drag_total_est']} {r['drag_unit']}"
            lines.append(f"{sid}: score {r['composite_score']}, RPP {r['rpp_total']}, drag {drag}, flags {r['qc_flags'] or '-'}")
        return "\n".join(lines)


def server(id: str, state) -> None:
    server_(id, state)
