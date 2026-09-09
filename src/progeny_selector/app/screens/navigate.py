"""Screen 3, Navigate: cross -> family -> generation tree with breadcrumbs.

Responsibility: derive the tree from the sample manifest (cross = RP x donor,
family_id, generation), let the user pick a node with mouse or arrow keys,
and write the selection into ``state.breadcrumb`` which filters the Rank screen.
Placeholder in M0: select inputs stand in for the tree control.

Interface:
    view(id) -> Tag
    server(id, state) -> None   (reads state.dataset; writes state.breadcrumb)
"""

from __future__ import annotations

from shiny import module, reactive, ui


@module.ui
def ui_(id: str = "") -> ui.Tag:
    return ui.card(
        ui.card_header("Navigate: cross > family > generation"),
        ui.input_select("family", "Family", choices=["(all)"]),
        ui.input_select("generation", "Generation", choices=["(all)"]),
        ui.output_ui("crumb"),
    )


def view(id: str) -> ui.Tag:
    return ui_(id)


@module.server
def server_(input, output, session, state) -> None:
    @reactive.effect
    def _fill_choices() -> None:
        dataset = state.dataset()
        if dataset is None:
            return
        fams = sorted({s.family_id for s in dataset.progeny if s.family_id})
        gens = sorted({s.generation for s in dataset.progeny if s.generation})
        ui.update_select("family", choices=["(all)", *fams])
        ui.update_select("generation", choices=["(all)", *gens])

    @reactive.effect
    def _publish() -> None:
        dataset = state.dataset()
        cross = f"{dataset.recurrent_parent.line_name} x {dataset.donor_parent.line_name}" if dataset else None
        fam = None if input.family() in (None, "(all)") else input.family()
        gen = None if input.generation() in (None, "(all)") else input.generation()
        state.breadcrumb.set({"cross": cross, "family": fam, "generation": gen})

    @__import__("shiny").render.ui
    def crumb():
        c = state.breadcrumb()
        parts = [p for p in (c.get("cross"), c.get("family"), c.get("generation")) if p]
        return ui.p(" > ".join(parts) if parts else "No data loaded")


def server(id: str, state) -> None:
    server_(id, state)
