"""Screen 3, Navigate: cross -> family -> generation tree with breadcrumbs.

Responsibility: show the tree from ``core.navigation.build_tree`` as an accordion
of families, each holding a radio group of its generations, so a node is picked
with mouse or keyboard (accordion headers are buttons, radios take arrow keys);
write the selection into ``state.breadcrumb``, which filters the Rank screen,
and render it as a Bootstrap breadcrumb whose ancestor crumbs step back up. The
radio offers every generation of the family and "(no generation)" when any
individual lacks one.

Interface:
    view(id) -> Tag
    server(id, state) -> None   (reads state.dataset; writes state.breadcrumb)
"""

from __future__ import annotations

from shiny import module, reactive, render, ui

from progeny_selector.app.present import ALL_GENERATIONS, decode_generation, generation_choices
from progeny_selector.core.navigation import UNASSIGNED, NavTree, build_tree, crumb_labels


@module.ui
def ui_(id: str = "") -> ui.Tag:
    return ui.card(
        ui.card_header("Navigate: cross > family > generation"),
        ui.output_ui("crumb"),
        ui.output_ui("tree"),
    )


def view(id: str) -> ui.Tag:
    return ui_(id)


@module.server
def server_(input, output, session, state) -> None:
    @reactive.calc
    def nav_tree() -> NavTree | None:
        dataset = state.dataset()
        return build_tree(dataset) if dataset else None

    def _open_index() -> int | None:
        """Index of the open family panel; None or an empty value means no family."""
        opened = input.families()
        if not opened:
            return None
        value = opened if isinstance(opened, str) else opened[0]
        return int(value.removeprefix("fam_"))

    # Count accordion reports from the browser. After a new tree is built, input.families
    # still holds the old DOM's value until the new accordion reports, so _publish ignores
    # the open panel (and its gen_i radio) until a report arrives for the new tree.
    families_reports = reactive.Value(0)
    tree_state: dict = {"tree": None, "reports_at_build": 0}

    @reactive.effect
    @reactive.event(input.families, ignore_none=False)
    def _count_families_report() -> None:
        with reactive.isolate():
            families_reports.set(families_reports() + 1)

    # Rendered while hidden too, so a new dataset rebuilds the accordion immediately.
    @output(suspend_when_hidden=False)
    @render.ui
    def tree() -> ui.TagChild:
        nav = nav_tree()
        if nav is None:
            return None
        panels = []
        for i, fam in enumerate(nav.families):
            label = fam.family_id if fam.family_id is not None else "(no family)"
            choices = generation_choices(fam)
            panels.append(
                ui.accordion_panel(
                    f"{label} ({fam.n})",
                    ui.input_radio_buttons(session.ns(f"gen_{i}"), "Generation", choices=choices, selected=ALL_GENERATIONS),
                    value=f"fam_{i}",
                )
            )
        return ui.accordion(*panels, id="families", multiple=False, open=False)

    @reactive.effect
    def _publish() -> None:
        nav = nav_tree()
        if nav is None:
            state.breadcrumb.set({"cross": None, "family": None, "generation": None})
            return
        reports = families_reports()
        if nav is not tree_state["tree"]:
            tree_state["tree"] = nav
            tree_state["reports_at_build"] = reports
        family = None
        generation = None
        i = _open_index() if reports > tree_state["reports_at_build"] else None
        if i is not None and i < len(nav.families):
            fam = nav.families[i]
            family = UNASSIGNED if fam.family_id is None else fam.family_id
            generation = decode_generation(input[f"gen_{i}"]())
        state.breadcrumb.set({"cross": nav.cross, "family": family, "generation": generation})

    @output(suspend_when_hidden=False)
    @render.ui
    def crumb() -> ui.TagChild:
        nav = nav_tree()
        if nav is None:
            return ui.p("No data loaded")
        c = state.breadcrumb()
        labels = crumb_labels(nav, c.get("family"), c.get("generation"))
        items = [
            ui.tags.li({"class": "breadcrumb-item"}, ui.input_action_link(session.ns(f"crumb_{depth}"), label))
            for depth, label in enumerate(labels[:-1])
        ]
        items.append(ui.tags.li({"class": "breadcrumb-item active", "aria-current": "page"}, labels[-1]))
        return ui.tags.nav({"aria-label": "breadcrumb"}, ui.tags.ol({"class": "breadcrumb"}, *items))

    @reactive.effect
    @reactive.event(input.crumb_0)
    def _to_cross() -> None:
        ui.update_accordion("families", show=False)

    @reactive.effect
    @reactive.event(input.crumb_1)
    def _to_family() -> None:
        i = _open_index()
        if i is not None:
            ui.update_radio_buttons(f"gen_{i}", selected=ALL_GENERATIONS)


def server(id: str, state) -> None:
    server_(id, state)
