"""Shiny application shell: top-level navigation across the seven screens.

Responsibility: assemble screen modules into a navset with breadcrumb state
(cross -> family -> generation -> individual) held in reactive values that
every screen reads. Status: M1 in progress; the Load screen runs the real
pipeline on uploaded files, the Navigate tree and breadcrumbs filter the Rank
grid, which shows per-locus status chips and recombinant flags; the remaining
screens are placeholders that state their interface.

Interface:
    app: shiny.App  (module attribute discovered by ``shiny run``)
    AppState: reactive container shared by screens (dataset, criteria, result, selection)
"""

from __future__ import annotations

from dataclasses import dataclass, field

from shiny import App, reactive, ui

from progeny_selector.app.screens import compare, export, load, navigate, qc, rank, select


@dataclass
class AppState:
    dataset: reactive.Value = field(default_factory=lambda: reactive.Value(None))
    criteria: reactive.Value = field(default_factory=lambda: reactive.Value(None))
    result: reactive.Value = field(default_factory=lambda: reactive.Value(None))
    breadcrumb: reactive.Value = field(default_factory=lambda: reactive.Value({"cross": None, "family": None, "generation": None}))
    selected_ids: reactive.Value = field(default_factory=lambda: reactive.Value([]))
    notes: reactive.Value = field(default_factory=lambda: reactive.Value({}))


SCREENS = (
    ("load", "1 Load", load),
    ("qc", "2 Validate and QC", qc),
    ("navigate", "3 Navigate", navigate),
    ("rank", "4 Rank", rank),
    ("compare", "5 Compare", compare),
    ("select", "6 Selection list", select),
    ("export", "7 Export", export),
)

app_ui = ui.page_navbar(
    *[ui.nav_panel(title, module.view(key), value=key) for key, title, module in SCREENS],
    title="progeny-selector",
    id="screen",
    fillable=True,
)


def server(input, output, session):
    state = AppState()
    for key, _title, module in SCREENS:
        module.server(key, state)


app = App(app_ui, server)
