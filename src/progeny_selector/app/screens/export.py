"""Screen 7, Export: results CSV, selected IDs, next-round sample manifest.

Responsibility: hand the three files produced by ``progeny_selector.io.export``
to the browser as downloads. Each handler yields the file's text, so nothing is
written to a temporary directory; under shinylive the files are generated in
the tab, under ``shiny run`` on the server side.

Interface:
    view(id) -> Tag
    server(id, state) -> None   (reads state.result, state.selected_ids, state.notes, state.dataset)
"""

from __future__ import annotations

from shiny import module, render, ui

from progeny_selector.io.export import MANIFEST_COLUMNS, next_round_manifest_text, results_csv_text, selection_csv_text


@module.ui
def ui_(id: str = "") -> ui.Tag:
    return ui.card(
        ui.card_header("Export"),
        ui.input_text("next_gen", "Next generation label", value="BC3F1"),
        ui.download_button("results", "results.csv (ranks and statuses)"),
        ui.download_button("selected", "selected.csv"),
        ui.download_button("manifest", "next_samples.csv (manifest for the next genotyping round)"),
    )


def view(id: str) -> ui.Tag:
    return ui_(id)


@module.server
def server_(input, output, session, state) -> None:
    # Each handler yields the file's text, byte-identical to the CLI writers; nothing touches the filesystem.
    @render.download_button(filename="results.csv")
    def results():
        result = state.result()
        yield results_csv_text(result.rows if result else [])

    @render.download_button(filename="selected.csv")
    def selected():
        result = state.result()
        rows = [result.row(s) for s in state.selected_ids()] if result else []
        yield selection_csv_text(rows, state.notes())

    @render.download_button(filename="next_samples.csv")
    def manifest():
        result, dataset = state.result(), state.dataset()
        rows = [result.row(s) for s in state.selected_ids()] if result else []
        if dataset:
            yield next_round_manifest_text(rows, input.next_gen(), dataset.recurrent_parent, dataset.donor_parent)
        else:
            # No parents yet: the header row alone, with the writers' CRLF line ending.
            yield ",".join(MANIFEST_COLUMNS) + "\r\n"


def server(id: str, state) -> None:
    server_(id, state)
