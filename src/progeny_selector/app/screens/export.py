"""Screen 7, Export: results CSV, selected IDs, next-round sample manifest.

Responsibility: hand the three files produced by ``progeny_selector.io.export``
to the browser as downloads. Under shinylive the files are generated in the
tab; under ``shiny run`` on the server side. Placeholder in M0.

Interface:
    view(id) -> Tag
    server(id, state) -> None   (reads state.result, state.selected_ids, state.notes, state.dataset)
"""

from __future__ import annotations

import tempfile
from pathlib import Path

from shiny import module, render, ui

from progeny_selector.io.export import write_next_round_manifest, write_results_csv, write_selection_csv


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
    def _tmp(name: str) -> Path:
        return Path(tempfile.mkdtemp()) / name

    @render.download(filename="results.csv")
    def results():
        path = _tmp("results.csv")
        write_results_csv(state.result().rows if state.result() else [], path)
        return str(path)

    @render.download(filename="selected.csv")
    def selected():
        path = _tmp("selected.csv")
        result = state.result()
        rows = [result.row(s) for s in state.selected_ids()] if result else []
        write_selection_csv(rows, path, state.notes())
        return str(path)

    @render.download(filename="next_samples.csv")
    def manifest():
        path = _tmp("next_samples.csv")
        result, dataset = state.result(), state.dataset()
        rows = [result.row(s) for s in state.selected_ids()] if result else []
        if dataset:
            write_next_round_manifest(rows, path, input.next_gen(), dataset.recurrent_parent, dataset.donor_parent)
        return str(path)


def server(id: str, state) -> None:
    server_(id, state)
