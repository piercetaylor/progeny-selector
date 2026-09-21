"""Screen 7, Export: results CSV, selected IDs, next-round sample manifest.

Responsibility: hand the three files produced by ``progeny_selector.io.export``
to the browser as downloads. Each handler yields the file's text, so nothing is
written to a temporary directory; under shinylive the files are generated in
the tab, under ``shiny run`` on the server side.

The manifest's placeholder-row count is checked by ``io.export``, not here: the screen passes the
field through, and on refusal shows the writer's own message in the status line and produces no
file, so a count the CLI rejects never becomes a silently parent-only manifest.

Interface:
    view(id) -> Tag
    server(id, state) -> None   (reads state.result, state.selected_ids, state.notes, state.dataset)
"""

from __future__ import annotations

from shiny import module, reactive, render, req, ui

from progeny_selector.io.export import (
    MANIFEST_COLUMNS,
    check_per_selected,
    next_round_manifest_text,
    results_csv_text,
    selection_csv_text,
)
from progeny_selector.model.dataset import DataContractError


@module.ui
def ui_(id: str = "") -> ui.Tag:
    return ui.card(
        ui.card_header("Export"),
        ui.input_text("next_gen", "Next generation label", value="BC3F1"),
        ui.input_numeric("per_selected", "Placeholder rows per selected individual", value=1, min=1, max=999),
        ui.download_button("results", "results.csv (ranks and statuses)"),
        ui.download_button("selected", "selected.csv"),
        ui.download_button("manifest", "next_samples.csv (manifest for the next genotyping round)"),
        ui.output_text_verbatim("status"),
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

    message = reactive.Value("")

    @reactive.effect
    def _check_count() -> None:
        # io/export owns the rule and its wording; the screen passes the field through and shows the
        # message. A reactive.Value set inside a download handler never reaches the client — that
        # request is outside the session's flush cycle — so the check runs here, on the input.
        per_selected = input.per_selected()
        if per_selected is None:
            message.set("")
            return
        try:
            check_per_selected(int(per_selected))
        except DataContractError as exc:
            message.set(f"error: {exc}")
        else:
            message.set("")

    @render.text
    def status() -> str:
        return message()

    @render.download_button(filename="next_samples.csv")
    def manifest():
        result, dataset = state.result(), state.dataset()
        rows = [result.row(s) for s in state.selected_ids()] if result else []
        if dataset:
            # A cleared numeric field reads as None — the min/max on the widget are client-side only.
            # The writer's own default then applies; no count is decided in the screen.
            per_selected = input.per_selected()
            rp, donor = dataset.recurrent_parent, dataset.donor_parent
            try:
                if per_selected is None:
                    text = next_round_manifest_text(rows, input.next_gen(), rp, donor)
                else:
                    text = next_round_manifest_text(rows, input.next_gen(), rp, donor, n_per_selected=int(per_selected))
            except DataContractError:
                # Refused by io/export: no file at all, rather than a manifest with no progeny in it.
                # The status line already carries the writer's message (the effect above).
                req(False)
                return
            yield text
        else:
            # No parents yet: the header row alone, with the writers' CRLF line ending.
            yield ",".join(MANIFEST_COLUMNS) + "\r\n"


def server(id: str, state) -> None:
    server_(id, state)
