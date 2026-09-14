"""Screen 1, Load: genotype file, samples.csv, optional markers.csv, criteria.yaml; criteria editor.

Responsibility: accept the four input files, parse them at the boundary via
``progeny_selector.io``, run the analysis, and publish dataset/criteria/result
into AppState. The criteria card is a YAML view over the applied ``Criteria``
in AppState: it is rewritten in canonical form after every successful Load or
Apply, Apply parses the text with the strict reader and re-analyses the loaded
dataset, and the download always serialises the applied criteria, never the
raw editor text. Errors from the data contract and the criteria reader are
shown verbatim; a failed Apply leaves AppState untouched.

Interface:
    view(id) -> Tag
    server(id, state: AppState) -> None   (writes state.dataset, state.criteria, state.result;
                                           Load resets state.breadcrumb and state.selected_ids on success;
                                           Apply clears state.selected_ids and keeps state.breadcrumb)
"""

from __future__ import annotations

from shiny import module, reactive, render, ui

from progeny_selector.core.navigation import build_tree
from progeny_selector.core.pipeline import run_analysis
from progeny_selector.io import load_dataset, read_criteria
from progeny_selector.io.criteria import dump_criteria_yaml, read_criteria_text
from progeny_selector.model.criteria import CriteriaError
from progeny_selector.model.dataset import DataContractError


@module.ui
def ui_(id: str = "") -> ui.Tag:
    return ui.div(
        ui.layout_columns(
            ui.card(
                ui.card_header("Input files"),
                ui.input_file("genotypes", "Genotypes (VCF, VCF.gz, HapMap, wide CSV)", accept=[".vcf", ".gz", ".txt", ".csv", ".hmp"]),
                ui.input_file("samples", "samples.csv", accept=[".csv"]),
                ui.input_file("markers", "markers.csv (optional)", accept=[".csv"]),
                ui.input_file("criteria", "criteria.yaml", accept=[".yaml", ".yml"]),
                ui.input_action_button("run", "Load and analyse", class_="btn-primary"),
            ),
            ui.card(ui.card_header("Status"), ui.output_text_verbatim("status")),
            col_widths=(5, 7),
        ),
        ui.card(
            ui.card_header("Criteria (editable)"),
            ui.input_text_area(
                "criteria_text",
                "criteria.yaml",
                rows=18,
                placeholder="Load a criteria.yaml above or paste one here",
                spellcheck="false",
            ),
            ui.input_action_button("apply", "Apply criteria and re-analyse"),
            ui.download_button("download_criteria", "Download criteria.yaml"),
            ui.output_text_verbatim("criteria_status"),
        ),
    )


def view(id: str) -> ui.Tag:
    return ui_(id)


@module.server
def server_(input, output, session, state) -> None:
    message = reactive.Value("Choose files and press Load.")
    criteria_message = reactive.Value("")

    @reactive.effect
    @reactive.event(input.run)
    def _run() -> None:
        try:
            g, s, c = input.genotypes(), input.samples(), input.criteria()
            if not (g and s and c):
                message.set("Genotypes, samples.csv and criteria.yaml are required.")
                return
            m = input.markers()
            dataset = load_dataset(g[0]["datapath"], s[0]["datapath"], m[0]["datapath"] if m else None)
            criteria = read_criteria(c[0]["datapath"])
            result = run_analysis(dataset, criteria)
            state.dataset.set(dataset)
            state.criteria.set(criteria)
            state.result.set(result)
            # A new dataset starts at the cross node with nothing selected.
            state.breadcrumb.set({"cross": build_tree(dataset).cross, "family": None, "generation": None})
            state.selected_ids.set([])
            ui.update_text_area("criteria_text", value=dump_criteria_yaml(criteria))
            criteria_message.set("criteria loaded; edit and Apply to re-analyse (YAML comments are not kept)")
            n_pass = sum(1 for r in result.rows if r["passes_filters"])
            lines = [f"{dataset.genotypes.n_markers} markers, {len(dataset.progeny)} progeny; {n_pass} pass hard filters"]
            lines += [f"warning: {w}" for w in result.warnings]
            message.set("\n".join(lines))
        except (DataContractError, CriteriaError) as exc:
            message.set(f"error: {exc}")

    @render.text
    def status() -> str:
        return message()

    @reactive.effect
    @reactive.event(input.apply)
    def _apply() -> None:
        dataset = state.dataset()
        if dataset is None:
            criteria_message.set("Load data first.")
            return
        try:
            criteria = read_criteria_text(input.criteria_text())
            result = run_analysis(dataset, criteria)
        except (DataContractError, CriteriaError) as exc:
            # Nothing in AppState has been written yet, so a failed Apply leaves the last good run in place.
            criteria_message.set(f"error: {exc}")
            return
        state.criteria.set(criteria)
        state.result.set(result)
        # The breadcrumb stays: the dataset, and so the tree, is unchanged.
        state.selected_ids.set([])
        ui.update_text_area("criteria_text", value=dump_criteria_yaml(criteria))
        n_pass = sum(1 for r in result.rows if r["passes_filters"])
        criteria_message.set(f"re-analysed: {n_pass} pass hard filters (YAML comments are not kept)")

    @render.text
    def criteria_status() -> str:
        return criteria_message()

    @render.download_button(filename="criteria.yaml", media_type="application/yaml")
    def download_criteria():
        # The applied Criteria, never the editor text, so the saved file is what the app ran.
        criteria = state.criteria()
        if criteria is None:
            yield "# No criteria applied yet. Load or apply criteria first.\n"
        else:
            yield dump_criteria_yaml(criteria)


def server(id: str, state) -> None:
    server_(id, state)
