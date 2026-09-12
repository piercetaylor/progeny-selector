"""Screen 1, Load: genotype file, samples.csv, optional markers.csv, criteria.yaml.

Responsibility: accept the four input files, parse them at the boundary via
``progeny_selector.io``, run the analysis, and publish dataset/criteria/result
into AppState. Errors from the data contract are shown verbatim.

Interface:
    view(id) -> Tag
    server(id, state: AppState) -> None   (writes state.dataset, state.criteria, state.result;
                                           resets state.breadcrumb and state.selected_ids on success)
"""

from __future__ import annotations

from shiny import module, reactive, render, ui

from progeny_selector.core.navigation import build_tree
from progeny_selector.core.pipeline import run_analysis
from progeny_selector.io import load_dataset, read_criteria
from progeny_selector.model.criteria import CriteriaError
from progeny_selector.model.dataset import DataContractError


@module.ui
def ui_(id: str = "") -> ui.Tag:
    return ui.layout_columns(
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
    )


def view(id: str) -> ui.Tag:
    return ui_(id)


@module.server
def server_(input, output, session, state) -> None:
    message = reactive.Value("Choose files and press Load.")

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
            n_pass = sum(1 for r in result.rows if r["passes_filters"])
            lines = [f"{dataset.genotypes.n_markers} markers, {len(dataset.progeny)} progeny; {n_pass} pass hard filters"]
            lines += [f"warning: {w}" for w in result.warnings]
            message.set("\n".join(lines))
        except (DataContractError, CriteriaError) as exc:
            message.set(f"error: {exc}")

    @render.text
    def status() -> str:
        return message()


def server(id: str, state) -> None:
    server_(id, state)
