"""Screen 1, Load: genotype file, samples.csv, optional markers.csv, criteria.yaml; criteria editor.

Responsibility: accept the four input files, parse them at the boundary via
``progeny_selector.io``, run the analysis, and publish dataset/criteria/result
into AppState. The criteria card is a YAML view over the applied ``Criteria``
in AppState: it is rewritten in canonical form after every successful Load or
Apply, Apply parses the text with the strict reader and re-analyses the loaded
dataset, and the download always serialises the applied criteria, never the
raw editor text. Errors from the data contract and the criteria reader are
shown verbatim; a failed Apply leaves AppState untouched. A "Token profile" select
(contract 1.4.0: the contract default or a built-in profile) and a "Custom token profile
(JSON, optional)" file input choose how HapMap and wide-CSV cells are read; a custom file,
when given, is used instead of the select, which is disabled while it is set (as backcross's
Upload screen); "Clear custom token profile" forgets it and resets the file input; an
unreadable or invalid file is a load error. A "Crop" select (contract 1.5.0: the nine built-in
chromosome schemes, soybean the default) chooses the scheme chromosome names are normalised and
ordered under; positions are not converted between assemblies.

Interface:
    read_profile_json(path) -> dict   (a custom token profile file; DataContractError when unreadable
                                       or not a JSON object)
    view(id) -> Tag
    server(id, state: AppState) -> None   (writes state.dataset, state.criteria, state.result;
                                           Load resets state.breadcrumb, state.selected_ids and state.notes
                                           on success; Apply clears state.selected_ids and keeps
                                           state.breadcrumb and state.notes)
"""

from __future__ import annotations

import json
import logging

from htmltools import Tag, TagChild
from shiny import module, reactive, render, ui

from progeny_selector.core.navigation import build_tree
from progeny_selector.core.pipeline import run_analysis
from progeny_selector.io import load_dataset, read_criteria
from progeny_selector.io.criteria import dump_criteria_yaml, read_criteria_text
from progeny_selector.io.crops import BUILTIN_CROPS, DEFAULT_CROP_ID
from progeny_selector.io.profiles import BUILTIN_PROFILES, DEFAULT_PROFILE_ID
from progeny_selector.model.criteria import CriteriaError
from progeny_selector.model.dataset import DataContractError


@module.ui
def ui_(id: str = "") -> ui.Tag:
    return ui.div(
        ui.layout_columns(
            ui.card(
                ui.card_header("Input files"),
                ui.input_file(
                    "genotypes",
                    "Genotypes (VCF, VCF.gz, HapMap, wide CSV)",
                    accept=[".vcf", ".gz", ".bgz", ".txt", ".csv", ".tsv", ".hmp", ".hapmap"],
                ),
                ui.input_file("samples", "samples.csv", accept=[".csv", ".tsv", ".txt"]),
                ui.input_file("markers", "markers.csv (optional)", accept=[".csv", ".tsv", ".txt"]),
                ui.output_ui("profile_select"),
                ui.output_ui("profile_file_input"),
                ui.input_action_button("profile_clear", "Clear custom token profile"),
                ui.input_select("crop", "Crop", CROP_CHOICES, selected=DEFAULT_CROP_ID),
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


PROFILE_CHOICES: dict[str, str] = {
    DEFAULT_PROFILE_ID: "Contract default (by format)",
    **{pid: p.name for pid, p in BUILTIN_PROFILES.items()},
}


# The crop chromosome schemes (contract 1.5.0), in file order with soybean first and selected.
CROP_CHOICES: dict[str, str] = {cid: scheme.name for cid, scheme in BUILTIN_CROPS.items()}


def _disable_select(node: TagChild) -> None:
    """Marks every <select> under ``node`` disabled."""
    if isinstance(node, Tag):
        if node.name == "select":
            node.attrs["disabled"] = "disabled"
        for child in node.children:
            _disable_select(child)


def profile_select_tag(selected: str, disabled: bool) -> Tag:
    """The Token profile select; disabled while a custom profile file is set."""
    tag = ui.input_select("profile", "Token profile", PROFILE_CHOICES, selected=selected)
    if disabled:
        _disable_select(tag)
    return tag


def read_profile_json(path: str) -> dict:
    """The custom token profile file as a dict; DataContractError when unreadable or not a JSON object."""
    try:
        with open(path, encoding="utf-8") as fh:
            loaded = json.load(fh)
    except (OSError, ValueError) as exc:
        raise DataContractError(f"token profile file: {exc}") from exc
    if not isinstance(loaded, dict):
        raise DataContractError("token profile: must be a JSON object")
    return loaded


def view(id: str) -> ui.Tag:
    return ui_(id)


@module.server
def server_(input, output, session, state) -> None:
    message = reactive.Value("Choose files and press Load.")
    criteria_message = reactive.Value("")
    custom_profile_path: reactive.Value[str | None] = reactive.Value(None)

    @reactive.effect
    @reactive.event(input.profile_file)
    def _profile_file() -> None:
        f = input.profile_file()
        custom_profile_path.set(f[0]["datapath"] if f else None)

    @reactive.effect
    @reactive.event(input.profile_clear)
    def _profile_clear() -> None:
        custom_profile_path.set(None)

    @render.ui
    def profile_select() -> Tag:
        with reactive.isolate():
            try:
                selected = input.profile() or DEFAULT_PROFILE_ID
            except Exception:
                selected = DEFAULT_PROFILE_ID
        return profile_select_tag(selected if selected in PROFILE_CHOICES else DEFAULT_PROFILE_ID, custom_profile_path() is not None)

    @render.ui
    def profile_file_input() -> Tag:
        input.profile_clear()  # a Clear re-renders the input, so the chosen file name disappears too
        return ui.input_file("profile_file", "Custom token profile (JSON, optional)", accept=[".json"])

    @reactive.effect
    @reactive.event(input.run)
    def _run() -> None:
        try:
            g, s, c = input.genotypes(), input.samples(), input.criteria()
            if not (g and s and c):
                message.set("Genotypes, samples.csv and criteria.yaml are required.")
                return
            m = input.markers()
            pf = custom_profile_path()
            profile: str | dict = input.profile() or DEFAULT_PROFILE_ID
            if pf is not None:
                profile = read_profile_json(pf)
            crop = input.crop() or DEFAULT_CROP_ID
            dataset = load_dataset(g[0]["datapath"], s[0]["datapath"], m[0]["datapath"] if m else None, profile=profile, crop=crop)
            criteria = read_criteria(c[0]["datapath"])
            result = run_analysis(dataset, criteria)
            state.dataset.set(dataset)
            state.criteria.set(criteria)
            state.result.set(result)
            # A new dataset starts at the cross node with nothing selected.
            state.breadcrumb.set({"cross": build_tree(dataset).cross, "family": None, "generation": None})
            state.selected_ids.set([])
            state.notes.set({})
            ui.update_text_area("criteria_text", value=dump_criteria_yaml(criteria))
            criteria_message.set("criteria loaded; edit and Apply to re-analyse (YAML comments are not kept)")
            n_pass = sum(1 for r in result.rows if r["passes_filters"])
            lines = [f"{dataset.genotypes.n_markers} markers, {len(dataset.progeny)} progeny; {n_pass} pass hard filters"]
            lines += [f"warning: {w}" for w in result.warnings]
            message.set("\n".join(lines))
        except (DataContractError, CriteriaError) as exc:
            message.set(f"error: {exc}")
        except Exception as exc:
            message.set(f"unexpected error: {type(exc).__name__}: {exc}")
            logging.getLogger(__name__).exception("load failed")

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
        except Exception as exc:
            # Nothing in AppState has been written yet, so a failed Apply leaves the last good run in place.
            criteria_message.set(f"unexpected error: {type(exc).__name__}: {exc}")
            logging.getLogger(__name__).exception("load failed")
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
