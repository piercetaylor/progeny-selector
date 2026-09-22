"""Criteria editor on the Load screen, criteria.yaml download, and the Export screen downloads."""

from __future__ import annotations

import csv
import dataclasses
import re
from pathlib import Path

from playwright.sync_api import Page, expect
from shiny.playwright import controller
from shiny.pytest import create_app_fixture
from shiny.run import ShinyAppProc

from progeny_selector.core.pipeline import run_analysis
from progeny_selector.io import load_dataset, read_criteria
from progeny_selector.io.criteria import read_criteria_text
from progeny_selector.io.export import FIXED_COLUMNS, results_csv_text
from tests.e2e.helpers import FIXTURE, load_fixture

app = create_app_fixture(["../../src/progeny_selector/app/app.py"])

# With filters.max_missing_rate raised from 0.2 to 0.5, BC2F1-F1-003 (the planted high-missing
# individual) also passes: 12, computed with run_analysis on the fixture before hard-coding.
N_PASS_RELAXED = 12
SAMPLES_HEADER = "sample_id,line_name,role,generation,family_id,notes"


def _open(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)


def _apply(page: Page, text: str) -> None:
    controller.InputTextArea(page, "load-criteria_text").set(text)
    controller.InputActionButton(page, "load-apply").click()


def _download(page: Page, output_id: str) -> tuple[str, bytes]:
    with page.expect_download() as dl:
        controller.DownloadButton(page, output_id).click()
    return dl.value.suggested_filename, Path(dl.value.path()).read_bytes()


def test_apply_edit_failed_apply_and_download(page: Page, app: ShinyAppProc) -> None:
    _open(page, app)
    editor = controller.InputTextArea(page, "load-criteria_text")
    editor.expect_value(re.compile(r"^name: synthetic BC2F1 fixture"))
    status = page.locator("#load-criteria_status")
    expect(status).to_contain_text("criteria loaded; edit and Apply to re-analyse (YAML comments are not kept)")
    text = editor.loc.input_value()
    assert "max_missing_rate: 0.2\n" in text
    _apply(page, text.replace("max_missing_rate: 0.2\n", "max_missing_rate: 0.5\n"))
    expect(status).to_contain_text(f"re-analysed: {N_PASS_RELAXED} pass hard filters (YAML comments are not kept)", timeout=60_000)
    editor.expect_value(re.compile(r"max_missing_rate: 0\.5\n"))

    navbar = controller.PageNavbar(page, "screen")
    navbar.set("rank")
    table = controller.OutputDataFrame(page, "rank-table")
    table.expect_nrow(N_PASS_RELAXED)

    navbar.set("load")
    _apply(page, "targets: []")
    expect(status).to_contain_text("criteria must define at least one target locus", timeout=60_000)
    navbar.set("rank")
    table.expect_nrow(N_PASS_RELAXED)

    # The download serialises the last applied criteria, not the rejected editor text.
    navbar.set("load")
    name, body = _download(page, "load-download_criteria")
    assert name == "criteria.yaml"
    expected = read_criteria(FIXTURE / "criteria.yaml")
    expected = dataclasses.replace(expected, filters=dataclasses.replace(expected.filters, max_missing_rate=0.5))
    assert read_criteria_text(body.decode("utf-8")) == expected


def test_region_shorthand_is_normalised(page: Page, app: ShinyAppProc) -> None:
    _open(page, app)
    text = (FIXTURE / "criteria.yaml").read_text(encoding="utf-8")
    assert "    marker_id: syn_Gm06_13\n" in text
    _apply(page, text.replace("    marker_id: syn_Gm06_13\n", '    region: "Gm06:25,000,000-26,000,000"\n'))
    expect(page.locator("#load-criteria_status")).to_contain_text("re-analysed:", timeout=60_000)
    editor = controller.InputTextArea(page, "load-criteria_text")
    editor.expect_value(re.compile(r"  chrom: Gm06\n  start_bp: 25000000\n  end_bp: 26000000\n"))
    assert "region:" not in editor.loc.input_value()


def test_apply_without_data(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    _apply(page, (FIXTURE / "criteria.yaml").read_text(encoding="utf-8"))
    expect(page.locator("#load-criteria_status")).to_have_text("Load data first.")
    name, body = _download(page, "load-download_criteria")
    assert name == "criteria.yaml"
    assert body == b"# No criteria applied yet. Load or apply criteria first.\n"


def test_manifest_before_load_is_header_only(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    controller.PageNavbar(page, "screen").set("export")
    name, body = _download(page, "export-manifest")
    assert name == "next_samples.csv"
    assert body == (SAMPLES_HEADER + "\r\n").encode("utf-8")


def test_results_before_load_is_header_only(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    controller.PageNavbar(page, "screen").set("export")
    name, body = _download(page, "export-results")
    assert name == "results.csv"
    # results.csv schema 1.1.0: with no results the header is exactly the fixed columns (docs/adr/0016).
    assert body == (",".join(FIXED_COLUMNS) + "\r\n").encode("utf-8")
    assert len(FIXED_COLUMNS) == 36


def test_export_downloads(page: Page, app: ShinyAppProc) -> None:
    _open(page, app)
    navbar = controller.PageNavbar(page, "screen")
    navbar.set("rank")
    table = controller.OutputDataFrame(page, "rank-table")
    table.expect_nrow(11)
    # Non-contiguous rows: the controller's shift-click path for contiguous ranges selected nothing (phase 3).
    table.select_rows([0, 2])
    table.expect_selected_num_rows(2)

    navbar.set("export")
    name, body = _download(page, "export-results")
    assert name == "results.csv"
    assert body.decode("utf-8").startswith("rank_overall,rank_in_family,sample_id")
    dataset = load_dataset(FIXTURE / "genotypes.vcf", FIXTURE / "samples.csv", FIXTURE / "markers.csv")
    result = run_analysis(dataset, read_criteria(FIXTURE / "criteria.yaml"))
    assert body == results_csv_text(result.rows).encode("utf-8")

    name, body = _download(page, "export-selected")
    assert name == "selected.csv"
    rows = list(csv.reader(body.decode("utf-8").splitlines()))
    assert len(rows) == 3

    name, body = _download(page, "export-manifest")
    assert name == "next_samples.csv"
    assert body.decode("utf-8").split("\r\n")[0] == SAMPLES_HEADER
