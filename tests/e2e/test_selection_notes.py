"""Selection list: editable notes reach selected.csv, keyed by sample_id and not by row."""

from __future__ import annotations

import csv
from pathlib import Path

from playwright.sync_api import Page
from shiny.playwright import controller
from shiny.pytest import create_app_fixture
from shiny.run import ShinyAppProc

from progeny_selector.app.screens.select import NOTES_COLUMNS
from progeny_selector.io.export import SELECTION_COLUMNS
from tests.e2e.helpers import load_fixture

app = create_app_fixture(["../../src/progeny_selector/app/app.py"])

NOTES_COL = NOTES_COLUMNS.index("notes")
SAMPLE_ID_COL = NOTES_COLUMNS.index("sample_id")


def _open(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)


def _download(page: Page, output_id: str) -> tuple[str, bytes]:
    with page.expect_download() as dl:
        controller.DownloadButton(page, output_id).click()
    return dl.value.suggested_filename, Path(dl.value.path()).read_bytes()


def test_notes_survive_sort_and_reach_selected_csv(page: Page, app: ShinyAppProc) -> None:
    _open(page, app)
    navbar = controller.PageNavbar(page, "screen")
    navbar.set("rank")
    rank_table = controller.OutputDataFrame(page, "rank-table")
    rank_table.expect_nrow(11)
    # rows[0, 2] unsorted are rank_overall 1 and 3: BC2F1-F1-001 and BC2F1-F1-010 (test_selection_follows_sorted_view).
    rank_table.select_rows([0, 2])
    rank_table.expect_selected_num_rows(2)

    navbar.set("select")
    select_table = controller.OutputDataFrame(page, "select-table")
    select_table.expect_nrow(2)
    # Descending by sample_id: "BC2F1-F1-010" > "BC2F1-F1-001", so view row 0 is BC2F1-F1-010.
    select_table.set_sort({"col": SAMPLE_ID_COL, "desc": True})
    select_table.expect_cell("BC2F1-F1-010", row=0, col=SAMPLE_ID_COL)
    select_table.set_cell("keep; vigorous", row=0, col=NOTES_COL, finish_key="Enter")

    navbar.set("export")
    name, body = _download(page, "export-selected")
    assert name == "selected.csv"
    rows = list(csv.DictReader(body.decode("utf-8").splitlines()))
    assert list(rows[0].keys()) == list(SELECTION_COLUMNS)
    by_id = {r["sample_id"]: r for r in rows}
    assert by_id["BC2F1-F1-010"]["notes"] == "keep; vigorous"
    assert by_id["BC2F1-F1-001"]["notes"] == ""

    # "Add top N per family" (N=1) adds BC2F1-F2-005 (its rank_in_family is 1) but must not disturb
    # the note already keyed to BC2F1-F1-010's sample_id (core.selection.select_top_n, unaffected by
    # this phase; confirmed against the fixture before writing this assertion).
    navbar.set("select")
    controller.InputNumeric(page, "select-top_n").set("1")
    controller.InputActionButton(page, "select-apply_top_n").click()
    select_table.expect_nrow(3)

    navbar.set("export")
    name, body = _download(page, "export-selected")
    rows = list(csv.DictReader(body.decode("utf-8").splitlines()))
    by_id = {r["sample_id"]: r for r in rows}
    assert by_id["BC2F1-F1-010"]["notes"] == "keep; vigorous"

    controller.InputNumeric(page, "export-per_selected").set("3")
    name, body = _download(page, "export-manifest")
    assert name == "next_samples.csv"
    manifest_rows = list(csv.reader(body.decode("utf-8").splitlines()))
    header, data = manifest_rows[0], manifest_rows[1:]
    assert header == ["sample_id", "line_name", "role", "generation", "family_id", "notes"]
    # 2 parents + 3 selected individuals x 3 placeholder rows each.
    assert len(data) == 2 + 3 * 3
    suffixes = sorted({r[0].rsplit("-", 1)[-1] for r in data if r[2] == "progeny"})
    assert suffixes == ["001", "002", "003"]
    assert any(r[0] == "BC2F1-F1-010-BC3F1-001" for r in data)
