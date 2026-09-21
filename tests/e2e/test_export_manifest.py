"""Export screen: a cleared 'rows per selected individual' field still produces a manifest.

The widget's min/max are client-side only, so the field can be empty when the download handler runs
and ``input.per_selected()`` is then None. Only a browser reaches this path: the handler is Shiny
reactive code, and ``int(None)`` raised a TypeError inside the download rather than anywhere pytest
could see it.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest
from playwright.sync_api import Error as PlaywrightError
from playwright.sync_api import Page
from shiny.playwright import controller
from shiny.pytest import create_app_fixture
from shiny.run import ShinyAppProc

from tests.e2e.helpers import load_fixture

app = create_app_fixture(["../../src/progeny_selector/app/app.py"])


def test_cleared_per_selected_still_downloads_one_row_per_selected(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)
    navbar = controller.PageNavbar(page, "screen")
    navbar.set("rank")
    rank_table = controller.OutputDataFrame(page, "rank-table")
    rank_table.expect_nrow(11)
    rank_table.select_rows([0])
    rank_table.expect_selected_num_rows(1)

    navbar.set("export")
    controller.InputNumeric(page, "export-per_selected").set("")
    with page.expect_download() as dl:
        controller.DownloadButton(page, "export-manifest").click()
    assert dl.value.suggested_filename == "next_samples.csv"
    rows = list(csv.reader(Path(dl.value.path()).read_bytes().decode("utf-8").splitlines()))
    header, data = rows[0], rows[1:]
    assert header == ["sample_id", "line_name", "role", "generation", "family_id", "notes"]
    # 2 parents plus the manifest writer's own default of one placeholder row for the one selection.
    assert len(data) == 3
    assert [r[0] for r in data if r[2] == "progeny"] == ["BC2F1-F1-001-BC3F1-001"]


def test_zero_per_selected_is_refused_and_no_file_is_downloaded(page: Page, app: ShinyAppProc) -> None:
    """0 passes the client-side min, so io/export refuses it; the screen shows that and downloads nothing.

    The silent failure this closes is a next_samples.csv holding the two parents and no progeny at
    all, which looks like a written file. The CLI already refuses --per-selected 0.
    """
    page.goto(app.url)
    load_fixture(page)
    navbar = controller.PageNavbar(page, "screen")
    navbar.set("rank")
    rank_table = controller.OutputDataFrame(page, "rank-table")
    rank_table.expect_nrow(11)
    rank_table.select_rows([0])

    navbar.set("export")
    controller.InputNumeric(page, "export-per_selected").set("0")
    controller.OutputTextVerbatim(page, "export-status").expect_value(
        "error: placeholder rows per selected individual must be 1 or more, got 0"
    )
    # No file at all: the download is cancelled, not saved as a parents-only manifest.
    with page.expect_download() as dl:
        controller.DownloadButton(page, "export-manifest").click()
    assert dl.value.failure() == "canceled"
    with pytest.raises(PlaywrightError):
        dl.value.path()
