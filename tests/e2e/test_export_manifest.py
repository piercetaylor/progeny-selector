"""Export screen: a cleared 'rows per selected individual' field still produces a manifest.

The widget's min/max are client-side only, so the field can be empty when the download handler runs
and ``input.per_selected()`` is then None. Only a browser reaches this path: the handler is Shiny
reactive code, and ``int(None)`` raised a TypeError inside the download rather than anywhere pytest
could see it.
"""

from __future__ import annotations

import csv
from pathlib import Path

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
