"""Load screen end to end: upload the synthetic fixture, run the analysis, read the status summary."""

from __future__ import annotations

from playwright.sync_api import Page, expect
from shiny.playwright import controller
from shiny.pytest import create_app_fixture
from shiny.run import ShinyAppProc

from tests.e2e.helpers import EXAMPLE_FIRST_LINE, LOAD_STATUS, load_fixture

app = create_app_fixture(["../../src/progeny_selector/app/app.py"])


def test_load_fixture(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)


def test_load_example_button_runs_the_shipped_example(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    controller.InputActionButton(page, "load-example").click()
    expect(page.locator("#load-status")).to_contain_text(EXAMPLE_FIRST_LINE)
    expect(page.locator("#load-status")).to_contain_text(LOAD_STATUS, timeout=60_000)
    page.locator("a.nav-link[data-value='rank']").click()
    expect(page.locator("#rank-table.html-fill-item div.shiny-data-grid table tbody tr")).to_have_count(11)
