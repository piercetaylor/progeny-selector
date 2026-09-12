"""Load screen end to end: upload the synthetic fixture, run the analysis, read the status summary."""

from __future__ import annotations

from playwright.sync_api import Page
from shiny.pytest import create_app_fixture
from shiny.run import ShinyAppProc

from tests.e2e.helpers import load_fixture

app = create_app_fixture(["../../src/progeny_selector/app/app.py"])


def test_load_fixture(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)
