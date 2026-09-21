"""Selection list: a cleared 'Top N per family' field adds nothing and leaves the screen working.

The widget's min is client-side only, so the field can be empty when the button is pressed and
``input.top_n()`` is then None. Only a browser reaches this path: the count is read inside a Shiny
reactive effect, where ``int(None)`` raised a TypeError that pytest never sees.
"""

from __future__ import annotations

from playwright.sync_api import Page
from shiny.playwright import controller
from shiny.pytest import create_app_fixture
from shiny.run import ShinyAppProc

from progeny_selector.app.screens.select import NOTES_COLUMNS
from tests.e2e.helpers import load_fixture

app = create_app_fixture(["../../src/progeny_selector/app/app.py"])

SAMPLE_ID_COL = NOTES_COLUMNS.index("sample_id")


def test_cleared_top_n_adds_nothing_and_the_screen_still_works(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)
    navbar = controller.PageNavbar(page, "screen")
    navbar.set("rank")
    rank_table = controller.OutputDataFrame(page, "rank-table")
    rank_table.expect_nrow(11)
    rank_table.select_rows([0])

    navbar.set("select")
    select_table = controller.OutputDataFrame(page, "select-table")
    select_table.expect_nrow(1)
    controller.InputNumeric(page, "select-top_n").set("")
    controller.InputActionButton(page, "select-apply_top_n").click()
    select_table.expect_nrow(1)

    # The effect is still alive after the cleared press: a real N still adds the rank-1 individual of
    # the other family. With int(None) raising, the session is gone and this second press does nothing.
    controller.InputNumeric(page, "select-top_n").set("1")
    controller.InputActionButton(page, "select-apply_top_n").click()
    select_table.expect_nrow(2)
