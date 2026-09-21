"""Executes docs/keyboard-walkthrough.md and asserts what it states."""

from __future__ import annotations

import csv

from playwright.sync_api import Page, expect
from shiny.playwright import controller
from shiny.pytest import create_app_fixture
from shiny.run import ShinyAppProc

from progeny_selector.app.screens.select import NOTES_COLUMNS
from tests.e2e.helpers import load_fixture

app = create_app_fixture(["../../src/progeny_selector/app/app.py"])


def test_navbar_arrow_keys_move_focus_enter_activates(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)
    controller.PageNavbar(page, "screen").set("rank")
    page.locator("a.nav-link", has_text="4 Rank").focus()
    page.keyboard.press("ArrowRight")
    page.keyboard.press("Enter")
    expect(page.locator(".nav-link.active")).to_have_text("5 Compare")


def test_rank_grid_needs_a_click_before_the_keyboard_works(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)
    controller.PageNavbar(page, "screen").set("rank")
    page.locator("#rank-only_pass").focus()
    for _ in range(6):
        page.keyboard.press("Tab")
    page.keyboard.press("ArrowDown")
    page.keyboard.press(" ")
    assert page.locator("#rank-table tr[aria-selected=true]").count() == 0

    table = controller.OutputDataFrame(page, "rank-table")
    table.cell_locator(0, 0).click()
    assert page.locator("#rank-table tr[aria-selected=true]").count() == 1
    page.keyboard.press("ArrowDown")
    page.keyboard.press("ArrowDown")
    assert page.locator("#rank-table tr[aria-selected=true]").count() == 1
    page.keyboard.press(" ")
    assert page.locator("#rank-table tr[aria-selected=true]").count() == 2


def test_selection_list_enter_commits(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)
    controller.PageNavbar(page, "screen").set("rank")
    controller.OutputDataFrame(page, "rank-table").select_rows([0, 2])
    controller.PageNavbar(page, "screen").set("select")
    select_table = controller.OutputDataFrame(page, "select-table")
    select_table.expect_nrow(2)
    notes_col = NOTES_COLUMNS.index("notes")

    select_table.set_cell("committed", row=1, col=notes_col, finish_key="Enter")
    with page.expect_download() as dl:
        controller.PageNavbar(page, "screen").set("export")
        controller.DownloadButton(page, "export-selected").click()
    rows = list(csv.DictReader(dl.value.path().read_text(encoding="utf-8").splitlines()))
    notes = {r["sample_id"]: r["notes"] for r in rows}
    assert "committed" in notes.values()


def test_selection_list_escape_cancels(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)
    controller.PageNavbar(page, "screen").set("rank")
    controller.OutputDataFrame(page, "rank-table").select_rows([0, 2])
    controller.PageNavbar(page, "screen").set("select")
    select_table = controller.OutputDataFrame(page, "select-table")
    select_table.expect_nrow(2)
    notes_col = NOTES_COLUMNS.index("notes")

    select_table.set_cell("should not stick", row=0, col=notes_col, finish_key="Escape")
    expect(select_table.cell_locator(0, notes_col)).to_have_text("")
