"""A second Load resets navigation and selection: no stale family filter, no stale selected ids."""

from __future__ import annotations

from playwright.sync_api import Page, expect
from shiny.playwright import controller
from shiny.pytest import create_app_fixture
from shiny.run import ShinyAppProc

from tests.e2e.helpers import load_fixture

app = create_app_fixture(["../../src/progeny_selector/app/app.py"])


def test_second_load_resets_state(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)
    navbar = controller.PageNavbar(page, "screen")
    navbar.set("navigate")
    families = controller.Accordion(page, "navigate-families")
    families.expect_panels(["fam_0", "fam_1"])
    families.set("fam_1")
    expect(page.locator("nav[aria-label=breadcrumb] li[aria-current=page]")).to_have_text("F2")

    navbar.set("rank")
    table = controller.OutputDataFrame(page, "rank-table")
    table.expect_nrow(7)
    table.select_rows([0, 2])
    table.expect_selected_num_rows(2)
    navbar.set("compare")
    expect(page.locator("#compare-panels .card")).to_have_count(2)
    expect(page.locator("#compare-panels")).to_contain_text("BC2F1-F2-")

    navbar.set("load")
    controller.InputActionButton(page, "load-run").click()
    expect(page.locator("#load-status")).to_contain_text("500 markers, 40 progeny; 11 pass hard filters", timeout=60_000)

    # Compare and Select first: the Rank grid is hidden and not re-rendered yet, so a stale
    # selection could only be cleared by the Load screen's reset.
    navbar.set("compare")
    panel = page.locator("#compare-panels")
    expect(panel).to_have_text("Select rows on the Rank screen.")
    expect(panel).not_to_contain_text("BC2F1-")
    expect(page.locator("#compare-panels.shiny-output-error")).to_have_count(0)

    navbar.set("select")
    summary = page.locator("#select-summary")
    expect(summary).to_have_text("No individuals selected.")
    expect(page.locator("#select-summary.shiny-output-error")).to_have_count(0)

    navbar.set("rank")
    controller.InputSwitch(page, "rank-only_pass").set(False)
    table.expect_nrow(40)

    navbar.set("navigate")
    families.expect_open([])
    expect(page.locator("nav[aria-label=breadcrumb] li")).to_have_count(1)
