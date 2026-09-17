"""Navigate tree and Rank grid end to end: breadcrumb filtering, status chips, keyboard use, view-aware selection."""

from __future__ import annotations

import csv
import shutil

from playwright.sync_api import Page, expect
from shiny.playwright import controller
from shiny.pytest import create_app_fixture
from shiny.run import ShinyAppProc

from progeny_selector.app.screens.rank import DISPLAY_COLUMNS
from tests.e2e.helpers import FIXTURE, LOAD_STATUS, load_fixture

app = create_app_fixture(["../../src/progeny_selector/app/app.py"])

# Frame columns are DISPLAY_COLUMNS, then status columns, then recombinant flags (rank.py).
STATUS_T1_COL = len(DISPLAY_COLUMNS)
RPP_TOTAL_COL = DISPLAY_COLUMNS.index("rpp_total")
SAMPLE_ID_COL = DISPLAY_COLUMNS.index("sample_id")
PASS_RGB = "rgb(0, 158, 115)"  # STATUS_COLORS["pass"], Okabe-Ito bluish green #009E73


def _open(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)


def test_family_filters_rank(page: Page, app: ShinyAppProc) -> None:
    _open(page, app)
    navbar = controller.PageNavbar(page, "screen")
    navbar.set("navigate")
    families = controller.Accordion(page, "navigate-families")
    families.expect_panels(["fam_0", "fam_1"])
    families.set("fam_1")
    expect(page.locator("nav[aria-label=breadcrumb]")).to_contain_text("F2")
    navbar.set("rank")
    controller.InputSwitch(page, "rank-only_pass").expect_checked(True)
    table = controller.OutputDataFrame(page, "rank-table")
    # expected_results.csv: 11 pass hard filters, 4 in F1 and 7 in F2.
    table.expect_nrow(7)
    controller.InputSwitch(page, "rank-only_pass").set(False)
    table.expect_nrow(20)


def test_status_columns_and_chip_colour(page: Page, app: ShinyAppProc) -> None:
    _open(page, app)
    controller.PageNavbar(page, "screen").set("rank")
    table = controller.OutputDataFrame(page, "rank-table")
    table.expect_nrow(11)
    expect(table.loc_column_label.nth(STATUS_T1_COL)).to_have_text("target_T1_status")
    labels = table.loc_column_label.all_inner_texts()
    for name in ("target_T1_status", "avoid_AV1_status", "recomb_T1_left", "recomb_T1_right"):
        assert name in [label.strip() for label in labels]
    table.expect_cell("pass", row=0, col=STATUS_T1_COL)
    expect(table.cell_locator(0, STATUS_T1_COL)).to_have_css("background-color", PASS_RGB)


def test_keyboard_navigation(page: Page, app: ShinyAppProc) -> None:
    _open(page, app)
    controller.PageNavbar(page, "screen").set("navigate")
    controller.Accordion(page, "navigate-families").expect_panels(["fam_0", "fam_1"])
    page.locator("#navigate-families .accordion-button").first.focus()
    page.keyboard.press("Enter")
    controller.Accordion(page, "navigate-families").expect_open(["fam_0"])
    page.keyboard.press("Tab")
    expect(page.locator("#navigate-gen_0 input[value='']")).to_be_focused()
    page.keyboard.press("ArrowDown")
    expect(page.locator("#navigate-gen_0 input[value='BC2F1']")).to_be_checked()
    crumb = page.locator("nav[aria-label=breadcrumb]")
    expect(crumb).to_contain_text("BC2F1")
    expect(crumb.locator("li[aria-current=page]")).to_have_text("BC2F1")


def test_selection_follows_sorted_view(page: Page, app: ShinyAppProc) -> None:
    _open(page, app)
    navbar = controller.PageNavbar(page, "screen")
    navbar.set("rank")
    table = controller.OutputDataFrame(page, "rank-table")
    table.expect_nrow(11)
    table.set_sort({"col": RPP_TOTAL_COL, "desc": True})
    # Highest rpp_total among the 11 passing individuals, per expected_results.csv.
    table.expect_cell("BC2F1-F1-001", row=0, col=SAMPLE_ID_COL)
    first = table.cell_locator(0, SAMPLE_ID_COL).inner_text().strip()
    third = table.cell_locator(2, SAMPLE_ID_COL).inner_text().strip()
    # Unsorted, view row 2 would be BC2F1-F1-010 (rank 3); sorted by rpp_total it is BC2F1-F2-015.
    assert third == "BC2F1-F2-015"
    table.select_rows([0, 2])
    table.expect_selected_num_rows(2)
    navbar.set("compare")
    cards = page.locator("#compare-panels .card")
    expect(cards).to_have_count(2)
    expect(cards.nth(0).locator(".card-header")).to_contain_text(first)
    expect(cards.nth(1).locator(".card-header")).to_contain_text(third)


def test_crumb_links_step_back(page: Page, app: ShinyAppProc) -> None:
    _open(page, app)
    navbar = controller.PageNavbar(page, "screen")
    navbar.set("navigate")
    families = controller.Accordion(page, "navigate-families")
    families.expect_panels(["fam_0", "fam_1"])
    families.set("fam_1")
    generations = controller.InputRadioButtons(page, "navigate-gen_1")
    generations.set("BC2F1")
    crumb = page.locator("nav[aria-label=breadcrumb]")
    current = crumb.locator("li[aria-current=page]")
    expect(current).to_have_text("BC2F1")

    controller.InputActionLink(page, "navigate-crumb_1").click()
    generations.expect_selected("")
    expect(current).to_have_text("F2")
    expect(crumb.locator("li")).to_have_count(2)

    controller.InputActionLink(page, "navigate-crumb_0").click()
    families.expect_open([])
    expect(crumb.locator("li")).to_have_count(1)
    expect(current).to_have_text("Williams 82 (synthetic) x PI synthetic donor")

    navbar.set("rank")
    controller.InputSwitch(page, "rank-only_pass").set(False)
    controller.OutputDataFrame(page, "rank-table").expect_nrow(40)


def test_no_family_node(page: Page, app: ShinyAppProc, tmp_path) -> None:
    fixture_copy = tmp_path / "synthetic_bc2f1"
    shutil.copytree(FIXTURE, fixture_copy)
    samples_path = fixture_copy / "samples.csv"
    with open(samples_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames
        rows = list(reader)
    for row in rows:
        if row["sample_id"] == "BC2F1-F1-004":
            row["family_id"] = ""
    with open(samples_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    page.goto(app.url)
    controller.InputFile(page, "load-genotypes").set(fixture_copy / "genotypes.vcf")
    controller.InputFile(page, "load-samples").set(samples_path)
    controller.InputFile(page, "load-markers").set(fixture_copy / "markers.csv")
    controller.InputFile(page, "load-criteria").set(fixture_copy / "criteria.yaml")
    controller.InputActionButton(page, "load-run").click()
    expect(page.locator("#load-status")).to_contain_text(LOAD_STATUS, timeout=60_000)

    navbar = controller.PageNavbar(page, "screen")
    navbar.set("navigate")
    families = controller.Accordion(page, "navigate-families")
    families.expect_panels(["fam_0", "fam_1", "fam_2"])
    families.set("fam_2")
    crumb = page.locator("nav[aria-label=breadcrumb]")
    expect(crumb.locator("li[aria-current=page]")).to_have_text("(no family)")

    navbar.set("rank")
    controller.InputSwitch(page, "rank-only_pass").set(False)
    table = controller.OutputDataFrame(page, "rank-table")
    table.expect_nrow(1)
    table.expect_cell("BC2F1-F1-004", row=0, col=DISPLAY_COLUMNS.index("sample_id"))
