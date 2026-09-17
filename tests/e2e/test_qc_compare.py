"""Validate/QC table and Compare screen end to end: highlighted QC rows, marker summary, side-by-side cards and strips."""

from __future__ import annotations

import csv

from playwright.sync_api import Page, expect
from shiny.playwright import controller
from shiny.pytest import create_app_fixture
from shiny.run import ShinyAppProc

from progeny_selector.app.screens.rank import DISPLAY_COLUMNS
from tests.e2e.helpers import FIXTURE, load_fixture

app = create_app_fixture(["../../src/progeny_selector/app/app.py"])

RPP_TOTAL_COL = DISPLAY_COLUMNS.index("rpp_total")
PASS_RGB = "rgb(0, 158, 115)"  # STATUS_COLORS["pass"], #009E73
EXCLUDED_RGBA = "rgba(213, 94, 0, 0.25)"  # STATUS_COLORS["fail"] #D55E00 at 25 % (present.qc_row_styles)
QC_COLUMNS = [
    "sample_id",
    "line_name",
    "family_id",
    "generation",
    "missing_rate",
    "het_rate",
    "expected_het",
    "hom_donor_rate",
    "nonparental_rate",
    "expected_rpp",
    "ibs_rp",
    "ibs_donor",
    "flags",
    "qc_excluded",
]


def _open(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)


def _progeny_index(sample_id: str) -> int:
    """QC rows follow the manifest's progeny order (core.qc.sample_qc)."""
    with open(FIXTURE / "samples.csv", newline="") as fh:
        progeny = [r["sample_id"] for r in csv.DictReader(fh) if r["role"] in ("progeny", "candidate")]
    return progeny.index(sample_id)


def test_qc_table_and_summary(page: Page, app: ShinyAppProc) -> None:
    _open(page, app)
    controller.PageNavbar(page, "screen").set("qc")
    table = controller.OutputDataFrame(page, "qc-table")
    table.expect_nrow(40)
    labels = [label.strip() for label in table.loc_column_label.all_inner_texts()]
    assert "flags" in labels and "qc_excluded" in labels
    assert labels == QC_COLUMNS
    row = _progeny_index("BC2F1-F2-002")
    excluded_col = QC_COLUMNS.index("qc_excluded")
    table.expect_cell("BC2F1-F2-002", row=row, col=0)
    # The DataGrid renders a pandas bool as lowercase text.
    table.expect_cell("true", row=row, col=excluded_col)
    expect(table.cell_locator(row, 0)).to_have_css("background-color", EXCLUDED_RGBA)
    summary = page.locator("#qc-summary dl")
    # expected: 500 markers, 475 informative; uninformative: 20 monomorphic, 3 donor missing, 2 RP heterozygous.
    expect(summary).to_contain_text("475")
    expect(summary).to_contain_text("parents identical (monomorphic): 20")


def _select_rank_rows(page: Page, rows: list[int]) -> None:
    table = controller.OutputDataFrame(page, "rank-table")
    table.expect_nrow(11)
    table.cell_locator(rows[0], 0).click()
    for r in rows[1:]:
        table.cell_locator(r, 0).click(modifiers=["ControlOrMeta"])
    table.expect_selected_num_rows(len(rows))


def test_compare_two_cards(page: Page, app: ShinyAppProc) -> None:
    _open(page, app)
    navbar = controller.PageNavbar(page, "screen")
    navbar.set("rank")
    table = controller.OutputDataFrame(page, "rank-table")
    table.expect_nrow(11)
    # Contiguous indices take the controller's shift-click path, which selected nothing in phase 3; [0, 2] ctrl-clicks.
    table.select_rows([0, 2])
    table.expect_selected_num_rows(2)
    navbar.set("compare")
    cards = page.locator("#compare-panels .card")
    expect(cards).to_have_count(2)
    first = cards.nth(0)
    expect(first).to_contain_text("BC2F1-F1-001")
    chip = first.locator(".status-chip").first
    expect(chip).to_have_text("pass")
    expect(chip).to_have_css("background-color", PASS_RGB)
    expect(first.locator("td", has_text="Gm06 (carrier)")).to_have_count(1)
    drag = first.locator("table[data-locus=T1]")
    expect(drag).to_contain_text("drag at T1")
    for label in ("left min", "left max", "right min", "right max", "total estimate", "total max", "recombinant left"):
        expect(drag).to_contain_text(label)
    # expected_results.csv: BC2F1-F1-001 recomb_left True, recomb_right False, drag_total_max_cm 15.284.
    expect(drag.locator("tr", has_text="recombinant left")).to_contain_text("yes")
    expect(drag.locator("tr", has_text="recombinant right")).to_contain_text("no")
    expect(drag.locator("tr", has_text="total max")).to_contain_text("15.28 cm")
    svg = first.locator("svg[role=img]")
    expect(svg).to_have_count(1)
    expect(svg.locator("text")).to_have_count(20)
    assert svg.locator("rect[fill='#0072B2']").count() >= 1
    # Gm06 for BC2F1-F1-001: A, U, A, H (donor segment), A.
    assert svg.locator("g[data-chrom=Gm06] rect").count() >= 3
    expect(first.locator("svg[role=img] rect title").first).not_to_be_empty()

    navbar.set("rank")
    table.set_sort({"col": RPP_TOTAL_COL, "desc": True})
    # Sorted by rpp_total descending, view row 2 is BC2F1-F2-015 (test_selection_follows_sorted_view).
    table.select_rows([0, 2])
    table.expect_selected_num_rows(2)
    navbar.set("compare")
    cards = page.locator("#compare-panels .card")
    expect(cards).to_have_count(2)
    second = cards.nth(1)
    expect(second).to_contain_text("BC2F1-F2-015")
    drag2 = second.locator("table[data-locus=T1]")
    # expected_results.csv: BC2F1-F2-015 drag_total_max_cm 109.534.
    expect(drag2.locator("tr", has_text="total max")).to_contain_text("109.53 cm")


def test_compare_caps_at_six(page: Page, app: ShinyAppProc) -> None:
    _open(page, app)
    navbar = controller.PageNavbar(page, "screen")
    navbar.set("rank")
    _select_rank_rows(page, list(range(7)))
    navbar.set("compare")
    expect(page.locator("#compare-panels .card")).to_have_count(6)
    expect(page.locator("#compare-panels")).to_contain_text("showing 6 of 7")


def test_six_cards_fit_1280(page: Page, app: ShinyAppProc) -> None:
    page.set_viewport_size({"width": 1280, "height": 800})
    _open(page, app)
    navbar = controller.PageNavbar(page, "screen")
    navbar.set("rank")
    _select_rank_rows(page, list(range(6)))
    navbar.set("compare")
    cards = page.locator("#compare-panels .card")
    expect(cards).to_have_count(6)
    boxes = [cards.nth(k).bounding_box() for k in range(6)]
    assert all(b is not None for b in boxes)
    tops = {round(b["y"]) for b in boxes if b is not None}
    print("card boxes at 1280 px:", [(round(b["x"]), round(b["width"])) for b in boxes if b is not None])
    assert len(tops) == 1, f"cards wrapped onto {len(tops)} rows"
    assert max(b["x"] + b["width"] for b in boxes if b is not None) <= 1280
    assert page.evaluate("document.documentElement.scrollWidth") <= 1280
