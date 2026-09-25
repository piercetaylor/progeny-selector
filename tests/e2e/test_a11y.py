"""Accessibility review: axe-core against every screen, plus the non-colour status check.

docs/accessibility.md records the method, tooling versions and known exceptions; this module is
the browser side of that review. ``KNOWN_A11Y_EXCEPTIONS`` holds rule ids for violations inside
markup this app does not author (the DataGrid's own internals, Bootstrap's navbar); every entry
here must also appear in docs/accessibility.md and in a dated line of docs/adr/0023.
"""

from __future__ import annotations

import pytest

pytest.importorskip("playwright")
pytest.importorskip("axe_playwright_python")

from axe_playwright_python.sync_playwright import Axe
from playwright.sync_api import Page, expect
from shiny.playwright import controller
from shiny.pytest import create_app_fixture
from shiny.run import ShinyAppProc

from tests.e2e.helpers import load_fixture

app = create_app_fixture(["../../src/progeny_selector/app/app.py"])

AXE_TAGS = ["wcag2a", "wcag2aa", "wcag21a", "wcag21aa", "wcag22aa"]

# axe rule id -> one-line reason. Every entry must also appear in docs/accessibility.md and in a
# dated line of docs/adr/0023-accessibility-baseline.md.
KNOWN_A11Y_EXCEPTIONS: dict[str, str] = {
    "aria-allowed-attr": (
        "shiny.render.DataGrid's own <table> markup (Validate, Rank and Selection-list grids); this app does not author it."
    ),
    "label": "shiny.render.DataGrid's own per-column filter <input>s (Validate and Rank grids); this app does not author them.",
    "scrollable-region-focusable": (
        "bslib's own fillable card-body container (ui.card under page_navbar(fillable=True)); "
        "this app does not author its scroll/tabindex behaviour."
    ),
}

SCREENS = ("load", "qc", "navigate", "rank", "compare", "select", "export")


def _violations(page: Page) -> list[dict]:
    results = Axe().run(page, options={"runOnly": {"type": "tag", "values": AXE_TAGS}})
    violations = results.response["violations"]
    return [v for v in violations if v["id"] not in KNOWN_A11Y_EXCEPTIONS]


def _report(violations: list[dict]) -> str:
    lines = []
    for v in violations:
        target = v["nodes"][0]["target"][0] if v["nodes"] else "?"
        lines.append(f"{v['id']} ({v['impact']}): {v['help']} [{target}]")
    return "\n".join(lines)


def test_load_screen_before_data(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    violations = _violations(page)
    assert not violations, _report(violations)


@pytest.mark.parametrize("screen", SCREENS)
def test_every_screen_after_load(page: Page, app: ShinyAppProc, screen: str) -> None:
    page.goto(app.url)
    load_fixture(page)
    navbar = controller.PageNavbar(page, "screen")
    navbar.set("rank")
    table = controller.OutputDataFrame(page, "rank-table")
    table.expect_nrow(11)
    table.select_rows([0, 2])
    table.expect_selected_num_rows(2)
    navbar.set(screen)
    violations = _violations(page)
    assert not violations, _report(violations)


def test_compare_empty_before_selection(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)
    controller.PageNavbar(page, "screen").set("compare")
    violations = _violations(page)
    assert not violations, _report(violations)


def test_exceptions_are_documented() -> None:
    from pathlib import Path

    doc = Path(__file__).parents[2] / "docs" / "accessibility.md"
    text = doc.read_text(encoding="utf-8")
    for rule_id in KNOWN_A11Y_EXCEPTIONS:
        assert rule_id in text, f"{rule_id} is excepted in test_a11y.py but missing from docs/accessibility.md"


def test_status_is_never_colour_alone(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)
    navbar = controller.PageNavbar(page, "screen")
    navbar.set("rank")
    table = controller.OutputDataFrame(page, "rank-table")
    table.expect_nrow(11)
    status_cols = [i for i, label in enumerate(table.loc_column_label.all_inner_texts()) if _is_status_column(label.strip())]
    assert status_cols
    for col in status_cols:
        for row in range(11):
            text = table.cell_locator(row, col).inner_text().strip()
            assert text, f"status cell (row {row}, col {col}) has no text, only colour"
    table.select_rows([0, 2])
    table.expect_selected_num_rows(2)
    navbar.set("compare")
    expect(page.locator("#compare-panels .card")).to_have_count(2)
    badges = page.locator("#compare-panels .status-chip")
    count = badges.count()
    assert count > 0
    for i in range(count):
        expect(badges.nth(i)).not_to_be_empty()
    svgs = page.locator("#compare-panels svg[role=img]")
    svg_count = svgs.count()
    assert svg_count > 0
    for i in range(svg_count):
        label = svgs.nth(i).get_attribute("aria-label")
        assert label


def _is_status_column(label: str) -> bool:
    return (label.startswith("target_") or label.startswith("avoid_")) and label.endswith("_status")
