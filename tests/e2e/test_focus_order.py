"""Focus walk: focus never lands on a hidden element, and the tab list and Load's button are reachable."""

from __future__ import annotations

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page
from shiny.pytest import create_app_fixture
from shiny.run import ShinyAppProc

from tests.e2e.helpers import load_fixture

app = create_app_fixture(["../../src/progeny_selector/app/app.py"])

MAX_TABS = 60


def _focused_info(page: Page) -> dict:
    return page.evaluate(
        """() => {
            const el = document.activeElement;
            if (!el) return {tagName: null, id: null, hidden: true};
            return {tagName: el.tagName, id: el.id, hidden: el.offsetParent === null};
        }"""
    )


def test_focus_order_after_load(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)
    page.locator("body").click()
    page.keyboard.press("Tab")  # away from wherever the click landed focus

    seen: list[dict] = []
    tab_role_focused = False
    load_run_focused = False
    for _ in range(MAX_TABS):
        info = _focused_info(page)
        seen.append(info)
        if info["tagName"] is not None and info["tagName"] != "BODY":
            assert not info["hidden"], f"focus landed on a hidden element: {info}"
        role = page.evaluate("document.activeElement ? document.activeElement.getAttribute('role') : null")
        if role == "tab":
            tab_role_focused = True
        if info["id"] == "load-run":
            load_run_focused = True
        page.keyboard.press("Tab")

    assert tab_role_focused, "no element with role='tab' was ever focused: the navbar is not reachable"
    assert load_run_focused, "#load-run was never focused: the Load screen's button is not reachable"
