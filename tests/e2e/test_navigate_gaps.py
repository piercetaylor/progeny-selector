"""Browser tests for the three screen-test gaps M1 deferred: the "(no generation)"
choice, the Rank caption naming an unassigned family or generation, and the
Navigate reload guard."""

from __future__ import annotations

import csv
import shutil

from playwright.sync_api import Page, expect
from shiny.playwright import controller
from shiny.pytest import create_app_fixture
from shiny.run import ShinyAppProc

from progeny_selector.app.screens.rank import DISPLAY_COLUMNS
from tests.e2e.helpers import FIXTURE, LOAD_STATUS, load_fixture

SAMPLE_ID_COL = DISPLAY_COLUMNS.index("sample_id")

app = create_app_fixture(["../../src/progeny_selector/app/app.py"])


def _blank_generation(fixture_copy, sample_ids: set[str]) -> None:
    samples_path = fixture_copy / "samples.csv"
    with open(samples_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames
        rows = list(reader)
    for row in rows:
        if row["sample_id"] in sample_ids:
            row["generation"] = ""
    with open(samples_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_no_generation_choice(page: Page, app: ShinyAppProc, tmp_path) -> None:
    fixture_copy = tmp_path / "synthetic_bc2f1"
    shutil.copytree(FIXTURE, fixture_copy)
    _blank_generation(fixture_copy, {"BC2F1-F1-004", "BC2F1-F1-007"})

    page.goto(app.url)
    controller.InputFile(page, "load-genotypes").set(fixture_copy / "genotypes.vcf")
    controller.InputFile(page, "load-samples").set(fixture_copy / "samples.csv")
    controller.InputFile(page, "load-markers").set(fixture_copy / "markers.csv")
    controller.InputFile(page, "load-criteria").set(fixture_copy / "criteria.yaml")
    controller.InputActionButton(page, "load-run").click()
    expect(page.locator("#load-status")).to_contain_text(LOAD_STATUS, timeout=60_000)

    navbar = controller.PageNavbar(page, "screen")
    navbar.set("navigate")
    families = controller.Accordion(page, "navigate-families")
    families.set("fam_0")
    radio = controller.InputRadioButtons(page, "navigate-gen_0")
    radio.expect_choices(["all", "g:BC2F1", "none"])
    radio.set("none")

    crumb = page.locator("nav[aria-label=breadcrumb]")
    expect(crumb.locator("li[aria-current=page]")).to_have_text("(no generation)")

    navbar.set("rank")
    controller.InputSwitch(page, "rank-only_pass").set(False)
    table = controller.OutputDataFrame(page, "rank-table")
    table.expect_nrow(2)
    sample_ids = {
        table.cell_locator(0, SAMPLE_ID_COL).inner_text().strip(),
        table.cell_locator(1, SAMPLE_ID_COL).inner_text().strip(),
    }
    assert sample_ids == {"BC2F1-F1-004", "BC2F1-F1-007"}
    caption = page.locator("#rank-caption")
    expect(caption).to_contain_text("2 individuals shown (")
    expect(caption).to_contain_text("F1 > (no generation))")


def test_stale_family_is_not_published_on_reload(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    load_fixture(page)

    navbar = controller.PageNavbar(page, "screen")
    navbar.set("navigate")
    families = controller.Accordion(page, "navigate-families")
    families.set("fam_1")
    crumb = page.locator("nav[aria-label=breadcrumb]")
    expect(crumb.locator("li[aria-current=page]")).to_have_text("F2")

    page.evaluate(
        "() => { window.__crumbs = []; const el = document.querySelector('#navigate-crumb'); "
        "new MutationObserver(() => window.__crumbs.push(el.innerText))"
        ".observe(el, {subtree: true, childList: true, characterData: true}); }"
    )

    navbar.set("load")
    controller.InputActionButton(page, "load-run").click()
    expect(page.locator("#load-status")).to_contain_text(LOAD_STATUS, timeout=60_000)

    crumbs = page.evaluate("window.__crumbs")
    assert all("F2" not in c for c in crumbs)

    navbar.set("navigate")
    expect(crumb.locator("li")).to_have_count(1)
    families.expect_open([])
