"""Load screen end to end: upload the synthetic fixture, run the analysis, read the status summary."""

from __future__ import annotations

from pathlib import Path

from playwright.sync_api import Page, expect
from shiny.playwright import controller
from shiny.pytest import create_app_fixture
from shiny.run import ShinyAppProc

FIXTURE = Path(__file__).parents[1] / "fixtures" / "synthetic_bc2f1"

app = create_app_fixture(["../../src/progeny_selector/app/app.py"])


def test_load_fixture(page: Page, app: ShinyAppProc) -> None:
    page.goto(app.url)
    controller.InputFile(page, "load-genotypes").set(FIXTURE / "genotypes.vcf")
    controller.InputFile(page, "load-samples").set(FIXTURE / "samples.csv")
    controller.InputFile(page, "load-markers").set(FIXTURE / "markers.csv")
    controller.InputFile(page, "load-criteria").set(FIXTURE / "criteria.yaml")
    controller.InputActionButton(page, "load-run").click()
    expect(page.locator("#load-status")).to_contain_text("500 markers, 40 progeny; 11 pass hard filters", timeout=60_000)
