"""Shared browser-test steps."""

from __future__ import annotations

from pathlib import Path

from playwright.sync_api import Page, expect
from shiny.playwright import controller

FIXTURE = Path(__file__).parents[1] / "fixtures" / "synthetic_bc2f1"

LOAD_STATUS = "500 markers, 40 progeny; 11 pass hard filters"


def load_fixture(page: Page) -> None:
    """Upload the synthetic fixture on the Load screen, run the analysis and wait for the status summary."""
    controller.InputFile(page, "load-genotypes").set(FIXTURE / "genotypes.vcf")
    controller.InputFile(page, "load-samples").set(FIXTURE / "samples.csv")
    controller.InputFile(page, "load-markers").set(FIXTURE / "markers.csv")
    controller.InputFile(page, "load-criteria").set(FIXTURE / "criteria.yaml")
    controller.InputActionButton(page, "load-run").click()
    expect(page.locator("#load-status")).to_contain_text(LOAD_STATUS, timeout=60_000)
