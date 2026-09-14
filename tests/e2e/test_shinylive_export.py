"""Shinylive export smoke test: serves the exported site and drives it in a real browser.

Skipped unless ``PS_SITE_DIR`` names a directory produced by
``scripts/build_shinylive.py`` (``$env:PS_SITE_DIR="site"; pytest -m e2e
tests/e2e/test_shinylive_export.py``). Every recorded request must start with the
page's own origin: a request to PyPI or a CDN would mean the staged package
(docs/adr/0009) did not import under Pyodide and the wheel-URL fallback
(resolution 5, docs/m1-phases.md) is needed instead.

The exported app renders inside a same-origin ``iframe`` (``/app_<id>/``), not in
the top-level document, so locators here are scoped through
``page.frame_locator("iframe")``. The ``shiny.playwright.controller`` classes used
by the rest of ``tests/e2e/`` assume no such frame and are not used here.
"""

from __future__ import annotations

import csv
import os
import re
import socket
import subprocess
import sys
import time
from pathlib import Path

import pytest
from playwright.sync_api import FrameLocator, Page, expect

from tests.e2e.helpers import FIXTURE, LOAD_STATUS

pytestmark = pytest.mark.skipif(
    not os.environ.get("PS_SITE_DIR"),
    reason="set PS_SITE_DIR to a directory built by scripts/build_shinylive.py",
)

REPO_ROOT = Path(__file__).resolve().parents[2]
HOST, PORT = "127.0.0.1", 8008
ORIGIN = f"http://{HOST}:{PORT}/"
LOAD_RUN_TIMEOUT_MS = 240_000  # hang guard; Pyodide + numpy + pandas is roughly 13 + 7.5 + 13 MB on first load.
STATUS_TIMEOUT_MS = 180_000


def _site_dir() -> Path:
    raw = Path(os.environ["PS_SITE_DIR"])
    return raw if raw.is_absolute() else REPO_ROOT / raw


@pytest.fixture(scope="module")
def site_server():
    site = _site_dir()
    assert (site / "app.json").exists(), f"{site} is not a shinylive export (no app.json); run scripts/build_shinylive.py"
    proc = subprocess.Popen(
        [sys.executable, "-m", "http.server", "--directory", str(site), "--bind", HOST, str(PORT)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    try:
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline:
            try:
                with socket.create_connection((HOST, PORT), timeout=0.5):
                    break
            except OSError:
                time.sleep(0.1)
        else:
            proc.terminate()
            raise RuntimeError(f"http.server did not start on {HOST}:{PORT}")
        yield
    finally:
        proc.terminate()
        proc.wait(timeout=10)


def _wait_upload_complete(frame: FrameLocator, field_id: str) -> None:
    # Shiny wraps each input_file in a div.shiny-input-container holding a
    # .progress-bar that animates to width: 100% once the upload is registered.
    # Clicking Run before that races the upload under Pyodide (the file input's
    # webworker round trip is slower than the DOM update): observed directly, the
    # status stays "... are required." if only the bar's presence is awaited.
    # shiny.playwright.controller.InputFile.expect_complete() checks the style
    # attribute directly rather than the computed style, which reports pixels for
    # a percentage width; the CSS transition also means "attached" alone is not
    # enough, since the bar exists before it reaches 100%.
    progress = (
        frame.locator(f"#{field_id}")
        .locator("xpath=ancestor::div[contains(@class,'shiny-input-container')][1]//div[contains(@class,'progress-bar')]")
        .first
    )
    expect(progress).to_have_attribute("style", re.compile(r"width:\s*100%"), timeout=30_000)


def _rank_cell(frame: FrameLocator, row: int, col: int = 0):
    return frame.locator(f"#rank-table.html-fill-item div.shiny-data-grid table tbody tr[data-index='{row}'] > td:not(.row-number)").nth(
        col
    )


def test_shinylive_export_runs_pipeline_in_browser(site_server: None, page: Page) -> None:
    requests: list[str] = []
    page.on("request", lambda req: requests.append(req.url))

    t_goto = time.monotonic()
    page.goto(ORIGIN)
    frame = page.frame_locator("iframe").first
    frame.locator("#load-run").wait_for(state="attached", timeout=LOAD_RUN_TIMEOUT_MS)
    print(f"time to #load-run: {time.monotonic() - t_goto:.1f}s")

    frame.locator("#load-genotypes").set_input_files(str(FIXTURE / "genotypes.vcf"))
    frame.locator("#load-samples").set_input_files(str(FIXTURE / "samples.csv"))
    frame.locator("#load-markers").set_input_files(str(FIXTURE / "markers.csv"))
    frame.locator("#load-criteria").set_input_files(str(FIXTURE / "criteria.yaml"))
    for field_id in ("load-genotypes", "load-samples", "load-markers", "load-criteria"):
        _wait_upload_complete(frame, field_id)

    t_click = time.monotonic()
    frame.locator("#load-run").click()
    expect(frame.locator("#load-status")).to_contain_text(LOAD_STATUS, timeout=STATUS_TIMEOUT_MS)
    print(f"time to #load-status: {time.monotonic() - t_click:.1f}s")

    assert requests, "no requests were recorded"
    for url in requests:
        assert url.startswith(ORIGIN), f"request left the page's origin: {url}"

    frame.locator("a.nav-link[data-value='rank']").click()
    expect(frame.locator("#rank-table.html-fill-item div.shiny-data-grid table tbody tr")).to_have_count(11)

    # [0, 2] is not contiguous: it ctrl-clicks rather than shift-clicking (tests/e2e/test_qc_compare.py).
    _rank_cell(frame, 0).click()
    _rank_cell(frame, 2).click(modifiers=["ControlOrMeta"])
    expect(frame.locator("#rank-table.html-fill-item div.shiny-data-grid table tbody tr[aria-selected='true']")).to_have_count(2)

    frame.locator("a.nav-link[data-value='compare']").click()
    expect(frame.locator("#compare-panels .card")).to_have_count(2)

    frame.locator("a.nav-link[data-value='export']").click()
    with page.expect_download() as download_info:
        frame.locator("#export-selected").click()
    download = download_info.value
    with open(download.path(), newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    assert len(rows) == 3, f"expected a header and 2 data rows in selected.csv, got {rows}"

    for url in requests:
        assert url.startswith(ORIGIN), f"request left the page's origin: {url}"
