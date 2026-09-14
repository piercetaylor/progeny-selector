"""Browser tests: every item under tests/e2e/ carries the e2e marker, which addopts excludes by default.

Run with ``pytest -m e2e`` after ``playwright install chromium``. Without Playwright the directory is skipped.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

pytest.importorskip("playwright")

from playwright.sync_api import Page, expect

E2E_DIR = Path(__file__).parent
REPO_ROOT = E2E_DIR.parents[1]

# `shiny run` in the app subprocess inherits this environment, so it imports the
# package even without an editable install.
os.environ.setdefault("PYTHONPATH", str(REPO_ROOT / "src"))

E2E_TIMEOUT_MS = int(os.environ.get("PS_E2E_TIMEOUT_MS", "30000"))


@pytest.fixture(autouse=True)
def _tolerant_timeouts(page: Page) -> None:
    # Playwright's 5 s default fails on a loaded laptop; shiny controllers pass timeout=None, which falls back to these.
    expect.set_options(timeout=E2E_TIMEOUT_MS)
    page.set_default_timeout(E2E_TIMEOUT_MS)


@pytest.hookimpl(tryfirst=True)
def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    # A module-level pytestmark does not cover a directory, so mark here before -m deselects.
    for item in items:
        if E2E_DIR in item.path.parents:
            item.add_marker(pytest.mark.e2e)
