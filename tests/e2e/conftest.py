"""Browser tests: every item under tests/e2e/ carries the e2e marker, which addopts excludes by default.

Run with ``pytest -m e2e`` after ``playwright install chromium``. Without Playwright the directory is skipped.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

pytest.importorskip("playwright")

E2E_DIR = Path(__file__).parent
REPO_ROOT = E2E_DIR.parents[1]

# `shiny run` in the app subprocess inherits this environment, so it imports the
# package even without an editable install.
os.environ.setdefault("PYTHONPATH", str(REPO_ROOT / "src"))


@pytest.hookimpl(tryfirst=True)
def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    # A module-level pytestmark does not cover a directory, so mark here before -m deselects.
    for item in items:
        if E2E_DIR in item.path.parents:
            item.add_marker(pytest.mark.e2e)
