"""Checks that the GitHub issue forms and dependabot config parse and keep their required content."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent
TEMPLATES = ROOT / ".github" / "ISSUE_TEMPLATE"
FORMS = ["bug_report.yml", "dataset_or_format_question.yml"]
WARNING = "Never attach real genotype data"


def _load(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def test_every_template_file_is_covered() -> None:
    found = {p.name for p in TEMPLATES.glob("*.yml")}
    assert found == {*FORMS, "config.yml"}


@pytest.mark.parametrize("name", FORMS)
def test_issue_form_shape(name: str) -> None:
    data = _load(TEMPLATES / name)
    assert data["name"]
    assert data["description"]
    assert isinstance(data["body"], list) and data["body"]


def test_bug_report_text() -> None:
    text = (TEMPLATES / "bug_report.yml").read_text(encoding="utf-8")
    assert "progeny-selector --version" in text
    body = _load(TEMPLATES / "bug_report.yml")["body"]
    assert any(WARNING in item.get("attributes", {}).get("value", "") for item in body if item["type"] == "markdown")


def test_config_blank_issues() -> None:
    assert "blank_issues_enabled" in _load(TEMPLATES / "config.yml")


def test_dependabot_single_actions_update() -> None:
    data = _load(ROOT / ".github" / "dependabot.yml")
    assert len(data["updates"]) == 1
    assert data["updates"][0]["package-ecosystem"] == "github-actions"
