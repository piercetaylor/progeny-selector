"""The version has one source, ``_version.py``, and the CLI and package metadata report it."""

from __future__ import annotations

import importlib.metadata
import re
import tomllib
from pathlib import Path

import pytest

import progeny_selector
from progeny_selector import __version__
from progeny_selector.cli import main

ROOT = Path(__file__).resolve().parent.parent


def test_pyproject_reads_the_version_from_the_package() -> None:
    with open(ROOT / "pyproject.toml", "rb") as fh:
        data = tomllib.load(fh)
    assert "version" not in data["project"]
    assert data["project"]["dynamic"] == ["version"]
    assert data["tool"]["hatch"]["version"]["path"] == "src/progeny_selector/_version.py"


def test_version_is_a_plain_semver() -> None:
    assert re.fullmatch(r"\d+\.\d+\.\d+", __version__)


def test_installed_metadata_matches_the_package() -> None:
    # A PackageNotFoundError is a failure, not a skip: after any bump, re-run `pip install -e .`.
    assert importlib.metadata.version("progeny-selector") == progeny_selector.__version__


def test_version_flag(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as e:
        main(["--version"])
    assert e.value.code == 0
    assert capsys.readouterr().out == f"progeny-selector {__version__}\n"


def test_no_other_version_literal() -> None:
    hits = [p.name for p in (ROOT / "src" / "progeny_selector").rglob("*.py") if '__version__ = "' in p.read_text(encoding="utf-8")]
    assert hits == ["_version.py"]
