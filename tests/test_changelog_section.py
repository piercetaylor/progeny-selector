"""scripts/changelog_section.py: section extraction and the CITATION.cff cross-check."""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location("changelog_section", ROOT / "scripts" / "changelog_section.py")
assert _spec is not None and _spec.loader is not None
cs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cs)

TEXT = """# Changelog

## [Unreleased]

- pending line

## [1.1.0] - 2026-10-02

### Added

- one-one line

## [1.0.0] - 2026-09-30

### Added

- one-zero line

## [0.9.0] - unreleased

- draft line

[Unreleased]: https://example.org/compare/v1.1.0...HEAD
[1.1.0]: https://example.org/compare/v1.0.0...v1.1.0
"""


def test_section_returns_date_and_only_its_lines() -> None:
    date, body = cs.section(TEXT, "1.0.0")
    assert date == "2026-09-30"
    assert "one-zero line" in body
    assert "one-one line" not in body and "draft line" not in body


def test_section_stops_before_older_section() -> None:
    date, body = cs.section(TEXT, "1.1.0")
    assert date == "2026-10-02"
    assert "one-one line" in body
    assert "one-zero line" not in body


def test_last_section_stops_at_link_block() -> None:
    text = TEXT.replace("## [0.9.0] - unreleased", "## [0.9.0] - 2026-09-01")
    _, body = cs.section(text, "0.9.0")
    assert "draft line" in body
    assert "example.org" not in body


def test_unreleased_date_refused(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as e:
        cs.section(TEXT, "0.9.0")
    assert e.value.code == 1
    assert "not YYYY-MM-DD" in capsys.readouterr().err


def test_missing_version_refused(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as e:
        cs.section(TEXT, "2.0.0")
    assert e.value.code == 1
    assert "no section" in capsys.readouterr().err


CFF_BASE = 'cff-version: 1.2.0\ntitle: "x"\nversion: 1.0.0\n'


def test_citation_without_date_released_accepted() -> None:
    cs.check_citation(CFF_BASE, "1.0.0", "2026-09-30")


def test_citation_date_equal_accepted() -> None:
    cs.check_citation(CFF_BASE + "date-released: 2026-09-30\n", "1.0.0", "2026-09-30")
    cs.check_citation(CFF_BASE + 'date-released: "2026-09-30"\n', "1.0.0", "2026-09-30")


def test_citation_date_different_refused(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as e:
        cs.check_citation(CFF_BASE + "date-released: 2026-09-29\n", "1.0.0", "2026-09-30")
    assert e.value.code == 1
    assert "date-released" in capsys.readouterr().err


def test_citation_version_mismatch_refused() -> None:
    with pytest.raises(SystemExit) as e:
        cs.check_citation(CFF_BASE, "1.1.0", "2026-10-02")
    assert e.value.code == 1


def test_main_writes_notes_and_checks_citation(tmp_path: Path) -> None:
    changelog = tmp_path / "CHANGELOG.md"
    cff = tmp_path / "CITATION.cff"
    out = tmp_path / "notes.md"
    changelog.write_text(TEXT, encoding="utf-8")
    cff.write_text(CFF_BASE + "date-released: 2026-09-30\n", encoding="utf-8")
    args = ["1.0.0", "--changelog", str(changelog), "--out", str(out), "--check-citation", str(cff)]
    assert cs.main(args) == 0
    assert "one-zero line" in out.read_text(encoding="utf-8")
    cff.write_text(CFF_BASE + "date-released: 2026-01-01\n", encoding="utf-8")
    with pytest.raises(SystemExit):
        cs.main(args)


def test_changelog_has_a_section_for_the_current_version() -> None:
    from progeny_selector import __version__

    text = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    # The date is not asserted: it is "unreleased" until the maintainer tags.
    assert re.search(rf"^## \[{re.escape(__version__)}\] - .+$", text, re.MULTILINE)


def test_citation_version_matches_the_package() -> None:
    import yaml

    from progeny_selector import __version__

    data = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    assert data["version"] == __version__
    assert data["cff-version"] == "1.2.0"
    assert data["authors"]
    # date-released may be absent before the tag; once present it must be a real date.
    released = data.get("date-released")
    if released is not None:
        assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", str(released))


def test_bracketed_body_line_stays_inside_section() -> None:
    text = TEXT.replace("- one-zero line", "- one-zero line\n[note] x")
    _, body = cs.section(text, "1.0.0")
    assert "[note] x" in body
