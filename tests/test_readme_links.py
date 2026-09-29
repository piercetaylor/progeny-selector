"""README links must be absolute (PyPI renders them) and point at files that exist."""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
DOCS = ("README.md", "docs/README.md")
BASE = re.compile(r"^https://github\.com/piercetaylor/progeny-selector/(?:blob|tree)/main/(.*)$")


def _targets(text: str) -> list[str]:
    return re.findall(r"\[[^\]]*\]\(([^)\s]+)\)", text)


def _slugs(text: str) -> set[str]:
    slugs = set()
    for line in text.splitlines():
        m = re.match(r"^#+\s+(.*?)\s*$", line)
        if m:
            slugs.add(re.sub(r"[^a-z0-9-]", "", m.group(1).lower().replace(" ", "-")))
    return slugs


@pytest.mark.parametrize("doc", DOCS)
def test_readme_links_are_absolute_and_resolve(doc: str) -> None:
    text = (ROOT / doc).read_text(encoding="utf-8")
    targets = _targets(text)
    assert targets
    slugs = _slugs(text)
    for t in targets:
        assert t.startswith(("https://", "#")), f"relative link: {t}"
        if t.startswith("#"):
            assert t[1:] in slugs, f"no heading for anchor: {t}"
            continue
        m = BASE.match(t)
        if m:
            path = m.group(1).split("#", 1)[0]
            assert (ROOT / path).exists(), f"missing path: {path}"
            assert not re.fullmatch(r"docs/[^/]*-phases\.md", path), f"phases doc is not in the sdist: {path}"
