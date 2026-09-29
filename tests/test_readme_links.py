"""README links must be absolute (PyPI renders them) and point at files that exist."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
README = (ROOT / "README.md").read_text(encoding="utf-8")
BASE = re.compile(r"^https://github\.com/piercetaylor/progeny-selector/(?:blob|tree)/main/(.*)$")


def _targets() -> list[str]:
    return re.findall(r"\[[^\]]*\]\(([^)\s]+)\)", README)


def _slugs() -> set[str]:
    slugs = set()
    for line in README.splitlines():
        m = re.match(r"^#+\s+(.*?)\s*$", line)
        if m:
            slugs.add(re.sub(r"[^a-z0-9-]", "", m.group(1).lower().replace(" ", "-")))
    return slugs


def test_readme_links_are_absolute_and_resolve() -> None:
    targets = _targets()
    assert targets
    slugs = _slugs()
    for t in targets:
        assert t.startswith(("https://", "#")), f"relative link: {t}"
        if t.startswith("#"):
            assert t[1:] in slugs, f"no heading for anchor: {t}"
            continue
        m = BASE.match(t)
        if m:
            path = m.group(1).split("#", 1)[0]
            assert (ROOT / path).exists(), f"missing path: {path}"
