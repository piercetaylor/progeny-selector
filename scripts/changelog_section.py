#!/usr/bin/env python3
"""Extract one released section from CHANGELOG.md and check CITATION.cff against it.

Responsibility: give the release workflow its notes and refuse a tag whose CHANGELOG
section is undated or whose CITATION.cff disagrees. Reads only; writes the notes file
when --out is given.

Usage: python scripts/changelog_section.py VERSION [--changelog CHANGELOG.md]
       [--out FILE] [--check-citation CITATION.cff]
Exit codes: 0 ok; 1 missing or undated section, or a CITATION.cff mismatch.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

DATE_RE = re.compile(r"\d{4}-\d{2}-\d{2}")
BOUNDARY_RE = re.compile(r"^(## \[|\[[^\]]+\]:\s)")


def section(text: str, version: str) -> tuple[str, str]:
    """Return (date, body) of the section `## [VERSION] - DATE`; SystemExit(1) if absent or undated."""
    heading = re.compile(rf"^## \[{re.escape(version)}\] - (?P<date>.+?)\s*$")
    lines = text.splitlines()
    start = None
    date = ""
    for i, line in enumerate(lines):
        m = heading.match(line)
        if m:
            start, date = i, m.group("date")
            break
    if start is None:
        print(f'CHANGELOG.md: no section "## [{version}] - YYYY-MM-DD"', file=sys.stderr)
        raise SystemExit(1)
    if not DATE_RE.fullmatch(date):
        print(f'CHANGELOG.md: section [{version}] has date "{date}", not YYYY-MM-DD', file=sys.stderr)
        raise SystemExit(1)
    end = len(lines)
    for j in range(start + 1, len(lines)):
        if BOUNDARY_RE.match(lines[j]):
            end = j
            break
    return date, "\n".join(lines[start + 1 : end]).strip()


def _field(text: str, key: str) -> str | None:
    m = re.search(rf"^{re.escape(key)}:\s*(.*?)\s*$", text, re.MULTILINE)
    if m is None:
        return None
    return m.group(1).strip("\"'")


def check_citation(cff_text: str, version: str, date: str) -> None:
    """SystemExit(1) unless version matches and date-released is absent or equals date."""
    cff_version = _field(cff_text, "version")
    if cff_version != version:
        print(f'CITATION.cff: version "{cff_version}" does not match "{version}"', file=sys.stderr)
        raise SystemExit(1)
    released = _field(cff_text, "date-released")
    if released is not None and released != date:
        print(f'CITATION.cff: date-released "{released}" does not match CHANGELOG date "{date}"', file=sys.stderr)
        raise SystemExit(1)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Extract a released CHANGELOG section.")
    p.add_argument("version")
    p.add_argument("--changelog", default="CHANGELOG.md")
    p.add_argument("--out")
    p.add_argument("--check-citation", metavar="CITATION.cff")
    args = p.parse_args(argv)
    date, body = section(Path(args.changelog).read_text(encoding="utf-8"), args.version)
    if args.check_citation:
        check_citation(Path(args.check_citation).read_text(encoding="utf-8"), args.version, date)
    if args.out:
        Path(args.out).write_text(body + "\n", encoding="utf-8", newline="\n")
    else:
        sys.stdout.write(body + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
