"""``select --top N`` is per family unless ``--overall``; the help text says so."""

from __future__ import annotations

import csv
from pathlib import Path

from progeny_selector.cli import build_parser, main


def _rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def test_top_is_per_family_unless_overall(tmp_path: Path, fixture_dir: Path) -> None:
    results = tmp_path / "results.csv"
    rank = [
        "rank",
        "--genotypes",
        str(fixture_dir / "genotypes.vcf"),
        "--samples",
        str(fixture_dir / "samples.csv"),
        "--markers",
        str(fixture_dir / "markers.csv"),
        "--criteria",
        str(fixture_dir / "criteria.yaml"),
        "--out",
        str(results),
    ]
    assert main(rank) == 0

    sel = tmp_path / "sel.csv"
    assert main(["select", "--results", str(results), "--top", "3", "--out", str(sel)]) == 0
    assert len(_rows(sel)) == 6

    assert main(["select", "--results", str(results), "--top", "3", "--overall", "--out", str(sel)]) == 0
    assert len(_rows(sel)) == 3

    assert main(["select", "--results", str(results), "--top", "1", "--out", str(sel)]) == 0
    rows = _rows(sel)
    assert len(rows) == 2
    assert {r["family_id"] for r in rows} == {"F1", "F2"}


def test_top_help_names_per_family() -> None:
    parser = build_parser()
    select = parser._subparsers._group_actions[0].choices["select"]  # type: ignore[union-attr]
    action = next(a for a in select._actions if "--top" in a.option_strings)
    assert action.help is not None
    assert "up to N per family" in action.help
    assert "--overall" in action.help
