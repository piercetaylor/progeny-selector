"""The results.csv ``assembly`` column records the table the run used (docs/adr/0015, amendment 2026-09-22).

The chromosome-length tables are Williams 82's, so ``chrom_length_bp`` returns its fallback under every
scheme but soybean: a maize run uses no length table at all, and must not record ``Wm82.a4``. Soybean is
the default, so a soybean-only test cannot tell the new behaviour from the old one; the maize cases are
the discriminators.
"""

from __future__ import annotations

import dataclasses
import re
from pathlib import Path

import pytest

from progeny_selector.core.chrom import SOYBEAN_SCHEME, chrom_length_bp, compile_scheme
from progeny_selector.core.pipeline import resolve_assembly, run_analysis
from progeny_selector.io import load_dataset, read_criteria
from progeny_selector.io.crops import resolve_crop
from progeny_selector.model.criteria import Criteria, CriteriaError, TargetSpec

ROOT = Path(__file__).resolve().parent.parent
CASES = ROOT / "contract" / "cases"

CRITERIA = """name: assembly case
targets:
  - locus_id: T1
    marker_id: r1
    required_state: either
background:
  model: count
  map_unit: bp
filters:
  max_missing_rate: 0.5
"""


def _run(tmp_path: Path, crop: str, assembly_line: str = ""):
    path = tmp_path / "criteria.yaml"
    path.write_text(CRITERIA + assembly_line, encoding="utf-8")
    case = CASES / f"crop-{crop}-spellings"
    dataset = load_dataset(case / "genotypes.csv", case / "samples.csv", crop=crop)
    return run_analysis(dataset, read_criteria(path))


def test_unset_assembly_under_maize_records_none(tmp_path: Path) -> None:
    result = _run(tmp_path, "maize")
    assert result.rows
    assert result.assembly == "none"
    assert all(row["assembly"] == "none" for row in result.rows)


def test_unset_assembly_under_soybean_still_records_the_default(tmp_path: Path) -> None:
    result = _run(tmp_path, "soybean")
    assert result.rows
    assert result.assembly == "Wm82.a4"
    assert all(row["assembly"] == "Wm82.a4" for row in result.rows)


def test_explicit_wm82_under_maize_is_refused(tmp_path: Path) -> None:
    with pytest.raises(CriteriaError, match=re.escape("assembly 'Wm82.a2'")):
        _run(tmp_path, "maize", "assembly: Wm82.a2\n")


def test_explicit_none_under_maize_is_accepted(tmp_path: Path) -> None:
    result = _run(tmp_path, "maize", "assembly: none\n")
    assert result.assembly == "none"
    assert all(row["assembly"] == "none" for row in result.rows)


def test_explicit_wm82_under_soybean_is_unchanged(tmp_path: Path) -> None:
    result = _run(tmp_path, "soybean", "assembly: Wm82.a2\n")
    assert result.assembly == "Wm82.a2"
    assert all(row["assembly"] == "Wm82.a2" for row in result.rows)


def test_the_maize_warning_names_the_resolved_table(tmp_path: Path) -> None:
    """The per-chromosome fallback warning quotes the same value the column records."""
    result = _run(tmp_path, "maize")
    assert any("the assembly none gives no length" in w for w in result.warnings)
    assert not any("Wm82" in w for w in result.warnings)


def test_resolution_follows_the_identity_test_chrom_length_bp_makes() -> None:
    """A scheme that only calls itself soybean gets no table from ``chrom_length_bp``, so it resolves to none."""
    unset = Criteria(targets=[TargetSpec("T1", marker_id="r1")])
    assert resolve_assembly(unset) == "Wm82.a4"
    assert resolve_assembly(unset, resolve_crop("maize")) == "none"

    lookalike = compile_scheme(dataclasses.replace(SOYBEAN_SCHEME, name="Soybean (copy)"))
    assert chrom_length_bp("Gm01", None, "Wm82.a4", lookalike) is None
    assert resolve_assembly(unset, lookalike) == "none"
