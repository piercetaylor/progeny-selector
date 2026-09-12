"""QC table rows and the uninformative-marker summary on the synthetic fixture."""

from __future__ import annotations

from pathlib import Path

import pytest

from progeny_selector.core.pipeline import run_analysis
from progeny_selector.core.qc import qc_table_rows, uninformative_summary
from progeny_selector.io import load_dataset, read_criteria


@pytest.fixture(scope="module")
def fixture_data(fixture_dir: Path):
    dataset = load_dataset(fixture_dir / "genotypes.vcf", fixture_dir / "samples.csv", fixture_dir / "markers.csv")
    result = run_analysis(dataset, read_criteria(fixture_dir / "criteria.yaml"))
    return dataset, result


def test_qc_table_rows(fixture_data):
    dataset, result = fixture_data
    rows = qc_table_rows(result.qc, dataset)
    assert len(rows) == 40
    assert list(rows[0]) == [
        "sample_id",
        "line_name",
        "family_id",
        "generation",
        "missing_rate",
        "het_rate",
        "expected_het",
        "hom_donor_rate",
        "nonparental_rate",
        "expected_rpp",
        "ibs_rp",
        "ibs_donor",
        "flags",
        "qc_excluded",
    ]
    by_id = {r["sample_id"]: r for r in rows}
    assert by_id["BC2F1-F2-002"]["qc_excluded"] is True
    assert "possible_self_or_outcross" in by_id["BC2F1-F2-002"]["flags"].split("|")
    # high_missing is advisory; the missing-rate hard filter, not the QC policy, excludes this plant.
    assert "high_missing" in by_id["BC2F1-F1-003"]["flags"].split("|")
    assert by_id["BC2F1-F1-003"]["qc_excluded"] is False
    assert by_id["BC2F1-F1-001"]["line_name"] == "LF1-001"


def test_uninformative_summary(fixture_data):
    _, result = fixture_data
    summary = uninformative_summary(result.classification)
    assert sum(n for _, n in summary) == 25
    assert summary[0] == ("parents identical (monomorphic)", 20)
    assert all(reason for reason, _ in summary)
    counts = [n for _, n in summary]
    assert counts == sorted(counts, reverse=True)
