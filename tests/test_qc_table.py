"""QC table rows and the uninformative-marker summary on the synthetic fixture."""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from progeny_selector.core.pipeline import run_analysis
from progeny_selector.core.qc import qc_table_rows, uninformative_summary
from progeny_selector.core.score import qc_excluding_hits
from progeny_selector.io import load_dataset, read_criteria
from progeny_selector.model.criteria import Filters


@pytest.fixture(scope="module")
def fixture_data(fixture_dir: Path):
    dataset = load_dataset(fixture_dir / "genotypes.vcf", fixture_dir / "samples.csv", fixture_dir / "markers.csv")
    criteria = read_criteria(fixture_dir / "criteria.yaml")
    result = run_analysis(dataset, criteria)
    return dataset, result, criteria


def test_qc_table_rows(fixture_data):
    dataset, result, criteria = fixture_data
    rows = qc_table_rows(result.qc, dataset, criteria.filters)
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


def test_qc_excluded_follows_filters(fixture_data):
    dataset, result, criteria = fixture_data
    sample = "BC2F1-F2-002"
    assert criteria.filters.exclude_qc_flagged is True
    assert result.row(sample)["passes_filters"] is False
    assert "qc:" in result.row(sample)["exclusion_reason"]
    assert {r["sample_id"]: r for r in qc_table_rows(result.qc, dataset, criteria.filters)}[sample]["qc_excluded"] is True

    lenient = dataclasses.replace(criteria, filters=dataclasses.replace(criteria.filters, exclude_qc_flagged=False))
    relaxed = run_analysis(dataset, lenient)
    assert relaxed.row(sample)["passes_filters"] is True
    assert relaxed.row(sample)["exclusion_reason"] == ""
    rows = {r["sample_id"]: r for r in qc_table_rows(relaxed.qc, dataset, lenient.filters)}
    assert rows[sample]["qc_excluded"] is False
    assert not any(r["qc_excluded"] for r in rows.values())


def test_uninformative_summary(fixture_data):
    _, result, _ = fixture_data
    summary = uninformative_summary(result.classification)
    assert sum(n for _, n in summary) == 25
    assert summary[0] == ("parents identical (monomorphic)", 20)
    assert all(reason for reason, _ in summary)
    counts = [n for _, n in summary]
    assert counts == sorted(counts, reverse=True)


def test_qc_excluding_hits_edge_cases():
    on, off = Filters(exclude_qc_flagged=True), Filters(exclude_qc_flagged=False)
    assert qc_excluding_hits([], on) == []
    assert qc_excluding_hits(["high_missing", "family_donor_outlier"], on) == []
    assert qc_excluding_hits(["high_missing", "possible_outcross"], on) == ["possible_outcross"]
    assert qc_excluding_hits(["possible_outcross"], off) == []
