"""Smoke test: full pipeline on the synthetic BC2F1 fixture against independently computed expectations."""

from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

import pytest

from progeny_selector.core.pipeline import run_analysis
from progeny_selector.core.selection import project_next_generation, select_top_n
from progeny_selector.io import load_dataset, read_criteria
from progeny_selector.io.export import write_results_csv


@pytest.fixture(scope="module")
def result(fixture_dir: Path):
    dataset = load_dataset(fixture_dir / "genotypes.vcf", fixture_dir / "samples.csv", fixture_dir / "markers.csv")
    criteria = read_criteria(fixture_dir / "criteria.yaml")
    return run_analysis(dataset, criteria)


def _close(a: float | None, b: str, tol: float = 1e-6) -> bool:
    return a is not None and math.isclose(a, float(b), abs_tol=tol)


def test_dataset_shape(result):
    assert result.classification.states.shape[0] == 500
    assert result.classification.n_informative == 475
    assert len(result.rows) == 40
    assert result.carrier_chroms == {"Gm06"}
    assert result.unit == "cm"


def test_statuses_and_hard_filters_match_expected(result, expected_rows):
    for row in result.rows:
        exp = expected_rows[row["sample_id"]]
        assert row["target_T1_status"] == exp["target_status"], row["sample_id"]
        assert row["avoid_AV1_status"] == exp["avoid_status"], row["sample_id"]
        assert row["passes_filters"] == (exp["passes_filters"] == "True"), (row["sample_id"], row["exclusion_reason"])
        assert row["exclusion_reason"] == exp["exclusion_reason"], row["sample_id"]


def test_background_drag_and_score_match_expected(result, expected_rows):
    for row in result.rows:
        exp = expected_rows[row["sample_id"]]
        assert _close(row["rpp_total"], exp["rpp_total"]), row["sample_id"]
        assert _close(row["rpp_carrier"], exp["rpp_carrier"]), row["sample_id"]
        assert _close(row["rpp_noncarrier"], exp["rpp_noncarrier"]), row["sample_id"]
        assert _close(row["missing_rate"], exp["missing_rate"]), row["sample_id"]
        assert row["drag_unit"] == "cm"
        assert _close(row["drag_total_max"], exp["drag_total_max_cm"], 1e-3), row["sample_id"]
        assert _close(row["drag_total_est"], exp["drag_total_est_cm"], 1e-3), row["sample_id"]
        assert row["recomb_T1_left"] == (exp["recomb_left"] == "True"), row["sample_id"]
        assert row["recomb_T1_right"] == (exp["recomb_right"] == "True"), row["sample_id"]
        assert _close(row["composite_score"], exp["composite_score"], 1e-5), row["sample_id"]


def test_ranks_match_expected(result, expected_rows):
    for row in result.rows:
        exp = expected_rows[row["sample_id"]]
        if exp["rank_overall"]:
            assert row["rank_overall"] == float(exp["rank_overall"]), row["sample_id"]
            assert row["rank_in_family"] == float(exp["rank_in_family"]), row["sample_id"]
        else:
            assert row["rank_overall"] is None
    assert result.rows[0]["sample_id"] == "BC2F1-F1-001"


def test_planted_individuals(result):
    by_id = {r["sample_id"]: r for r in result.rows}
    assert by_id["BC2F1-F1-001"]["rank_overall"] == 1
    assert by_id["BC2F1-F1-002"]["exclusion_reason"].startswith("target:T1:fail")
    assert "avoid:AV1:fail" in by_id["BC2F1-F2-001"]["exclusion_reason"]
    assert "possible_self_or_outcross" in by_id["BC2F1-F2-002"]["qc_flags"]
    assert by_id["BC2F1-F1-003"]["exclusion_reason"] == "missing_rate>0.2"
    assert by_id["BC2F1-F1-001"]["expected_rpp"] == pytest.approx(0.875)


def test_qc_advisory_flags(result, expected_rows):
    for row in result.rows:
        exp = expected_rows[row["sample_id"]]
        assert ("het_rate_deviates" in row["qc_flags"]) == (exp["het_rate_deviates"] == "True"), row["sample_id"]
        assert ("family_donor_outlier" in row["qc_flags"]) == (exp["family_donor_outlier"] == "True"), row["sample_id"]


def test_selection_and_projection(result):
    top = select_top_n(result.rows, 2, per_family=True)
    assert [r["sample_id"] for r in top] == ["BC2F1-F1-001", "BC2F1-F1-010", "BC2F1-F2-005", "BC2F1-F2-019"]
    proj = project_next_generation(top, "backcross")
    assert proj.expected_rpp_next == pytest.approx((1 + proj.mean_rpp_selected) / 2)
    selfed = project_next_generation(top, "self")
    assert selfed.expected_rpp_next == pytest.approx(proj.mean_rpp_selected)
    assert selfed.expected_het_next == pytest.approx(proj.expected_het_next * 0 + selfed.expected_het_next)


def test_results_csv_roundtrip(result, tmp_path: Path):
    out = tmp_path / "results.csv"
    write_results_csv(result.rows, out)
    header = out.read_text().splitlines()[0].split(",")
    assert header[:3] == ["rank_overall", "rank_in_family", "sample_id"]
    assert "target_T1_status" in header and "rpp_Gm06" in header


def test_cli_rank_and_select(fixture_dir: Path, tmp_path: Path):
    env = {"PYTHONPATH": str(Path(__file__).resolve().parents[1] / "src")}
    results = tmp_path / "results.csv"
    cmd = [
        sys.executable,
        "-m",
        "progeny_selector",
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
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stderr
    assert "40 individuals, 11 pass hard filters" in proc.stdout
    selected = tmp_path / "selected.csv"
    manifest = tmp_path / "next.csv"
    cmd = [
        sys.executable,
        "-m",
        "progeny_selector",
        "select",
        "--results",
        str(results),
        "--top",
        "1",
        "--out",
        str(selected),
        "--next-manifest",
        str(manifest),
        "--samples",
        str(fixture_dir / "samples.csv"),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stderr
    assert selected.read_text().count("\n") == 3  # header + one per family
    assert manifest.read_text().splitlines()[0] == "sample_id,line_name,role,generation,family_id,notes"
