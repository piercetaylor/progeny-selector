"""Two-generation round trip: the BC2F1 selection's next_samples.csv is the BC3F1 samples.csv.

The two sides are independent. ``io/export.py`` writes the manifest; ``scripts/make_fixture.py``
builds ``tests/fixtures/synthetic_bc3f1/samples.csv`` from its own simulation without importing the
writer, and generates the BC3F1 genotypes and expectations from the parents' true states. A test that
compares them therefore pins the manifest format with two implementations (docs/adr/0018).
"""

from __future__ import annotations

import csv
import dataclasses
import math
import subprocess
import sys
from pathlib import Path

import pytest

from progeny_selector.core.navigation import build_tree
from progeny_selector.core.pipeline import run_analysis
from progeny_selector.core.selection import select_top_n
from progeny_selector.io import load_dataset, read_criteria, read_samples
from progeny_selector.io.export import MANIFEST_COLUMNS, next_round_manifest_text
from progeny_selector.model.criteria import RankingOptions

BC3F1_DIR = Path(__file__).parent / "fixtures" / "synthetic_bc3f1"
PER_SELECTED = 10
PARENT_ROW = "RP_Williams82,Williams 82 (synthetic),recurrent_parent,,,synthetic recurrent parent"
DONOR_ROW = "DONOR_PI_synthetic,PI synthetic donor,donor_parent,,,synthetic donor"


@pytest.fixture(scope="module")
def bc3f1_dir() -> Path:
    return BC3F1_DIR


@pytest.fixture(scope="module")
def manifest_fixture_text() -> str:
    """The committed BC3F1 samples.csv, read as bytes so a stray CR would show up rather than vanish."""
    return (BC3F1_DIR / "samples.csv").read_bytes().decode("utf-8")


@pytest.fixture(scope="module")
def gen1(fixture_dir: Path):
    dataset = load_dataset(fixture_dir / "genotypes.vcf", fixture_dir / "samples.csv", fixture_dir / "markers.csv")
    criteria = read_criteria(fixture_dir / "criteria.yaml")
    return dataset, criteria, run_analysis(dataset, criteria)


@pytest.fixture(scope="module")
def gen2(fixture_dir: Path, bc3f1_dir: Path):
    """The BC3F1 fixture read with the BC2F1 markers.csv and criteria.yaml, which it does not duplicate."""
    dataset = load_dataset(bc3f1_dir / "genotypes.vcf", bc3f1_dir / "samples.csv", fixture_dir / "markers.csv")
    criteria = read_criteria(fixture_dir / "criteria.yaml")
    return dataset, criteria, run_analysis(dataset, criteria)


@pytest.fixture(scope="module")
def gen2_staged(fixture_dir: Path, bc3f1_dir: Path):
    """The same BC3F1 dataset under ``ranking.mode = staged`` (docs/adr/0007 amendment)."""
    dataset = load_dataset(bc3f1_dir / "genotypes.vcf", bc3f1_dir / "samples.csv", fixture_dir / "markers.csv")
    criteria = read_criteria(fixture_dir / "criteria.yaml")
    return run_analysis(dataset, dataclasses.replace(criteria, ranking=RankingOptions(mode="staged")))


@pytest.fixture(scope="module")
def bc3f1_expected(bc3f1_dir: Path) -> dict[str, dict]:
    with open(bc3f1_dir / "expected_results.csv", newline="") as fh:
        return {r["sample_id"]: r for r in csv.DictReader(fh)}


def _close(a: float | None, b: str, tol: float = 1e-6) -> bool:
    return a is not None and math.isclose(a, float(b), abs_tol=tol)


# (1) the writer against the generator's independent copy


def test_manifest_text_is_the_bc3f1_samples_file(gen1, fixture_dir: Path, manifest_fixture_text: str):
    _dataset, _criteria, result = gen1
    chosen = select_top_n(result.rows, 2, per_family=True)
    assert [r["sample_id"] for r in chosen] == ["BC2F1-F1-001", "BC2F1-F1-010", "BC2F1-F2-005", "BC2F1-F2-019"]
    samples = read_samples(fixture_dir / "samples.csv")
    rp = next(s for s in samples if s.role == "recurrent_parent")
    donor = next(s for s in samples if s.role == "donor_parent")
    text = next_round_manifest_text(chosen, "BC3F1", rp, donor, n_per_selected=PER_SELECTED)
    assert text.replace("\r\n", "\n") == manifest_fixture_text


def test_manifest_keeps_empty_cells_not_na(manifest_fixture_text: str):
    """next_samples.csv is the input contract's samples.csv: a missing cell is empty, never ``NA``.

    backcross's manifest reader would take ``NA`` for a family named "NA" (docs/adr/0016, decision 4),
    so this pins the literal parent rows as well as the absence of the token.
    """
    lines = manifest_fixture_text.rstrip("\n").split("\n")
    assert lines[0] == ",".join(MANIFEST_COLUMNS)
    assert lines[1] == PARENT_ROW
    assert lines[2] == DONOR_ROW
    assert len(lines) == 3 + 4 * PER_SELECTED
    for line in lines:
        assert ",NA," not in line and not line.endswith(",NA"), line
    assert lines[3] == "BC2F1-F1-001-BC3F1-001,LF1-001,progeny,BC3F1,BC2F1-F1-001,derived from BC2F1-F1-001"


# (2) the next generation against its own expectations


def test_bc3f1_shape_and_tree(gen2):
    dataset, _criteria, result = gen2
    assert result.classification.states.shape[0] == 500
    assert result.classification.n_informative == 475
    assert len(result.rows) == 40
    tree = build_tree(dataset)
    assert len(tree.families) == 4
    for fam in tree.families:
        assert fam.n == 10
        assert [g.generation for g in fam.generations] == ["BC3F1"]


def test_bc3f1_statuses_and_filters_match_expected(gen2, bc3f1_expected):
    _dataset, _criteria, result = gen2
    for row in result.rows:
        exp = bc3f1_expected[row["sample_id"]]
        assert row["target_T1_status"] == exp["target_status"], row["sample_id"]
        assert row["avoid_AV1_status"] == exp["avoid_status"], row["sample_id"]
        assert row["passes_filters"] == (exp["passes_filters"] == "True"), (row["sample_id"], row["exclusion_reason"])
        assert row["exclusion_reason"] == exp["exclusion_reason"], row["sample_id"]


def test_bc3f1_metrics_and_ranks_match_expected(gen2, bc3f1_expected):
    _dataset, _criteria, result = gen2
    for row in result.rows:
        exp = bc3f1_expected[row["sample_id"]]
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
        if exp["rank_overall"]:
            assert row["rank_overall"] == float(exp["rank_overall"]), row["sample_id"]
            assert row["rank_in_family"] == float(exp["rank_in_family"]), row["sample_id"]
        else:
            assert row["rank_overall"] is None, row["sample_id"]


def test_bc3f1_staged_ranks_match_expected(gen2, gen2_staged, bc3f1_expected):
    """The staged order of the new generation, against the generator's independent lexicographic sort."""
    _dataset, _criteria, weighted = gen2
    default_scores = {r["sample_id"]: r["composite_score"] for r in weighted.rows}
    ranked = 0
    for row in gen2_staged.rows:
        exp = bc3f1_expected[row["sample_id"]]
        assert row["rank_mode"] == "staged", row["sample_id"]
        if exp["rank_overall_staged"]:
            assert row["rank_overall"] == float(exp["rank_overall_staged"]), row["sample_id"]
            assert row["rank_in_family"] == float(exp["rank_in_family_staged"]), row["sample_id"]
            ranked += 1
        else:
            assert row["rank_overall"] is None, row["sample_id"]
        assert row["composite_score"] == default_scores[row["sample_id"]], row["sample_id"]
    assert ranked == 21  # the passing individuals, so the comparison is not vacuous
    staged_order = [r["sample_id"] for r in gen2_staged.rows if r["rank_overall"] is not None]
    weighted_order = [r["sample_id"] for r in weighted.rows if r["rank_overall"] is not None]
    assert staged_order != weighted_order  # the two modes really do disagree on this generation


def test_bc3f1_generation_expectations_and_flags(gen2, bc3f1_expected):
    _dataset, _criteria, result = gen2
    for row in result.rows:
        exp = bc3f1_expected[row["sample_id"]]
        assert row["generation"] == "BC3F1", row["sample_id"]
        assert row["expected_rpp"] == 0.9375, row["sample_id"]
        assert row["expected_het"] == 0.125, row["sample_id"]
        assert ("het_rate_deviates" in row["qc_flags"]) == (exp["het_rate_deviates"] == "True"), row["sample_id"]
        assert ("family_donor_outlier" in row["qc_flags"]) == (exp["family_donor_outlier"] == "True"), row["sample_id"]
        assert ("possible_duplicate" in row["qc_flags"]) == (exp["possible_duplicate"] == "True"), row["sample_id"]


# (3) the same round trip through the CLI


def _cli(*args: str) -> subprocess.CompletedProcess[str]:
    env = {"PYTHONPATH": str(Path(__file__).resolve().parents[1] / "src")}
    return subprocess.run([sys.executable, "-m", "progeny_selector", *args], capture_output=True, text=True, env=env)


def test_cli_round_trip(fixture_dir: Path, bc3f1_dir: Path, tmp_path: Path, manifest_fixture_text: str):
    results = tmp_path / "results.csv"
    proc = _cli(
        "rank",
        "--genotypes", str(fixture_dir / "genotypes.vcf"),
        "--samples", str(fixture_dir / "samples.csv"),
        "--markers", str(fixture_dir / "markers.csv"),
        "--criteria", str(fixture_dir / "criteria.yaml"),
        "--out", str(results),
    )  # fmt: skip
    assert proc.returncode == 0, proc.stderr
    manifest = tmp_path / "m.csv"
    proc = _cli(
        "select",
        "--results", str(results),
        "--top", "2",
        "--per-selected", str(PER_SELECTED),
        "--out", str(tmp_path / "selected.csv"),
        "--next-manifest", str(manifest),
        "--samples", str(fixture_dir / "samples.csv"),
    )  # fmt: skip
    assert proc.returncode == 0, proc.stderr
    assert manifest.read_bytes().decode("utf-8").replace("\r\n", "\n") == manifest_fixture_text
    proc = _cli(
        "rank",
        "--genotypes", str(bc3f1_dir / "genotypes.vcf"),
        "--samples", str(manifest),
        "--markers", str(fixture_dir / "markers.csv"),
        "--criteria", str(fixture_dir / "criteria.yaml"),
        "--out", str(tmp_path / "results_bc3f1.csv"),
    )  # fmt: skip
    assert proc.returncode == 0, proc.stderr
    assert "40 individuals" in proc.stdout


def test_cli_per_selected_zero_is_a_usage_error(fixture_dir: Path, tmp_path: Path):
    proc = _cli(
        "select",
        "--results", str(tmp_path / "missing.csv"),
        "--top", "2",
        "--per-selected", "0",
        "--out", str(tmp_path / "selected.csv"),
    )  # fmt: skip
    assert proc.returncode == 2
    assert "1 or more" in proc.stderr
