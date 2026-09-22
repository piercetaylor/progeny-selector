"""results.csv schema 1.1.0: the fixed column list, NA for missing values, a header with no rows (docs/adr/0016).

The column list is written out here rather than imported, so a change to ``FIXED_COLUMNS`` has to be
made twice, deliberately, and the schema version bumped with it.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from progeny_selector.cli import _read_results
from progeny_selector.core.pipeline import run_analysis
from progeny_selector.core.selection import select_top_n
from progeny_selector.io import load_dataset, read_criteria
from progeny_selector.io.export import (
    SELECTION_COLUMNS,
    next_round_manifest_text,
    results_csv_text,
    selection_csv_text,
    write_results_csv,
)
from progeny_selector.model.criteria import BackgroundOptions, Criteria, TargetSpec
from progeny_selector.model.dataset import DataContractError, GenotypeMatrix, Marker, Sample
from tests.conftest import make_dataset, make_matrix
from tests.test_map_warnings import without_cm

# The fixed prefix of results.csv: 35 names, ``crop`` (contract 1.5.0) last of them, then the
# dynamic per-locus and per-chromosome columns, then ``token_profile`` (contract 1.4.0), which is
# the 36th name of the empty header.
FIXED_PREFIX = [
    "rank_overall",
    "rank_in_family",
    "sample_id",
    "line_name",
    "family_id",
    "generation",
    "passes_filters",
    "exclusion_reason",
    "composite_score",
    "foreground_all_pass",
    "avoid_all_pass",
    "rpp_total",
    "rpp_carrier",
    "rpp_noncarrier",
    "expected_rpp",
    "drag_total_est",
    "drag_total_max",
    "drag_unit",
    "ibs_rp",
    "ibs_donor",
    "missing_rate",
    "het_rate",
    "expected_het",
    "qc_flags",
    "role",
    "n_informative_called",
    "frac_a",
    "frac_h",
    "frac_b",
    "background_model",
    "background_unit",
    "rank_mode",
    "assembly",
    "results_schema",
    "crop",
]
EMPTY_HEADER = [*FIXED_PREFIX, "token_profile"]


@pytest.fixture(scope="module")
def result(fixture_dir: Path):
    dataset = load_dataset(fixture_dir / "genotypes.vcf", fixture_dir / "samples.csv", fixture_dir / "markers.csv")
    return run_analysis(dataset, read_criteria(fixture_dir / "criteria.yaml"))


def test_fixed_prefix_is_35_names() -> None:
    assert len(FIXED_PREFIX) == 35
    assert len(set(EMPTY_HEADER)) == len(EMPTY_HEADER)


def test_header_only_with_no_rows() -> None:
    assert results_csv_text([]) == ",".join(EMPTY_HEADER) + "\r\n"


def test_fixture_header_prefix_and_trailing_token_profile(result) -> None:
    header = results_csv_text(result.rows).split("\r\n")[0].split(",")
    assert header[: len(FIXED_PREFIX)] == FIXED_PREFIX
    assert header[-1] == "token_profile"


def test_metadata_columns_on_every_row(result) -> None:
    for row in result.rows:
        assert row["results_schema"] == "1.1.0"
        assert row["background_model"] == "count"
        assert row["background_unit"] == "cm"
        assert row["rank_mode"] == "weighted"
        assert row["assembly"] == "Wm82.a4"


def test_no_literal_nan_and_na_for_excluded_rank(result) -> None:
    text = results_csv_text(result.rows)
    assert ",nan," not in text
    assert ",nan\r\n" not in text
    rows = {r["sample_id"]: r for r in csv.DictReader(text.splitlines())}
    excluded = rows["BC2F1-F1-002"]
    assert excluded["passes_filters"] == "FALSE"
    assert excluded["rank_overall"] == "NA"
    assert excluded["rank_in_family"] == "NA"
    # empty text is an empty cell, never NA
    assert all(r["qc_flags"] != "NA" for r in rows.values())
    assert any(r["qc_flags"] == "" for r in rows.values())
    assert rows["BC2F1-F1-001"]["exclusion_reason"] == ""


def test_selection_csv_header_and_na(result) -> None:
    row = dict(result.row("BC2F1-F1-001"))
    row["family_id"] = None
    row["generation"] = None
    lines = selection_csv_text([row]).split("\r\n")
    assert lines[0] == ",".join(SELECTION_COLUMNS)
    cells = next(csv.reader([lines[1]]))
    assert cells[SELECTION_COLUMNS.index("family_id")] == "NA"
    assert cells[SELECTION_COLUMNS.index("generation")] == "NA"
    assert cells[SELECTION_COLUMNS.index("notes")] == ""
    assert cells[SELECTION_COLUMNS.index("results_schema")] == "1.1.0"


def test_missing_qc_numbers_are_na_not_nan() -> None:
    """The pipeline routes het_rate, expected_het and expected_rpp through _num, so a NaN or None writes NA."""
    gm = make_matrix(["A", "N", "N", "N", "N", "N"])
    dataset = make_dataset(gm, generation="")  # no generation: expected_het and expected_rpp are unknown
    criteria = Criteria(targets=[TargetSpec("T1", marker_id="m1", required_state="either")])
    result = run_analysis(dataset, criteria)
    row = result.row("P1")
    assert row["expected_rpp"] is None and row["expected_het"] is None
    text = results_csv_text(result.rows)
    header = text.split("\r\n")[0].split(",")
    cells = text.split("\r\n")[1].split(",")
    cell = dict(zip(header, cells, strict=False))
    assert cell["expected_rpp"] == "NA" and cell["expected_het"] == "NA"
    assert "nan" not in text

    # het_rate itself is NaN when no informative marker is called.
    gm_all_missing = make_matrix(["N"] * 6)
    result = run_analysis(make_dataset(gm_all_missing), criteria)
    assert np.isnan(result.qc[0].het_rate)
    assert result.row("P1")["het_rate"] is None
    text = results_csv_text(result.rows)
    assert ",NA," in text and "nan" not in text


def test_background_unit_records_the_bp_fallback() -> None:
    """A cM request the map cannot honour is recorded on the row, not only in the warnings."""
    criteria = Criteria(
        targets=[TargetSpec("T1", marker_id="m2", required_state="either")],
        background=BackgroundOptions(model="count", map_unit="cm"),
    )
    result = run_analysis(make_dataset(without_cm(make_matrix(["A", "H", "B", "A", "A", "H"]))), criteria)
    assert result.unit == "bp"
    assert all(r["background_unit"] == "bp" for r in result.rows)


def test_chromosome_colliding_with_a_fixed_column_is_an_error() -> None:
    """normalize_chrom keeps an unknown name verbatim, so rpp_total could be overwritten by a chromosome."""
    gm = make_matrix(["A", "H", "B", "A", "A", "H"])
    markers = [Marker(m.marker_id, "total", m.pos_bp, m.cm) for m in gm.markers]
    gm = GenotypeMatrix(markers=markers, sample_ids=list(gm.sample_ids), alleles=gm.alleles, calls=gm.calls)
    criteria = Criteria(targets=[TargetSpec("T1", marker_id="m2", required_state="either")])
    with pytest.raises(DataContractError, match="rpp_total collides"):
        run_analysis(make_dataset(gm), criteria)


def test_a_line_named_na_survives_rank_to_select(tmp_path: Path) -> None:
    """NA means missing only where a missing value is meaningful: an id or a line name is read as written."""
    rows = [
        {
            "sample_id": "NA",
            "line_name": "NA",
            "family_id": None,
            "generation": "BC2F1",
            "rank_overall": 1.0,
            "rank_in_family": None,
            "passes_filters": True,
            "composite_score": 0.5,
            "rpp_total": 0.9,
            "frac_a": 0.9,
            "frac_h": 0.1,
            "frac_b": 0.0,
        }
    ]
    out = tmp_path / "results.csv"
    write_results_csv(rows, out)
    back = _read_results(str(out))
    assert back[0]["sample_id"] == "NA"
    assert back[0]["line_name"] == "NA"
    assert back[0]["family_id"] is None
    assert back[0]["rank_in_family"] is None
    rp = Sample("RP", "RP", "recurrent_parent")
    donor = Sample("DONOR", "DONOR", "donor_parent")
    manifest = next_round_manifest_text(select_top_n(back, 1, per_family=False), "BC3F1", rp, donor)
    assert "NA-BC3F1-001" in manifest
    assert "None-BC3F1-001" not in manifest
