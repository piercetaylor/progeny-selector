"""The in-memory export wrappers produce the same bytes as the path writers, CRLF line endings included."""

from __future__ import annotations

from pathlib import Path

import pytest

from progeny_selector.core.pipeline import run_analysis
from progeny_selector.io import load_dataset, read_criteria
from progeny_selector.io.export import (
    next_round_manifest_text,
    results_csv_text,
    selection_csv_text,
    write_next_round_manifest,
    write_results_csv,
    write_selection_csv,
)


@pytest.fixture(scope="module")
def loaded(fixture_dir: Path):
    dataset = load_dataset(fixture_dir / "genotypes.vcf", fixture_dir / "samples.csv", fixture_dir / "markers.csv")
    return dataset, run_analysis(dataset, read_criteria(fixture_dir / "criteria.yaml"))


def test_results_csv_text_matches_file(loaded, tmp_path: Path) -> None:
    _dataset, result = loaded
    out = tmp_path / "results.csv"
    write_results_csv(result.rows, out)
    text = results_csv_text(result.rows)
    assert text.encode("utf-8") == out.read_bytes()
    assert "\r\n" in text


def test_results_csv_text_empty_matches_file(tmp_path: Path) -> None:
    out = tmp_path / "results.csv"
    write_results_csv([], out)
    assert results_csv_text([]).encode("utf-8") == out.read_bytes()


def test_selection_csv_text_matches_file(loaded, tmp_path: Path) -> None:
    _dataset, result = loaded
    rows = [result.row("BC2F1-F1-001"), result.row("BC2F1-F2-005")]
    notes = {"BC2F1-F1-001": "keep, vigorous; énergie", "BC2F1-F2-005": 'tall "lodging" risk\nrecheck at R2'}
    out = tmp_path / "selected.csv"
    write_selection_csv(rows, out, notes)
    text = selection_csv_text(rows, notes)
    assert text.encode("utf-8") == out.read_bytes()
    assert '"tall ""lodging"" risk\nrecheck at R2"' in text


def test_selection_csv_text_empty_matches_file(tmp_path: Path) -> None:
    out = tmp_path / "selected.csv"
    write_selection_csv([], out)
    text = selection_csv_text([])
    assert text.encode("utf-8") == out.read_bytes()
    assert text.count("\r\n") == 1


def test_next_round_manifest_text_matches_file(loaded, tmp_path: Path) -> None:
    dataset, result = loaded
    rows = [result.row("BC2F1-F1-001"), result.row("BC2F1-F2-005")]
    out = tmp_path / "next_samples.csv"
    write_next_round_manifest(rows, out, "BC3F1", dataset.recurrent_parent, dataset.donor_parent)
    text = next_round_manifest_text(rows, "BC3F1", dataset.recurrent_parent, dataset.donor_parent)
    assert text.encode("utf-8") == out.read_bytes()
    assert text.split("\r\n")[0] == "sample_id,line_name,role,generation,family_id,notes"


def test_next_round_manifest_text_empty_matches_file(loaded, tmp_path: Path) -> None:
    dataset, _result = loaded
    out = tmp_path / "next_samples.csv"
    write_next_round_manifest([], out, "BC3F1", dataset.recurrent_parent, dataset.donor_parent)
    text = next_round_manifest_text([], "BC3F1", dataset.recurrent_parent, dataset.donor_parent)
    assert text.encode("utf-8") == out.read_bytes()
    assert text.count("\r\n") == 3  # header and both parents
