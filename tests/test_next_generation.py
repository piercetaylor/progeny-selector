"""The next-round manifest refuses a blank next-generation label and trims a padded one.

A blank label used to be stamped straight into placeholder ids (``BC2F1-F1-001--001``) and into an
empty ``generation`` column, in a file that looked written. The rule lives in
``io.export.check_next_generation``, where both the CLI and the Export screen arrive; the CLI also
reports it at the argparse level, which is what gives ``--next-generation ""`` a usage exit code.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from progeny_selector.cli import main
from progeny_selector.core.pipeline import run_analysis
from progeny_selector.io import load_dataset, read_criteria
from progeny_selector.io.export import (
    next_round_manifest_text,
    write_next_round_manifest,
    write_results_csv,
)
from progeny_selector.model.dataset import DataContractError


@pytest.fixture(scope="module")
def loaded(fixture_dir: Path):
    dataset = load_dataset(fixture_dir / "genotypes.vcf", fixture_dir / "samples.csv", fixture_dir / "markers.csv")
    return dataset, run_analysis(dataset, read_criteria(fixture_dir / "criteria.yaml"))


@pytest.fixture(scope="module")
def parents(loaded):
    dataset, _result = loaded
    return dataset.recurrent_parent, dataset.donor_parent


def _progeny_rows(text: str) -> list[list[str]]:
    return [r for r in csv.reader(text.splitlines()) if len(r) > 2 and r[2] == "progeny"]


def test_empty_label_is_refused(loaded, parents) -> None:
    _dataset, result = loaded
    rp, donor = parents
    with pytest.raises(DataContractError, match="must not be empty"):
        next_round_manifest_text(result.rows[:1], "", rp, donor)


def test_whitespace_label_is_refused(loaded, parents) -> None:
    _dataset, result = loaded
    rp, donor = parents
    with pytest.raises(DataContractError, match="must not be empty"):
        next_round_manifest_text(result.rows[:1], "  \t", rp, donor)


def test_none_label_is_refused(loaded, parents) -> None:
    """A cleared Shiny field can hand the writer ``None``; without the guard it becomes the label "None"."""
    _dataset, result = loaded
    rp, donor = parents
    with pytest.raises(DataContractError, match="must not be empty"):
        next_round_manifest_text(result.rows[:1], None, rp, donor)  # type: ignore[arg-type]


def test_padded_label_is_trimmed_in_the_id_and_the_generation_cell(loaded, parents) -> None:
    """The discriminating input: without the strip the placeholder id carries the spaces."""
    _dataset, result = loaded
    rp, donor = parents
    row = result.rows[0]
    sid = row["sample_id"]
    text = next_round_manifest_text([row], " BC3F1 ", rp, donor)
    progeny = _progeny_rows(text)
    assert len(progeny) == 1
    assert progeny[0][0] == f"{sid}-BC3F1-001"
    assert progeny[0][3] == "BC3F1"


def test_zero_is_a_valid_label(loaded, parents) -> None:
    """``"0"`` is non-blank, so it is a label like any other: the check is on blankness, not on truthiness."""
    _dataset, result = loaded
    rp, donor = parents
    row = result.rows[0]
    progeny = _progeny_rows(next_round_manifest_text([row], "0", rp, donor))
    assert progeny[0][0] == f"{row['sample_id']}-0-001"
    assert progeny[0][3] == "0"


def test_inner_space_survives_trimming(loaded, parents) -> None:
    """Q2: only the ends are trimmed. A label with an inner space is written verbatim, unparseable or not."""
    _dataset, result = loaded
    rp, donor = parents
    row = result.rows[0]
    progeny = _progeny_rows(next_round_manifest_text([row], "  BC3 F1  ", rp, donor))
    assert progeny[0][0] == f"{row['sample_id']}-BC3 F1-001"
    assert progeny[0][3] == "BC3 F1"


def test_written_file_trims_the_label_like_the_text_writer(loaded, parents, tmp_path: Path) -> None:
    """The path writer and the text wrapper must agree: both strip, neither writes the padding."""
    _dataset, result = loaded
    rp, donor = parents
    row = result.rows[0]
    out = tmp_path / "next_samples.csv"
    write_next_round_manifest([row], out, " BC3F1 ", rp, donor)
    text = out.read_bytes().decode("utf-8")
    progeny = _progeny_rows(text)
    assert progeny[0][0] == f"{row['sample_id']}-BC3F1-001"
    assert progeny[0][3] == "BC3F1"
    assert text == next_round_manifest_text([row], " BC3F1 ", rp, donor)


def test_refused_write_leaves_no_file(loaded, parents, tmp_path: Path) -> None:
    _dataset, result = loaded
    rp, donor = parents
    with pytest.raises(DataContractError, match="must not be empty"):
        write_next_round_manifest(result.rows, tmp_path / "m.csv", "", rp, donor)
    assert not (tmp_path / "m.csv").exists()


def test_refused_count_also_leaves_no_file(loaded, parents, tmp_path: Path) -> None:
    _dataset, result = loaded
    rp, donor = parents
    with pytest.raises(DataContractError, match="1 or more"):
        write_next_round_manifest(result.rows, tmp_path / "m2.csv", "BC3F1", rp, donor, 0)
    assert not (tmp_path / "m2.csv").exists()


def test_cli_blank_label_is_a_usage_error(loaded, fixture_dir: Path, tmp_path: Path) -> None:
    _dataset, result = loaded
    results_csv = tmp_path / "results.csv"
    write_results_csv(result.rows, results_csv)
    out = tmp_path / "selected.csv"
    m = tmp_path / "next_samples.csv"
    with pytest.raises(SystemExit) as excinfo:
        main(
            [
                "select",
                "--results",
                str(results_csv),
                "--top",
                "1",
                "--out",
                str(out),
                "--next-manifest",
                str(m),
                "--next-generation",
                "  ",
                "--samples",
                str(fixture_dir / "samples.csv"),
            ]
        )
    assert excinfo.value.code == 2
    assert not m.exists()
