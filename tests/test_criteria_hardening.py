"""criteria.yaml hardening: aliases, document size, non-finite numbers, numeric ids and whole-number positions."""

from __future__ import annotations

import pytest

from progeny_selector.io.criteria import read_criteria_text
from progeny_selector.model.criteria import CriteriaError


def _region_target(**extra: object) -> str:
    fields = {"locus_id": "T1", "chrom": "Gm06", "start_bp": 100, "end_bp": 500, **extra}
    body = "\n".join(f"    {k}: {v}" for k, v in fields.items())
    return "targets:\n  -\n" + body + "\n"


def test_aliases_rejected():
    with pytest.raises(CriteriaError, match="aliases are not accepted"):
        read_criteria_text("a: &x [1]\nb: *x")


def test_anchor_without_alias_loads():
    assert read_criteria_text("name: &x trial\n" + _region_target()).name == "trial"


def test_inline_merge_without_alias_loads():
    crit = read_criteria_text("targets: [{locus_id: T1, marker_id: M1, <<: {min_markers: 1}}]")
    assert crit.targets[0].min_markers == 1


def test_deep_nesting_rejected():
    with pytest.raises(CriteriaError, match="nested too deeply"):
        read_criteria_text("[" * 1500 + "]" * 1500)


def test_document_over_1_mb_rejected():
    text = "name: big\n" + "# " + "x" * 1_150_000 + "\n"
    with pytest.raises(CriteriaError, match="larger than 1 MB"):
        read_criteria_text(text)


def _padded(n_bytes: int) -> str:
    head = "name: big\n" + _region_target() + "# "
    text = head + "x" * (n_bytes - len(head) - 1) + "\n"
    assert len(text.encode("utf-8")) == n_bytes
    return text


def test_document_of_exactly_1_mb_accepted():
    assert read_criteria_text(_padded(1_048_576)).name == "big"


def test_document_one_byte_over_1_mb_rejected():
    with pytest.raises(CriteriaError, match="larger than 1 MB"):
        read_criteria_text(_padded(1_048_577))


@pytest.mark.parametrize("value", [".nan", ".inf", "-.inf"])
def test_non_finite_weight_rejected(value: str):
    with pytest.raises(CriteriaError, match="finite"):
        read_criteria_text(f"weights: {{rpp_noncarrier: {value}}}")


@pytest.mark.parametrize("value", [".nan", ".inf"])
def test_non_finite_filter_rejected(value: str):
    with pytest.raises(CriteriaError, match="finite"):
        read_criteria_text(f"filters: {{max_missing_rate: {value}}}")


@pytest.mark.parametrize("value", [".nan", ".inf"])
def test_non_finite_flank_window_rejected(value: str):
    with pytest.raises(CriteriaError, match="finite"):
        read_criteria_text(f"flank_window: {value}")


@pytest.mark.parametrize("value", [".nan", ".inf"])
def test_non_finite_start_bp_rejected(value: str):
    with pytest.raises(CriteriaError, match="finite"):
        read_criteria_text(_region_target(start_bp=value))


@pytest.mark.parametrize(
    "locus",
    [
        "{locus_id: [T1], marker_id: M1}",
        "{locus_id: T1, marker_id: {a: 1}}",
        "{locus_id: true, marker_id: M1}",
        "{locus_id: null, marker_id: M1}",
        "{locus_id: T1, chrom: [Gm06], start_bp: 1, end_bp: 2}",
    ],
)
def test_non_text_locus_keys_rejected(locus: str):
    with pytest.raises(CriteriaError, match="must be text"):
        read_criteria_text(f"targets: [{locus}]")


def test_numeric_ids_read_as_text():
    crit = read_criteria_text("targets: [{locus_id: 1, marker_id: 12345}]")
    assert crit.targets[0].marker_id == "12345"
    assert crit.targets[0].locus_id == "1"


def test_numeric_flanking_markers_and_chrom_read_as_text():
    crit = read_criteria_text("avoid: [{locus_id: A1, left_marker: 11, right_marker: 12.5}]\n" + _region_target(chrom=6))
    assert crit.avoid[0].left_marker == "11"
    assert crit.avoid[0].right_marker == "12.5"
    assert crit.targets[0].chrom == "6"


def test_fractional_start_bp_rejected():
    with pytest.raises(CriteriaError, match="must be an integer"):
        read_criteria_text(_region_target(start_bp=100.5))


def test_whole_float_start_bp_becomes_int():
    spec = read_criteria_text(_region_target(start_bp="100.0", end_bp="500.0")).targets[0]
    assert spec.start_bp == 100 and type(spec.start_bp) is int
    assert spec.end_bp == 500 and type(spec.end_bp) is int


def test_whole_float_min_markers_accepted():
    spec = read_criteria_text(_region_target(min_markers="2.0")).targets[0]
    assert spec.min_markers == 2 and type(spec.min_markers) is int


@pytest.mark.parametrize("value", ["2.5", "true", ".nan"])
def test_bad_min_markers_rejected(value: str):
    with pytest.raises(CriteriaError, match="must be an integer"):
        read_criteria_text(_region_target(min_markers=value))


def test_whole_float_anchor_bp_accepted():
    spec = read_criteria_text(_region_target(rule="run", anchor_bp="300.0")).targets[0]
    assert spec.anchor_bp == 300 and type(spec.anchor_bp) is int


@pytest.mark.parametrize("value", ["300.5", ".inf"])
def test_bad_anchor_bp_rejected(value: str):
    with pytest.raises(CriteriaError, match="must be an integer"):
        read_criteria_text(_region_target(rule="run", anchor_bp=value))


def test_whole_float_min_run_accepted():
    spec = read_criteria_text(_region_target(rule="run", min_run="2.0")).targets[0]
    assert spec.min_run == 2 and type(spec.min_run) is int


@pytest.mark.parametrize("value", ["2.5", ".nan", ".inf"])
def test_bad_min_run_rejected(value: str):
    with pytest.raises(CriteriaError, match="must be an integer"):
        read_criteria_text(_region_target(rule="run", min_run=value))
