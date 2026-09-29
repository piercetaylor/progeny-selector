"""Single-cross Purdy generation labels (docs/adr/0034): dose, filial guard and the QC checks it skips."""

from __future__ import annotations

import pytest

from progeny_selector.core.classify import classify
from progeny_selector.core.generation import Generation, expected_rpp, parse_generation
from progeny_selector.core.qc import sample_qc
from progeny_selector.model.criteria import Filters
from tests.conftest import make_dataset, make_matrix


@pytest.mark.parametrize("label", ["RP*3/DON", "3*RP/DON", "DON/RP*3", "DON/3*RP", "  RP*3/DON  ", "Williams82*3/PI88788"])
def test_every_placement_of_the_dose_reads_the_same(label):
    gen = parse_generation(label)
    assert gen == Generation(n_backcross=2, n_filial=1, filial_known=False)
    assert gen.label() == "BC2"


@pytest.mark.parametrize(("label", "n_backcross"), [("RP*2/DON", 1), ("RP*4/DON", 3), ("DON/7*RP", 6), ("RP*1/DON", 0)])
def test_the_dose_counts_the_initial_cross(label, n_backcross):
    # n >= 2 is BC(n-1); n = 1 is the F1 itself (n_backcross 0).
    gen = parse_generation(label)
    assert gen is not None and gen.n_backcross == n_backcross and not gen.filial_known


def test_label_has_no_filial_part():
    # BC{n}, except that a dose of 1 is the F1 itself.
    assert parse_generation("RP*1/DON").label() == "F1"
    assert parse_generation("RP*2/DON").label() == "BC1"


@pytest.mark.parametrize(
    "label",
    ["A/B", "A*2/B*3", "A*2/B/C", "A*2//B", "A**2/B", "A*/B", "*2/B", "A * 2/B", "RP*0/DON", "A*2", "A/B*", ""],
)
def test_anything_else_stays_unparsed(label):
    assert parse_generation(label) is None


def test_the_bcnfm_grammar_is_unchanged_and_filial_known():
    gen = parse_generation("BC2F1")
    assert gen == Generation(n_backcross=2, n_filial=1) and gen.filial_known
    assert gen.label() == "BC2F1"
    assert parse_generation("F2").label() == "F2"


def _qc(generation: str):
    # Ten markers: 5 A, 5 B, no H. A BCnF1 cannot carry B, so a filial-known label fires both checks.
    gm = make_matrix(["A", "A", "A", "A", "A", "B", "B", "B", "B", "B"])
    dataset = make_dataset(gm, generation=generation)
    cls = classify(gm, "RP", "DONOR")
    return sample_qc(gm, dataset, cls, Filters())[0]


def test_a_purdy_line_sets_expected_rpp_only():
    q = _qc("RP*2/DON")
    assert q.expected_rpp == pytest.approx(expected_rpp(1)) == pytest.approx(0.75)
    assert q.expected_het is None
    assert "generation_unparsed" not in q.flags


def test_a_selfed_purdy_line_never_fires_het_rate_deviates():
    """A BC1F2 written as RP*2/DON: het 0 against a BC1F1's 0.5 would fire under a filial-known label."""
    control = _qc("BC1F1")
    assert "het_rate_deviates" in control.flags and "possible_self_or_outcross" in control.flags
    purdy = _qc("RP*2/DON")
    assert "het_rate_deviates" not in purdy.flags
    assert "possible_self_or_outcross" not in purdy.flags


def test_a_plain_cross_is_flagged_unparsed():
    q = _qc("RP/DON")
    assert "generation_unparsed" in q.flags
    assert q.expected_rpp is None and q.expected_het is None
