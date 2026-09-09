"""Parent-of-origin classification and generation expectations on hand-built cases."""

from __future__ import annotations

import numpy as np

from progeny_selector.constants import STATE_A, STATE_B, STATE_H, STATE_N, STATE_U, STATE_X
from progeny_selector.core.chrom import chrom_sort_key, normalize_chrom
from progeny_selector.core.classify import classify, state_labels
from progeny_selector.core.generation import expected_fractions, expected_rpp, parse_generation
from progeny_selector.model.dataset import GenotypeMatrix, Marker
from tests.conftest import make_matrix


def test_six_states():
    gm = make_matrix(["A", "H", "B", "N", "X"])
    cls = classify(gm, "RP", "DONOR")
    assert list(cls.states[:, 2]) == [STATE_A, STATE_H, STATE_B, STATE_N, STATE_X]
    assert cls.informative.all()
    assert list(state_labels(cls.states[:, 2])) == ["A", "H", "B", "N", "X"]


def test_uninformative_markers_reasons():
    markers = [Marker(f"m{i}", "Gm01", i * 1000) for i in range(1, 5)]
    calls = np.zeros((4, 3, 2), dtype=np.int8)
    calls[0, 0], calls[0, 1], calls[0, 2] = (0, 0), (0, 0), (0, 0)  # monomorphic
    calls[1, 0], calls[1, 1], calls[1, 2] = (0, 1), (1, 1), (0, 1)  # RP het
    calls[2, 0], calls[2, 1], calls[2, 2] = (0, 0), (-1, -1), (1, 1)  # donor missing
    calls[3, 0], calls[3, 1], calls[3, 2] = (0, 0), (0, 1), (0, 0)  # donor het
    gm = GenotypeMatrix(markers=markers, sample_ids=["RP", "DONOR", "P1"], alleles=[["A", "T"]] * 4, calls=calls)
    cls = classify(gm, "RP", "DONOR")
    assert not cls.informative.any()
    assert (cls.states == STATE_U).all()
    assert list(cls.uninformative_reason) == [
        "parents identical (monomorphic)",
        "recurrent parent heterozygous",
        "donor parent missing",
        "donor parent heterozygous",
    ]


def test_generation_parsing_and_expectations():
    assert parse_generation("BC2F1").n_backcross == 2
    assert parse_generation("bc3f2").n_filial == 2
    assert parse_generation("F2").n_backcross == 0 and parse_generation("F2").n_filial == 2
    assert parse_generation("BC1").n_filial == 1
    assert parse_generation("BC2S1").n_filial == 2
    assert parse_generation("RIL") is None and parse_generation("") is None
    assert expected_rpp(1) == 0.75 and expected_rpp(2) == 0.875 and expected_rpp(3) == 0.9375
    e = expected_fractions(parse_generation("BC2F1"))
    assert (e.het, e.hom_donor, e.rpp) == (0.25, 0.0, 0.875)
    e2 = expected_fractions(parse_generation("BC2F2"))
    assert e2.het == 0.125 and e2.hom_donor == 0.0625 and e2.rpp == 0.875


def test_chromosome_normalisation():
    for name in ("Gm06", "gm6", "chr6", "Chr06", "6", "06", "chromosome6"):
        assert normalize_chrom(name) == "Gm06"
    assert normalize_chrom("scaffold_12") == "scaffold_12"
    assert chrom_sort_key("Gm20") < chrom_sort_key("scaffold_12")
    assert chrom_sort_key("Gm02") < chrom_sort_key("Gm10")
