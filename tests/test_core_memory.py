"""Analysis-memory bound for ``classify``, ``rpp`` and ``ibs_to_sample`` (docs/adr/0025).

The three reference implementations below are copied verbatim from the pre-phase
bodies of ``core.background.rpp``, ``core.similarity.ibs_to_sample`` and
``core.classify.classify``. They are the second implementation the chunked code is
measured against: the equality tests compare two implementations, and the memory
tests state the peak the pre-phase code could not meet.
"""

from __future__ import annotations

import tracemalloc
from collections.abc import Callable

import numpy as np
import pytest

from progeny_selector.constants import STATE_A, STATE_B, STATE_H, STATE_N, STATE_U, STATE_X
from progeny_selector.core.background import rpp
from progeny_selector.core.classify import Classification, classify, is_called_informative, rpp_contribution
from progeny_selector.core.similarity import _shared_alleles, ibs_to_sample
from progeny_selector.model.dataset import GenotypeMatrix, Marker

# --- pre-phase reference implementations, copied verbatim ------------------------------------


def _rpp_reference(states: np.ndarray, weights: np.ndarray | None = None, marker_mask: np.ndarray | None = None) -> np.ndarray:
    contrib = rpp_contribution(states)
    counted = is_called_informative(states)
    if marker_mask is not None:
        counted = counted & marker_mask[:, None]
    w = np.ones(states.shape[0], dtype=float) if weights is None else weights
    w2 = np.where(counted, w[:, None], 0.0)
    denom = w2.sum(axis=0)
    numer = np.nansum(w2 * np.nan_to_num(contrib), axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = numer / denom
    out[denom == 0] = np.nan
    return out


def _ibs_reference(gm: GenotypeMatrix, ref_sample_id: str) -> np.ndarray:
    ref = gm.calls[:, gm.sample_index(ref_sample_id), :]
    called = np.asarray((gm.calls >= 0).all(axis=2)) & np.asarray((ref >= 0).all(axis=1))[:, None]
    shared = _shared_alleles(gm.calls, ref[:, None, :]) / 2.0
    denom = called.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(called, shared, 0.0).sum(axis=0) / denom
    out[denom == 0] = np.nan
    return out


def _classify_reference(gm: GenotypeMatrix, rp_id: str, donor_id: str) -> Classification:
    rp = gm.calls[:, gm.sample_index(rp_id), :]
    donor = gm.calls[:, gm.sample_index(donor_id), :]

    rp_called = (rp >= 0).all(axis=1)
    donor_called = (donor >= 0).all(axis=1)
    rp_hom = rp[:, 0] == rp[:, 1]
    donor_hom = donor[:, 0] == donor[:, 1]
    differ = rp[:, 0] != donor[:, 0]
    informative = rp_called & donor_called & rp_hom & donor_hom & differ

    reason = np.full(gm.n_markers, "", dtype=object)
    reason[~rp_called] = "recurrent parent missing"
    reason[rp_called & ~donor_called] = "donor parent missing"
    reason[rp_called & donor_called & ~rp_hom] = "recurrent parent heterozygous"
    reason[rp_called & donor_called & rp_hom & ~donor_hom] = "donor parent heterozygous"
    reason[rp_called & donor_called & rp_hom & donor_hom & ~differ] = "parents identical (monomorphic)"

    rp_allele = np.where(informative, rp[:, 0], -1).astype(np.int8)
    donor_allele = np.where(informative, donor[:, 0], -1).astype(np.int8)

    a1 = gm.calls[:, :, 0]
    a2 = gm.calls[:, :, 1]
    r = rp_allele[:, None]
    d = donor_allele[:, None]

    missing = (a1 < 0) | (a2 < 0)
    is_a = (a1 == r) & (a2 == r)
    is_b = (a1 == d) & (a2 == d)
    is_h = ((a1 == r) & (a2 == d)) | ((a1 == d) & (a2 == r))

    states = np.full((gm.n_markers, gm.n_samples), STATE_X, dtype=np.int8)
    states[is_h] = STATE_H
    states[is_b] = STATE_B
    states[is_a] = STATE_A
    states[missing] = STATE_N
    states[~informative, :] = STATE_U
    return Classification(
        states=states,
        informative=informative,
        uninformative_reason=reason,
        rp_allele=rp_allele,
        donor_allele=donor_allele,
    )


# --- inputs ------------------------------------------------------------------------------------


def _random_matrix(n_markers: int, n_samples: int, seed: int = 7) -> GenotypeMatrix:
    """Random calls in -1..2 with RP and DONOR planted as in ``tests/conftest.py::make_matrix``.

    Column order is RP, DONOR, then progeny. RP is homozygous allele 0 and DONOR homozygous
    allele 2 at every marker, so every marker is informative.
    """
    rng = np.random.default_rng(seed)
    calls = rng.integers(-1, 3, size=(n_markers, n_samples, 2), dtype=np.int8)
    calls[:, 0, :] = 0
    calls[:, 1, :] = 2
    markers = [Marker(f"m{i + 1}", "Gm01", (i + 1) * 1000, (i + 1) * 0.0025) for i in range(n_markers)]
    sample_ids = ["RP", "DONOR"] + [f"P{i}" for i in range(n_samples - 2)]
    return GenotypeMatrix(markers=markers, sample_ids=sample_ids, alleles=[["A", "G", "T"]] * n_markers, calls=calls)


def _random_states(n_markers: int, n_samples: int, seed: int = 11) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.integers(0, 6, size=(n_markers, n_samples)).astype(np.int8)


def _peak_bytes(fn: Callable[[], object]) -> int:
    """tracemalloc peak of one call, with the inputs already built and a warm-up discarded."""
    fn()
    tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        fn()
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


# --- bit-identity -------------------------------------------------------------------------------
# 700 is not a multiple of 64, so the last block of every loop is partial. 65, 129, 193 and 2049
# leave a one-column tail, which numpy reduces pairwise where a wider block accumulates rows in
# order; ``sample_blocks`` folds that tail into the block before it so weighted RPP stays exact.

BLOCK_SHAPES = [64, 65, 128, 129, 193, 700, 2049]


@pytest.mark.parametrize("n_s", BLOCK_SHAPES)
def test_rpp_bit_identical_at_block_boundaries(n_s: int) -> None:
    states = _random_states(3_000, n_s, seed=100 + n_s)
    rng = np.random.default_rng(200 + n_s)
    w = rng.uniform(0.1, 5.0, size=3_000)
    mask = rng.random(3_000) < 0.4
    assert np.array_equal(rpp(states, w), _rpp_reference(states, w), equal_nan=True)
    assert np.array_equal(rpp(states, w, mask), _rpp_reference(states, w, mask), equal_nan=True)
    assert np.array_equal(rpp(states, None), _rpp_reference(states, None), equal_nan=True)


@pytest.mark.parametrize("n_s", [64, 65, 128, 129, 193])
def test_ibs_and_classify_bit_identical_at_block_boundaries(n_s: int) -> None:
    gm = _random_matrix(600, n_s, seed=300 + n_s)
    assert np.array_equal(ibs_to_sample(gm, "RP"), _ibs_reference(gm, "RP"), equal_nan=True)
    assert np.array_equal(classify(gm, "RP", "DONOR").states, _classify_reference(gm, "RP", "DONOR").states)


def test_rpp_with_no_samples() -> None:
    states = np.zeros((100, 0), dtype=np.int8)
    w = np.random.default_rng(21).uniform(0.1, 5.0, size=100)
    out = rpp(states, w)
    assert out.shape == (0,)
    assert np.array_equal(out, _rpp_reference(states, w), equal_nan=True)


def test_rpp_with_an_infinite_weight() -> None:
    """An infinite weight makes the columns it enters NaN and leaves every other column finite."""
    states = _random_states(500, 130, seed=22)
    w = np.random.default_rng(23).uniform(0.1, 5.0, size=500)
    w[7] = np.inf
    got = rpp(states, w)
    want = _rpp_reference(states, w)
    assert np.array_equal(got, want, equal_nan=True)
    assert not np.all(np.isnan(got)), "an infinite weight must not poison every column"
    assert np.isnan(got[np.flatnonzero(is_called_informative(states[7]))[0]])


def test_rpp_bit_identical_to_reference() -> None:
    states = _random_states(3_000, 700)
    rng = np.random.default_rng(3)
    w = rng.uniform(0.1, 5.0, size=3_000)
    mask = rng.random(3_000) < 0.4
    assert np.array_equal(rpp(states, w), _rpp_reference(states, w), equal_nan=True)
    assert np.array_equal(rpp(states, w, mask), _rpp_reference(states, w, mask), equal_nan=True)
    assert np.array_equal(rpp(states, None), _rpp_reference(states, None), equal_nan=True)
    assert np.array_equal(rpp(states, None, mask), _rpp_reference(states, None, mask), equal_nan=True)


def test_rpp_edge_masks_bit_identical() -> None:
    states = _random_states(3_000, 700)
    w = np.random.default_rng(4).uniform(0.1, 5.0, size=3_000)
    none_mask = np.zeros(3_000, dtype=bool)
    zero_w = np.zeros(3_000)
    assert np.array_equal(rpp(states, w, none_mask), _rpp_reference(states, w, none_mask), equal_nan=True)
    assert np.all(np.isnan(rpp(states, w, none_mask)))
    assert np.array_equal(rpp(states, zero_w), _rpp_reference(states, zero_w), equal_nan=True)


def test_ibs_to_sample_bit_identical_to_reference() -> None:
    gm = _random_matrix(3_000, 700)
    assert np.array_equal(ibs_to_sample(gm, "RP"), _ibs_reference(gm, "RP"), equal_nan=True)
    assert np.array_equal(ibs_to_sample(gm, "DONOR"), _ibs_reference(gm, "DONOR"), equal_nan=True)


def test_classify_bit_identical_to_reference() -> None:
    gm = _random_matrix(3_000, 700)
    got = classify(gm, "RP", "DONOR")
    want = _classify_reference(gm, "RP", "DONOR")
    assert np.array_equal(got.states, want.states)
    assert np.array_equal(got.informative, want.informative)
    assert np.array_equal(got.uninformative_reason, want.uninformative_reason)
    assert np.array_equal(got.rp_allele, want.rp_allele)
    assert np.array_equal(got.donor_allele, want.donor_allele)


def test_classify_uninformative_markers_bit_identical() -> None:
    """Every block path is overwritten with U after the loop, including a block-boundary marker."""
    gm = _random_matrix(3_000, 700, seed=31)
    gm.calls[5, 0, :] = -1  # recurrent parent missing
    gm.calls[64, 0, :] = (0, 2)  # recurrent parent heterozygous
    gm.calls[129, 1, :] = (0, 2)  # donor parent heterozygous
    gm.calls[699, 1, :] = 0  # parents identical
    got = classify(gm, "RP", "DONOR")
    want = _classify_reference(gm, "RP", "DONOR")
    assert not got.informative.all()
    assert (got.states[~got.informative, :] == STATE_U).all()
    assert np.array_equal(got.states, want.states)
    assert np.array_equal(got.informative, want.informative)
    assert np.array_equal(got.uninformative_reason, want.uninformative_reason)
    assert np.array_equal(got.rp_allele, want.rp_allele)
    assert np.array_equal(got.donor_allele, want.donor_allele)


def test_ibs_all_missing_column_is_nan() -> None:
    """A sample with no call anywhere has denom 0 and reads NaN, in the first block and the last."""
    gm = _random_matrix(600, 129, seed=32)
    gm.calls[:, 3, :] = -1
    gm.calls[:, 128, :] = -1
    got = ibs_to_sample(gm, "RP")
    assert np.isnan(got[3]) and np.isnan(got[128])
    assert np.array_equal(got, _ibs_reference(gm, "RP"), equal_nan=True)


def test_fewer_samples_than_one_block() -> None:
    gm = _random_matrix(120, 40)
    states = _random_states(120, 40, seed=5)
    w = np.random.default_rng(6).uniform(0.1, 5.0, size=120)
    assert np.array_equal(rpp(states, w), _rpp_reference(states, w), equal_nan=True)
    assert np.array_equal(ibs_to_sample(gm, "RP"), _ibs_reference(gm, "RP"), equal_nan=True)
    assert np.array_equal(classify(gm, "RP", "DONOR").states, _classify_reference(gm, "RP", "DONOR").states)


# --- memory -------------------------------------------------------------------------------------
# Pre-phase peaks on the same inputs, measured against the reference implementations above:
# rpp 136,017,736 B (34.0 x states.nbytes), ibs_to_sample 76,002,696 B (9.5 x calls.nbytes),
# classify 28,066,408 B (3.5 x calls.nbytes).


def test_rpp_memory_within_four_times_states() -> None:
    states = _random_states(4_000, 1_000, seed=13)
    w = np.random.default_rng(14).uniform(0.1, 5.0, size=4_000)
    peak = _peak_bytes(lambda: rpp(states, w))
    assert peak < 4 * states.nbytes, f"rpp peak {peak / states.nbytes:.1f} x states.nbytes"


def test_ibs_memory_within_four_times_calls() -> None:
    gm = _random_matrix(4_000, 1_000, seed=15)
    peak = _peak_bytes(lambda: ibs_to_sample(gm, "RP"))
    assert peak < 4 * gm.calls.nbytes, f"ibs_to_sample peak {peak / gm.calls.nbytes:.1f} x calls.nbytes"


def test_classify_memory_within_three_times_calls() -> None:
    gm = _random_matrix(4_000, 1_000, seed=16)
    peak = _peak_bytes(lambda: classify(gm, "RP", "DONOR"))
    assert peak < 3 * gm.calls.nbytes, f"classify peak {peak / gm.calls.nbytes:.1f} x calls.nbytes"
