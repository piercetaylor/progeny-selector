"""Duplicate detection in the pipeline: the advisory possible_duplicate flag (docs/adr/0017)."""

from __future__ import annotations

import numpy as np
import pytest

from progeny_selector.core.pipeline import run_analysis
from progeny_selector.core.similarity import MAX_IBS_MARKERS, pairwise_ibs
from progeny_selector.model.criteria import Criteria, TargetSpec
from progeny_selector.model.dataset import Dataset, GenotypeMatrix, Marker, Sample

# 30 markers on Gm01, every one informative: RP homozygous 'A', donor homozygous 'T'.
# P1: H at the first eight markers, A elsewhere (het rate 8/30, near the BC2F1 expectation 0.25).
# P2: P1 with one call missing.  P3: P1 with ten A markers turned H, so it shares far less.
N_MARKERS = 30
P1_STATES = ["H"] * 8 + ["A"] * 22
P2_STATES = [*P1_STATES[:-1], "N"]
P3_STATES = [*P1_STATES[:10], *(["H"] * 10), *P1_STATES[20:]]
CALLS = {"A": (0, 0), "H": (0, 2), "B": (2, 2), "N": (-1, -1)}


def build(progeny: dict[str, list[str]], n_informative: int | None = None, n_monomorphic: int = 0) -> Dataset:
    """Dataset from per-progeny state columns, optionally preceded by monomorphic (uninformative) markers.

    A monomorphic marker gives RP, donor and every progeny the same homozygous call, so
    ``classify`` leaves it uninformative while IBS over all markers still counts it.
    """
    n_inf = n_informative if n_informative is not None else len(next(iter(progeny.values())))
    n = n_monomorphic + n_inf
    markers = [Marker(f"m{i + 1}", "Gm01", (i + 1) * 1_000_000, (i + 1) * 2.5) for i in range(n)]
    sample_ids = ["RP", "DONOR", *progeny]
    calls = np.zeros((n, len(sample_ids), 2), dtype=np.int8)
    calls[:n_monomorphic, :, :] = 0  # everybody homozygous for the same allele
    calls[n_monomorphic:, 0] = CALLS["A"]
    calls[n_monomorphic:, 1] = CALLS["B"]
    for j, states in enumerate(progeny.values(), start=2):
        for i, state in enumerate(states):
            calls[n_monomorphic + i, j] = CALLS[state]
    gm = GenotypeMatrix(markers=markers, sample_ids=sample_ids, alleles=[["A", "G", "T"] for _ in markers], calls=calls)
    samples = [Sample("RP", "RP", "recurrent_parent"), Sample("DONOR", "DONOR", "donor_parent")]
    samples += [Sample(s, s, "progeny", "BC2F1", "F1") for s in progeny]
    return Dataset(genotypes=gm, samples=samples)


def _criteria(marker_id: str = "m1") -> Criteria:
    return Criteria(targets=[TargetSpec("T1", marker_id=marker_id)])


def test_duplicate_pair_flagged_advisory():
    dataset = build({"P1": P1_STATES, "P2": P2_STATES, "P3": P3_STATES})
    result = run_analysis(dataset, _criteria())
    assert len(result.duplicates) == 1
    a, b, ibs = result.duplicates[0]
    assert (a, b) == ("P1", "P2")
    assert ibs >= 0.995
    rows = {r["sample_id"]: r for r in result.rows}
    assert "possible_duplicate" in rows["P1"]["qc_flags"]
    assert "possible_duplicate" in rows["P2"]["qc_flags"]
    assert "possible_duplicate" not in rows["P3"]["qc_flags"]
    # Advisory only: the flag never excludes, and the target passes for all three.
    assert all(rows[s]["passes_filters"] for s in ("P1", "P2", "P3"))
    assert all(rows[s]["exclusion_reason"] == "" for s in ("P1", "P2", "P3"))


def test_uninformative_markers_do_not_make_a_duplicate():
    """IBS is restricted to the informative markers, so a shared monomorphic bulk cannot flag a pair.

    200 monomorphic markers plus 10 informative, P1 and P2 differing at one informative marker:
    over all 210 markers IBS is 209.5 / 210 = 0.9976 and the pair would be flagged; over the 10
    informative markers it is 9.5 / 10 = 0.95 and must not be (docs/adr/0017, decision 5).
    """
    p1 = ["H", "H", "H", "A", "A", "A", "A", "A", "A", "A"]
    p2 = [*p1[:3], "H", *p1[4:]]  # one informative marker A -> H
    dataset = build({"P1": p1, "P2": p2}, n_informative=10, n_monomorphic=200)
    gm = dataset.genotypes
    all_markers = pairwise_ibs(gm, np.array([2, 3]))
    informative_only = pairwise_ibs(gm, np.array([2, 3]), marker_idx=np.arange(200, 210))
    assert all_markers[0, 1] == pytest.approx(209.5 / 210) and all_markers[0, 1] >= 0.995
    assert informative_only[0, 1] == pytest.approx(0.95)
    result = run_analysis(dataset, _criteria("m201"))
    assert result.duplicates == []
    for row in result.rows:
        assert "possible_duplicate" not in row["qc_flags"], row["sample_id"]


def test_high_missing_member_is_still_flagged():
    """IBS uses the markers called in both, so an excluded high-missing sample is still flagged."""
    p2 = [*P1_STATES[:8], *(["N"] * 10), *P1_STATES[18:]]  # 10 of 30 calls missing -> high_missing
    result = run_analysis(build({"P1": P1_STATES, "P2": p2}), _criteria())
    assert [(a, b) for a, b, _ in result.duplicates] == [("P1", "P2")]
    rows = {r["sample_id"]: r for r in result.rows}
    assert "possible_duplicate" in rows["P1"]["qc_flags"] and rows["P1"]["passes_filters"] is True
    assert "high_missing" in rows["P2"]["qc_flags"] and "possible_duplicate" in rows["P2"]["qc_flags"]
    assert rows["P2"]["passes_filters"] is False and "missing_rate" in rows["P2"]["exclusion_reason"]


def test_skipped_above_two_thousand_individuals():
    """The pairwise scan is quadratic, so above 2,000 individuals it is skipped with a warning."""
    states = ["H", *(["A"] * 9)]
    result = run_analysis(build({f"P{i:04d}": list(states) for i in range(2001)}), _criteria())
    assert result.duplicates == []
    assert "duplicate detection skipped above 2000 individuals" in result.warnings
    assert not any("possible_duplicate" in r["qc_flags"] for r in result.rows)


def test_sparse_sample_is_not_a_duplicate_of_everyone():
    """A near-empty sample agrees with whatever it overlaps, so a minimum overlap is required.

    P2 is called at 2 of 30 markers. Its IBS is 1.0 against both P1 and P3, which differ from each
    other at 16 markers and are plainly different plants; the overlap floor (half the 30 markers
    used) rejects both pairs (docs/adr/0017, amendment 2026-09-21).
    """
    p3 = [*(["H"] * 2), *(["A"] * 8), *(["H"] * 10), *(["A"] * 10)]
    p2 = [*(["H"] * 2), *(["N"] * 28)]
    dataset = build({"P1": P1_STATES, "P2": p2, "P3": p3})
    ibs, counts = pairwise_ibs(dataset.genotypes, np.array([2, 3, 4]), return_counts=True)
    assert ibs[0, 1] == pytest.approx(1.0) and ibs[1, 2] == pytest.approx(1.0)
    assert counts[0, 1] == 2 and counts[0, 2] == 30
    assert sum(1 for a, b in zip(P1_STATES, p3, strict=True) if a != b) == 16
    result = run_analysis(dataset, _criteria())
    assert result.duplicates == []
    for row in result.rows:
        assert "possible_duplicate" not in row["qc_flags"], row["sample_id"]


def test_overlap_just_above_the_floor_is_still_reported():
    """The floor is ceil(0.5 * 30) = 15: a pair overlapping on 15 markers is still a duplicate."""
    p2 = [*P1_STATES[:15], *(["N"] * 15)]  # 15 calls, all agreeing with P1
    dataset = build({"P1": P1_STATES, "P2": p2, "P3": P3_STATES})
    _, counts = pairwise_ibs(dataset.genotypes, np.array([2, 3]), return_counts=True)
    assert counts[0, 1] == 15
    result = run_analysis(dataset, _criteria())
    assert [(a, b) for a, b, _ in result.duplicates] == [("P1", "P2")]
    assert result.duplicates[0][2] == pytest.approx(1.0)
    # One call fewer and the same pair drops out.
    thinner = build({"P1": P1_STATES, "P2": [*P1_STATES[:14], *(["N"] * 16)], "P3": P3_STATES})
    assert run_analysis(thinner, _criteria()).duplicates == []


def test_boolean_marker_mask_is_rejected():
    """A boolean mask would be cast to 0/1 and index marker 1 repeatedly, giving IBS 1.0 for every pair."""
    gm = build({"P1": P1_STATES, "P2": P2_STATES}).genotypes
    mask = np.zeros(N_MARKERS, dtype=bool)
    mask[:10] = True
    with pytest.raises(TypeError, match="boolean mask"):
        pairwise_ibs(gm, np.array([2, 3]), marker_idx=mask)


def test_repeated_marker_indices_are_deduplicated():
    """A marker named twice must not weigh twice in the mean.

    P1 vs P3 over markers 0 (H/H), 9 (A/A) and 10 (A/H, so 0.5): 2.5 / 3 = 0.8333. Counting
    marker 0 a second time would give 3.5 / 4 = 0.875, which is why the indices are deduplicated.
    """
    gm = build({"P1": P1_STATES, "P2": P2_STATES, "P3": P3_STATES}).genotypes
    p1_p3 = np.array([2, 4])
    repeated = pairwise_ibs(gm, p1_p3, marker_idx=np.array([9, 0, 10, 0]))
    once = pairwise_ibs(gm, p1_p3, marker_idx=np.array([0, 9, 10]))
    assert once[0, 1] == pytest.approx(2.5 / 3)
    assert repeated[0, 1] == pytest.approx(2.5 / 3)  # not 3.5 / 4 = 0.875
    assert repeated == pytest.approx(once, nan_ok=True)


def test_marker_index_order_does_not_change_the_subsample():
    """Above the subsampling limit the seed draws from the sorted pool, so caller order cannot matter."""
    n = MAX_IBS_MARKERS + 500
    p1 = ["H" if i % 4 == 0 else "A" for i in range(n)]
    p2 = ["H" if i % 3 == 0 else "A" for i in range(n)]
    gm = build({"P1": p1, "P2": p2}).genotypes
    pool = np.arange(n)
    forward = pairwise_ibs(gm, np.array([2, 3]), marker_idx=pool)
    reversed_ = pairwise_ibs(gm, np.array([2, 3]), marker_idx=pool[::-1])
    assert n > MAX_IBS_MARKERS  # subsampling is in play, so the draw itself is being compared
    assert 0.0 < forward[0, 1] < 1.0  # a degenerate all-identical pair would make this test null
    assert reversed_ == pytest.approx(forward, nan_ok=True)
