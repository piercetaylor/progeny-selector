"""Parent-of-origin classification of every call relative to the two parents.

Responsibility: turn allele-index calls into the six-state code defined in
``progeny_selector.constants`` (A, H, B, X, N, U). Pure function of the
genotype matrix and the two parent sample ids. This is the definition shared
with the sibling backcross project (PLAN.md, algorithm 1).

Interface:
    classify(gm: GenotypeMatrix, rp_id: str, donor_id: str) -> Classification
    Classification.states  int8 (n_markers, n_samples)
    Classification.informative  bool (n_markers,)
    Classification.uninformative_reason  str (n_markers,) ('' when informative)
    state_labels(states) -> np.ndarray of 'A'/'H'/'B'/'X'/'N'/'U'
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from progeny_selector.constants import (
    STATE_A,
    STATE_B,
    STATE_H,
    STATE_LABELS,
    STATE_N,
    STATE_U,
    STATE_X,
)
from progeny_selector.model.dataset import GenotypeMatrix


@dataclass
class Classification:
    states: np.ndarray
    informative: np.ndarray
    uninformative_reason: np.ndarray
    rp_allele: np.ndarray
    donor_allele: np.ndarray

    @property
    def n_informative(self) -> int:
        return int(self.informative.sum())


def classify(gm: GenotypeMatrix, rp_id: str, donor_id: str) -> Classification:
    """Classify each call as A/H/B/X/N/U.

    A marker is informative only when both parents are called, both are homozygous,
    and their alleles differ. All calls at an uninformative marker are coded U,
    including the parents' own calls.
    """
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


def state_labels(states: np.ndarray) -> np.ndarray:
    """Map an int8 state array to single-letter labels."""
    lookup = np.array([STATE_LABELS[i] for i in range(len(STATE_LABELS))], dtype=object)
    return lookup[states]


def rpp_contribution(states: np.ndarray) -> np.ndarray:
    """Per-call recurrent-parent contribution: A=1, H=0.5, B=0, otherwise NaN (excluded)."""
    out = np.full(states.shape, np.nan, dtype=float)
    out[states == STATE_A] = 1.0
    out[states == STATE_H] = 0.5
    out[states == STATE_B] = 0.0
    return out


def is_called_informative(states: np.ndarray) -> np.ndarray:
    """True for calls that enter RPP denominators (A, H or B)."""
    return (states == STATE_A) | (states == STATE_H) | (states == STATE_B)


__all__ = ["STATE_N", "STATE_U", "Classification", "classify", "is_called_informative", "rpp_contribution", "state_labels"]
