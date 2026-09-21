"""Identity-by-state similarity to each parent and between individuals.

Responsibility: IBS proportions computed on raw allele calls. ``ibs_to_sample``
uses every marker, not only the informative ones; ``pairwise_ibs`` uses every
marker unless the caller restricts it to a marker subset (the pipeline restricts
duplicate detection to the informative markers, docs/adr/0017). Definition: per
marker, the size of the multiset
intersection of the two diploid genotypes divided by 2 (identical homozygotes
1, one shared allele 0.5, none 0), averaged over markers called in both. For
biallelic markers this equals SNPRelate's ``1 - |g1 - g2| / 2`` averaged over
SNPs [web] https://rdrr.io/bioc/SNPRelate/man/snpgdsIBS.html (PLAN.md, algorithm 6).

Interface:
    ibs_to_sample(gm, ref_sample_id) -> float (n_samples,)
    pairwise_ibs(gm, sample_idx, max_markers=2000, seed=0, marker_idx=None, return_counts=False)
        -> float (k, k), or (float (k, k), int (k, k)) with return_counts
"""

from __future__ import annotations

from typing import Literal, overload

import numpy as np

from progeny_selector.model.dataset import GenotypeMatrix

MAX_IBS_MARKERS = 2000  # pairwise IBS subsamples down to this many markers


def _shared_alleles(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Multiset intersection size (0, 1, 2) of diploid calls; a, b have trailing dimension 2."""
    a1, a2 = a[..., 0], a[..., 1]
    b1, b2 = b[..., 0], b[..., 1]
    same_pair = ((a1 == b1) & (a2 == b2)) | ((a1 == b2) & (a2 == b1))
    any_shared = (a1 == b1) | (a1 == b2) | (a2 == b1) | (a2 == b2)
    return np.where(same_pair, 2, np.where(any_shared, 1, 0))


def ibs_to_sample(gm: GenotypeMatrix, ref_sample_id: str) -> np.ndarray:
    """Mean shared-allele proportion between every sample and ``ref_sample_id``; NaN when no overlap."""
    ref = gm.calls[:, gm.sample_index(ref_sample_id), :]
    called = np.asarray((gm.calls >= 0).all(axis=2)) & np.asarray((ref >= 0).all(axis=1))[:, None]
    shared = _shared_alleles(gm.calls, ref[:, None, :]) / 2.0
    denom = called.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(called, shared, 0.0).sum(axis=0) / denom
    out[denom == 0] = np.nan
    return out


@overload
def pairwise_ibs(
    gm: GenotypeMatrix,
    sample_idx: np.ndarray,
    max_markers: int = ...,
    seed: int = ...,
    marker_idx: np.ndarray | None = ...,
    return_counts: Literal[False] = ...,
) -> np.ndarray: ...


@overload
def pairwise_ibs(
    gm: GenotypeMatrix,
    sample_idx: np.ndarray,
    max_markers: int = ...,
    seed: int = ...,
    marker_idx: np.ndarray | None = ...,
    *,
    return_counts: Literal[True],
) -> tuple[np.ndarray, np.ndarray]: ...


def pairwise_ibs(
    gm: GenotypeMatrix,
    sample_idx: np.ndarray,
    max_markers: int = MAX_IBS_MARKERS,
    seed: int = 0,
    marker_idx: np.ndarray | None = None,
    return_counts: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Pairwise IBS among the given samples on a random subset of at most ``max_markers`` markers.

    With ``return_counts`` the return is ``(ibs, counts)``, where ``counts`` is the integer matrix
    of markers called in both members of each pair — the denominator of the IBS average, which a
    caller needs to tell a confident match from an agreement over a handful of calls. The default
    returns the IBS matrix alone.

    ``marker_idx`` restricts the pool the subset is drawn from; the default is every marker. It is
    positional indices, never a boolean mask: a mask would be cast to 0/1 and index marker 1 over
    and over, giving IBS 1.0 for every pair, so a boolean array is rejected. The indices are
    deduplicated and sorted before subsampling, so the subset a seed draws does not depend on the
    order the caller passed (docs/adr/0017).
    """
    rng = np.random.default_rng(seed)
    if marker_idx is None:
        m_idx = np.arange(gm.n_markers)
    else:
        given = np.asarray(marker_idx)
        if given.dtype == bool:
            raise TypeError("pairwise_ibs: marker_idx must be positional indices, not a boolean mask; use np.flatnonzero")
        m_idx = np.unique(given.astype(int))
    if m_idx.size > max_markers:
        m_idx = np.sort(rng.choice(m_idx, size=max_markers, replace=False))
    calls = gm.calls[m_idx][:, sample_idx, :]
    k = len(sample_idx)
    out = np.full((k, k), np.nan)
    counts = np.zeros((k, k), dtype=int)
    called = np.asarray((calls >= 0).all(axis=2))
    for i in range(k):
        shared = _shared_alleles(calls, calls[:, i : i + 1, :]) / 2.0
        both = called & called[:, i : i + 1]
        denom = both.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            row = np.where(both, shared, 0.0).sum(axis=0) / denom
        row[denom == 0] = np.nan
        out[i, :] = row
        counts[i, :] = denom
    return (out, counts) if return_counts else out
