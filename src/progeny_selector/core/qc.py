"""Quality control of parents and progeny before ranking.

Responsibility: per-individual missing, heterozygosity, homozygous-donor and
non-parental rates; comparison with the generation expectation; parent
heterozygosity; heuristic flags for possible selfs, outcrosses, sample swaps
and duplicates (PLAN.md, algorithm 8). Flags are advisory strings; the
ranking step decides whether flagged individuals are excluded.

Interface:
    sample_qc(gm, dataset, classification, filters) -> list[SampleQC]
    parent_qc(gm, dataset, classification, filters) -> list[str] warnings
    duplicate_pairs(gm, sample_ids, threshold=0.995) -> list[tuple[str, str, float]]
    uninformative_summary(classification) -> list[tuple[str, int]]   (reason, count), count desc then reason
    qc_table_rows(qc, dataset) -> list[dict]   one flat row per SampleQC for the Validate screen
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import numpy as np

from progeny_selector.constants import STATE_A, STATE_B, STATE_H, STATE_N, STATE_X
from progeny_selector.core.classify import Classification
from progeny_selector.core.generation import expected_fractions, parse_generation
from progeny_selector.core.score import QC_EXCLUDING_FLAGS
from progeny_selector.core.similarity import ibs_to_sample, pairwise_ibs
from progeny_selector.model.criteria import Filters
from progeny_selector.model.dataset import Dataset, GenotypeMatrix


@dataclass
class SampleQC:
    sample_id: str
    missing_rate: float
    het_rate: float
    hom_donor_rate: float
    nonparental_rate: float
    expected_het: float | None
    expected_rpp: float | None
    ibs_rp: float
    ibs_donor: float
    flags: list[str] = field(default_factory=list)


def sample_qc(
    gm: GenotypeMatrix,
    dataset: Dataset,
    classification: Classification,
    filters: Filters,
) -> list[SampleQC]:
    """QC metrics and flags for every progeny/candidate in the manifest."""
    states = classification.states
    inf = classification.informative
    n_all = gm.n_markers
    ibs_rp = ibs_to_sample(gm, dataset.recurrent_parent.sample_id)
    ibs_dp = ibs_to_sample(gm, dataset.donor_parent.sample_id)
    out: list[SampleQC] = []
    for s in dataset.progeny:
        j = gm.sample_index(s.sample_id)
        col = states[:, j]
        missing_rate = float(((gm.calls[:, j, :] < 0).any(axis=1)).sum() / n_all) if n_all else float("nan")
        n_a = int((col == STATE_A).sum())
        n_h = int((col == STATE_H).sum())
        n_b = int((col == STATE_B).sum())
        n_x = int((col == STATE_X).sum())
        n_called_inf = n_a + n_h + n_b
        n_inf_nonmissing = int((inf & (col != STATE_N)).sum())
        het_rate = n_h / n_called_inf if n_called_inf else float("nan")
        hom_donor_rate = n_b / n_called_inf if n_called_inf else float("nan")
        nonparental_rate = n_x / n_inf_nonmissing if n_inf_nonmissing else float("nan")
        gen = parse_generation(s.generation)
        exp = expected_fractions(gen) if gen else None
        qc = SampleQC(
            sample_id=s.sample_id,
            missing_rate=missing_rate,
            het_rate=het_rate,
            hom_donor_rate=hom_donor_rate,
            nonparental_rate=nonparental_rate,
            expected_het=exp.het if exp else None,
            expected_rpp=exp.rpp if exp else None,
            ibs_rp=float(ibs_rp[j]),
            ibs_donor=float(ibs_dp[j]),
        )
        if missing_rate > filters.max_missing_rate:
            qc.flags.append("high_missing")
        if gen is None and s.generation:
            qc.flags.append("generation_unparsed")
        if n_called_inf:
            if gen is not None and gen.n_filial == 1 and gen.n_backcross >= 1 and hom_donor_rate > filters.max_hom_donor_rate_bcf1:
                qc.flags.append("possible_self_or_outcross")
            if exp is not None and abs(het_rate - exp.het) > filters.het_rate_tolerance:
                qc.flags.append("het_rate_deviates")
            if nonparental_rate > filters.max_nonparental_rate:
                qc.flags.append("possible_outcross")
            if ibs_rp[j] > 0.995 and het_rate < 0.005:
                qc.flags.append("possible_rp_sample")
            if ibs_dp[j] > 0.995 and het_rate < 0.005:
                qc.flags.append("possible_donor_sample")
        out.append(qc)
    return out


def parent_qc(gm: GenotypeMatrix, dataset: Dataset, classification: Classification, filters: Filters) -> list[str]:
    """Warnings about parent calls: heterozygosity, missingness, informativeness."""
    warnings: list[str] = []
    for parent in (dataset.recurrent_parent, dataset.donor_parent):
        j = gm.sample_index(parent.sample_id)
        calls = gm.calls[:, j, :]
        called = (calls >= 0).all(axis=1)
        het = (calls[:, 0] != calls[:, 1]) & called
        het_rate = het.sum() / called.sum() if called.any() else float("nan")
        if het_rate > filters.parent_max_het_rate:
            warnings.append(
                f"{parent.role} {parent.sample_id}: heterozygous call rate {het_rate:.3f} exceeds {filters.parent_max_het_rate}"
            )
        miss = 1.0 - called.mean() if gm.n_markers else float("nan")
        if miss > 0.1:
            warnings.append(f"{parent.role} {parent.sample_id}: missing rate {miss:.3f}")
    n_inf = classification.n_informative
    if n_inf == 0:
        warnings.append("no informative markers: parents identical, heterozygous or missing everywhere")
    elif n_inf < 0.2 * gm.n_markers:
        warnings.append(f"only {n_inf} of {gm.n_markers} markers are informative between the parents")
    return warnings


def duplicate_pairs(
    gm: GenotypeMatrix, sample_ids: list[str], threshold: float = 0.995, max_samples: int = 2000
) -> list[tuple[str, str, float]]:
    """Pairs of samples with IBS above ``threshold`` (marker-subsampled); empty when too many samples."""
    if len(sample_ids) < 2 or len(sample_ids) > max_samples:
        return []
    idx = np.array([gm.sample_index(s) for s in sample_ids])
    ibs = pairwise_ibs(gm, idx)
    pairs: list[tuple[str, str, float]] = []
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            if ibs[a, b] >= threshold:
                pairs.append((sample_ids[a], sample_ids[b], float(ibs[a, b])))
    return pairs


def uninformative_summary(classification: Classification) -> list[tuple[str, int]]:
    """Count of uninformative markers per reason, most frequent first, ties by reason; informative markers excluded."""
    counts = Counter(str(r) for r in classification.uninformative_reason.tolist() if r)
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))


def qc_table_rows(qc: list[SampleQC], dataset: Dataset) -> list[dict]:
    """Flat rows for the QC table, in ``qc`` order.

    ``qc_excluded`` is true when a QC flag alone excludes the individual (any flag in
    ``QC_EXCLUDING_FLAGS``); the missing-rate and locus hard filters are reported on the Rank
    screen. NaN rates become None, as in ``AnalysisResult.rows``.
    """
    rows: list[dict] = []
    for q in qc:
        s = dataset.sample(q.sample_id)
        rows.append(
            {
                "sample_id": q.sample_id,
                "line_name": s.line_name,
                "family_id": s.family_id,
                "generation": s.generation,
                "missing_rate": _num(q.missing_rate),
                "het_rate": _num(q.het_rate),
                "expected_het": q.expected_het,
                "hom_donor_rate": _num(q.hom_donor_rate),
                "nonparental_rate": _num(q.nonparental_rate),
                "expected_rpp": q.expected_rpp,
                "ibs_rp": _num(q.ibs_rp),
                "ibs_donor": _num(q.ibs_donor),
                "flags": "|".join(q.flags),
                "qc_excluded": any(f in QC_EXCLUDING_FLAGS for f in q.flags),
            }
        )
    return rows


def _num(x: float) -> float | None:
    """NaN -> None so rows serialise cleanly (core.pipeline has the same helper; importing it here would be circular)."""
    return None if np.isnan(x) else float(x)
