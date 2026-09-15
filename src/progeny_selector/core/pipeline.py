"""End-to-end analysis: dataset + criteria -> per-individual result rows.

Responsibility: orchestrate classification, foreground, background, drag,
avoid, similarity, QC, hard filters, composite score and ranking into one
list of flat dict rows (one per progeny/candidate) plus per-run metadata.
No file I/O; the UI and CLI call ``run_analysis`` and hand the rows to
``progeny_selector.io.export``.

Interface:
    run_analysis(dataset: Dataset, criteria: Criteria) -> AnalysisResult
    AnalysisResult.rows: list[dict]   (column names listed in docs/data-formats.md, "results.csv")
    AnalysisResult.states, .classification, .warnings, .carrier_chroms, .unit
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from progeny_selector.constants import STATE_A, STATE_B, STATE_H, STATE_N, STATUS_LABELS
from progeny_selector.core.avoid import avoid_status
from progeny_selector.core.background import carrier_mask, marker_weights, n_called_informative, rpp, rpp_per_chromosome
from progeny_selector.core.chrom import chrom_length_bp
from progeny_selector.core.classify import Classification, classify
from progeny_selector.core.drag import DragResult, donor_segment
from progeny_selector.core.foreground import ResolvedLocus, foreground_status, resolve_locus
from progeny_selector.core.qc import SampleQC, parent_qc, sample_qc
from progeny_selector.core.score import composite_score, hard_filters, rank_rows
from progeny_selector.core.similarity import ibs_to_sample
from progeny_selector.model.criteria import Criteria
from progeny_selector.model.dataset import Dataset


@dataclass
class AnalysisResult:
    rows: list[dict]
    classification: Classification
    resolved_targets: dict[str, ResolvedLocus]
    resolved_avoid: dict[str, ResolvedLocus]
    drag: dict[str, DragResult]
    qc: list[SampleQC]
    carrier_chroms: set[str]
    unit: str
    warnings: list[str] = field(default_factory=list)

    def row(self, sample_id: str) -> dict:
        return next(r for r in self.rows if r["sample_id"] == sample_id)


def _pick_unit(dataset: Dataset, requested: str) -> str:
    if requested == "cm" and dataset.genotypes.has_cm():
        return "cm"
    if requested == "auto" and dataset.genotypes.has_cm():
        return "cm"
    return "bp"


def run_analysis(dataset: Dataset, criteria: Criteria) -> AnalysisResult:
    criteria.validate()
    gm = dataset.genotypes.sorted_by_position()
    dataset = Dataset(
        genotypes=gm, samples=dataset.samples, warnings=list(dataset.warnings), synthetic_sample_ids=dataset.synthetic_sample_ids
    )
    warnings = list(dataset.warnings)
    rp_id, donor_id = dataset.recurrent_parent.sample_id, dataset.donor_parent.sample_id
    cls = classify(gm, rp_id, donor_id)
    warnings += parent_qc(gm, dataset, cls, criteria.filters)

    progeny = dataset.progeny
    pidx = np.array([gm.sample_index(s.sample_id) for s in progeny], dtype=int)
    states = cls.states[:, pidx]
    sample_ids = [s.sample_id for s in progeny]
    family_ids = [s.family_id for s in progeny]

    unit = _pick_unit(dataset, criteria.background.map_unit)
    flank_unit = _pick_unit(dataset, criteria.flank_unit)
    if criteria.flank_unit == "cm" and flank_unit == "bp":
        warnings.append("flank_unit is cm but the map has no cM; windows interpreted in bp")
    chrom_lengths = {c: float(chrom_length_bp(c, None) or gm.positions("bp")[gm.chroms() == c].max()) for c in set(gm.chroms().tolist())}

    # Foreground
    resolved_t = {t.locus_id: resolve_locus(t, gm) for t in criteria.targets}
    positions_bp = gm.positions("bp")
    target_status = {
        t.locus_id: foreground_status(
            states,
            resolved_t[t.locus_id],
            t.required_state,
            t.rule,
            t.min_markers,
            positions=positions_bp,
            min_run=t.min_run,
            anchor_bp=t.anchor_bp,
            tolerate_isolated=t.tolerate_isolated,
        )
        for t in criteria.targets
    }
    carrier = {r.chrom for r in resolved_t.values()}

    # Background
    weights = None
    if criteria.background.model == "weighted":
        weights = marker_weights(
            gm, cls.informative, unit, criteria.background.max_marker_coverage, chrom_lengths if unit == "bp" else None
        )
    cmask = carrier_mask(gm, carrier)
    rpp_total = rpp(states, weights)
    rpp_carrier = rpp(states, weights, cmask)
    rpp_noncarrier = rpp(states, weights, ~cmask)
    rpp_chrom = rpp_per_chromosome(states, gm, weights)
    n_inf_called = n_called_informative(states)

    # Linkage drag and recombinants
    drag: dict[str, DragResult] = {}
    for t in criteria.targets:
        r = resolved_t[t.locus_id]
        wl = criteria.flank_window if t.flank_left is None else t.flank_left
        wr = criteria.flank_window if t.flank_right is None else t.flank_right
        clen = chrom_lengths.get(r.chrom) if flank_unit == "bp" else None
        drag[t.locus_id] = donor_segment(states, gm, r, flank_unit, clen, wl, wr)

    # Avoid
    resolved_a = {a.locus_id: resolve_locus(a, gm) for a in criteria.avoid}
    avoid_st = {a.locus_id: avoid_status(states, resolved_a[a.locus_id], a.allow_het, a.rule, a.min_markers) for a in criteria.avoid}

    # Similarity and QC
    ibs_rp_all = ibs_to_sample(gm, rp_id)[pidx]
    ibs_dp_all = ibs_to_sample(gm, donor_id)[pidx]
    qc = sample_qc(gm, dataset, cls, criteria.filters)
    missing_rate = np.array([q.missing_rate for q in qc])
    qc_flags = [q.flags for q in qc]

    # Components
    drag_total_est = np.zeros(len(progeny))
    drag_total_max = np.zeros(len(progeny))
    drag_norm = np.zeros(len(progeny))
    rec_frac = np.zeros(len(progeny))
    for t in criteria.targets:
        d = drag[t.locus_id]
        drag_total_est += d.total_est
        drag_total_max += d.total_max
        r = resolved_t[t.locus_id]
        span = chrom_lengths[r.chrom] if flank_unit == "bp" else float(np.nanmax(gm.positions("cm")[gm.chroms() == r.chrom]))
        drag_norm += np.clip(d.total_est / span, 0.0, 1.0) if span > 0 else 1.0
        rec_frac += d.recombinant_left.astype(float) + d.recombinant_right.astype(float)
    n_t = max(len(criteria.targets), 1)
    components = {
        "rpp_noncarrier": np.where(np.isnan(rpp_noncarrier), rpp_total, rpp_noncarrier),
        "rpp_carrier": rpp_carrier,
        "drag": 1.0 - drag_norm / n_t,
        "recombinant": rec_frac / (2 * n_t),
        "similarity_rp": ibs_rp_all,
        "completeness": 1.0 - missing_rate,
    }
    score, _ = composite_score(components, criteria.weights)
    passes, reasons = hard_filters(target_status, avoid_st, missing_rate, qc_flags, criteria.filters)
    rank_overall, rank_in_family = rank_rows(score, rpp_total, drag_total_est, missing_rate, sample_ids, passes, family_ids)

    # Rows
    counted = (states == STATE_A) | (states == STATE_H) | (states == STATE_B)
    denom = np.maximum(counted.sum(axis=0), 1)
    frac_a = (states == STATE_A).sum(axis=0) / denom
    frac_h = (states == STATE_H).sum(axis=0) / denom
    frac_b = (states == STATE_B).sum(axis=0) / denom
    rows: list[dict] = []
    for i, s in enumerate(progeny):
        row: dict = {
            "sample_id": s.sample_id,
            "line_name": s.line_name,
            "family_id": s.family_id,
            "generation": s.generation,
            "role": s.role,
            "rank_overall": _num(rank_overall[i]),
            "rank_in_family": _num(rank_in_family[i]),
            "passes_filters": bool(passes[i]),
            "exclusion_reason": reasons[i],
            "composite_score": _num(score[i]),
            "foreground_all_pass": all(target_status[t][i] == 1 for t in target_status),
            "avoid_all_pass": all(avoid_st[a][i] == 1 for a in avoid_st) if avoid_st else True,
            "rpp_total": _num(rpp_total[i]),
            "rpp_carrier": _num(rpp_carrier[i]),
            "rpp_noncarrier": _num(rpp_noncarrier[i]),
            "expected_rpp": qc[i].expected_rpp,
            "drag_total_est": _num(drag_total_est[i]),
            "drag_total_max": _num(drag_total_max[i]),
            "drag_unit": flank_unit,
            "ibs_rp": _num(ibs_rp_all[i]),
            "ibs_donor": _num(ibs_dp_all[i]),
            "missing_rate": _num(missing_rate[i]),
            "het_rate": qc[i].het_rate,
            "expected_het": qc[i].expected_het,
            "n_informative_called": int(n_inf_called[i]),
            "frac_a": float(frac_a[i]),
            "frac_h": float(frac_h[i]),
            "frac_b": float(frac_b[i]),
            "qc_flags": "|".join(qc[i].flags),
        }
        for t in criteria.targets:
            d = drag[t.locus_id]
            row[f"target_{t.locus_id}_status"] = STATUS_LABELS[int(target_status[t.locus_id][i])]
            row[f"drag_{t.locus_id}_left_max"] = _num(d.left_max[i])
            row[f"drag_{t.locus_id}_right_max"] = _num(d.right_max[i])
            row[f"recomb_{t.locus_id}_left"] = bool(d.recombinant_left[i])
            row[f"recomb_{t.locus_id}_right"] = bool(d.recombinant_right[i])
        for a in criteria.avoid:
            row[f"avoid_{a.locus_id}_status"] = STATUS_LABELS[int(avoid_st[a.locus_id][i])]
        for chrom, values in rpp_chrom.items():
            row[f"rpp_{chrom}"] = _num(values[i])
        rows.append(row)
    rows.sort(key=lambda r: (r["rank_overall"] if r["rank_overall"] is not None else float("inf"), r["sample_id"]))
    return AnalysisResult(
        rows=rows,
        classification=cls,
        resolved_targets=resolved_t,
        resolved_avoid=resolved_a,
        drag=drag,
        qc=qc,
        carrier_chroms=carrier,
        unit=unit,
        warnings=warnings,
    )


def _num(x: float) -> float | None:
    """NaN -> None so rows serialise cleanly."""
    x = float(x)
    return None if np.isnan(x) else x


def n_missing(states: np.ndarray) -> np.ndarray:
    return (states == STATE_N).sum(axis=0)
