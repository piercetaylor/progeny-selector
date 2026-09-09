"""Selection lists and next-generation projections.

Responsibility: pick the top N ranked individuals overall or per family from a
results table, and project the expected genome composition of their offspring
under one more backcross or one selfing (PLAN.md, algorithm 7). Assumes
Mendelian segregation, unlinked loci, no selection within the next
generation, no distortion, no genotyping error.

Interface:
    select_top_n(rows, n, per_family=True) -> list[dict]
    project_next_generation(rows, step) -> Projection
"""

from __future__ import annotations

from dataclasses import dataclass
from statistics import fmean


@dataclass(frozen=True)
class Projection:
    step: str  # 'backcross' or 'self'
    n_selected: int
    mean_rpp_selected: float
    expected_rpp_next: float
    expected_het_next: float
    expected_target_carrier_fraction: float
    note: str


def select_top_n(rows: list[dict], n: int, per_family: bool = True) -> list[dict]:
    """Rows with rank_in_family <= n (per family) or rank_overall <= n; excluded rows never qualify."""
    key = "rank_in_family" if per_family else "rank_overall"
    chosen = [r for r in rows if r.get("passes_filters") and r.get(key) is not None and r[key] == r[key] and r[key] <= n]
    return sorted(chosen, key=lambda r: (r.get("family_id") or "", r[key]))


def project_next_generation(rows: list[dict], step: str = "backcross") -> Projection:
    """Expected RPP, heterozygosity and target-carrier fraction one generation on.

    Computed from the mean A/H/B fractions of the selected individuals:
      backcross: A' = A + H/2, H' = H/2 + B, B' = 0   (so RPP' = (1 + RPP) / 2)
      self:      A' = A + H/4, H' = H/2,     B' = B + H/4   (RPP' = RPP)
    Target carrier fraction for a locus heterozygous in the parent: backcross 0.5 (all
    carriers heterozygous); self 0.75 (0.25 homozygous donor, 0.5 heterozygous).
    """
    if not rows:
        raise ValueError("no selected rows to project")
    a = fmean(r["frac_a"] for r in rows)
    h = fmean(r["frac_h"] for r in rows)
    b = fmean(r["frac_b"] for r in rows)
    if step == "backcross":
        a2, h2, b2 = a + h / 2, h / 2 + b, 0.0
        carrier = 0.5
        note = "next generation is BC(n+1)F1: donor content halves; every locus with donor is heterozygous"
    elif step == "self":
        a2, h2, b2 = a + h / 4, h / 2, b + h / 4
        carrier = 0.75
        note = "next generation is BCnF(k+1): heterozygosity halves; RPP unchanged in expectation"
    else:
        raise ValueError("step must be 'backcross' or 'self'")
    rpp_now = a + h / 2
    return Projection(
        step=step,
        n_selected=len(rows),
        mean_rpp_selected=rpp_now,
        expected_rpp_next=a2 + h2 / 2,
        expected_het_next=h2,
        expected_target_carrier_fraction=carrier,
        note=note + f"; assumes Mendelian segregation and unlinked loci (b'={b2:.3f})",
    )
