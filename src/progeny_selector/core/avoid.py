"""Negative selection: status at user-defined "avoid" loci or regions.

Responsibility: decide per individual whether an avoid locus is recurrent-parent
homozygous (pass), carries donor alleles (fail) or has no data (unknown).
Shares locus resolution with foreground (PLAN.md, algorithm 5).

Interface:
    avoid_status(states, resolved, allow_het=False, rule="all", min_markers=1) -> int8 per sample
"""

from __future__ import annotations

import numpy as np

from progeny_selector.core.foreground import ResolvedLocus, locus_status, marker_predicate


def avoid_status(
    states: np.ndarray,
    resolved: ResolvedLocus,
    allow_het: bool = False,
    rule: str = "all",
    min_markers: int = 1,
) -> np.ndarray:
    """Pass when every counted marker in the locus is A (or A/H when allow_het)."""
    required = "rp_or_het" if allow_het else "hom_rp"
    return locus_status(states, resolved, marker_predicate(states, required), rule, min_markers)
