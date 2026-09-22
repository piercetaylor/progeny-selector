"""Per-chromosome state segments and locus ticks for one individual's genome strip.

Responsibility: collapse one individual's per-marker parent-of-origin states into
runs along each chromosome, expressed as fractions of the chromosome length, and
place target and avoid loci on the same scale. A run boundary sits midway between
the two markers whose states differ; the first run starts at 0 and the last ends
at 1. Markers are ordered by position within each chromosome (stable, so ties
keep input order). Missing (N) and uninformative (U) calls keep their own runs,
so the strip never fills a gap with a neighbouring state. Markers at one
position with different states give zero-width runs, which are kept so every
call stays visible. Pixel geometry lives in
``progeny_selector.app.present``.

Interface:
    Segment, ChromStrip, LocusTick (frozen dataclasses)
    chromosome_strips(states: int8 (n_markers,), gm, assembly, scheme) -> list[ChromStrip]
    locus_ticks(targets, avoid, strips) -> list[LocusTick]
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from progeny_selector.constants import DEFAULT_ASSEMBLY
from progeny_selector.core.chrom import SOYBEAN, CompiledScheme, chrom_length_bp, chrom_sort_key
from progeny_selector.core.foreground import ResolvedLocus
from progeny_selector.model.dataset import GenotypeMatrix


@dataclass(frozen=True)
class Segment:
    start: float  # fraction of the chromosome length, 0 <= start <= end <= 1; start == end for co-located markers
    end: float
    state: int


@dataclass(frozen=True)
class ChromStrip:
    chrom: str
    length_bp: float
    n_markers: int
    segments: tuple[Segment, ...]


@dataclass(frozen=True)
class LocusTick:
    chrom: str
    locus_id: str
    kind: str  # 'target' | 'avoid'
    start: float
    end: float


def chromosome_strips(
    states: np.ndarray, gm: GenotypeMatrix, assembly: str = DEFAULT_ASSEMBLY, scheme: CompiledScheme = SOYBEAN
) -> list[ChromStrip]:
    """One strip per chromosome, in ``chrom_sort_key`` order under ``scheme``; ``assembly`` names the chromosome-length table."""
    chroms = gm.chroms()
    pos = gm.positions("bp")
    strips: list[ChromStrip] = []
    for chrom in sorted(set(chroms.tolist()), key=lambda c: chrom_sort_key(c, scheme)):
        idx = np.where(chroms == chrom)[0]
        idx = idx[np.argsort(pos[idx], kind="stable")]
        p = pos[idx]
        s = states[idx]
        # A marker beyond the assembly length (or an unknown chromosome) stretches the strip to that marker.
        length = max(float(chrom_length_bp(chrom, None, assembly, scheme) or 0), float(p.max()))
        if length <= 0:
            length = 1.0
        segments: list[Segment] = []
        start = 0.0
        for k in range(1, len(idx)):
            if s[k] == s[k - 1]:
                continue
            boundary = (p[k - 1] + p[k]) / 2.0 / length
            segments.append(Segment(start, boundary, int(s[k - 1])))
            start = boundary
        segments.append(Segment(start, 1.0, int(s[-1])))
        strips.append(ChromStrip(chrom=chrom, length_bp=length, n_markers=len(idx), segments=tuple(segments)))
    return strips


def locus_ticks(targets: dict[str, ResolvedLocus], avoid: dict[str, ResolvedLocus], strips: list[ChromStrip]) -> list[LocusTick]:
    """Ticks spanning start_bp..end_bp of each locus on its strip; loci on chromosomes without a strip are skipped."""
    lengths = {s.chrom: s.length_bp for s in strips}
    ticks: list[LocusTick] = []
    for kind, loci in (("target", targets), ("avoid", avoid)):
        for r in loci.values():
            length = lengths.get(r.chrom)
            if length is None:
                continue
            ticks.append(LocusTick(r.chrom, r.locus_id, kind, r.start_bp / length, r.end_bp / length))
    return ticks
