"""Generation labels and Mendelian expectations.

Responsibility: parse generation strings such as "BC2F1", "BC3F2", "F2", "BC1"
into (number of backcrosses, filial generation) and derive the expected
genome fractions under Mendelian segregation with no selection, no
segregation distortion, unlinked loci and no genotyping error.

Interface:
    parse_generation(label: str | None) -> Generation | None
    expected_fractions(gen: Generation) -> ExpectedFractions
    expected_rpp(n_backcross: int) -> float
"""

from __future__ import annotations

import re
from dataclasses import dataclass

_PATTERN = re.compile(r"^\s*(?:BC(?P<bc>\d+))?\s*(?:[FS](?P<f>\d+))?\s*$", re.IGNORECASE)


@dataclass(frozen=True)
class Generation:
    """n_backcross: number of backcrosses to the recurrent parent after the F1.

    n_filial: filial generation counted so that F1 == 1 (no selfing), F2 == 2 (one selfing).
    "BC2S1" is read as one selfing after BC2, i.e. n_filial == 2.
    """

    n_backcross: int
    n_filial: int

    @property
    def n_selfings(self) -> int:
        return max(self.n_filial - 1, 0)

    def label(self) -> str:
        return f"BC{self.n_backcross}F{self.n_filial}" if self.n_backcross else f"F{self.n_filial}"


@dataclass(frozen=True)
class ExpectedFractions:
    """Expected genome fractions among informative markers."""

    rpp: float
    hom_rp: float
    het: float
    hom_donor: float


def parse_generation(label: str | None) -> Generation | None:
    """Parse "BC2F1", "bc2f1", "BC2", "F2", "BC1S1"; return None when not parseable."""
    if label is None:
        return None
    match = _PATTERN.match(str(label))
    if not match or (match.group("bc") is None and match.group("f") is None):
        return None
    n_bc = int(match.group("bc")) if match.group("bc") else 0
    if match.group("f") is None:
        n_f = 1
    else:
        n_f = int(match.group("f"))
        # "BC2S1" convention: S counts selfings, F counts filial generations (F1 = no selfing).
        if str(label).upper().replace(" ", "").find("S") != -1 and "F" not in str(label).upper():
            n_f += 1
    if n_f < 1:
        return None
    return Generation(n_backcross=n_bc, n_filial=n_f)


def expected_rpp(n_backcross: int) -> float:
    """Expected recurrent-parent genome proportion after n backcrosses (unlinked loci, no selection).

    1 - (1/2)^(n+1): BC1 = 0.75, BC2 = 0.875, BC3 = 0.9375.
    """
    return 1.0 - 0.5 ** (n_backcross + 1)


def expected_fractions(gen: Generation) -> ExpectedFractions:
    """Expected A/H/B fractions among informative markers for a BCnFk individual.

    After n backcrosses every locus is H with probability (1/2)^n, otherwise A.
    Each selfing halves the heterozygous fraction and splits the remainder equally into A and B.
    RPP = A + H/2 stays 1 - (1/2)^(n+1) regardless of selfing.
    """
    het_bc = 0.5**gen.n_backcross
    het = het_bc * 0.5**gen.n_selfings
    hom_donor = (het_bc - het) / 2.0
    hom_rp = 1.0 - het - hom_donor
    return ExpectedFractions(rpp=hom_rp + het / 2.0, hom_rp=hom_rp, het=het, hom_donor=hom_donor)
