"""Shared constants: parent-of-origin state codes, soybean chromosome tables, palette.

Responsibility: single source of truth for small tables used across io, core and app.
Interface: module-level constants only; no functions with side effects.
"""

from __future__ import annotations

import numpy as np

# Parent-of-origin state codes (int8). See PLAN.md "Core algorithms and metrics", 1.
STATE_A: int = 0  # homozygous recurrent-parent allele
STATE_H: int = 1  # heterozygous: one recurrent-parent allele, one donor allele
STATE_B: int = 2  # homozygous donor allele
STATE_X: int = 3  # carries an allele found in neither parent (non-parental)
STATE_N: int = 4  # missing call in the progeny
STATE_U: int = 5  # marker uninformative (parents identical, heterozygous, or missing)

STATE_LABELS: dict[int, str] = {
    STATE_A: "A",
    STATE_H: "H",
    STATE_B: "B",
    STATE_X: "X",
    STATE_N: "N",
    STATE_U: "U",
}
LABEL_TO_STATE: dict[str, int] = {v: k for k, v in STATE_LABELS.items()}

# Recurrent-parent contribution per state for RPP; NaN means "excluded from numerator and denominator".
RPP_CONTRIBUTION: np.ndarray = np.array([1.0, 0.5, 0.0, np.nan, np.nan, np.nan], dtype=float)

# Locus status codes.
STATUS_FAIL: int = 0
STATUS_PASS: int = 1
STATUS_UNKNOWN: int = 2
STATUS_LABELS: dict[int, str] = {STATUS_FAIL: "fail", STATUS_PASS: "pass", STATUS_UNKNOWN: "unknown"}

# Soybean (Glycine max) chromosomes. Lengths in bp per assembly; the criteria.yaml key ``assembly``
# picks the table (docs/adr/0015). a1, a2 and a4 are as listed by SoyBase,
# https://www.soybase.org/resources/genome_info/ (a2 and a4 fetched 2026-09-04, a1 2026-09-16).
# Checked 2026-09-19 against NCBI's assembly sequence reports (datasets v2alpha) for GCF_000004515.3
# (v1.1 = a1), .4 (v2.0 = a2) and .6 (v4.0 = a4): every a1 and a4 chromosome is exactly one greater
# than the NCBI length, which is SoyBase's end convention, not an off-by-one here. Ten a2 chromosomes
# are 300 bp to 3.6 kb longer than NCBI's, because SoyBase's Wm82.a2.v1 is not NCBI's copy of it;
# SoyBase governs, since SoySNP50K positions are published against it.
SOYBEAN_CHROMOSOMES: tuple[str, ...] = tuple(f"Gm{i:02d}" for i in range(1, 21))

# Wm82.a1.v1, SoyBase https://www.soybase.org/resources/genome_info/ (fetched 2026-09-16).
# The values are the 1-based chromosome ends as published, so Gm01 ends at 55,915,596 and Gm18 at 62,308,141.
SOYBEAN_CHROM_LENGTHS_BP_WM82A1: dict[str, int] = {
    "Gm01": 55_915_596,
    "Gm02": 51_656_714,
    "Gm03": 47_781_077,
    "Gm04": 49_243_853,
    "Gm05": 41_936_505,
    "Gm06": 50_722_822,
    "Gm07": 44_683_158,
    "Gm08": 46_995_533,
    "Gm09": 46_843_751,
    "Gm10": 50_969_636,
    "Gm11": 39_172_791,
    "Gm12": 40_113_141,
    "Gm13": 44_408_972,
    "Gm14": 49_711_205,
    "Gm15": 50_939_161,
    "Gm16": 37_397_386,
    "Gm17": 41_906_775,
    "Gm18": 62_308_141,
    "Gm19": 50_589_442,
    "Gm20": 46_773_168,
}

SOYBEAN_CHROM_LENGTHS_BP_WM82A4: dict[str, int] = {
    "Gm01": 57_932_356,
    "Gm02": 50_400_359,
    "Gm03": 46_951_867,
    "Gm04": 51_203_390,
    "Gm05": 42_274_531,
    "Gm06": 50_945_865,
    "Gm07": 44_949_257,
    "Gm08": 47_227_185,
    "Gm09": 50_572_669,
    "Gm10": 51_638_688,
    "Gm11": 39_643_746,
    "Gm12": 41_531_200,
    "Gm13": 45_225_049,
    "Gm14": 49_893_279,
    "Gm15": 53_754_296,
    "Gm16": 38_112_071,
    "Gm17": 41_740_657,
    "Gm18": 58_286_271,
    "Gm19": 51_272_881,
    "Gm20": 47_846_027,
}

SOYBEAN_CHROM_LENGTHS_BP_WM82A2: dict[str, int] = {
    "Gm01": 56_831_625,
    "Gm02": 48_577_506,
    "Gm03": 45_779_782,
    "Gm04": 52_389_147,
    "Gm05": 42_234_499,
    "Gm06": 51_416_487,
    "Gm07": 44_630_647,
    "Gm08": 47_837_941,
    "Gm09": 50_189_765,
    "Gm10": 51_566_899,
    "Gm11": 34_766_868,
    "Gm12": 40_091_315,
    "Gm13": 45_874_163,
    "Gm14": 49_042_193,
    "Gm15": 51_756_344,
    "Gm16": 37_887_015,
    "Gm17": 41_641_367,
    "Gm18": 58_018_743,
    "Gm19": 50_746_917,
    "Gm20": 47_904_182,
}


# The criteria.yaml ``assembly`` default; ``model.criteria`` re-exports it, and ``core`` reads it from here
# because ``model`` imports ``core.chrom``.
DEFAULT_ASSEMBLY: str = "Wm82.a4"

# Chromosome lengths keyed by the criteria.yaml ``assembly`` value; "none" and unknown names have no table
# and fall back to the last marker on the chromosome (docs/adr/0015).
SOYBEAN_CHROM_LENGTHS_BP: dict[str, dict[str, int]] = {
    "Wm82.a1": SOYBEAN_CHROM_LENGTHS_BP_WM82A1,
    "Wm82.a2": SOYBEAN_CHROM_LENGTHS_BP_WM82A2,
    "Wm82.a4": SOYBEAN_CHROM_LENGTHS_BP_WM82A4,
}

# Okabe-Ito colorblind-safe palette, hex values as shipped in R's grDevices::palette.colors
# ("Okabe-Ito"), https://raw.githubusercontent.com/wch/r-source/trunk/src/library/grDevices/R/colorstuff.R
PALETTE_OKABE_ITO: dict[str, str] = {
    "black": "#000000",
    "orange": "#E69F00",
    "skyblue": "#56B4E9",
    "bluishgreen": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddishpurple": "#CC79A7",
    "gray": "#999999",
}

# Colour per parent-of-origin state used by every strip and chip in the UI.
STATE_COLORS: dict[str, str] = {
    "A": PALETTE_OKABE_ITO["blue"],
    "H": PALETTE_OKABE_ITO["bluishgreen"],
    "B": PALETTE_OKABE_ITO["vermillion"],
    "X": PALETTE_OKABE_ITO["reddishpurple"],
    "N": PALETTE_OKABE_ITO["gray"],
    "U": PALETTE_OKABE_ITO["yellow"],
}
STATUS_COLORS: dict[str, str] = {
    "pass": PALETTE_OKABE_ITO["bluishgreen"],
    "fail": PALETTE_OKABE_ITO["vermillion"],
    "unknown": PALETTE_OKABE_ITO["gray"],
}

SAMPLE_ROLES: tuple[str, ...] = ("recurrent_parent", "donor_parent", "candidate", "progeny")
# Wide-CSV missing tokens per mode (contract/data-contract.md 1.1.0, "Wide CSV"). `?` is not missing.
WIDE_NUCLEOTIDE_MISSING: frozenset[str] = frozenset({"", "N", "NN", "NA", "-", "--", ".", "./.", ".|."})
WIDE_CODED_MISSING: frozenset[str] = frozenset({"", "N", "NA"})
# HapMap missing tokens (contract 1.1.0, "HapMap"): the wide-CSV nucleotide list plus TASSEL's X/XX.
HAPMAP_MISSING: frozenset[str] = WIDE_NUCLEOTIDE_MISSING | {"X", "XX"}
