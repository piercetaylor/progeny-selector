"""Generate the synthetic BC2F1 and BC3F1 fixtures under tests/fixtures/.

Deterministic (seed 20260904). Design:
  - 20 chromosomes Gm01..Gm20, 25 evenly spaced markers each (500 markers), positions
    inside the Wm82.a4.v1 chromosome lengths; a synthetic genetic map at 2.5 cM/Mb.
  - Recurrent parent (RP) homozygous for a random allele; donor homozygous for the
    other, except 20 monomorphic markers, 3 markers with a missing donor call and
    2 with a heterozygous RP call (25 uninformative, 475 informative).
  - Two families of 20 BC2F1 progeny. Each family descends from one simulated BC1F1
    plant (selected as a target carrier); BC2F1 gametes are Haldane recombinants of
    that plant's two haplotypes; the other gamete is RP.
  - Target locus: marker syn_Gm06_13 (donor allele required, state 'either').
    Avoid locus: marker syn_Gm13_10 (must be RP). Flank windows 6 cM.
  - Planted individuals: F1-001 the best possible passer (target H with a two-marker
    donor segment, everything else A); F1-002 fails foreground (A at target);
    F2-001 fails avoid (H at the avoid marker); F2-002 is a selfed (BC2F2) contaminant
    (B calls present); F1-003 has 30 % missing calls; 1 % random missing elsewhere.
  - expected_results.csv holds the count-model metrics computed here with the
    formulas in PLAN.md, independently of the package (only the chromosome-length
    table is imported), including the advisory family_donor_outlier flag (docs/adr/0012).
    It also holds the weighted-model RPP columns (rpp_*_weighted, cM, 10 cM cap),
    the staged rank columns (rank_*_staged, docs/adr/0007 amendment) and the advisory
    possible_duplicate flag (docs/adr/0017), each computed here by the same independent
    route, so a test compares two implementations rather than a function against itself.
  - tests/fixtures/synthetic_bc3f1/ is the next generation of the same cross (docs/adr/0018):
    the generator's own top 2 per family by rank_in_family are backcrossed once more, 10 BC3F1
    progeny each, from the true (pre-missing) BC2F1 states. Its samples.csv is written in the
    layout write_next_round_manifest produces, so the round-trip test compares the writer's
    output with this independently built file. criteria.yaml and markers.csv are not duplicated;
    the BC3F1 fixture reuses the BC2F1 ones.

Usage: python scripts/make_fixture.py
"""

from __future__ import annotations

import csv
import json
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from progeny_selector.constants import SOYBEAN_CHROM_LENGTHS_BP_WM82A4

FIXTURES = Path(__file__).resolve().parents[1] / "tests" / "fixtures"
OUT = FIXTURES / "synthetic_bc2f1"
OUT_BC3F1 = FIXTURES / "synthetic_bc3f1"
OUT_BRAPI = FIXTURES / "brapi"
SEED = 20260904
PER_CHROM = 25
CM_PER_BP = 2.5 / 1_000_000
RP_ID, DONOR_ID = "RP_Williams82", "DONOR_PI_synthetic"
TARGET = "syn_Gm06_13"
AVOID = "syn_Gm13_10"
FLANK_CM = 6.0
MAX_COVERAGE_CM = 10.0  # weighted model: each marker covers at most 10 cM, half to each side
DUPLICATE_THRESHOLD = 0.995  # docs/adr/0017
MIN_DUPLICATE_OVERLAP_FRAC = 0.5  # a pair needs calls in common at this fraction of the markers used
MAX_IBS_MARKERS = 2000  # core/similarity.py subsamples above this; the generator does not (see duplicate_flags)
WEIGHTS = {"rpp_noncarrier": 0.5, "rpp_carrier": 0.2, "drag": 0.2, "recombinant": 0.1}
TOP_PER_FAMILY = 2  # BC3F1 parents: the generator's own top 2 per family by rank_in_family
BC3F1_PER_PARENT = 10  # placeholder progeny per selected parent, matching --per-selected 10
BC3F1_MISSING_RATE = 0.01
BC3F1_PARENTS = ("BC2F1-F1-001", "BC2F1-F1-010", "BC2F1-F2-005", "BC2F1-F2-019")
NUC = "ACGT"
BRAPI_VARIANT_SET = "vs1"
BRAPI_PAGE_VARIANTS = 13  # 25 Gm06 variants over two pages, 13 then 12
BRAPI_PAGE_CALL_SETS = 5  # 8 call sets over two pages, 5 then 3
BRAPI_PROGENY_PER_FAMILY = 3
BRAPI_MAX_BYTES = 64 * 1024

rng = random.Random(SEED)


def build_markers() -> list[dict]:
    markers = []
    for c in range(1, 21):
        chrom = f"Gm{c:02d}"
        length = SOYBEAN_CHROM_LENGTHS_BP_WM82A4[chrom]
        for k in range(PER_CHROM):
            pos = round((k + 0.5) * length / PER_CHROM)
            rp = rng.choice(NUC)
            donor = rng.choice([n for n in NUC if n != rp])
            markers.append(
                {
                    "id": f"syn_{chrom}_{k + 1:02d}",
                    "chrom": chrom,
                    "pos": pos,
                    "cm": round(pos * CM_PER_BP, 3),
                    "rp": rp,
                    "donor": donor,
                    "rp_call": (rp, rp),
                    "donor_call": (donor, donor),
                    "informative": True,
                }
            )
    # uninformative markers: monomorphic, donor missing, RP heterozygous (never the target/avoid markers)
    pool = [m for m in markers if m["id"] not in (TARGET, AVOID)]
    for m in rng.sample(pool, 20):
        m["donor"] = m["rp"]
        m["donor_call"] = (m["rp"], m["rp"])
        m["informative"] = False
    rest = [m for m in pool if m["informative"]]
    for m in rng.sample(rest, 3):
        m["donor_call"] = None
        m["informative"] = False
    rest = [m for m in pool if m["informative"]]
    for m in rng.sample(rest, 2):
        m["rp_call"] = (m["rp"], m["donor"])
        m["informative"] = False
    return markers


def haldane_gamete(hap_a: list[int], hap_b: list[int], cm: list[float]) -> list[int]:
    """Recombine two haplotypes (per-marker 0 = RP allele, 1 = donor allele) along one chromosome."""
    current = rng.randrange(2)
    out = []
    for i, _ in enumerate(cm):
        if i > 0:
            d = (cm[i] - cm[i - 1]) / 100.0
            p_rec = 0.5 * (1 - math.exp(-2 * d))
            if rng.random() < p_rec:
                current ^= 1
        out.append((hap_a if current == 0 else hap_b)[i])
    return out


def index_by_chrom(markers: list[dict]) -> dict[str, list[int]]:
    """Marker indices grouped by chromosome, in file order."""
    by_chrom: dict[str, list[int]] = {}
    for i, m in enumerate(markers):
        by_chrom.setdefault(m["chrom"], []).append(i)
    return by_chrom


def bc_gamete(markers: list[dict], by_chrom: dict[str, list[int]], hap_donor_side: list[int]) -> list[int]:
    """One gamete from a plant whose haplotypes are (RP all 0, ``hap_donor_side``), chromosome by chromosome."""
    gam = [0] * len(markers)
    for idx in by_chrom.values():
        cm = [markers[i]["cm"] for i in idx]
        g = haldane_gamete([0] * len(idx), [hap_donor_side[i] for i in idx], cm)
        for j, i in enumerate(idx):
            gam[i] = g[j]
    return gam


def simulate(markers: list[dict]) -> tuple[list[str], dict[str, list[int]], dict[str, list[int]], dict[str, str], dict[str, str]]:
    """Return progeny ids, states (0=A,1=H,2=B,4=N), the true pre-missing states, family map, generation map."""
    by_chrom = index_by_chrom(markers)
    n = len(markers)
    t_idx = next(i for i, m in enumerate(markers) if m["id"] == TARGET)

    def bc1f1_carrier() -> list[int]:
        while True:
            g = bc_gamete(markers, by_chrom, [1] * n)  # F1 gamete: recombinant of RP and donor haplotypes
            if g[t_idx] == 1:
                return g  # BC1F1 haplotype pair is (RP all 0, g)

    ids: list[str] = []
    states: dict[str, list[int]] = {}
    family: dict[str, str] = {}
    generation: dict[str, str] = {}
    for fam in ("F1", "F2"):
        parent_hap = bc1f1_carrier()
        for k in range(1, 21):
            sid = f"BC2F1-{fam}-{k:03d}"
            g = bc_gamete(markers, by_chrom, parent_hap)  # BC1F1 gamete from haplotypes (0..., parent_hap)
            st = [1 if a == 1 else 0 for a in g]  # other gamete is RP -> H where donor, else A
            ids.append(sid)
            states[sid] = st
            family[sid] = fam
            generation[sid] = "BC2F1"
    # planted individuals
    best = [0] * n
    best[t_idx] = 1
    best[t_idx + 1] = 1
    for chrom, first, last in (("Gm02", 1, 10), ("Gm09", 5, 13), ("Gm17", 16, 24)):  # 28 H markers off the carrier chromosome
        for i in by_chrom[chrom][first - 1 : last]:
            best[i] = 1
    states["BC2F1-F1-001"] = best
    fail_fg = list(states["BC2F1-F1-002"])
    fail_fg[t_idx] = 0
    fail_fg[t_idx - 1] = 0
    fail_fg[t_idx + 1] = 0
    states["BC2F1-F1-002"] = fail_fg
    a_idx = next(i for i, m in enumerate(markers) if m["id"] == AVOID)
    fail_avoid = list(states["BC2F1-F2-001"])
    fail_avoid[t_idx] = 1
    fail_avoid[a_idx] = 1
    fail_avoid[a_idx + 1] = 1
    states["BC2F1-F2-001"] = fail_avoid
    # self contaminant: a BC2F2 plant (selfed progeny of BC2F1-F1-010) mislabelled as BC2F1 in family F2
    parent_hap = [1 if s == 1 else 0 for s in states["BC2F1-F1-010"]]
    g1, g2 = bc_gamete(markers, by_chrom, parent_hap), bc_gamete(markers, by_chrom, parent_hap)
    selfed = [g1[i] + g2[i] for i in range(n)]  # 0=A,1=H,2=B
    selfed[t_idx] = max(selfed[t_idx], 1)
    states["BC2F1-F2-002"] = selfed
    states["BC2F1-F1-003"][t_idx] = 1  # the high-missing individual still carries the target
    # the true states, before missing calls are punched in: the BC3F1 fixture breeds from these
    true_states = {sid: list(states[sid]) for sid in ids}
    # missing calls (never at the target or avoid marker, never in the planted best individual)
    a_idx = next(i for i, m in enumerate(markers) if m["id"] == AVOID)
    for sid in ids:
        if sid == "BC2F1-F1-001":
            continue
        rate = 0.30 if sid == "BC2F1-F1-003" else 0.01
        st = states[sid]
        for i in range(n):
            if i not in (t_idx, a_idx) and rng.random() < rate:
                st[i] = 4
    return ids, states, true_states, family, generation


def marker_weights_cm(markers: list[dict]) -> list[float]:
    """Map-interval weight in cM per informative marker; 0 for uninformative markers.

    Half the gap to each neighbouring informative marker on the same chromosome, each side capped
    at half of MAX_COVERAGE_CM; the outer side of the first and last informative marker is the cap's
    half, because no genetic length is known for a chromosome (the pipeline passes no cM lengths).
    """
    half = MAX_COVERAGE_CM / 2
    by_chrom: dict[str, list[int]] = {}
    for i, m in enumerate(markers):
        if m["informative"]:
            by_chrom.setdefault(m["chrom"], []).append(i)
    weights = [0.0] * len(markers)
    for idx in by_chrom.values():
        idx = sorted(idx, key=lambda i: markers[i]["cm"])
        for k, i in enumerate(idx):
            left = half if k == 0 else min((markers[i]["cm"] - markers[idx[k - 1]]["cm"]) / 2, half)
            right = half if k == len(idx) - 1 else min((markers[idx[k + 1]]["cm"] - markers[i]["cm"]) / 2, half)
            weights[i] = left + right
    return weights


def progeny_call(marker: dict, state: int) -> tuple[str, str] | None:
    """The diploid call written to the VCF for this state: 0=A, 1=H, 2=B, 4=missing."""
    rp, donor = marker["rp"], marker["donor"]
    return None if state == 4 else {0: (rp, rp), 1: (rp, donor), 2: (donor, donor)}[state]


def shared_alleles(a: tuple[str, str], b: tuple[str, str]) -> int:
    """Size of the multiset intersection of two diploid calls (0, 1 or 2)."""
    rest = list(b)
    n = 0
    for allele in a:
        if allele in rest:
            rest.remove(allele)
            n += 1
    return n


def duplicate_flags(markers: list[dict], ids: list[str], states: dict[str, list[int]]) -> set[str]:
    """Ids in any progeny pair whose IBS over informative markers called in both is >= the threshold.

    IBS per pair: mean over those markers of shared alleles / 2 (docs/adr/0017). The fixture has
    475 informative markers, below the 2,000-marker subsampling limit, so no subsampling applies.
    A pair whose members share calls at fewer than MIN_DUPLICATE_OVERLAP_FRAC of those markers is
    not reported, however high its IBS (docs/adr/0017, amendment 2026-09-21).
    """
    inf = [i for i, m in enumerate(markers) if m["informative"]]
    assert len(inf) <= MAX_IBS_MARKERS, (
        f"{len(inf)} informative markers exceeds MAX_IBS_MARKERS={MAX_IBS_MARKERS}: core/qc.py computes the overlap "
        "floor as ceil(frac * min(n_pool, MAX_IBS_MARKERS)) over a seed-0 subsample of that size, so this "
        "unsubsampled floor over the whole pool would no longer be the production rule. Subsample here the same "
        "way, or keep the fixture below the limit."
    )
    floor = max(1, math.ceil(MIN_DUPLICATE_OVERLAP_FRAC * len(inf)))
    flagged: set[str] = set()
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            sa, sb = states[ids[a]], states[ids[b]]
            total, n = 0.0, 0
            for i in inf:
                ca, cb = progeny_call(markers[i], sa[i]), progeny_call(markers[i], sb[i])
                if ca is None or cb is None:
                    continue
                total += shared_alleles(ca, cb) / 2
                n += 1
            if n >= floor and total / n >= DUPLICATE_THRESHOLD:
                flagged.update((ids[a], ids[b]))
    return flagged


def _desc(value: float) -> float:
    """Sort key for a descending numeric key, with NaN last."""
    return float("inf") if math.isnan(value) else -value


def _asc(value: float) -> float:
    """Sort key for an ascending numeric key, with NaN last."""
    return float("inf") if math.isnan(value) else value


def expected_metrics(
    markers: list[dict],
    ids: list[str],
    states: dict[str, list[int]],
    family: dict[str, str],
    expected_het: float = 0.25,
    design_checks: bool = True,
) -> list[dict]:
    """Expected per-individual metrics. ``expected_het`` is the generation's Mendelian heterozygosity
    (BC2F1 0.25, BC3F1 0.125); ``design_checks`` asserts the planted BC2F1 individuals came out as designed."""
    n = len(markers)
    t_idx = next(i for i, m in enumerate(markers) if m["id"] == TARGET)
    a_idx = next(i for i, m in enumerate(markers) if m["id"] == AVOID)
    carrier_chrom = markers[t_idx]["chrom"]
    chrom_len = SOYBEAN_CHROM_LENGTHS_BP_WM82A4[carrier_chrom]
    t_cm = markers[t_idx]["cm"]
    t_bp = markers[t_idx]["pos"]
    left = [i for i in range(n) if markers[i]["chrom"] == carrier_chrom and i < t_idx and markers[i]["informative"]][::-1]
    right = [i for i in range(n) if markers[i]["chrom"] == carrier_chrom and i > t_idx and markers[i]["informative"]]
    weights_cm = marker_weights_cm(markers)
    duplicates = duplicate_flags(markers, ids, states)
    rows = []
    for sid in ids:
        st = states[sid]
        inf = [i for i in range(n) if markers[i]["informative"]]
        called = [i for i in inf if st[i] != 4]
        n_a = sum(1 for i in called if st[i] == 0)
        n_h = sum(1 for i in called if st[i] == 1)
        n_b = sum(1 for i in called if st[i] == 2)
        missing_rate = sum(1 for i in range(n) if st[i] == 4) / n

        def rpp(sub: list[int], st: list[int] = st) -> float:
            den = len(sub)
            return float("nan") if den == 0 else sum(1.0 if st[i] == 0 else 0.5 if st[i] == 1 else 0.0 for i in sub) / den

        rpp_total = rpp(called)
        car = [i for i in called if markers[i]["chrom"] == carrier_chrom]
        non = [i for i in called if markers[i]["chrom"] != carrier_chrom]
        rpp_car, rpp_non = rpp(car), rpp(non)

        def wrpp(sub: list[int], st: list[int] = st) -> float:
            den = sum(weights_cm[i] for i in sub)
            if den == 0:
                return float("nan")
            return sum(weights_cm[i] * (1.0 if st[i] == 0 else 0.5 if st[i] == 1 else 0.0) for i in sub) / den

        target_state = st[t_idx]
        target_status = "unknown" if target_state == 4 else ("pass" if target_state in (1, 2) else "fail")
        avoid_state = st[a_idx]
        avoid_status = "unknown" if avoid_state == 4 else ("pass" if avoid_state == 0 else "fail")

        def walk(idx: list[int], end_bp: float, end_cm: float, st: list[int] = st) -> tuple[float, float, float, float]:
            first_a = next((i for i in idx if st[i] == 0), None)
            if first_a is None:
                max_bp, max_cm = end_bp, end_cm
                donor_idx = [i for i in idx if st[i] in (1, 2)]
            else:
                max_bp, max_cm = abs(markers[first_a]["pos"] - t_bp), abs(markers[first_a]["cm"] - t_cm)
                donor_idx = [i for i in idx[: idx.index(first_a)] if st[i] in (1, 2)]
            if donor_idx:
                last = donor_idx[-1]
                min_bp, min_cm = abs(markers[last]["pos"] - t_bp), abs(markers[last]["cm"] - t_cm)
            else:
                min_bp = min_cm = 0.0
            return min_bp, max_bp, min_cm, max_cm

        _lmin_bp, _lmax_bp, lmin_cm, lmax_cm = walk(left, t_bp, t_cm)
        span_cm = max(markers[i]["cm"] for i in range(n) if markers[i]["chrom"] == carrier_chrom)
        # in cM the chromosome end is the last mapped marker (no genetic length is known)
        _rmin_bp, _rmax_bp, rmin_cm, rmax_cm = walk(right, chrom_len - t_bp, span_cm - t_cm)
        drag_max_cm = lmax_cm + rmax_cm
        drag_est_cm = (lmin_cm + rmin_cm + drag_max_cm) / 2
        rec_left, rec_right = lmax_cm <= FLANK_CM, rmax_cm <= FLANK_CM
        comp = {
            "rpp_noncarrier": rpp_non,
            "rpp_carrier": rpp_car,
            "drag": 1 - min(1.0, drag_est_cm / span_cm),
            "recombinant": (rec_left + rec_right) / 2,
        }
        score = sum(WEIGHTS[k] * v for k, v in comp.items()) / sum(WEIGHTS.values())
        hom_donor_rate = n_b / (n_a + n_h + n_b)
        het_rate = n_h / (n_a + n_h + n_b)
        qc_self = hom_donor_rate > 0.02  # B calls are impossible in a true BCnF1 (excluding flag)
        qc_het = abs(het_rate - expected_het) > 0.15  # advisory flag only, never excludes
        reasons = []
        if target_status != "pass":
            reasons.append(f"target:T1:{target_status}")
        if avoid_status == "fail":
            reasons.append("avoid:AV1:fail")
        if missing_rate > 0.2:
            reasons.append("missing_rate>0.2")
        if qc_self:
            reasons.append("qc:possible_self_or_outcross")
        rows.append(
            {
                "sample_id": sid,
                "family_id": family[sid],
                "target_status": target_status,
                "avoid_status": avoid_status,
                "rpp_total": rpp_total,
                "rpp_carrier": rpp_car,
                "rpp_noncarrier": rpp_non,
                "rpp_total_weighted": wrpp(called),
                "rpp_carrier_weighted": wrpp(car),
                "rpp_noncarrier_weighted": wrpp(non),
                "drag_total_max_cm": drag_max_cm,
                "drag_total_est_cm": drag_est_cm,
                "recomb_left": rec_left,
                "recomb_right": rec_right,
                "composite_score": score,
                "missing_rate": missing_rate,
                "het_rate": het_rate,
                "het_rate_deviates": qc_het,
                "family_donor_outlier": False,  # set below, once every family member's fraction is known
                "possible_duplicate": sid in duplicates,
                "donor_fraction": 1 - rpp_total,  # count model: the non-recurrent share of called informative markers
                "passes_filters": not reasons,
                "exclusion_reason": ";".join(reasons),
            }
        )
    # family_donor_outlier (advisory, docs/adr/0012): donor fraction above family median + 2.5 * max(1.4826 * MAD, 0.01);
    # plants over the 0.2 missing-rate limit (high_missing) neither count nor get the flag
    by_family: dict[str, list[dict]] = {}
    for r in rows:
        if r["missing_rate"] <= 0.2:
            by_family.setdefault(r["family_id"], []).append(r)
    for members in by_family.values():
        if len(members) < 6:
            continue
        fractions = [r["donor_fraction"] for r in members]
        centre = middle_value(fractions)
        spread = max(1.4826 * middle_value([abs(f - centre) for f in fractions]), 0.01)
        for r in members:
            r["family_donor_outlier"] = r["donor_fraction"] > centre + 2.5 * spread
    for r in rows:
        del r["donor_fraction"]
    passing = [r for r in rows if r["passes_filters"]]
    passing.sort(key=lambda r: (-r["composite_score"], -r["rpp_total"], r["drag_total_est_cm"], r["missing_rate"], r["sample_id"]))
    fam_counter: dict[str, int] = {}
    for rank, r in enumerate(passing, start=1):
        r["rank_overall"] = rank
        fam_counter[r["family_id"]] = fam_counter.get(r["family_id"], 0) + 1
        r["rank_in_family"] = fam_counter[r["family_id"]]
    # staged ranking (docs/adr/0007 amendment): hard filters first, then recombinant flanks (count)
    # descending, rpp_carrier descending, rpp_noncarrier descending, drag estimate ascending,
    # missing rate ascending, sample_id ascending. sample_id breaks every tie, so ranks are 1..k.
    staged = sorted(
        passing,
        key=lambda r: (
            -(r["recomb_left"] + r["recomb_right"]),
            _desc(r["rpp_carrier"]),
            _desc(r["rpp_noncarrier"]),
            _asc(r["drag_total_est_cm"]),
            _asc(r["missing_rate"]),
            r["sample_id"],
        ),
    )
    fam_counter = {}
    for rank, r in enumerate(staged, start=1):
        r["rank_overall_staged"] = rank
        fam_counter[r["family_id"]] = fam_counter.get(r["family_id"], 0) + 1
        r["rank_in_family_staged"] = fam_counter[r["family_id"]]
    for r in rows:
        r.setdefault("rank_overall", "")
        r.setdefault("rank_in_family", "")
        r.setdefault("rank_overall_staged", "")
        r.setdefault("rank_in_family_staged", "")
    if design_checks:
        design = {r["sample_id"]: r for r in rows}
        assert design["BC2F1-F1-001"]["rank_overall"] == 1, "planted best individual must rank first"
        assert design["BC2F1-F1-002"]["exclusion_reason"].startswith("target:"), "planted foreground failure"
        assert "avoid:" in design["BC2F1-F2-001"]["exclusion_reason"], "planted avoid failure"
        assert "qc:possible_self_or_outcross" in design["BC2F1-F2-002"]["exclusion_reason"], "planted self contaminant"
        assert "missing_rate" in design["BC2F1-F1-003"]["exclusion_reason"], "planted high-missing individual"
        outliers = {r["sample_id"] for r in rows if r["family_donor_outlier"]}
        assert outliers == {"BC2F1-F2-002"}, f"family donor outliers {sorted(outliers)}"
        assert not duplicates, f"fixture progeny at or above IBS {DUPLICATE_THRESHOLD}: {sorted(duplicates)}"
    return rows


def middle_value(values: list[float]) -> float:
    """Median: the middle value, or the mean of the two middle values for an even count."""
    s = sorted(values)
    k = len(s)
    return s[k // 2] if k % 2 == 1 else (s[k // 2 - 1] + s[k // 2]) / 2


def write_vcf(path: Path, markers: list[dict], ids: list[str], states: dict[str, list[int]]) -> None:
    """The genotype matrix: both parents then the progeny, one row per marker, GT only."""
    samples = [RP_ID, DONOR_ID, *ids]
    with open(path, "w", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n##source=progeny-selector make_fixture.py (synthetic)\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for i, m in enumerate(markers):
            ref, alt = m["rp"], m["donor"]
            if alt == ref:
                alt = next(x for x in NUC if x != ref)  # monomorphic: ALT is an unused allele
            alleles = {ref: "0", alt: "1"}

            def gt(call: tuple[str, str] | None, alleles: dict[str, str] = alleles) -> str:
                return "./." if call is None else "/".join(alleles[a] for a in call)

            calls = [gt(m["rp_call"]), gt(m["donor_call"])]
            for sid in ids:
                s = states[sid][i]
                calls.append({0: "0/0", 1: "0/1", 2: "1/1", 4: "./."}[s])
            fh.write(f"{m['chrom']}\t{m['pos']}\t{m['id']}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(calls) + "\n")


def write_expected_results(path: Path, rows: list[dict]) -> None:
    """The generator's independently computed expectations, one row per progeny."""
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.10g}" if isinstance(v, float) else v) for k, v in r.items()})


def write_outputs(
    markers: list[dict], ids: list[str], states: dict[str, list[int]], family: dict[str, str], generation: dict[str, str], rows: list[dict]
) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    write_vcf(OUT / "genotypes.vcf", markers, ids, states)
    with open(OUT / "samples.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["sample_id", "line_name", "role", "generation", "family_id", "notes"])
        w.writerow([RP_ID, "Williams 82 (synthetic)", "recurrent_parent", "", "", "synthetic recurrent parent"])
        w.writerow([DONOR_ID, "PI synthetic donor", "donor_parent", "", "", "synthetic donor"])
        notes = {
            "BC2F1-F1-001": "planted: best possible (target H, two-marker donor segment, all else A)",
            "BC2F1-F1-002": "planted: fails foreground (A at target)",
            "BC2F1-F1-003": "planted: 30% missing calls",
            "BC2F1-F2-001": "planted: fails avoid (H at avoid marker)",
            "BC2F1-F2-002": "planted: selfed BC2F2 contaminant mislabelled as BC2F1 (B calls present)",
        }
        for sid in ids:
            w.writerow([sid, sid.replace("BC2F1-", "L"), "progeny", generation[sid], family[sid], notes.get(sid, "")])
    with open(OUT / "markers.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["marker_id", "chrom", "pos_bp", "cm"])
        for m in markers:
            w.writerow([m["id"], m["chrom"], m["pos"], m["cm"]])
    (OUT / "criteria.yaml").write_text(
        f"""name: synthetic BC2F1 fixture
targets:
  - locus_id: T1
    marker_id: {TARGET}
    required_state: either
avoid:
  - locus_id: AV1
    marker_id: {AVOID}
flank_window: {FLANK_CM}
flank_unit: cm
background:
  model: count
  map_unit: cm
weights:
  rpp_noncarrier: {WEIGHTS["rpp_noncarrier"]}
  rpp_carrier: {WEIGHTS["rpp_carrier"]}
  drag: {WEIGHTS["drag"]}
  recombinant: {WEIGHTS["recombinant"]}
filters:
  max_missing_rate: 0.2
  unknown_target_is: fail
  unknown_avoid_is: pass
  exclude_qc_flagged: true
""",
        newline="\n",
    )
    write_expected_results(OUT / "expected_results.csv", rows)


def line_name_bc2f1(sid: str) -> str:
    """The line_name samples.csv carries for a BC2F1 individual; write_next_round_manifest copies it to the progeny."""
    return sid.replace("BC2F1-", "L")


def select_bc3f1_parents(rows: list[dict]) -> list[str]:
    """The generator's own top TOP_PER_FAMILY per family by rank_in_family, in select_top_n order (family, rank)."""
    chosen = [r for r in rows if r["passes_filters"] and r["rank_in_family"] != "" and r["rank_in_family"] <= TOP_PER_FAMILY]
    chosen.sort(key=lambda r: (r["family_id"], r["rank_in_family"]))
    parents = [r["sample_id"] for r in chosen]
    assert tuple(parents) == BC3F1_PARENTS, f"BC3F1 parents {parents}"
    return parents


def simulate_bc3f1(
    markers: list[dict], parents: list[str], true_states: dict[str, list[int]]
) -> tuple[list[str], dict[str, list[int]], dict[str, str]]:
    """One more backcross to the RP from each selected BC2F1 plant's true (pre-missing) states.

    The parent's haplotypes are (RP all 0, 1 where its true state is H); one gamete is a Haldane
    recombinant of those, the other is RP, so the progeny is H where the gamete carries the donor
    allele and A otherwise. Missing calls are punched in at BC3F1_MISSING_RATE, never at the target
    or avoid marker. Returns ids in samples.csv order, states and the family map (family = parent id).
    """
    by_chrom = index_by_chrom(markers)
    n = len(markers)
    t_idx = next(i for i, m in enumerate(markers) if m["id"] == TARGET)
    a_idx = next(i for i, m in enumerate(markers) if m["id"] == AVOID)
    ids: list[str] = []
    states: dict[str, list[int]] = {}
    family: dict[str, str] = {}
    for parent in parents:
        hap = [1 if state == 1 else 0 for state in true_states[parent]]
        for k in range(1, BC3F1_PER_PARENT + 1):
            sid = f"{parent}-BC3F1-{k:03d}"
            gamete = bc_gamete(markers, by_chrom, hap)
            st = [1 if a == 1 else 0 for a in gamete]
            for i in range(n):
                if i not in (t_idx, a_idx) and rng.random() < BC3F1_MISSING_RATE:
                    st[i] = 4
            ids.append(sid)
            states[sid] = st
            family[sid] = parent
    return ids, states, family


BC3F1_README = f"""# Synthetic BC3F1 fixture

Generated by `python scripts/make_fixture.py` (seed 20260904); do not edit by hand. The next
generation of `../synthetic_bc2f1`: the {TOP_PER_FAMILY} best individuals per family by
`rank_in_family` ({{parents}}) are each backcrossed once more to the recurrent parent, giving
{BC3F1_PER_PARENT} BC3F1 progeny apiece, bred from the parents' true (pre-missing) states. The same
500 markers and the same two parents; `criteria.yaml` and `markers.csv` are not duplicated, the tests
read the BC2F1 copies. {BC3F1_MISSING_RATE:.0%} of calls are missing, never at the target or avoid marker.

`samples.csv` is exactly what `write_next_round_manifest(..., n_per_selected={BC3F1_PER_PARENT})` writes
from the BC2F1 selection, built here independently of that writer, so `tests/test_round_trip.py` compares
two implementations of the manifest format (docs/adr/0018). It follows the writer in copying each selected
parent's `line_name` to all {BC3F1_PER_PARENT} of its progeny, so a line_name is shared within a family;
ids and line names are placeholders to be edited after planting. The parent rows' `generation` and
`family_id` are empty cells, never `NA`, because this file is the input contract's samples.csv.

`expected_results.csv` holds statuses, RPP, drag bounds (cM), recombinant flags, composite scores,
exclusion reasons, ranks and advisory QC flags computed by the generator's independent implementation,
with the BC3F1 expectations (expected RPP 0.9375, expected heterozygosity 0.125).

## The `possible_duplicate` flags are an artefact of one planted parent

Eight of these forty individuals carry the advisory `possible_duplicate` flag (docs/adr/0017), and all
eight are in the family of `BC2F1-F1-001`. **Do not read that as typical of a BC3F1 generation.**
`BC2F1-F1-001` is the deliberately extreme individual planted in the BC2F1 fixture as the best possible
passer: 28 heterozygous markers in four blocks, RPP 0.9705 and heterozygosity 0.059, which is
BC4/BC5-equivalent genome recovery wearing a BC2F1 label. Its progeny inherit half of what little donor
genome it has, so they differ from one another at a handful of markers and reach IBS 0.995 honestly.

The three families descending from realistic parents show nothing of the kind. Maximum within-family IBS,
by parent: `BC2F1-F1-001` 0.99892 (five pairs flagged), `BC2F1-F1-010` (RPP 0.783) 0.93629,
`BC2F1-F2-005` (RPP 0.902) 0.98191, `BC2F1-F2-019` (RPP 0.851) 0.96269 — not one flag between the three.
So the fixture shows 80 % of one deliberately extreme family flagged and 0 % of the three normal ones.
Every pair overlaps at more than 460 of the 475 informative markers, far above the 238-marker overlap
floor, so the threshold alone produces this and the floor plays no part. docs/adr/0017 records the same.
"""


def write_bc3f1_outputs(
    markers: list[dict], parents: list[str], ids: list[str], states: dict[str, list[int]], family: dict[str, str], rows: list[dict]
) -> None:
    OUT_BC3F1.mkdir(parents=True, exist_ok=True)
    write_vcf(OUT_BC3F1 / "genotypes.vcf", markers, ids, states)
    with open(OUT_BC3F1 / "samples.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["sample_id", "line_name", "role", "generation", "family_id", "notes"])
        w.writerow([RP_ID, "Williams 82 (synthetic)", "recurrent_parent", "", "", "synthetic recurrent parent"])
        w.writerow([DONOR_ID, "PI synthetic donor", "donor_parent", "", "", "synthetic donor"])
        for sid in ids:
            parent = family[sid]
            w.writerow([sid, line_name_bc2f1(parent), "progeny", "BC3F1", parent, f"derived from {parent}"])
    write_expected_results(OUT_BC3F1 / "expected_results.csv", rows)
    (OUT_BC3F1 / "README.md").write_text(BC3F1_README.format(parents=", ".join(parents)), encoding="utf-8", newline="\n")


def vcf_alleles(marker: dict) -> tuple[str, str]:
    """REF and ALT exactly as ``write_vcf`` writes them, so the BrAPI pages carry the same symbols."""
    ref, alt = marker["rp"], marker["donor"]
    if alt == ref:
        alt = next(x for x in NUC if x != ref)
    return ref, alt


def vcf_token(marker: dict, sample_id: str, states: dict[str, list[int]], index: int) -> str:
    """The GT the VCF holds for one cell, as a slash-separated diploid token."""
    ref, alt = vcf_alleles(marker)
    codes = {ref: "0", alt: "1"}
    if sample_id in (RP_ID, DONOR_ID):
        call = marker["rp_call"] if sample_id == RP_ID else marker["donor_call"]
        return "./." if call is None else "/".join(codes[a] for a in call)
    return {0: "0/0", 1: "0/1", 2: "1/1", 4: "./."}[states[sample_id][index]]


def brapi_token(token: str, v: int, c: int, phased_het_used: list[bool]) -> str:
    """The same call spelled as one of the forms a BrAPI server may send.

    Page 0,0 collapses homozygotes to a single index and phases the first heterozygote it holds;
    page 1,1 writes a missing call as the bare unknown string. The allele indices never change,
    so the loaded matrix equals the VCF's.
    """
    if token == "./.":
        return "." if (v == 1 and c == 1) else "./."
    left, right = token.split("/")
    if left == right:
        return left if (v == 0 and c == 0) else token
    if v == 0 and c == 0 and not phased_het_used[0]:
        phased_het_used[0] = True
        return token.replace("/", "|")
    return token


def brapi_call_set_ids(ids: list[str], family: dict[str, str]) -> list[str]:
    """Both parents and the first BRAPI_PROGENY_PER_FAMILY progeny of each family, in manifest order."""
    chosen = [RP_ID, DONOR_ID]
    counts: dict[str, int] = {}
    for sid in ids:
        fam = family[sid]
        counts[fam] = counts.get(fam, 0) + 1
        if counts[fam] <= BRAPI_PROGENY_PER_FAMILY:
            chosen.append(sid)
    return chosen


def _write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


def _pages(items: list, size: int) -> list[list]:
    return [items[i : i + size] for i in range(0, len(items), size)]


BRAPI_README = """# BrAPI fixture (generated)

`scripts/make_fixture.py` writes these files; do not edit them by hand. They are recorded
BrAPI v2.1 pages for one variant set, `{variant_set}`, holding exactly the calls the BC2F1
fixture's `genotypes.vcf` holds for {n_variants} variants on {chrom} and {n_call_sets} call sets
(both parents and the first {per_family} progeny of each family).

Paging is deliberately small so every path is exercised with no network. `/callsets` returns
{n_call_sets} call sets over two pages of {page_call_sets}, `/variants` returns {n_variants}
variants over two pages of {page_variants}, and `/allelematrix` therefore has four pages,
`allelematrix.v0.c0.json` through `allelematrix.v1.c1.json`, read variant page outermost.

GT tokens are spelled differently per page while denoting the same alleles, which is what makes
`tests/test_brapi.py` a test of the parser and not of the fixture. Page 0,0 collapses homozygotes
to a single index (`0`, `1`) and phases one heterozygote (`0|1`); the other pages expand
homozygotes (`0/0`, `1/1`). A missing call is `.` on page 1,1 and `./.` elsewhere.

`variants-nopos.p0.json` repeats the same {n_variants} variants in one page with `referenceName`,
`start`, `end` and `referenceBases` null and `alternateBases` empty. It is the D1 fallback case:
positions then come from markers.csv, and a load without markers.csv fails naming the first marker.

`samples.csv` declares the roles of the {n_call_sets} sample ids the server's `callSetName` values
produce. `criteria.yaml` is the BC2F1 criteria with the avoid locus dropped, because that locus is
on Gm13 and this variant set holds {chrom} alone.
"""


def write_brapi_fixture(
    markers: list[dict], ids: list[str], states: dict[str, list[int]], family: dict[str, str], generation: dict[str, str]
) -> tuple[int, int]:
    """Recorded BrAPI v2.1 pages for one variant set carrying the same calls as the BC2F1 VCF.

    Returns the variant and call-set counts for the print line. Nothing here draws on ``rng``, so
    the BC3F1 simulation that follows in ``main`` is unaffected.
    """
    OUT_BRAPI.mkdir(parents=True, exist_ok=True)
    chrom = next(m["chrom"] for m in markers if m["id"] == TARGET)
    chosen = [(i, m) for i, m in enumerate(markers) if m["chrom"] == chrom]
    sample_ids = brapi_call_set_ids(ids, family)
    variant_db_ids = [f"var{k}" for k in range(len(chosen))]
    call_set_db_ids = [f"cs{k}" for k in range(len(sample_ids))]
    variant_pages = _pages(list(range(len(chosen))), BRAPI_PAGE_VARIANTS)
    call_set_pages = _pages(list(range(len(sample_ids))), BRAPI_PAGE_CALL_SETS)

    for p, page in enumerate(call_set_pages):
        _write_json(
            OUT_BRAPI / f"callsets.p{p}.json",
            {
                "metadata": {
                    "pagination": {
                        "currentPage": p,
                        "pageSize": BRAPI_PAGE_CALL_SETS,
                        "totalCount": len(sample_ids),
                        "totalPages": len(call_set_pages),
                    }
                },
                "result": {
                    "data": [
                        {
                            "callSetDbId": call_set_db_ids[k],
                            "callSetName": sample_ids[k],
                            "sampleDbId": f"smp{k}",
                            "studyDbId": "study1",
                            "variantSetDbIds": [BRAPI_VARIANT_SET],
                        }
                        for k in page
                    ]
                },
            },
        )

    def variant_record(k: int, with_position: bool) -> dict:
        marker = chosen[k][1]
        ref, alt = vcf_alleles(marker)
        if not with_position:
            return {
                "alternateBases": [],
                "end": None,
                "referenceBases": None,
                "referenceName": None,
                "start": None,
                "variantDbId": variant_db_ids[k],
                "variantNames": [marker["id"]],
                "variantSetDbId": BRAPI_VARIANT_SET,
            }
        return {
            "alternateBases": [alt],
            "end": marker["pos"],
            "referenceBases": ref,
            "referenceName": marker["chrom"],
            # BrAPI start is 0-based with end exclusive; VCF POS is 1-based (D1).
            "start": marker["pos"] - 1,
            "variantDbId": variant_db_ids[k],
            "variantNames": [marker["id"]],
            "variantSetDbId": BRAPI_VARIANT_SET,
        }

    for p, page in enumerate(variant_pages):
        _write_json(
            OUT_BRAPI / f"variants.p{p}.json",
            {
                "metadata": {
                    "pagination": {
                        "currentPage": p,
                        "nextPageToken": "",
                        "pageSize": BRAPI_PAGE_VARIANTS,
                        "totalCount": len(chosen),
                        "totalPages": len(variant_pages),
                    }
                },
                "result": {"data": [variant_record(k, True) for k in page]},
            },
        )
    _write_json(
        OUT_BRAPI / "variants-nopos.p0.json",
        {
            "metadata": {
                "pagination": {
                    "currentPage": 0,
                    "nextPageToken": "",
                    "pageSize": len(chosen),
                    "totalCount": len(chosen),
                    "totalPages": 1,
                }
            },
            "result": {"data": [variant_record(k, False) for k in range(len(chosen))]},
        },
    )

    phased_het_used = [False]
    for v, vpage in enumerate(variant_pages):
        for c, cpage in enumerate(call_set_pages):
            matrix = [
                [brapi_token(vcf_token(chosen[k][1], sample_ids[j], states, chosen[k][0]), v, c, phased_het_used) for j in cpage]
                for k in vpage
            ]
            _write_json(
                OUT_BRAPI / f"allelematrix.v{v}.c{c}.json",
                {
                    "result": {
                        "callSetDbIds": [call_set_db_ids[j] for j in cpage],
                        "dataMatrices": [
                            {
                                "dataMatrix": matrix,
                                "dataMatrixAbbreviation": "GT",
                                "dataMatrixName": "Genotype",
                                "dataType": "string",
                            }
                        ],
                        "expandHomozygotes": False,
                        "pagination": [
                            {
                                "dimension": "VARIANTS",
                                "page": v,
                                "pageSize": BRAPI_PAGE_VARIANTS,
                                "totalCount": len(chosen),
                                "totalPages": len(variant_pages),
                            },
                            {
                                "dimension": "CALLSETS",
                                "page": c,
                                "pageSize": BRAPI_PAGE_CALL_SETS,
                                "totalCount": len(sample_ids),
                                "totalPages": len(call_set_pages),
                            },
                        ],
                        "sepPhased": "|",
                        "sepUnphased": "/",
                        "unknownString": ".",
                        "variantDbIds": [variant_db_ids[k] for k in vpage],
                        "variantSetDbIds": [BRAPI_VARIANT_SET],
                    }
                },
            )

    with open(OUT_BRAPI / "samples.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["sample_id", "line_name", "role", "generation", "family_id", "notes"])
        w.writerow([RP_ID, "Williams 82 (synthetic)", "recurrent_parent", "", "", "synthetic recurrent parent"])
        w.writerow([DONOR_ID, "PI synthetic donor", "donor_parent", "", "", "synthetic donor"])
        for sid in sample_ids[2:]:
            w.writerow([sid, line_name_bc2f1(sid), "progeny", generation[sid], family[sid], ""])

    # The BC2F1 criteria minus the avoid locus: syn_Gm13_10 is not in this variant set.
    (OUT_BRAPI / "criteria.yaml").write_text(
        f"""name: synthetic BrAPI fixture ({chrom} only)
targets:
  - locus_id: T1
    marker_id: {TARGET}
    required_state: either
flank_window: {FLANK_CM}
flank_unit: cm
background:
  model: count
  map_unit: cm
weights:
  rpp_noncarrier: {WEIGHTS["rpp_noncarrier"]}
  rpp_carrier: {WEIGHTS["rpp_carrier"]}
  drag: {WEIGHTS["drag"]}
  recombinant: {WEIGHTS["recombinant"]}
filters:
  max_missing_rate: 0.2
  unknown_target_is: fail
  unknown_avoid_is: pass
  exclude_qc_flagged: true
""",
        newline="\n",
    )
    (OUT_BRAPI / "README.md").write_text(
        BRAPI_README.format(
            variant_set=BRAPI_VARIANT_SET,
            n_variants=len(chosen),
            chrom=chrom,
            n_call_sets=len(sample_ids),
            per_family=BRAPI_PROGENY_PER_FAMILY,
            page_call_sets=BRAPI_PAGE_CALL_SETS,
            page_variants=BRAPI_PAGE_VARIANTS,
        ),
        encoding="utf-8",
        newline="\n",
    )
    total = sum(f.stat().st_size for f in OUT_BRAPI.iterdir() if f.is_file())
    if total >= BRAPI_MAX_BYTES:
        raise AssertionError(f"{OUT_BRAPI} is {total} bytes, over the {BRAPI_MAX_BYTES}-byte cap")
    return len(chosen), len(sample_ids)


def main() -> None:
    markers = build_markers()
    ids, states, true_states, family, generation = simulate(markers)
    rows = expected_metrics(markers, ids, states, family)
    write_outputs(markers, ids, states, family, generation, rows)
    n_pass = sum(1 for r in rows if r["passes_filters"])
    n_inf = sum(m["informative"] for m in markers)
    print(f"wrote fixture to {OUT}: {len(markers)} markers, {len(ids)} progeny, {n_pass} pass; informative {n_inf}")
    n_variants, n_call_sets = write_brapi_fixture(markers, ids, states, family, generation)
    print(f"wrote fixture to {OUT_BRAPI}: {n_variants} variants x {n_call_sets} call sets")
    parents = select_bc3f1_parents(rows)
    ids3, states3, family3 = simulate_bc3f1(markers, parents, true_states)
    rows3 = expected_metrics(markers, ids3, states3, family3, expected_het=0.125, design_checks=False)
    write_bc3f1_outputs(markers, parents, ids3, states3, family3, rows3)
    n_pass3 = sum(1 for r in rows3 if r["passes_filters"])
    print(f"wrote fixture to {OUT_BC3F1}: {len(parents)} selected parents, {len(ids3)} BC3F1 progeny, {n_pass3} pass")


if __name__ == "__main__":
    main()
