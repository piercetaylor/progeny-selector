"""Generate the synthetic BC2F1 fixture under tests/fixtures/synthetic_bc2f1/.

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

Usage: python scripts/make_fixture.py
"""

from __future__ import annotations

import csv
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from progeny_selector.constants import SOYBEAN_CHROM_LENGTHS_BP_WM82A4

OUT = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "synthetic_bc2f1"
SEED = 20260904
PER_CHROM = 25
CM_PER_BP = 2.5 / 1_000_000
RP_ID, DONOR_ID = "RP_Williams82", "DONOR_PI_synthetic"
TARGET = "syn_Gm06_13"
AVOID = "syn_Gm13_10"
FLANK_CM = 6.0
WEIGHTS = {"rpp_noncarrier": 0.5, "rpp_carrier": 0.2, "drag": 0.2, "recombinant": 0.1}
NUC = "ACGT"

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


def simulate(markers: list[dict]) -> tuple[list[str], dict[str, list[int]], dict[str, str], dict[str, str]]:
    """Return progeny ids, per-progeny states (0=A,1=H,2=B,4=N), family map, generation map."""
    by_chrom: dict[str, list[int]] = {}
    for i, m in enumerate(markers):
        by_chrom.setdefault(m["chrom"], []).append(i)
    n = len(markers)
    t_idx = next(i for i, m in enumerate(markers) if m["id"] == TARGET)

    def bc_gamete(hap_donor_side: list[int]) -> list[int]:
        gam = [0] * n
        for idx in by_chrom.values():
            cm = [markers[i]["cm"] for i in idx]
            g = haldane_gamete([0] * len(idx), [hap_donor_side[i] for i in idx], cm)
            for j, i in enumerate(idx):
                gam[i] = g[j]
        return gam

    def bc1f1_carrier() -> list[int]:
        while True:
            g = bc_gamete([1] * n)  # F1 gamete: recombinant of RP and donor haplotypes
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
            g = bc_gamete(parent_hap)  # BC1F1 gamete from haplotypes (0..., parent_hap)
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
    g1, g2 = bc_gamete(parent_hap), bc_gamete(parent_hap)
    selfed = [g1[i] + g2[i] for i in range(n)]  # 0=A,1=H,2=B
    selfed[t_idx] = max(selfed[t_idx], 1)
    states["BC2F1-F2-002"] = selfed
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
    states["BC2F1-F1-003"][t_idx] = 1
    return ids, states, family, generation


def expected_metrics(markers: list[dict], ids: list[str], states: dict[str, list[int]], family: dict[str, str]) -> list[dict]:
    n = len(markers)
    t_idx = next(i for i, m in enumerate(markers) if m["id"] == TARGET)
    a_idx = next(i for i, m in enumerate(markers) if m["id"] == AVOID)
    carrier_chrom = markers[t_idx]["chrom"]
    chrom_len = SOYBEAN_CHROM_LENGTHS_BP_WM82A4[carrier_chrom]
    t_cm = markers[t_idx]["cm"]
    t_bp = markers[t_idx]["pos"]
    left = [i for i in range(n) if markers[i]["chrom"] == carrier_chrom and i < t_idx and markers[i]["informative"]][::-1]
    right = [i for i in range(n) if markers[i]["chrom"] == carrier_chrom and i > t_idx and markers[i]["informative"]]
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
        qc_het = abs(het_rate - 0.25) > 0.15  # advisory flag only, never excludes
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
                "drag_total_max_cm": drag_max_cm,
                "drag_total_est_cm": drag_est_cm,
                "recomb_left": rec_left,
                "recomb_right": rec_right,
                "composite_score": score,
                "missing_rate": missing_rate,
                "het_rate": het_rate,
                "het_rate_deviates": qc_het,
                "family_donor_outlier": False,  # set below, once every family member's fraction is known
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
    for r in rows:
        r.setdefault("rank_overall", "")
        r.setdefault("rank_in_family", "")
    design = {r["sample_id"]: r for r in rows}
    assert design["BC2F1-F1-001"]["rank_overall"] == 1, "planted best individual must rank first"
    assert design["BC2F1-F1-002"]["exclusion_reason"].startswith("target:"), "planted foreground failure"
    assert "avoid:" in design["BC2F1-F2-001"]["exclusion_reason"], "planted avoid failure"
    assert "qc:possible_self_or_outcross" in design["BC2F1-F2-002"]["exclusion_reason"], "planted self contaminant"
    assert "missing_rate" in design["BC2F1-F1-003"]["exclusion_reason"], "planted high-missing individual"
    outliers = {r["sample_id"] for r in rows if r["family_donor_outlier"]}
    assert outliers == {"BC2F1-F2-002"}, f"family donor outliers {sorted(outliers)}"
    return rows


def middle_value(values: list[float]) -> float:
    """Median: the middle value, or the mean of the two middle values for an even count."""
    s = sorted(values)
    k = len(s)
    return s[k // 2] if k % 2 == 1 else (s[k // 2 - 1] + s[k // 2]) / 2


def write_outputs(
    markers: list[dict], ids: list[str], states: dict[str, list[int]], family: dict[str, str], generation: dict[str, str], rows: list[dict]
) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    samples = [RP_ID, DONOR_ID, *ids]
    with open(OUT / "genotypes.vcf", "w", newline="\n") as fh:
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
    with open(OUT / "expected_results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.10g}" if isinstance(v, float) else v) for k, v in r.items()})


def main() -> None:
    markers = build_markers()
    ids, states, family, generation = simulate(markers)
    rows = expected_metrics(markers, ids, states, family)
    write_outputs(markers, ids, states, family, generation, rows)
    n_pass = sum(1 for r in rows if r["passes_filters"])
    n_inf = sum(m["informative"] for m in markers)
    print(f"wrote fixture to {OUT}: {len(markers)} markers, {len(ids)} progeny, {n_pass} pass; informative {n_inf}")


if __name__ == "__main__":
    main()
