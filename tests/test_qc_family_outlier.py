"""QC flag family_donor_outlier (docs/adr/0012): hand-built cases for the within-family robust rule."""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
import pytest

from progeny_selector.core.pipeline import run_analysis
from progeny_selector.core.qc import FAMILY_OUTLIER_MIN_SCALE, FAMILY_OUTLIER_Z, SampleQC, _family_donor_outliers, qc_table_rows
from progeny_selector.model.criteria import Criteria, Filters, TargetSpec
from progeny_selector.model.dataset import Dataset, GenotypeMatrix, Marker, Sample

FLAG = "family_donor_outlier"


def criteria(exclude_qc_flagged: bool = True) -> Criteria:
    """One marker target at m1 (donor allele, either state); the target only makes criteria valid."""
    return Criteria(targets=[TargetSpec("T", marker_id="m1")], filters=Filters(exclude_qc_flagged=exclude_qc_flagged))


def matrix(progeny: dict[str, str]) -> GenotypeMatrix:
    """Nucleotide matrix with RP, DONOR and one column per progeny; each progeny is a string of states over
    the same markers ('A','H','B','N'). RP allele 'A', donor allele 'T'; every marker is informative."""
    n = len(next(iter(progeny.values())))
    ids = ["RP", "DONOR", *progeny]
    calls = np.zeros((n, len(ids), 2), dtype=np.int8)
    calls[:, 1] = (1, 1)
    for j, states in enumerate(progeny.values(), start=2):
        for i, s in enumerate(states):
            calls[i, j] = {"A": (0, 0), "H": (0, 1), "B": (1, 1), "N": (-1, -1)}[s]
    markers = [Marker(f"m{i + 1}", "Gm01", (i + 1) * 1_000_000, (i + 1) * 2.5) for i in range(n)]
    return GenotypeMatrix(markers=markers, sample_ids=ids, alleles=[["A", "T"] for _ in range(n)], calls=calls)


def dataset_for(gm: GenotypeMatrix, family: dict[str, str | None], generation: str | None = None) -> Dataset:
    samples = [Sample("RP", "RP", "recurrent_parent"), Sample("DONOR", "DONOR", "donor_parent")]
    samples += [Sample(s, s, "progeny", generation, family[s]) for s in gm.sample_ids[2:]]
    return Dataset(genotypes=gm, samples=samples)


def qc_from_fractions(groups: dict[str | None, list[float]]) -> tuple[list[SampleQC], Dataset]:
    """SampleQC rows whose donor fraction is the given value (all of it as hom_donor_rate, het_rate 0)."""
    qc: list[SampleQC] = []
    family: dict[str, str | None] = {}
    for fam, fractions in groups.items():
        for k, f in enumerate(fractions):
            sid = f"{fam}-{k:02d}"
            het = float("nan") if math.isnan(f) else 0.0
            qc.append(SampleQC(sid, 0.0, het, f, 0.0, None, None, 0.0, 0.0))
            family[sid] = fam
    ids = [q.sample_id for q in qc]
    gm = GenotypeMatrix(
        markers=[Marker("m1", "Gm01", 1_000_000)],
        sample_ids=["RP", "DONOR", *ids],
        alleles=[["A", "T"]],
        calls=np.zeros((1, len(ids) + 2, 2), dtype=np.int8),
    )
    return qc, dataset_for(gm, family)


def test_one_high_individual_in_family_of_eight():
    qc, ds = qc_from_fractions({"F": [0.03] * 7 + [0.14]})
    assert _family_donor_outliers(qc, ds) == {"F-07"}


def test_family_below_minimum_size_flags_nobody():
    qc, ds = qc_from_fractions({"F": [0.03] * 4 + [0.14]})
    assert _family_donor_outliers(qc, ds) == set()


def test_uniform_family_uses_scale_floor():
    qc, ds = qc_from_fractions({"F": [0.02] * 10 + [0.03]})
    assert _family_donor_outliers(qc, ds) == set()  # 0.03 < 0.02 + 2.5 * 0.01
    qc, ds = qc_from_fractions({"F": [0.02] * 10 + [0.05]})
    assert _family_donor_outliers(qc, ds) == {"F-10"}


def test_cut_is_strict():
    # median 0.02, MAD 0, floor scale 0.01: cut 0.02 + 2.5 * 0.01, and a value equal to it is not flagged.
    cut = 0.02 + FAMILY_OUTLIER_Z * FAMILY_OUTLIER_MIN_SCALE
    qc, ds = qc_from_fractions({"F": [0.02] * 10 + [cut]})
    assert _family_donor_outliers(qc, ds) == set()
    qc, ds = qc_from_fractions({"F": [0.02] * 10 + [0.045]})
    assert _family_donor_outliers(qc, ds) == set()


def test_progeny_and_candidates_in_a_family_are_grouped_together():
    qc, ds = qc_from_fractions({"F": [0.03] * 5 + [0.14]})
    samples = [s if s.sample_id not in ("F-00", "F-01", "F-02") else replace(s, role="candidate") for s in ds.samples]
    ds = Dataset(genotypes=ds.genotypes, samples=samples)
    assert _family_donor_outliers(qc, ds) == {"F-05"}


def test_even_length_median_is_mean_of_middle_values():
    # median (0.02 + 0.04) / 2 = 0.03; deviations 0.01 x5 and 0.07 -> MAD 0.01; scale 0.014826; cut 0.067065.
    # A lower-middle median (0.02, MAD 0, floor scale) would give cut 0.045 and flag 0.066; the mean of the middle values does not.
    qc, ds = qc_from_fractions({"F": [0.02, 0.02, 0.02, 0.04, 0.04, 0.066]})
    assert _family_donor_outliers(qc, ds) == set()
    qc, ds = qc_from_fractions({"F": [0.02, 0.02, 0.02, 0.04, 0.04, 0.10]})
    assert _family_donor_outliers(qc, ds) == {"F-05"}


def test_low_individual_is_never_flagged():
    qc, ds = qc_from_fractions({"F": [0.10] * 7 + [0.0]})
    assert _family_donor_outliers(qc, ds) == set()


def test_families_are_assessed_separately():
    qc, ds = qc_from_fractions({"F1": [0.03] * 7 + [0.14], "F2": [0.14] * 7 + [0.15]})
    assert _family_donor_outliers(qc, ds) == {"F1-07"}


def test_none_and_empty_family_form_one_group():
    qc, ds = qc_from_fractions({None: [0.03] * 3, "": [0.03] * 2 + [0.14]})
    assert _family_donor_outliers(qc, ds) == {"-02"}


def test_nan_fractions_neither_count_nor_flag():
    nan = float("nan")
    qc, ds = qc_from_fractions({"F": [0.03] * 4 + [0.14] + [nan] * 3})
    assert _family_donor_outliers(qc, ds) == set()  # only 5 assessable individuals
    qc, ds = qc_from_fractions({"F": [0.03] * 7 + [0.14] + [nan] * 3})
    assert _family_donor_outliers(qc, ds) == {"F-07"}


def test_het_counts_half_through_sample_qc():
    progeny = {f"P{k}": "H" + "A" * 19 for k in range(6)}
    progeny["ALLH"] = "H" * 20
    progeny["NOCALL"] = "N" * 20
    gm = matrix(progeny)
    ds = dataset_for(gm, dict.fromkeys(progeny, "F"))
    result = run_analysis(ds, criteria())
    by_id = {q.sample_id: q for q in result.qc}
    q = by_id["ALLH"]
    assert q.hom_donor_rate + q.het_rate / 2 == pytest.approx(0.5)
    assert FLAG in q.flags
    assert all(FLAG not in by_id[s].flags for s in progeny if s != "ALLH")
    assert math.isnan(by_id["NOCALL"].het_rate)


def test_high_missing_individual_neither_counts_nor_is_flagged():
    # "H" + 19 N: one call, donor fraction 0.5, missing 0.95 > max_missing_rate 0.2 -> high_missing.
    progeny = {f"P{k}": "H" + "A" * 19 for k in range(7)}
    progeny["FEW"] = "H" + "N" * 19
    ds = dataset_for(matrix(progeny), dict.fromkeys(progeny, "F"))
    by_id = {q.sample_id: q for q in run_analysis(ds, criteria()).qc}
    assert "high_missing" in by_id["FEW"].flags
    assert by_id["FEW"].hom_donor_rate + by_id["FEW"].het_rate / 2 == pytest.approx(0.5)
    assert all(FLAG not in q.flags for q in by_id.values())
    # Not counted: four normal (0.025) + HIGH (0.2) is five assessable individuals, below the minimum, so nobody
    # is flagged; counting FEW (one A call, 0.0) would make six, median 0.025, MAD 0, cut 0.05, and flag HIGH.
    progeny = {f"P{k}": "H" + "A" * 19 for k in range(4)}
    progeny["HIGH"] = "H" * 8 + "A" * 12
    progeny["FEW"] = "A" + "N" * 19
    ds = dataset_for(matrix(progeny), dict.fromkeys(progeny, "F"))
    assert all(FLAG not in q.flags for q in run_analysis(ds, criteria()).qc)
    progeny["P4"] = "H" + "A" * 19
    ds = dataset_for(matrix(progeny), dict.fromkeys(progeny, "F"))
    assert [q.sample_id for q in run_analysis(ds, criteria()).qc if FLAG in q.flags] == ["HIGH"]


def test_flag_does_not_exclude():
    progeny = {f"P{k}": "H" + "A" * 19 for k in range(6)}
    progeny["OUT"] = "HHHHHHHHAAAAAAAAAAAA"
    gm = matrix(progeny)
    ds = dataset_for(gm, dict.fromkeys(progeny, "F"))
    result = run_analysis(ds, criteria(exclude_qc_flagged=True))
    out = next(q for q in result.qc if q.sample_id == "OUT")
    assert out.flags == [FLAG]
    assert result.row("OUT")["passes_filters"] is True
    rows = {r["sample_id"]: r for r in qc_table_rows(result.qc, ds)}
    assert rows["OUT"]["flags"] == FLAG
    assert rows["OUT"]["qc_excluded"] is False


def test_flag_is_appended_last():
    progeny = {f"P{k}": "H" + "A" * 19 for k in range(6)}
    progeny["SELF"] = "BBBBHHAAAAAAAAAAAAAA"
    gm = matrix(progeny)
    ds = dataset_for(gm, dict.fromkeys(progeny, "F"), generation="BC2F1")
    first = run_analysis(ds, criteria())
    again = run_analysis(ds, criteria())
    flags = next(q for q in first.qc if q.sample_id == "SELF").flags
    assert "possible_self_or_outcross" in flags
    assert flags[-1] == FLAG
    assert flags.count(FLAG) == 1
    assert flags == next(q for q in again.qc if q.sample_id == "SELF").flags
