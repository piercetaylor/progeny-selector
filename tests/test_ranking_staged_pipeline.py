"""``ranking.mode: staged`` through ``run_analysis`` on a two-chromosome hand-built dataset."""

from __future__ import annotations

import dataclasses

import numpy as np

from progeny_selector.core.pipeline import run_analysis
from progeny_selector.model.criteria import Criteria, Filters, RankingOptions, TargetSpec
from progeny_selector.model.dataset import Dataset, GenotypeMatrix, Marker, Sample

CODE = {"A": (0, 0), "H": (0, 2), "B": (2, 2)}
# Four markers on the carrier chromosome (the target is c2) and four on a non-carrier chromosome, so
# rpp_noncarrier is a real number rather than NaN and the staged order can use it.
LAYOUT = (("Gm06", "c", 4), ("Gm01", "n", 4))


def build(states: dict[str, str]) -> Dataset:
    markers = [Marker(f"{prefix}{i + 1}", chrom, (i + 1) * 1_000_000, (i + 1) * 2.5) for chrom, prefix, k in LAYOUT for i in range(k)]
    ids = ["RP", "DONOR", *states]
    calls = np.zeros((len(markers), len(ids), 2), dtype=np.int8)
    for m in range(len(markers)):
        calls[m, 0] = (0, 0)
        calls[m, 1] = (2, 2)
        for j, sid in enumerate(states, start=2):
            calls[m, j] = CODE[states[sid][m]]
    gm = GenotypeMatrix(markers=markers, sample_ids=ids, alleles=[["A", "G", "T"]] * len(markers), calls=calls)
    samples = [Sample("RP", "RP", "recurrent_parent"), Sample("DONOR", "DONOR", "donor_parent")]
    samples += [Sample(s, s, "progeny", "BC2F1", "F1") for s in states]
    return Dataset(genotypes=gm, samples=samples)


# P1 recovers the carrier chromosome, P2 the rest of the genome; both flanks of the target recombine in both.
DATASET = build({"P1": "AHAA" + "BBBA", "P2": "BHBB" + "AAAA"})
CRITERIA = Criteria(
    targets=[TargetSpec("T1", marker_id="c2", required_state="either")],
    filters=Filters(max_hom_donor_rate_bcf1=1.0, het_rate_tolerance=1.0, max_nonparental_rate=1.0),
)


def ranks(mode: str) -> dict[str, float]:
    result = run_analysis(DATASET, dataclasses.replace(CRITERIA, ranking=RankingOptions(mode=mode)))
    assert all(row["passes_filters"] for row in result.rows)
    return {row["sample_id"]: row["rank_overall"] for row in result.rows}


def test_weighted_ranks_by_the_composite_score():
    # P2 wins on the composite: rpp_noncarrier carries weight 0.5 against rpp_carrier's 0.2.
    assert ranks("weighted") == {"P2": 1.0, "P1": 2.0}


def test_staged_puts_the_carrier_chromosome_first():
    # Equal recombinant flank counts, so rpp_carrier decides and the order reverses.
    assert ranks("staged") == {"P1": 1.0, "P2": 2.0}


def test_staged_keeps_the_composite_score_and_the_rest_of_the_row():
    weighted = run_analysis(DATASET, dataclasses.replace(CRITERIA, ranking=RankingOptions(mode="weighted")))
    staged = run_analysis(DATASET, dataclasses.replace(CRITERIA, ranking=RankingOptions(mode="staged")))
    by_id = {row["sample_id"]: row for row in staged.rows}
    for row in weighted.rows:
        other = by_id[row["sample_id"]]
        assert other["composite_score"] == row["composite_score"]
        assert other["rpp_carrier"] == row["rpp_carrier"]
        assert other["rank_in_family"] == other["rank_overall"]  # one family
