"""Warnings for a cM request the map cannot honour and for assembly lengths (docs/adr/0015)."""

from __future__ import annotations

import dataclasses

import numpy as np

from progeny_selector import run_analysis
from progeny_selector.core.background import rpp
from progeny_selector.model.criteria import BackgroundOptions, Criteria, Filters, TargetSpec
from progeny_selector.model.dataset import Dataset, GenotypeMatrix, Marker, Sample
from tests.conftest import make_dataset, make_matrix

STATES = ["A", "H", "B", "A", "A", "H"]
NO_LENGTH = "chromosome {c}: the assembly {a} gives no length; the last marker is the chromosome end for drag bounds"


def without_cm(gm: GenotypeMatrix) -> GenotypeMatrix:
    """The same matrix with no cM on any marker, which is what disables cM mode for the dataset."""
    markers = [Marker(m.marker_id, m.chrom, m.pos_bp, None) for m in gm.markers]
    return GenotypeMatrix(markers=markers, sample_ids=list(gm.sample_ids), alleles=gm.alleles, calls=gm.calls)


def criteria_cm() -> Criteria:
    return Criteria(
        targets=[TargetSpec("T1", marker_id="m2", required_state="either")],
        flank_unit="cm",
        background=BackgroundOptions(model="count", map_unit="cm"),
    )


def test_cm_requested_without_a_map_warns_for_both_settings():
    result = run_analysis(make_dataset(without_cm(make_matrix(STATES))), criteria_cm())
    assert "background.map_unit is cm but the map has no cM; RPP weights computed in bp" in result.warnings
    assert "flank_unit is cm but the map has no cM; windows interpreted in bp" in result.warnings
    assert result.unit == "bp"


def test_weighted_model_in_bp_reports_the_cap_only_when_cm_was_wanted():
    cap = "weighted RPP in bp: no cM map, weights capped at 4000000 bp per marker"
    for map_unit in ("cm", "auto"):
        criteria = dataclasses.replace(criteria_cm(), background=BackgroundOptions(model="weighted", map_unit=map_unit))
        result = run_analysis(make_dataset(without_cm(make_matrix(STATES))), criteria)
        assert cap in result.warnings
    # bp asked for explicitly on a dataset that has cM: the unit is what the user chose, so no warning.
    criteria = dataclasses.replace(criteria_cm(), flank_unit="bp", background=BackgroundOptions(model="weighted", map_unit="bp"))
    result = run_analysis(make_dataset(make_matrix(STATES)), criteria)
    assert result.unit == "bp"
    assert result.warnings == []


def test_marker_beyond_the_assembly_length_warns_and_keeps_drag_bounds_non_negative():
    gm = without_cm(make_matrix(STATES))
    markers = list(gm.markers)
    markers[-1] = Marker(markers[-1].marker_id, "Gm01", 60_000_000, None)
    gm = GenotypeMatrix(markers=markers, sample_ids=list(gm.sample_ids), alleles=gm.alleles, calls=gm.calls)
    result = run_analysis(make_dataset(gm), criteria_cm())
    assert any("beyond the Wm82.a4 length 57932356; check the assembly setting" in w for w in result.warnings)
    assert any(w.startswith("chromosome Gm01: marker positions reach 60000000") for w in result.warnings)
    drag = result.drag["T1"]
    assert np.all(np.nan_to_num(drag.right_max) >= 0)


def test_chromosome_without_an_assembly_length_warns():
    gm = without_cm(make_matrix(STATES, chrom="scaffold_1"))
    result = run_analysis(make_dataset(gm), criteria_cm())
    assert NO_LENGTH.format(c="scaffold_1", a="Wm82.a4") in result.warnings


def many_chromosome_dataset(n: int, prefix: str, pos_bp: int | None = None) -> Dataset:
    """One informative marker on each of ``n`` chromosomes named ``prefix``1..n, plus a target marker."""
    markers = [Marker("m1", "Gm01", 1_000_000, None)]
    markers += [Marker(f"x{i + 1}", f"{prefix}{i + 1}", 1_000_000 if pos_bp is None else pos_bp, None) for i in range(n)]
    calls = np.zeros((len(markers), 3, 2), dtype=np.int8)
    for m in range(len(markers)):
        calls[m, 0] = (0, 0)
        calls[m, 1] = (2, 2)
        calls[m, 2] = (0, 2)
    gm = GenotypeMatrix(markers=markers, sample_ids=["RP", "DONOR", "P1"], alleles=[["A", "G", "T"]] * len(markers), calls=calls)
    samples = [
        Sample("RP", "RP", "recurrent_parent"),
        Sample("DONOR", "DONOR", "donor_parent"),
        Sample("P1", "P1", "progeny", "BC2F1", "F1"),
    ]
    return Dataset(genotypes=gm, samples=samples)


def criteria_marker(marker_id: str) -> Criteria:
    return dataclasses.replace(criteria_cm(), targets=[TargetSpec("T1", marker_id=marker_id, required_state="either")])


def test_missing_length_warnings_are_capped_and_the_continuation_names_itself():
    result = run_analysis(many_chromosome_dataset(7, "scaffold_"), criteria_marker("m1"))
    listed = [w for w in result.warnings if w.startswith("chromosome scaffold_")]
    assert len(listed) == 5
    assert "... and 2 more chromosomes to which the assembly Wm82.a4 gives no length" in result.warnings


def test_beyond_length_warnings_are_capped_the_same_way():
    # 70 Mb is past every Wm82.a4 chromosome, so each of the six named chromosomes is beyond its length.
    result = run_analysis(many_chromosome_dataset(6, "Gm1", pos_bp=70_000_000), criteria_marker("m1"))
    listed = [w for w in result.warnings if "marker positions reach" in w]
    assert len(listed) == 5
    assert "... and 1 more chromosomes with marker positions beyond the Wm82.a4 length" in result.warnings


def test_terminal_marker_keeps_the_cap_weight_on_a_chromosome_without_a_length():
    """Two markers on an unplaced chromosome: outer weights are the cap's half, not the distance to the last marker."""
    gm = without_cm(make_matrix(["A", "B"], chrom="scaffold_1"))
    criteria = dataclasses.replace(
        criteria_cm(),
        targets=[TargetSpec("T1", marker_id="m1", required_state="either")],
        flank_unit="bp",
        background=BackgroundOptions(model="weighted", map_unit="bp"),
        filters=Filters(max_hom_donor_rate_bcf1=1.0, het_rate_tolerance=1.0, max_nonparental_rate=1.0),
    )
    result = run_analysis(make_dataset(gm), criteria)
    # cap/2 = 2,000,000 on each outer side and 500,000 inside, so both markers weigh 2,500,000 and RPP is 0.5.
    # Ending the chromosome at the last marker would give weights 1,500,000 and 500,000, and RPP 0.75.
    assert result.row("P1")["rpp_total"] == 0.5
    assert rpp(result.classification.states[:, 2:]) == 0.5  # the count model, for comparison
