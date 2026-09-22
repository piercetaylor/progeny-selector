"""criteria.yaml serialisation: canonical dump and round trip through the reader."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from progeny_selector.io import read_criteria
from progeny_selector.io.criteria import criteria_from_dict, criteria_to_dict, dump_criteria_yaml, read_criteria_text
from progeny_selector.model.criteria import AvoidSpec, BackgroundOptions, Criteria, CriteriaError, TargetSpec


@pytest.fixture(scope="module")
def fixture_criteria(fixture_dir: Path) -> Criteria:
    return read_criteria(fixture_dir / "criteria.yaml")


def hand_built() -> Criteria:
    return Criteria(
        targets=[
            TargetSpec("T1", chrom="Gm06", start_bp=12_000_000, end_bp=14_000_000, required_state="het", flank_left=2.5),
            TargetSpec("RUN", chrom="Gm08", start_bp=1_000_000, end_bp=3_000_000, rule="run", min_run=4, anchor_bp=1_500_000),
        ],
        avoid=[AvoidSpec("AV1", left_marker="m1", right_marker="m2", allow_het=True, notes="flanking avoid")],
        background=BackgroundOptions(max_marker_coverage=8),
    )


def test_dict_roundtrip(fixture_criteria):
    for c in (fixture_criteria, hand_built()):
        assert criteria_from_dict(criteria_to_dict(c)) == c


def test_yaml_roundtrip(fixture_criteria):
    for c in (fixture_criteria, hand_built()):
        text = dump_criteria_yaml(c)
        assert "\r" not in text
        assert read_criteria_text(text) == c


def test_dump_starts_with_name(fixture_criteria):
    assert dump_criteria_yaml(fixture_criteria).startswith("name: synthetic BC2F1 fixture")
    assert dump_criteria_yaml(hand_built()).startswith("name: null")  # an unset optional key is explicit null


def test_canonical_form():
    doc = criteria_to_dict(hand_built())
    assert list(doc) == [
        "name",
        "targets",
        "avoid",
        "flank_window",
        "flank_unit",
        "assembly",
        "background",
        "ranking",
        "weights",
        "filters",
    ]
    target = doc["targets"][0]
    assert list(target) == [
        "locus_id",
        "chrom",
        "start_bp",
        "end_bp",
        "rule",
        "min_markers",
        "notes",
        "required_state",
        "flank_left",
        "flank_right",
    ]
    assert (doc["name"], target["notes"], target["flank_right"]) == (None, None, None)
    assert list(doc["avoid"][0]) == ["locus_id", "left_marker", "right_marker", "rule", "min_markers", "notes", "allow_het"]
    assert doc["background"] == {"model": "weighted", "map_unit": "auto", "max_marker_coverage": 8}
    assert doc["ranking"] == {"mode": "weighted"}
    assert doc["assembly"] is None  # unset: the crop decides (core.pipeline.resolve_assembly)
    assert criteria_to_dict(Criteria(targets=[TargetSpec("T", marker_id="m")]))["background"]["max_marker_coverage"] is None
    text = dump_criteria_yaml(hand_built())
    assert "region" not in text
    loaded = yaml.safe_load(text)["targets"][0]
    assert (loaded["chrom"], loaded["start_bp"], loaded["end_bp"]) == ("Gm06", 12_000_000, 14_000_000)


def test_fixture_reads_unchanged(fixture_dir: Path, fixture_criteria):
    c = read_criteria(fixture_dir / "criteria.yaml")
    assert c == fixture_criteria
    assert (c.name, c.flank_window, c.flank_unit, c.background.model) == ("synthetic BC2F1 fixture", 6.0, "cm", "count")
    assert c.targets[0].min_markers == 1 and c.filters.exclude_qc_flagged is True


def test_malformed_yaml_raises_criteria_error():
    with pytest.raises(CriteriaError, match="not valid YAML"):
        read_criteria_text("targets:\n  - locus_id: T1\n marker_id: m1\n")


@pytest.mark.parametrize(
    ("text", "key"),
    [
        ("targets:\n  - locus_id: T1\n    marker_id: m1\n    min_markers: null\n", "min_markers"),
        ("targets:\n  - locus_id: T1\n    marker_id: m1\nweights: {rpp_noncarrier: 'a'}\n", "rpp_noncarrier"),
        ("targets:\n  - locus_id: T1\n    marker_id: m1\n    min_markers: true\n", "min_markers"),
        ("targets:\n  - locus_id: T1\n    marker_id: m1\n    flank_left: '2'\n", "flank_left"),
        ("targets:\n  - locus_id: T1\n    marker_id: m1\navoid:\n  - locus_id: A\n    marker_id: m2\n    allow_het: 'yes'\n", "allow_het"),
        ("targets:\n  - locus_id: T1\n    marker_id: m1\nfilters: {exclude_qc_flagged: 1}\n", "exclude_qc_flagged"),
        ("targets:\n  - locus_id: T1\n    marker_id: m1\nfilters: {max_missing_rate: null}\n", "max_missing_rate"),
        ("targets:\n  - locus_id: T1\n    marker_id: m1\nbackground: {max_marker_coverage: false}\n", "max_marker_coverage"),
        ("targets:\n  - locus_id: T1\n    marker_id: m1\nflank_window: six\n", "flank_window"),
    ],
)
def test_wrong_typed_fields_raise_criteria_error(text: str, key: str):
    with pytest.raises(CriteriaError, match=key):
        read_criteria_text(text)


T1 = "targets:\n  - locus_id: T1\n    marker_id: m1\n"


@pytest.mark.parametrize(
    ("text", "key"),
    [
        ("targets:\n  - locus_id: T1\n    chrom: Gm06\n    start_bp: twelve\n    end_bp: 14000000\n", "start_bp"),
        ("targets:\n  - locus_id: T1\n    chrom: Gm06\n    start_bp: 12000000\n    end_bp: [1]\n", "end_bp"),
        ("targets:\n  - locus_id: T1\n    region: Gm06-12000000-14000000\n", "region"),
        ("targets:\n  - locus_id: T1\n    region: 'Gm06:12000000'\n", "region"),
        ("targets:\n  - locus_id: T1\n    region: 'Gm06:a-b'\n", "region"),
        ("targets:\n  - locus_id: T1\n    region: 6\n", "region"),
        (T1 + "weights: 5\n", "weights"),
        (T1 + "filters: [max_missing_rate]\n", "filters"),
        (T1 + "background: count\n", "background"),
        ("targets: T1\n", "targets"),
        ("targets:\n  - T1\n", r"targets\[0\]"),
        (T1 + "avoid: {locus_id: AV1}\n", "avoid"),
        (T1 + "avoid:\n  - 7\n", r"avoid\[0\]"),
    ],
)
def test_malformed_shapes_raise_criteria_error(text: str, key: str):
    with pytest.raises(CriteriaError, match=key):
        read_criteria_text(text)


def test_region_shorthand_still_expands():
    c = read_criteria_text("targets:\n  - locus_id: T1\n    region: 'Gm06:12,000,000-14,000,000'\n")
    assert (c.targets[0].chrom, c.targets[0].start_bp, c.targets[0].end_bp) == ("Gm06", 12_000_000, 14_000_000)


def test_read_criteria_text_rejects_non_mappings():
    with pytest.raises(CriteriaError, match="must be a mapping"):
        read_criteria_text("")
    with pytest.raises(CriteriaError, match="must be a mapping"):
        read_criteria_text("- a")


def test_ranking_and_assembly_round_trip(fixture_criteria):
    text = dump_criteria_yaml(fixture_criteria)
    assert "ranking:\n  mode: weighted" in text
    assert "assembly: null" in text  # unset in the fixture, so it round-trips as unset and not as the resolved value
    assert read_criteria_text(text).assembly is None
    assert read_criteria_text(T1 + "assembly: Wm82.a2\n").assembly == "Wm82.a2"
    assert read_criteria_text(text).ranking.mode == "weighted"
    assert read_criteria_text(T1 + "ranking: {mode: staged}\nassembly: none\n").ranking.mode == "staged"


@pytest.mark.parametrize(
    ("text", "key"),
    [
        (T1 + "ranking: {mode: lexicographic}\n", "ranking.mode"),
        (T1 + "assembly: Wm82.a9\n", "assembly"),
        (T1 + "assembly: 4\n", "assembly must be text"),
        (T1 + "ranking: [mode]\n", "ranking must be a mapping"),
    ],
)
def test_bad_ranking_or_assembly_raises(text: str, key: str):
    with pytest.raises(CriteriaError, match=key):
        read_criteria_text(text)
