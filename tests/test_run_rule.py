"""Foreground rule: run (docs/adr/0011): hand-built cases, a randomised reference comparison, validation, pipeline."""

from __future__ import annotations

import random

import numpy as np
import pytest

from progeny_selector.constants import LABEL_TO_STATE, STATUS_FAIL, STATUS_LABELS, STATUS_PASS, STATUS_UNKNOWN
from progeny_selector.core.foreground import ResolvedLocus, foreground_status, locus_status, marker_predicate
from progeny_selector.core.pipeline import run_analysis
from progeny_selector.io.criteria import criteria_to_dict, dump_criteria_yaml, read_criteria_text
from progeny_selector.model.criteria import AvoidSpec, Criteria, CriteriaError, TargetSpec
from tests.conftest import make_dataset, make_matrix

MB = 1_000_000


def status_of(calls: str, required_state: str = "hom_donor", min_run: int = 3, tolerate: bool = True, min_markers: int = 1) -> int:
    """One sample; markers at 1..N Mb in string order. '|' marks an anchor between two markers, '^' before a
    character puts the anchor exactly on that marker."""
    anchor = None
    labels: list[str] = []
    k = 0
    while k < len(calls):
        ch = calls[k]
        if ch == "|":
            anchor = len(labels) * MB + MB // 2
        elif ch == "^":
            anchor = (len(labels) + 1) * MB
        else:
            labels.append(ch)
        k += 1
    assert anchor is not None
    n = len(labels)
    states = np.array([[LABEL_TO_STATE[c]] for c in labels], dtype=np.int8).reshape(n, 1)
    positions = np.arange(1, n + 1, dtype=np.int64) * MB
    resolved = ResolvedLocus("T", "region", "Gm01", 0, (n + 1) * MB, np.arange(n))
    out = locus_status(
        states,
        resolved,
        marker_predicate(states, required_state),
        "run",
        min_markers,
        positions=positions,
        min_run=min_run,
        anchor_bp=anchor,
        tolerate_isolated=tolerate,
    )
    return int(out[0])


def test_contiguous_block_min_run():
    assert status_of("BB|BBB", min_run=5) == STATUS_PASS
    assert status_of("BB|BBB", min_run=6) == STATUS_FAIL


def test_run_across_anchor_and_alternating_calls():
    assert status_of("AABB|BAA", min_run=3) == STATUS_PASS
    assert status_of("ABA|BAB", min_run=2, tolerate=False) == STATUS_FAIL
    # Tolerance on: both inner A calls bridge, the run spans calls 2-6 and counts its 3 B calls.
    assert status_of("ABA|BAB", min_run=3, tolerate=True) == STATUS_PASS
    assert status_of("ABA|BAB", min_run=4, tolerate=True) == STATUS_FAIL


def test_isolated_mismatch_at_anchor_neighbour():
    assert status_of("BBB|ABBB", min_run=5, tolerate=True) == STATUS_PASS
    assert status_of("BBB|ABBB", min_run=6, tolerate=True) == STATUS_PASS
    assert status_of("BBB|ABBB", min_run=7, tolerate=True) == STATUS_FAIL
    assert status_of("BBB|ABBB", min_run=1, tolerate=False) == STATUS_FAIL


def test_two_mismatches_do_not_bridge():
    assert status_of("BBAA|BBB", min_run=1) == STATUS_FAIL


def test_run_entirely_right_of_anchor():
    assert status_of("AAA|ABBBBBB", min_run=1) == STATUS_FAIL


def test_missing_uninformative_and_nonparental_are_ignored():
    for gap in "NUX":
        assert status_of(f"BB{gap}|BB", min_run=4, tolerate=False) == STATUS_PASS
        assert status_of(f"BB{gap}|BB", min_run=5, tolerate=False) == STATUS_FAIL
        assert status_of(f"B{gap}B^{gap}B", min_run=3, tolerate=False) == STATUS_PASS


def test_anchor_on_marker():
    assert status_of("BB^BBB", min_run=5, tolerate=False) == STATUS_PASS
    assert status_of("BB^ABB", min_run=4, tolerate=True) == STATUS_PASS
    assert status_of("BB^ABB", min_run=5, tolerate=True) == STATUS_FAIL
    assert status_of("BB^ABB", min_run=1, tolerate=False) == STATUS_FAIL


def test_no_counted_marker_on_one_side_and_min_markers():
    assert status_of("|BBBB", min_run=1) == STATUS_FAIL
    assert status_of("BBBB|", min_run=1) == STATUS_FAIL
    assert status_of("NN|BBBB", min_run=1) == STATUS_FAIL
    assert status_of("BB|BB", min_run=1, min_markers=5) == STATUS_UNKNOWN
    assert status_of("BBN|NBB", min_run=1, min_markers=5) == STATUS_UNKNOWN
    assert status_of("NN|NN", min_run=1) == STATUS_UNKNOWN


def test_required_states():
    assert status_of("AHH|HA", required_state="het", min_run=3) == STATUS_PASS
    assert status_of("ABB|BA", required_state="het", min_run=1) == STATUS_FAIL
    assert status_of("AHB|HBA", required_state="either", min_run=4) == STATUS_PASS
    assert status_of("AHB|HBA", required_state="hom_donor", min_run=1, tolerate=False) == STATUS_FAIL


def test_positions_required_and_unsorted_marker_idx():
    states = np.array([[2], [2], [0]], dtype=np.int8)  # markers 0, 1, 2 = B, B, A
    resolved = ResolvedLocus("T", "region", "Gm01", 0, 4 * MB, np.array([2, 0, 1]))
    pred = marker_predicate(states, "hom_donor")
    with pytest.raises(ValueError, match="positions"):
        locus_status(states, resolved, pred, "run")
    positions = np.array([1, 2, 3]) * MB  # indexed by marker: position order is B@1, B@2, A@3

    def run(anchor: int) -> int:
        out = locus_status(states, resolved, pred, "run", positions=positions, min_run=2, anchor_bp=anchor, tolerate_isolated=False)
        return int(out[0])

    # Reading marker_idx order (A, B, B) as position order would flip both results.
    assert run(MB + MB // 2) == STATUS_PASS
    assert run(2 * MB + MB // 2) == STATUS_FAIL


def tie_status(calls: str, positions_mb: list[int], anchor_mb: int, min_run: int, tolerate: bool) -> int:
    n = len(calls)
    states = np.array([LABEL_TO_STATE[c] for c in calls], dtype=np.int8).reshape(n, 1)
    resolved = ResolvedLocus("T", "region", "Gm01", 0, 10 * MB, np.arange(n))
    out = locus_status(
        states,
        resolved,
        marker_predicate(states, "hom_donor"),
        "run",
        positions=np.array(positions_mb, dtype=np.int64) * MB,
        min_run=min_run,
        anchor_bp=anchor_mb * MB,
        tolerate_isolated=tolerate,
    )
    return int(out[0])


def test_ties_on_the_anchor_must_all_join():
    assert tie_status("BBAABB", [1, 2, 2, 2, 2, 3], 2, min_run=2, tolerate=False) == STATUS_FAIL
    assert tie_status("BBAABB", [1, 2, 2, 2, 2, 3], 2, min_run=1, tolerate=True) == STATUS_FAIL
    assert tie_status("BBABB", [1, 2, 2, 2, 3], 2, min_run=1, tolerate=False) == STATUS_FAIL
    assert tie_status("BBABB", [1, 2, 2, 2, 3], 2, min_run=4, tolerate=True) == STATUS_PASS
    assert tie_status("BBABB", [1, 2, 2, 2, 3], 2, min_run=5, tolerate=True) == STATUS_FAIL
    assert tie_status("ABBBA", [1, 2, 2, 2, 3], 2, min_run=3, tolerate=False) == STATUS_PASS


def test_terminal_mismatch_is_not_bridged():
    assert status_of("BB|A", min_run=1, tolerate=True) == STATUS_FAIL
    assert status_of("A|BB", min_run=1, tolerate=True) == STATUS_FAIL


def test_flanking_locus_ignores_run_kwargs():
    states = np.array([[2, 0], [2, 2]], dtype=np.int8)
    resolved = ResolvedLocus("F", "flanking", "Gm01", MB, 2 * MB, np.array([0, 1]))
    out = foreground_status(states, resolved, "hom_donor", "run", 1, positions=None, min_run=99, anchor_bp=0, tolerate_isolated=False)
    assert list(out) == [STATUS_PASS, STATUS_FAIL]


def reference_status(labels: list[str], positions: list[int], anchor: int, required: str, min_run: int, tolerate: bool, min_markers: int):
    """Independent per-sample loop over the written semantics."""
    good = {"hom_donor": {"B"}, "het": {"H"}, "either": {"H", "B"}}[required]
    c: list[str] = []
    p: list[int] = []
    for lab, pos in zip(labels, positions, strict=True):
        if lab in ("A", "H", "B"):
            c.append(lab)
            p.append(pos)
    n = len(c)
    if n < min_markers:
        return "unknown"
    sat = [x in good for x in c]
    joined = list(sat)
    if tolerate:
        for k in range(1, n - 1):
            if not sat[k] and sat[k - 1] and sat[k + 1]:
                joined[k] = True
    left = [k for k in range(n) if p[k] <= anchor]
    right = [k for k in range(n) if p[k] >= anchor]
    if not left or not right:
        return "fail"
    # Required calls: the nearest call on each side plus every call sitting exactly on the anchor.
    required = {max(left), min(right)} | {k for k in range(n) if p[k] == anchor}
    stretches: list[tuple[int, int]] = []
    k = 0
    while k < n:
        if joined[k]:
            start = k
            while k + 1 < n and joined[k + 1]:
                k += 1
            stretches.append((start, k))
        k += 1
    for start, end in stretches:
        if all(start <= r <= end for r in required):
            length = sum(1 for m in range(start, end + 1) if sat[m])
            return "pass" if length >= min_run else "fail"
    return "fail"


def test_randomised_against_reference_loop():
    rng = random.Random(20260915)
    for _ in range(2000):
        n = rng.randint(0, 40)
        labels = [rng.choice("AHBNUX") for _ in range(n)]
        positions = sorted(rng.choices(range(1, 25), k=n))  # ties are common
        anchor = rng.choice(positions) if positions and rng.random() < 0.5 else rng.randint(0, 25)
        required = rng.choice(["hom_donor", "het", "either"])
        min_run = rng.randint(1, 6)
        tolerate = rng.choice([True, False])
        min_markers = rng.randint(1, 4)
        states = np.array([LABEL_TO_STATE[x] for x in labels], dtype=np.int8).reshape(n, 1)
        resolved = ResolvedLocus("T", "region", "Gm01", 0, 100, np.arange(n, dtype=int))
        got = locus_status(
            states,
            resolved,
            marker_predicate(states, required),
            "run",
            min_markers,
            positions=np.array(positions, dtype=np.int64),
            min_run=min_run,
            anchor_bp=anchor,
            tolerate_isolated=tolerate,
        )
        expected = reference_status(labels, positions, anchor, required, min_run, tolerate, min_markers)
        assert STATUS_LABELS[int(got[0])] == expected, (labels, positions, anchor, required, min_run, tolerate, min_markers)


REGION_TARGET = "targets:\n  - locus_id: T1\n    chrom: Gm01\n    start_bp: 1000000\n    end_bp: 9000000\n"


@pytest.mark.parametrize(
    ("text", "match"),
    [
        ("targets:\n  - locus_id: T1\n    marker_id: m1\n    rule: run\n", "rule run needs a region locus"),
        ("targets:\n  - locus_id: T1\n    left_marker: m1\n    right_marker: m2\n    rule: run\n", "rule run needs a region locus"),
        (REGION_TARGET + "avoid:\n  - locus_id: A1\n    marker_id: m2\n    rule: run\n", "'A1'.*run is for targets only"),
        (REGION_TARGET + "avoid:\n  - locus_id: A1\n    marker_id: m2\n    min_run: 3\n", "unknown keys"),
        (REGION_TARGET + "    rule: run\n    min_run: 0\n", "min_run must be >= 1"),
        (REGION_TARGET + "    rule: run\n    anchor_bp: 9000001\n", "anchor_bp"),
        (REGION_TARGET + "    rule: run\n    anchor_bp: 999999\n", "anchor_bp"),
        (REGION_TARGET + "    rule: any\n    min_run: 3\n", "only with rule: run"),
        (REGION_TARGET + "    tolerate_isolated: false\n", "only with rule: run"),
        (REGION_TARGET + "    rule: run\n    min_run: 2.5\n", "min_run must be an integer"),
        (REGION_TARGET + "    rule: run\n    anchor_bp: '5,000,000'\n", "anchor_bp must be an integer"),
        (REGION_TARGET + "    rule: run\n    tolerate_isolated: 'yes'\n", "tolerate_isolated must be true or false"),
    ],
)
def test_run_validation(text, match):
    with pytest.raises(CriteriaError, match=match):
        read_criteria_text(text)


def test_run_defaults_and_model_validation():
    c = read_criteria_text(REGION_TARGET + "    rule: run\n")
    t = c.targets[0]
    assert (t.rule, t.min_run, t.anchor_bp, t.tolerate_isolated) == ("run", 3, None, True)
    with pytest.raises(CriteriaError, match="run is for targets only"):
        AvoidSpec("A", marker_id="m1", rule="run").validate()
    with pytest.raises(CriteriaError, match="min_run"):
        TargetSpec("T", chrom="Gm01", start_bp=1, end_bp=9, rule="run", min_run=0).validate()
    for kw in ({"min_run": 4}, {"anchor_bp": 5}, {"tolerate_isolated": False}):
        with pytest.raises(CriteriaError, match="apply only with rule run"):
            TargetSpec("T", chrom="Gm01", start_bp=1, end_bp=9, rule="any", **kw).validate()  # type: ignore[arg-type]
    TargetSpec("T", chrom="Gm01", start_bp=1, end_bp=9, rule="all").validate()


def test_region_shorthand_with_run():
    c = read_criteria_text('targets:\n  - locus_id: T1\n    region: "Gm08:1,000,000-3,000,000"\n    rule: run\n    anchor_bp: 2500000\n')
    t = c.targets[0]
    assert (t.chrom, t.start_bp, t.end_bp, t.rule, t.anchor_bp) == ("Gm08", 1_000_000, 3_000_000, "run", 2_500_000)
    assert read_criteria_text(dump_criteria_yaml(c)) == c


def test_run_roundtrip_and_canonical_dump():
    c = Criteria(
        targets=[
            TargetSpec("R", chrom="Gm01", start_bp=1_000_001, end_bp=9_000_000, rule="run", min_run=4, anchor_bp=2_000_000),
            TargetSpec("D", chrom="Gm02", start_bp=1_000_001, end_bp=9_000_000, rule="run", tolerate_isolated=False),
            TargetSpec("ALL", chrom="Gm03", start_bp=1, end_bp=9),
            TargetSpec("ANY", marker_id="m1", rule="any"),
        ],
        avoid=[AvoidSpec("AV", marker_id="m2")],
    )
    doc = criteria_to_dict(c)
    r, d = doc["targets"][0], doc["targets"][1]
    assert list(r) == [
        "locus_id",
        "chrom",
        "start_bp",
        "end_bp",
        "rule",
        "min_markers",
        "min_run",
        "anchor_bp",
        "tolerate_isolated",
        "required_state",
    ]
    assert (r["min_run"], r["anchor_bp"], r["tolerate_isolated"]) == (4, 2_000_000, True)
    assert (d["min_run"], d["anchor_bp"], d["tolerate_isolated"]) == (3, 5_000_000, False)  # (1_000_001 + 9_000_000) // 2
    for locus in doc["targets"][2:] + doc["avoid"]:
        assert not {"min_run", "anchor_bp", "tolerate_isolated"} & set(locus)
    text = dump_criteria_yaml(c)
    back = read_criteria_text(text)
    # A None anchor is written resolved, so the reloaded spec carries the integer; everything else is equal.
    c.targets[1].anchor_bp = 5_000_000
    assert back == c
    assert dump_criteria_yaml(back) == text


def test_pipeline_run_target_status_column():
    # P1 at 1..10 Mb: donor segment 3..7 Mb with one recurrent call at 5 Mb.
    gm = make_matrix(["A", "A", "B", "B", "A", "B", "B", "A", "A", "A"])
    dataset = make_dataset(gm)

    def status(**kw) -> str:
        crit = Criteria(
            targets=[TargetSpec("T", chrom="Gm01", start_bp=2 * MB, end_bp=8 * MB, rule="run", required_state="hom_donor", **kw)]
        )
        return run_analysis(dataset, crit).row("P1")["target_T_status"]

    assert status() == "pass"  # anchor 5 Mb sits on the bridged A; run of 4 B
    assert status(min_run=5) == "fail"
    assert status(tolerate_isolated=False) == "fail"
    assert status(anchor_bp=3 * MB, tolerate_isolated=False, min_run=2) == "pass"
