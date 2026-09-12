"""Navigation tree, row filtering and breadcrumb labels."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from progeny_selector.core.navigation import FamilyNode, GenerationNode, build_tree, crumb_labels, filter_rows
from progeny_selector.core.pipeline import run_analysis
from progeny_selector.io import load_dataset, read_criteria
from progeny_selector.model.dataset import Dataset, Sample
from tests.conftest import make_matrix


@pytest.fixture(scope="module")
def fixture_data(fixture_dir: Path):
    dataset = load_dataset(fixture_dir / "genotypes.vcf", fixture_dir / "samples.csv", fixture_dir / "markers.csv")
    result = run_analysis(dataset, read_criteria(fixture_dir / "criteria.yaml"))
    return dataset, result


def test_fixture_tree(fixture_data):
    dataset, _ = fixture_data
    tree = build_tree(dataset)
    assert tree.cross == "Williams 82 (synthetic) x PI synthetic donor"
    assert tree.n == 40
    assert tree.families == (
        FamilyNode("F1", 20, (GenerationNode("BC2F1", 20),)),
        FamilyNode("F2", 20, (GenerationNode("BC2F1", 20),)),
    )


def test_missing_family_sorts_last():
    gm = make_matrix(["A", "H"])
    gm = gm.with_samples_added(["P0", "P2"], np.zeros((2, 2, 2), dtype=np.int8))
    samples = [
        Sample("RP", "RP", "recurrent_parent"),
        Sample("DONOR", "DONOR", "donor_parent"),
        Sample("P0", "P0", "progeny", "BC1F1", None),
        Sample("P1", "P1", "progeny", "BC1F1", "F1"),
        Sample("P2", "P2", "progeny", None, "F1"),
    ]
    tree = build_tree(Dataset(genotypes=gm, samples=samples))
    assert tree.cross == "RP x DONOR"
    assert tree.n == 3
    assert tree.families[0] == FamilyNode("F1", 2, (GenerationNode("BC1F1", 1), GenerationNode(None, 1)))
    assert tree.families[-1] == FamilyNode(None, 1, (GenerationNode("BC1F1", 1),))


def test_filter_rows(fixture_data):
    _, result = fixture_data
    rows = result.rows
    f1 = filter_rows(rows, "F1", None)
    assert len(f1) == 20
    assert f1 == [r for r in rows if r["family_id"] == "F1"]
    assert filter_rows(rows, None, "BC3F1") == []
    assert filter_rows(rows, None, None) == rows


def test_crumb_labels(fixture_data):
    tree = build_tree(fixture_data[0])
    cross = "Williams 82 (synthetic) x PI synthetic donor"
    assert crumb_labels(tree, None, None) == [cross]
    assert crumb_labels(tree, "F1", None) == [cross, "F1"]
    assert crumb_labels(tree, "F1", "BC2F1") == [cross, "F1", "BC2F1"]
