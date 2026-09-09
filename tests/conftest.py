"""Shared fixtures: paths to the synthetic BC2F1 dataset and small hand-built matrices."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from progeny_selector.model.dataset import Dataset, GenotypeMatrix, Marker, Sample

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "synthetic_bc2f1"


@pytest.fixture(scope="session")
def fixture_dir() -> Path:
    return FIXTURE_DIR


@pytest.fixture(scope="session")
def expected_rows() -> dict[str, dict]:
    with open(FIXTURE_DIR / "expected_results.csv", newline="") as fh:
        return {r["sample_id"]: r for r in csv.DictReader(fh)}


def make_matrix(states: list[str], chrom: str = "Gm01", spacing_bp: int = 1_000_000, cm_per_mb: float = 2.5) -> GenotypeMatrix:
    """Build a tiny nucleotide matrix from per-marker progeny states ('A','H','B','N','X') for one progeny.

    Column order: RP, DONOR, P1. Alleles: RP = 'A', donor = 'T', non-parental = 'G'.
    """
    markers = [Marker(f"m{i + 1}", chrom, (i + 1) * spacing_bp, (i + 1) * spacing_bp / 1e6 * cm_per_mb) for i in range(len(states))]
    calls = np.zeros((len(states), 3, 2), dtype=np.int8)
    alleles = [["A", "G", "T"] for _ in states]
    for i, s in enumerate(states):
        calls[i, 0] = (0, 0)
        calls[i, 1] = (2, 2)
        calls[i, 2] = {"A": (0, 0), "H": (0, 2), "B": (2, 2), "N": (-1, -1), "X": (0, 1)}[s]
    return GenotypeMatrix(markers=markers, sample_ids=["RP", "DONOR", "P1"], alleles=alleles, calls=calls)


def make_dataset(gm: GenotypeMatrix, generation: str = "BC2F1") -> Dataset:
    samples = [Sample("RP", "RP", "recurrent_parent"), Sample("DONOR", "DONOR", "donor_parent")]
    samples += [Sample(s, s, "progeny", generation, "F1") for s in gm.sample_ids[2:]]
    return Dataset(genotypes=gm, samples=samples)
