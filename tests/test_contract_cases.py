"""The shared data contract (contract/README.md), checked against this repository's loaders.

Every directory under contract/cases/ is loaded through load_dataset (genotypes.<ext>,
samples.csv, optional markers.csv), normalised to the language-neutral shape of
expected.json and compared; an error case must raise DataContractError matching the
pattern its kind maps to below. The kind-to-pattern table lives here, not in the
contract, so rewording a message is a change to this test. The manifest and the
version string are checked too. Mirrors isoline-browser tests/contract-cases.test.ts.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pytest

from progeny_selector.core.chrom import chrom_sort_key
from progeny_selector.io import load_dataset
from progeny_selector.model.dataset import DataContractError, Dataset

CONTRACT = Path(__file__).resolve().parent.parent / "contract"
CASES = CONTRACT / "cases"
VERSION = (CONTRACT / "VERSION").read_bytes().decode("utf-8").strip()

ERROR_KIND_PATTERNS: dict[str, str] = {
    "manifest.roles": r"exactly one (recurrent_parent|donor_parent) required|no progeny or candidate samples",
    "manifest.unknown_role": r"has role .* expected one of",
    "manifest.duplicate_sample": r"duplicate sample_id",
    "dataset.sample_missing": r"samples in manifest but not in genotype file",
    "genotypes.duplicate_marker": r"duplicate marker_id in genotype file",
    "genotypes.no_gt": r"FORMAT has no GT",
    "genotypes.no_header": (
        r"VCF data line before #CHROM header|VCF contains no variant records|"
        r"HapMap header must start with rs#|wide CSV must start with columns"
    ),
    "genotypes.unknown_cell": r"unrecognised (coded|nucleotide) call",
    "genotypes.invalid_position": r"invalid position",
}

CASE_NAMES = sorted(p.name for p in CASES.iterdir() if p.is_dir())


def normalise_dataset(ds: Dataset, version: str) -> dict:
    """expected.json shape: markers in (chromosome order, position), cm None when absent; sampleIds in
    manifest order excluding synthetic parents; calls as sorted allele symbols or None."""
    gm = ds.genotypes.sorted_by_position()
    sample_ids = [s for s in gm.sample_ids if s not in ds.synthetic_sample_ids]
    calls: dict[str, list[list[str] | None]] = {}
    for sid in sample_ids:
        col = gm.sample_index(sid)
        per_marker: list[list[str] | None] = []
        for m in range(gm.n_markers):
            a, b = int(gm.calls[m, col, 0]), int(gm.calls[m, col, 1])
            per_marker.append(None if a < 0 or b < 0 else sorted([gm.alleles[m][a], gm.alleles[m][b]]))
        calls[sid] = per_marker
    return {
        "contractVersion": version,
        "coded": gm.coded,
        "chromosomeOrder": sorted({m.chrom for m in gm.markers}, key=chrom_sort_key),
        "markers": [{"id": m.marker_id, "chrom": m.chrom, "posBp": m.pos_bp, "cm": m.cm} for m in gm.markers],
        "sampleIds": sample_ids,
        "calls": calls,
    }


def _genotype_file(case_dir: Path) -> Path:
    names = [p for p in case_dir.iterdir() if p.name.startswith("genotypes.")]
    assert len(names) == 1, case_dir.name
    return names[0]


def _load(case_dir: Path) -> Dataset:
    markers = case_dir / "markers.csv"
    return load_dataset(_genotype_file(case_dir), case_dir / "samples.csv", markers if markers.exists() else None)


def test_cases_exist_with_one_expectation_each() -> None:
    assert CASE_NAMES
    for name in CASE_NAMES:
        files = {p.name for p in (CASES / name).iterdir()}
        assert len({"expected.json", "expected-error.json"} & files) == 1, name


@pytest.mark.parametrize("name", CASE_NAMES)
def test_case(name: str) -> None:
    case_dir = CASES / name
    expected_path = case_dir / "expected.json"
    if expected_path.exists():
        expected = json.loads(expected_path.read_bytes().decode("utf-8"))
        assert normalise_dataset(_load(case_dir), VERSION) == expected
    else:
        err = json.loads((case_dir / "expected-error.json").read_bytes().decode("utf-8"))
        assert err["contractVersion"] == VERSION
        assert err["kind"] in ERROR_KIND_PATTERNS, f"unknown error kind {err['kind']}"
        with pytest.raises(DataContractError, match=ERROR_KIND_PATTERNS[err["kind"]]):
            _load(case_dir)


def test_manifest_matches_files() -> None:
    paths = sorted(p.relative_to(CONTRACT).as_posix() for p in CONTRACT.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256")
    recomputed = "".join(f"{hashlib.sha256((CONTRACT / p).read_bytes()).hexdigest()}  {p}\n" for p in paths)
    assert (CONTRACT / "MANIFEST.sha256").read_bytes().decode("utf-8") == recomputed


def test_version_in_docs() -> None:
    assert re.fullmatch(r"\d+\.\d+\.\d+", VERSION)
    assert f"Contract version: {VERSION}\n" in (CONTRACT / "data-contract.md").read_bytes().decode("utf-8")
    formats = (CONTRACT.parent / "docs" / "data-formats.md").read_bytes().decode("utf-8")
    assert f"version {VERSION}" in formats
