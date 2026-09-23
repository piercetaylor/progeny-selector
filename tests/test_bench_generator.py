"""Tests for scripts/bench_pipeline.py::generate_case at fixture scale (cheap; no browser)."""

from __future__ import annotations

import filecmp
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from progeny_selector.core.pipeline import run_analysis
from progeny_selector.io import load_dataset, read_criteria

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "bench_pipeline.py"


def _load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location("bench_pipeline", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


bench_pipeline = _load_module()


def test_generate_case_writes_five_files(tmp_path: Path) -> None:
    out = tmp_path / "case"
    bench_pipeline.generate_case(out, 40, 6)
    for name in ("genotypes.vcf", "genotypes.vcf.gz", "samples.csv", "markers.csv", "criteria.yaml"):
        assert (out / name).exists(), name


def test_generated_case_loads_and_analyses(tmp_path: Path) -> None:
    out = tmp_path / "case"
    bench_pipeline.generate_case(out, 40, 6)
    dataset = load_dataset(out / "genotypes.vcf", out / "samples.csv", out / "markers.csv")
    assert dataset.genotypes.n_markers == 40
    chroms = {m.chrom for m in dataset.genotypes.markers}
    assert chroms == {f"Gm{i:02d}" for i in range(1, 21)}
    assert len(dataset.progeny) == 6

    criteria = read_criteria(out / "criteria.yaml")
    result = run_analysis(dataset, criteria)
    assert all(row["background_model"] == "weighted" for row in result.rows)


def test_generate_case_is_deterministic(tmp_path: Path) -> None:
    out1, out2 = tmp_path / "one", tmp_path / "two"
    bench_pipeline.generate_case(out1, 40, 6)
    bench_pipeline.generate_case(out2, 40, 6)
    for name in ("genotypes.vcf", "samples.csv", "markers.csv", "criteria.yaml"):
        assert filecmp.cmp(out1 / name, out2 / name, shallow=False), name
