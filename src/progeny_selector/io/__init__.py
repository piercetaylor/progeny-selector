"""Boundary parsers and writers. Format detection by extension lives in ``load_dataset``.

Interface:
    load_genotypes(path, coding='auto') -> GenotypeMatrix
    load_dataset(genotype_path, samples_path, markers_path=None, coding='auto') -> Dataset
"""

from __future__ import annotations

from pathlib import Path

from progeny_selector.io.criteria import read_criteria
from progeny_selector.io.hapmap import read_hapmap
from progeny_selector.io.manifest import apply_marker_map, build_dataset, read_markers, read_samples
from progeny_selector.io.vcf import read_vcf
from progeny_selector.io.wide_csv import read_wide_csv
from progeny_selector.model.dataset import DataContractError, Dataset, GenotypeMatrix


def load_genotypes(path: str | Path, coding: str = "auto") -> GenotypeMatrix:
    name = str(path).lower()
    if name.endswith(".gz"):
        base = name[:-3]
    elif name.endswith(".bgz"):
        base = name[:-4]
    else:
        base = name
    if base.endswith(".vcf"):
        return read_vcf(path)
    if base.endswith((".hmp.txt", ".hmp", ".hapmap")):
        return read_hapmap(path)
    if base.endswith((".csv", ".tsv", ".txt")):
        return read_wide_csv(path, coding=coding)
    raise DataContractError(f"cannot infer genotype format from extension: {path}")


def load_dataset(
    genotype_path: str | Path, samples_path: str | Path, markers_path: str | Path | None = None, coding: str = "auto"
) -> Dataset:
    gm = load_genotypes(genotype_path, coding=coding)
    warnings: list[str] = []
    if markers_path is not None:
        gm, warnings = apply_marker_map(gm, read_markers(markers_path))
    dataset = build_dataset(gm, read_samples(samples_path))
    dataset.warnings = warnings + dataset.warnings
    return dataset


__all__ = ["load_dataset", "load_genotypes", "read_criteria", "read_hapmap", "read_markers", "read_samples", "read_vcf", "read_wide_csv"]
