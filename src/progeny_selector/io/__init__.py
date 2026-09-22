"""Boundary parsers and writers. Format detection by extension lives in ``load_dataset``.

A token profile (profiles.py, contract 1.4.0) is a built-in id, a dict of the profile JSON, or a
TokenProfile; it applies to HapMap and wide CSV, and a profile with a VCF is an error.
``load_dataset`` records its label ("default", a built-in id or "custom:<id>") in Dataset.token_profile.

A crop (crops.py, contract 1.5.0) is a built-in scheme id; it chooses the chromosome naming scheme
every reader normalises and orders names under. "soybean" is the default and reproduces contract
1.2.0 exactly; ``load_dataset`` resolves it and puts the compiled scheme on Dataset.scheme, from
which Dataset.crop is derived.

Interface:
    load_genotypes(path, coding='auto', profile=None, crop='soybean') -> GenotypeMatrix
    load_dataset(genotype_path, samples_path, markers_path=None, coding='auto', profile=None, crop='soybean') -> Dataset
"""

from __future__ import annotations

from pathlib import Path

from progeny_selector.io.criteria import read_criteria
from progeny_selector.io.crops import DEFAULT_CROP_ID, resolve_crop
from progeny_selector.io.hapmap import read_hapmap
from progeny_selector.io.manifest import apply_marker_map, build_dataset, read_markers, read_samples
from progeny_selector.io.profiles import TokenProfile, profile_label, resolve_profile
from progeny_selector.io.vcf import read_vcf
from progeny_selector.io.wide_csv import read_wide_csv
from progeny_selector.model.dataset import DataContractError, Dataset, GenotypeMatrix


def load_genotypes(
    path: str | Path, coding: str = "auto", profile: str | dict | TokenProfile | None = None, crop: str | None = DEFAULT_CROP_ID
) -> GenotypeMatrix:
    resolved = resolve_profile(profile)
    scheme = resolve_crop(crop)
    name = str(path).lower()
    if name.endswith(".gz"):
        base = name[:-3]
    elif name.endswith(".bgz"):
        base = name[:-4]
    else:
        base = name
    if base.endswith(".vcf"):
        if resolved is not None:
            raise DataContractError(f'token profile "{resolved.id}" applies to HapMap and wide CSV; the genotype file is VCF')
        return read_vcf(path, scheme=scheme)
    if base.endswith((".hmp.txt", ".hmp", ".hapmap")):
        return read_hapmap(path, coding=coding, profile=resolved, scheme=scheme)
    if base.endswith((".csv", ".tsv", ".txt")):
        return read_wide_csv(path, coding=coding, profile=resolved, scheme=scheme)
    raise DataContractError(f"cannot infer genotype format from extension: {path}")


def load_dataset(
    genotype_path: str | Path,
    samples_path: str | Path,
    markers_path: str | Path | None = None,
    coding: str = "auto",
    profile: str | dict | TokenProfile | None = None,
    crop: str | None = DEFAULT_CROP_ID,
) -> Dataset:
    resolved = resolve_profile(profile)
    scheme = resolve_crop(crop)
    gm = load_genotypes(genotype_path, coding=coding, profile=resolved, crop=crop)
    warnings: list[str] = []
    if markers_path is not None:
        gm, warnings = apply_marker_map(gm, read_markers(markers_path, scheme=scheme))
    dataset = build_dataset(gm, read_samples(samples_path))
    dataset.warnings = warnings + dataset.warnings
    dataset.token_profile = profile_label(resolved)
    dataset.scheme = scheme
    return dataset


__all__ = ["load_dataset", "load_genotypes", "read_criteria", "read_hapmap", "read_markers", "read_samples", "read_vcf", "read_wide_csv"]
