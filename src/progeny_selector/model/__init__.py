"""Typed data model shared by parsers, compute core and UI."""

from progeny_selector.model.criteria import (
    AvoidSpec,
    BackgroundOptions,
    Criteria,
    CriteriaError,
    Filters,
    LocusSpec,
    TargetSpec,
    Weights,
)
from progeny_selector.model.dataset import DataContractError, Dataset, GenotypeMatrix, Marker, Sample

__all__ = [
    "AvoidSpec",
    "BackgroundOptions",
    "Criteria",
    "CriteriaError",
    "DataContractError",
    "Dataset",
    "Filters",
    "GenotypeMatrix",
    "LocusSpec",
    "Marker",
    "Sample",
    "TargetSpec",
    "Weights",
]
