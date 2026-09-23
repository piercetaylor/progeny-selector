"""In-memory data model for genotypes, markers and the sample manifest.

Responsibility: hold validated, boundary-checked data as plain dataclasses and
numpy arrays. Parsers in ``progeny_selector.io`` construct these; the compute
core in ``progeny_selector.core`` only reads them.

Interface:
    Marker, Sample, CallSetRef, GenotypeMatrix, Dataset (dataclasses)
    GenotypeMatrix.sample_index(sample_id) -> int
    GenotypeMatrix.marker_index(marker_id) -> int
    GenotypeMatrix.sorted_by_position(scheme=SOYBEAN) -> GenotypeMatrix
    Dataset.scheme: CompiledScheme (the crop scheme; Dataset.crop is its id)
    GenotypeMatrix.positions(unit) -> np.ndarray
    GenotypeMatrix.select_samples(sample_ids) -> GenotypeMatrix
    Dataset.recurrent_parent / donor_parent / progeny -> Sample(s)
    Dataset.synthetic_sample_ids
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from progeny_selector.core.chrom import SOYBEAN, CompiledScheme, chrom_sort_key


class DataContractError(ValueError):
    """Raised when an input file violates contract/data-contract.md or docs/data-formats.md."""


@dataclass(frozen=True)
class Marker:
    marker_id: str
    chrom: str
    pos_bp: int
    cm: float | None = None


@dataclass(frozen=True)
class Sample:
    sample_id: str
    line_name: str
    role: str
    generation: str | None = None
    family_id: str | None = None
    notes: str | None = None


@dataclass(frozen=True)
class CallSetRef:
    """One BrAPI call set and the ``sample_id`` it was loaded under (docs/adr/0024)."""

    sample_id: str
    call_set_name: str
    call_set_db_id: str
    sample_db_id: str


@dataclass
class GenotypeMatrix:
    """Diploid calls as allele indices.

    calls: int8 array of shape (n_markers, n_samples, 2); -1 marks a missing allele.
    alleles: per-marker list of allele strings; ``calls[m, s, k]`` indexes ``alleles[m]``.
    coded: True when the input used A/B/H coding (alleles are exactly ["A", "B"]).
    """

    markers: list[Marker]
    sample_ids: list[str]
    alleles: list[list[str]]
    calls: np.ndarray
    coded: bool = False
    _sample_pos: dict[str, int] = field(default_factory=dict, repr=False)
    _marker_pos: dict[str, int] = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        n_m, n_s = len(self.markers), len(self.sample_ids)
        if self.calls.shape != (n_m, n_s, 2):
            raise DataContractError(f"calls shape {self.calls.shape} != ({n_m}, {n_s}, 2)")
        if len(self.alleles) != n_m:
            raise DataContractError("alleles list length differs from marker count")
        if len(set(self.sample_ids)) != n_s:
            raise DataContractError("duplicate sample_id in genotype file")
        ids = [m.marker_id for m in self.markers]
        if len(set(ids)) != n_m:
            raise DataContractError("duplicate marker_id in genotype file")
        self._sample_pos = {s: i for i, s in enumerate(self.sample_ids)}
        self._marker_pos = {m: i for i, m in enumerate(ids)}

    @property
    def n_markers(self) -> int:
        return len(self.markers)

    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)

    def sample_index(self, sample_id: str) -> int:
        try:
            return self._sample_pos[sample_id]
        except KeyError as exc:
            raise DataContractError(f"sample_id {sample_id!r} not in genotype file") from exc

    def marker_index(self, marker_id: str) -> int:
        try:
            return self._marker_pos[marker_id]
        except KeyError as exc:
            raise DataContractError(f"marker_id {marker_id!r} not in genotype file") from exc

    def chroms(self) -> np.ndarray:
        return np.array([m.chrom for m in self.markers], dtype=object)

    def positions(self, unit: str = "bp") -> np.ndarray:
        """Marker positions in 'bp' or 'cm' (NaN where cM is absent)."""
        if unit == "bp":
            return np.array([m.pos_bp for m in self.markers], dtype=float)
        if unit == "cm":
            return np.array([np.nan if m.cm is None else m.cm for m in self.markers], dtype=float)
        raise ValueError(f"unknown unit {unit!r}")

    def has_cm(self) -> bool:
        return all(m.cm is not None for m in self.markers)

    def sorted_by_position(self, scheme: CompiledScheme = SOYBEAN) -> GenotypeMatrix:
        order = sorted(range(self.n_markers), key=lambda i: (chrom_sort_key(self.markers[i].chrom, scheme), self.markers[i].pos_bp))
        if order == list(range(self.n_markers)):
            return self
        return GenotypeMatrix(
            markers=[self.markers[i] for i in order],
            sample_ids=list(self.sample_ids),
            alleles=[self.alleles[i] for i in order],
            calls=self.calls[order],
            coded=self.coded,
        )

    def with_samples_added(self, new_ids: list[str], new_calls: np.ndarray) -> GenotypeMatrix:
        return GenotypeMatrix(
            markers=self.markers,
            sample_ids=list(self.sample_ids) + list(new_ids),
            alleles=self.alleles,
            calls=np.concatenate([self.calls, new_calls], axis=1),
            coded=self.coded,
        )

    def select_samples(self, sample_ids: list[str]) -> GenotypeMatrix:
        """Columns for ``sample_ids`` in that order; every id must be present."""
        idx = [self.sample_index(s) for s in sample_ids]
        if idx == list(range(self.n_samples)):
            return self
        return GenotypeMatrix(
            markers=self.markers,
            sample_ids=list(sample_ids),
            alleles=self.alleles,
            calls=self.calls[:, idx],
            coded=self.coded,
        )


@dataclass
class Dataset:
    """A genotype matrix joined with its sample manifest.

    synthetic_sample_ids: parents of a coded file that had no genotype column and were synthesised by
    build_dataset; they are not part of the loaded sample list the contract describes.
    token_profile: the token profile the genotype file was read with (contract 1.4.0): "default",
    a built-in id, or "custom:<id>".
    scheme: the compiled crop chromosome scheme the names were read under (contract 1.5.0),
    which every later normalisation, ordering and chromosome-length lookup uses; ``SOYBEAN`` is
    the default and reproduces contract 1.2.0. ``crop`` is its id, derived so the two cannot
    disagree.
    call_sets: every call set of the BrAPI variant set the dataset came from (docs/adr/0024), in
    server order; empty for a file. It is not aligned with ``genotypes.sample_ids``, which
    ``build_dataset`` has already reordered to the manifest and stripped of unlisted call sets;
    join the two on ``CallSetRef.sample_id``.
    """

    genotypes: GenotypeMatrix
    samples: list[Sample]
    warnings: list[str] = field(default_factory=list)
    synthetic_sample_ids: tuple[str, ...] = ()
    token_profile: str = "default"
    scheme: CompiledScheme = SOYBEAN
    call_sets: tuple[CallSetRef, ...] = ()

    @property
    def crop(self) -> str:
        """The id of the crop chromosome scheme (contract 1.5.0)."""
        return self.scheme.id

    def __post_init__(self) -> None:
        roles = [s.role for s in self.samples]
        if roles.count("recurrent_parent") != 1 or roles.count("donor_parent") != 1:
            raise DataContractError("samples.csv must contain exactly one recurrent_parent and one donor_parent")
        for s in self.samples:
            self.genotypes.sample_index(s.sample_id)

    @property
    def recurrent_parent(self) -> Sample:
        return next(s for s in self.samples if s.role == "recurrent_parent")

    @property
    def donor_parent(self) -> Sample:
        return next(s for s in self.samples if s.role == "donor_parent")

    @property
    def progeny(self) -> list[Sample]:
        return [s for s in self.samples if s.role in ("progeny", "candidate")]

    def sample(self, sample_id: str) -> Sample:
        return next(s for s in self.samples if s.sample_id == sample_id)
