"""Selection criteria model (targets, avoid loci, windows, weights, filters).

Responsibility: typed representation of criteria.yaml (docs/data-formats.md,
section "criteria.yaml"). Validation of field values happens here so that
``progeny_selector.io.criteria`` only maps YAML keys onto these classes.

Interface:
    LocusSpec, TargetSpec, AvoidSpec, Weights, Filters, BackgroundOptions, Criteria
    Criteria.validate() -> None (raises CriteriaError)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import ClassVar

REQUIRED_STATES: tuple[str, ...] = ("hom_donor", "het", "either")
LOCUS_RULES: tuple[str, ...] = ("all", "any")
TARGET_RULES: tuple[str, ...] = ("all", "any", "run")
BACKGROUND_MODELS: tuple[str, ...] = ("count", "weighted")
MAP_UNITS: tuple[str, ...] = ("auto", "bp", "cm")
UNKNOWN_POLICIES: tuple[str, ...] = ("fail", "pass")
COMPONENT_NAMES: tuple[str, ...] = (
    "rpp_noncarrier",
    "rpp_carrier",
    "drag",
    "recombinant",
    "similarity_rp",
    "completeness",
)


class CriteriaError(ValueError):
    """Raised when criteria.yaml is inconsistent with the data contract."""


@dataclass
class LocusSpec:
    """A locus is defined by exactly one of: marker_id; chrom+start_bp+end_bp; left_marker+right_marker."""

    locus_id: str
    marker_id: str | None = None
    chrom: str | None = None
    start_bp: int | None = None
    end_bp: int | None = None
    left_marker: str | None = None
    right_marker: str | None = None
    rule: str = "all"
    min_markers: int = 1
    notes: str | None = None

    allowed_rules: ClassVar[tuple[str, ...]] = LOCUS_RULES

    def kind(self) -> str:
        if self.marker_id:
            return "marker"
        if self.left_marker and self.right_marker:
            return "flanking"
        if self.chrom and self.start_bp is not None and self.end_bp is not None:
            return "region"
        raise CriteriaError(f"locus {self.locus_id!r}: give marker_id, or chrom+start_bp+end_bp, or left_marker+right_marker")

    def validate(self) -> None:
        self.kind()
        if self.rule not in self.allowed_rules:
            raise CriteriaError(f"locus {self.locus_id!r}: rule must be one of {self.allowed_rules}")
        if self.min_markers < 1:
            raise CriteriaError(f"locus {self.locus_id!r}: min_markers must be >= 1")
        if self.kind() == "region" and self.start_bp > self.end_bp:  # type: ignore[operator]
            raise CriteriaError(f"locus {self.locus_id!r}: start_bp > end_bp")


@dataclass
class TargetSpec(LocusSpec):
    required_state: str = "either"
    flank_left: float | None = None  # window size in flank_unit; None -> Criteria.flank_window
    flank_right: float | None = None
    min_run: int = 3  # rule run: minimum contiguous predicate calls through the anchor
    anchor_bp: int | None = None  # rule run: None -> (start_bp + end_bp) // 2
    tolerate_isolated: bool = True  # rule run: bridge a single non-predicate call flanked by predicate calls

    allowed_rules: ClassVar[tuple[str, ...]] = TARGET_RULES

    def resolved_anchor_bp(self) -> int:
        """The run anchor: ``anchor_bp`` when set, otherwise the region midpoint rounded down."""
        if self.anchor_bp is not None:
            return int(self.anchor_bp)
        return (int(self.start_bp) + int(self.end_bp)) // 2  # type: ignore[arg-type]

    def validate(self) -> None:
        super().validate()
        if self.required_state not in REQUIRED_STATES:
            raise CriteriaError(f"target {self.locus_id!r}: required_state must be one of {REQUIRED_STATES}")
        if self.min_run < 1:
            raise CriteriaError(f"target {self.locus_id!r}: min_run must be >= 1")
        if self.rule == "run":
            if self.kind() != "region":
                raise CriteriaError(f"target {self.locus_id!r}: rule run needs a region locus")
            if self.anchor_bp is not None and not self.start_bp <= self.anchor_bp <= self.end_bp:  # type: ignore[operator]
                raise CriteriaError(f"target {self.locus_id!r}: anchor_bp {self.anchor_bp} is outside [{self.start_bp}, {self.end_bp}]")
        elif self.min_run != 3 or self.anchor_bp is not None or self.tolerate_isolated is not True:
            raise CriteriaError(f"target {self.locus_id!r}: min_run, anchor_bp and tolerate_isolated apply only with rule run")


@dataclass
class AvoidSpec(LocusSpec):
    allow_het: bool = False

    def validate(self) -> None:
        if self.rule == "run":
            raise CriteriaError(f"avoid locus {self.locus_id!r}: rule run is for targets only")
        super().validate()


@dataclass
class Weights:
    rpp_noncarrier: float = 0.5
    rpp_carrier: float = 0.2
    drag: float = 0.2
    recombinant: float = 0.1
    similarity_rp: float = 0.0
    completeness: float = 0.0

    def as_dict(self) -> dict[str, float]:
        return {name: float(getattr(self, name)) for name in COMPONENT_NAMES}

    def validate(self) -> None:
        values = self.as_dict()
        if any(v < 0 for v in values.values()):
            raise CriteriaError("weights must be >= 0")
        if sum(values.values()) <= 0:
            raise CriteriaError("at least one weight must be > 0")


@dataclass
class Filters:
    max_missing_rate: float = 0.2
    unknown_target_is: str = "fail"
    unknown_avoid_is: str = "pass"
    exclude_qc_flagged: bool = True
    max_hom_donor_rate_bcf1: float = 0.02  # QC: B calls in a BCnF1 above this rate flag a possible self/outcross
    max_nonparental_rate: float = 0.02  # QC: X calls above this rate flag a possible outcross
    het_rate_tolerance: float = 0.15  # QC: |observed - expected| het above this flags the individual
    parent_max_het_rate: float = 0.02

    def validate(self) -> None:
        if not 0 <= self.max_missing_rate <= 1:
            raise CriteriaError("filters.max_missing_rate must be in [0, 1]")
        if self.unknown_target_is not in UNKNOWN_POLICIES or self.unknown_avoid_is not in UNKNOWN_POLICIES:
            raise CriteriaError(f"unknown_*_is must be one of {UNKNOWN_POLICIES}")


@dataclass
class BackgroundOptions:
    model: str = "weighted"
    map_unit: str = "auto"
    max_marker_coverage: float | None = None  # in map_unit; default 10 cM or 4,000,000 bp

    def validate(self) -> None:
        if self.model not in BACKGROUND_MODELS:
            raise CriteriaError(f"background.model must be one of {BACKGROUND_MODELS}")
        if self.map_unit not in MAP_UNITS:
            raise CriteriaError(f"background.map_unit must be one of {MAP_UNITS}")


@dataclass
class Criteria:
    targets: list[TargetSpec] = field(default_factory=list)
    avoid: list[AvoidSpec] = field(default_factory=list)
    weights: Weights = field(default_factory=Weights)
    filters: Filters = field(default_factory=Filters)
    background: BackgroundOptions = field(default_factory=BackgroundOptions)
    flank_window: float = 5.0  # default recombinant window on each side of a target
    flank_unit: str = "cm"  # "cm" or "bp"; falls back to bp when the map has no cM
    name: str | None = None

    def validate(self) -> None:
        if not self.targets:
            raise CriteriaError("criteria must define at least one target locus")
        ids = [t.locus_id for t in self.targets] + [a.locus_id for a in self.avoid]
        if len(ids) != len(set(ids)):
            raise CriteriaError("locus_id values must be unique across targets and avoid loci")
        for t in self.targets:
            t.validate()
        for a in self.avoid:
            a.validate()
        self.weights.validate()
        self.filters.validate()
        self.background.validate()
        if self.flank_unit not in ("cm", "bp"):
            raise CriteriaError("flank_unit must be 'cm' or 'bp'")
        if self.flank_window <= 0:
            raise CriteriaError("flank_window must be > 0")
