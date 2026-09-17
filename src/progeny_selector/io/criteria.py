"""criteria.yaml reader and writer.

Responsibility: map the YAML document described in docs/data-formats.md onto
``progeny_selector.model.criteria`` dataclasses, rejecting unknown keys so
typos fail at the boundary, and serialise a Criteria back to that document in
one canonical form. PyYAML is the only dependency (available in Pyodide, so
the reader also runs under shinylive). Numeric and boolean fields are
type-checked before the dataclasses are built, and malformed YAML is reported
as CriteriaError, so those mistakes never surface as native exceptions.

Interface:
    read_criteria(path) -> Criteria
    read_criteria_text(text: str) -> Criteria
    criteria_from_dict(doc: dict) -> Criteria
    criteria_to_dict(criteria: Criteria) -> dict   canonical key order; criteria_from_dict inverts it
    dump_criteria_yaml(criteria: Criteria) -> str  LF line endings
"""

from __future__ import annotations

import math
import re
from dataclasses import fields
from pathlib import Path

import yaml

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

TOP_KEYS = {"name", "targets", "avoid", "weights", "filters", "background", "flank_window", "flank_unit"}

# Field types checked at the boundary so a wrong-typed value raises CriteriaError here instead of
# a TypeError inside Criteria.validate(). bool is a subclass of int, so it is rejected explicitly.
_NUMBER_KEYS: dict[str, tuple[str, ...]] = {
    "weights": tuple(f.name for f in fields(Weights)),
    "filters": ("max_missing_rate", "max_hom_donor_rate_bcf1", "max_nonparental_rate", "het_rate_tolerance", "parent_max_het_rate"),
}
_BOOL_KEYS: dict[str, tuple[str, ...]] = {"filters": ("exclude_qc_flagged",)}
_REGION_RE = re.compile(r"^\s*([^:]+):\s*([\d,]+)\s*-\s*([\d,]+)\s*$")
_MAX_CRITERIA_BYTES = 1_048_576
_TEXT_LOCUS_KEYS = ("marker_id", "left_marker", "right_marker", "chrom", "locus_id")
_POSITION_KEYS = ("start_bp", "end_bp", "anchor_bp", "min_markers", "min_run")


class _NoAliasLoader(yaml.SafeLoader):
    """SafeLoader that refuses YAML aliases, so a small document cannot expand into a huge one; a bare anchor loads."""

    def compose_node(self, parent, index):
        if self.check_event(yaml.AliasEvent):
            raise CriteriaError("criteria.yaml: YAML aliases are not accepted")
        return super().compose_node(parent, index)


def _check_number(where: str, key: str, value: object, optional: bool = False) -> None:
    if optional and value is None:
        return
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise CriteriaError(f"{where}{key} must be a number, got {value!r}")
    if not math.isfinite(value):
        raise CriteriaError(f"{where}{key} must be a finite number, got {value!r}")


def _check_int(where: str, key: str, value: object) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise CriteriaError(f"{where}{key} must be an integer, got {value!r}")


def _is_whole_float(value: object) -> bool:
    """A finite float with no fractional part, such as 2.0; converted to int in ``_normalise_locus``."""
    return isinstance(value, float) and math.isfinite(value) and value.is_integer()


def _check_bool(where: str, key: str, value: object) -> None:
    if not isinstance(value, bool):
        raise CriteriaError(f"{where}{key} must be true or false, got {value!r}")


def _check_shapes(doc: dict) -> None:
    """Sections must be mappings and locus lists must be lists of mappings; absent or null means default."""
    for section in ("weights", "filters", "background"):
        if doc.get(section) is not None and not isinstance(doc[section], dict):
            raise CriteriaError(f"{section} must be a mapping, got {doc[section]!r}")
    for section in ("targets", "avoid"):
        loci = doc.get(section)
        if loci is None:
            continue
        if not isinstance(loci, list):
            raise CriteriaError(f"{section} must be a list of loci, got {loci!r}")
        for i, locus in enumerate(loci):
            if not isinstance(locus, dict):
                raise CriteriaError(f"{section}[{i}] must be a mapping, got {locus!r}")


def _check_types(doc: dict) -> None:
    """Reject wrong-typed numeric and boolean fields, naming the key and the value. Call after _check_shapes."""
    if "flank_window" in doc:
        _check_number("", "flank_window", doc["flank_window"])
    for section in ("weights", "filters"):
        sub = doc.get(section) or {}
        for key in _NUMBER_KEYS.get(section, ()):
            if key in sub:
                _check_number(f"{section}.", key, sub[key])
        for key in _BOOL_KEYS.get(section, ()):
            if key in sub:
                _check_bool(f"{section}.", key, sub[key])
    background = doc.get("background") or {}
    if "max_marker_coverage" in background:
        _check_number("background.", "max_marker_coverage", background["max_marker_coverage"], optional=True)
    for section in ("targets", "avoid"):
        for i, locus in enumerate(doc.get(section) or []):
            where = f"{section}[{i}]."
            for key in ("start_bp", "end_bp"):
                if key in locus:
                    _check_number(where, key, locus[key], optional=True)
            if "region" in locus and not isinstance(locus["region"], str):
                raise CriteriaError(f"{where}region must be text like 'Gm06:12000000-14000000', got {locus['region']!r}")
            if "min_markers" in locus and not _is_whole_float(locus["min_markers"]):
                _check_int(where, "min_markers", locus["min_markers"])
            for key in ("flank_left", "flank_right"):
                if key in locus:
                    _check_number(where, key, locus[key], optional=True)
            if "allow_het" in locus:
                _check_bool(where, "allow_het", locus["allow_het"])
            if section == "targets":
                _check_run_keys(where, locus)


_RUN_KEYS: tuple[str, ...] = ("min_run", "anchor_bp", "tolerate_isolated")


def _check_run_keys(where: str, locus: dict) -> None:
    """Type-check the rule-run keys of a target and reject them unless the rule is run, so a typo in rule never hides them.

    On avoid loci these keys are left to the unknown-key check in ``_build``.
    """
    for key in ("min_run", "anchor_bp"):
        if _is_whole_float(locus.get(key)):
            continue
        if key in locus:
            _check_int(where, key, locus[key])
    if "tolerate_isolated" in locus:
        _check_bool(where, "tolerate_isolated", locus["tolerate_isolated"])
    present = [key for key in _RUN_KEYS if key in locus]
    if present and locus.get("rule", "all") != "run":
        raise CriteriaError(f"{where[:-1]}: {', '.join(present)} apply only with rule: run, got rule {locus.get('rule', 'all')!r}")


def _build(cls, doc: dict, where: str):
    allowed = {f.name for f in fields(cls)}
    unknown = set(doc) - allowed
    if unknown:
        raise CriteriaError(f"{where}: unknown keys {sorted(unknown)}; allowed {sorted(allowed)}")
    try:
        return cls(**doc)
    except TypeError as exc:
        raise CriteriaError(f"{where}: {exc}") from exc


def criteria_from_dict(doc: dict) -> Criteria:
    if not isinstance(doc, dict):
        raise CriteriaError("criteria document must be a mapping")
    unknown = set(doc) - TOP_KEYS
    if unknown:
        raise CriteriaError(f"unknown top-level keys {sorted(unknown)}; allowed {sorted(TOP_KEYS)}")
    _check_shapes(doc)
    _check_types(doc)
    targets = [_build(TargetSpec, _normalise_locus(t, f"targets[{i}]"), f"targets[{i}]") for i, t in enumerate(doc.get("targets") or [])]
    avoid = [_build(AvoidSpec, _normalise_locus(a, f"avoid[{i}]"), f"avoid[{i}]") for i, a in enumerate(doc.get("avoid") or [])]
    crit = Criteria(
        name=doc.get("name"),
        targets=targets,
        avoid=avoid,
        weights=_build(Weights, doc.get("weights") or {}, "weights"),
        filters=_build(Filters, doc.get("filters") or {}, "filters"),
        background=_build(BackgroundOptions, doc.get("background") or {}, "background"),
        flank_window=float(doc.get("flank_window", 5.0)),
        flank_unit=str(doc.get("flank_unit", "cm")),
    )
    crit.validate()
    return crit


def _normalise_locus(doc: dict, where: str) -> dict:
    """Rename ``id`` and expand the region shorthand; shapes and numeric types are checked beforehand."""
    out = dict(doc)
    if "id" in out and "locus_id" not in out:
        out["locus_id"] = out.pop("id")
    if "region" in out:  # "Gm06:12000000-14000000", thousands separators allowed
        region = out.pop("region")
        match = _REGION_RE.match(region)
        if match is None:
            raise CriteriaError(f"{where}.region must look like 'Gm06:12000000-14000000', got {region!r}")
        chrom, start, end = match.groups()
        out.update(chrom=chrom.strip(), start_bp=int(start.replace(",", "")), end_bp=int(end.replace(",", "")))
    for key in _TEXT_LOCUS_KEYS:
        if key not in out:
            continue
        value = out[key]
        if isinstance(value, bool):  # before int: YAML true is an int in Python
            raise CriteriaError(f"{where}.{key} must be text, got {value!r}")
        if isinstance(value, int | float):
            out[key] = value = str(value)
        if not isinstance(value, str):
            raise CriteriaError(f"{where}.{key} must be text, got {value!r}")
    for key in _POSITION_KEYS:
        value = out.get(key)
        if isinstance(value, float):
            if not math.isfinite(value) or value != int(value):
                raise CriteriaError(f"{where}.{key} must be an integer, got {value!r}")
            out[key] = int(value)
    return out


def read_criteria(path: str | Path) -> Criteria:
    return read_criteria_text(Path(path).read_text(encoding="utf-8"))


def read_criteria_text(text: str) -> Criteria:
    """Parse a criteria.yaml document held in memory.

    Malformed YAML, empty text, documents over 1 MB, YAML aliases and nesting deep enough to exhaust
    the recursion limit raise CriteriaError.
    """
    if len(text.encode("utf-8")) > _MAX_CRITERIA_BYTES:
        raise CriteriaError("criteria.yaml is larger than 1 MB")
    try:
        doc = yaml.load(text, Loader=_NoAliasLoader)
    except yaml.YAMLError as exc:
        raise CriteriaError(f"criteria.yaml is not valid YAML: {exc}") from exc
    except RecursionError as exc:
        raise CriteriaError("criteria.yaml is nested too deeply") from exc
    return criteria_from_dict(doc)


# Keys that define each locus kind (LocusSpec.kind); a dumped locus carries only its own kind's keys.
_KIND_KEYS: dict[str, tuple[str, ...]] = {
    "marker": ("marker_id",),
    "region": ("chrom", "start_bp", "end_bp"),
    "flanking": ("left_marker", "right_marker"),
}


def _locus_to_dict(spec: LocusSpec) -> dict:
    out: dict = {"locus_id": spec.locus_id}
    for key in _KIND_KEYS[spec.kind()]:  # never the "region" shorthand
        out[key] = getattr(spec, key)
    out["rule"] = spec.rule
    out["min_markers"] = spec.min_markers
    if isinstance(spec, TargetSpec) and spec.rule == "run":
        out["min_run"] = spec.min_run
        out["anchor_bp"] = spec.resolved_anchor_bp()
        out["tolerate_isolated"] = spec.tolerate_isolated
    if spec.notes is not None:
        out["notes"] = spec.notes
    if isinstance(spec, TargetSpec):
        out["required_state"] = spec.required_state
        if spec.flank_left is not None:
            out["flank_left"] = spec.flank_left
        if spec.flank_right is not None:
            out["flank_right"] = spec.flank_right
    elif isinstance(spec, AvoidSpec):
        out["allow_het"] = spec.allow_het
    return out


def criteria_to_dict(criteria: Criteria) -> dict:
    """Canonical document: every weight and filter explicit, optional values only when set."""
    doc: dict = {}
    if criteria.name is not None:
        doc["name"] = criteria.name
    doc["targets"] = [_locus_to_dict(t) for t in criteria.targets]
    doc["avoid"] = [_locus_to_dict(a) for a in criteria.avoid]
    doc["flank_window"] = criteria.flank_window
    doc["flank_unit"] = criteria.flank_unit
    background: dict = {"model": criteria.background.model, "map_unit": criteria.background.map_unit}
    if criteria.background.max_marker_coverage is not None:
        background["max_marker_coverage"] = criteria.background.max_marker_coverage
    doc["background"] = background
    doc["weights"] = {f.name: getattr(criteria.weights, f.name) for f in fields(Weights)}
    doc["filters"] = {f.name: getattr(criteria.filters, f.name) for f in fields(Filters)}
    return doc


def dump_criteria_yaml(criteria: Criteria) -> str:
    """criteria.yaml text in the canonical form of ``criteria_to_dict``; PyYAML emits LF line endings."""
    return yaml.safe_dump(criteria_to_dict(criteria), sort_keys=False, allow_unicode=True)
