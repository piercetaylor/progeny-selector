"""criteria.yaml reader.

Responsibility: map the YAML document described in docs/data-formats.md onto
``progeny_selector.model.criteria`` dataclasses, rejecting unknown keys so
typos fail at the boundary. PyYAML is the only dependency (available in
Pyodide, so the reader also runs under shinylive).

Interface:
    read_criteria(path) -> Criteria
    criteria_from_dict(doc: dict) -> Criteria
"""

from __future__ import annotations

from dataclasses import fields
from pathlib import Path

import yaml

from progeny_selector.model.criteria import AvoidSpec, BackgroundOptions, Criteria, CriteriaError, Filters, TargetSpec, Weights

TOP_KEYS = {"name", "targets", "avoid", "weights", "filters", "background", "flank_window", "flank_unit"}


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
    targets = [_build(TargetSpec, _normalise_locus(t), f"targets[{i}]") for i, t in enumerate(doc.get("targets") or [])]
    avoid = [_build(AvoidSpec, _normalise_locus(a), f"avoid[{i}]") for i, a in enumerate(doc.get("avoid") or [])]
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


def _normalise_locus(doc: dict) -> dict:
    if not isinstance(doc, dict):
        raise CriteriaError("each locus must be a mapping")
    out = dict(doc)
    if "id" in out and "locus_id" not in out:
        out["locus_id"] = out.pop("id")
    if "region" in out:  # "Gm06:12000000-14000000"
        region = str(out.pop("region"))
        chrom, span = region.split(":")
        start, end = span.replace(",", "").split("-")
        out.update(chrom=chrom, start_bp=int(start), end_bp=int(end))
    for key in ("start_bp", "end_bp", "min_markers"):
        if key in out and out[key] is not None:
            out[key] = int(out[key])
    return out


def read_criteria(path: str | Path) -> Criteria:
    with open(path, encoding="utf-8") as fh:
        doc = yaml.safe_load(fh)
    return criteria_from_dict(doc)
