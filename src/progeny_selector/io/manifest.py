"""samples.csv and markers.csv readers with boundary validation.

Responsibility: read the sample manifest (sample_id, line_name, role,
generation, family_id, notes) and the optional marker map (marker_id, chrom,
pos_bp, cm), validate them against docs/data-formats.md, and join the map onto
a GenotypeMatrix. Validation errors raise DataContractError; recoverable
issues are returned as warning strings.

Interface:
    read_samples(path) -> list[Sample]
    read_markers(path) -> dict[marker_id, Marker]
    apply_marker_map(gm, marker_map) -> (GenotypeMatrix, warnings)
    build_dataset(gm, samples) -> Dataset
"""

from __future__ import annotations

import csv
from pathlib import Path

from progeny_selector.constants import SAMPLE_ROLES
from progeny_selector.core.chrom import normalize_chrom
from progeny_selector.io.wide_csv import add_synthetic_parents
from progeny_selector.model.dataset import DataContractError, Dataset, GenotypeMatrix, Marker, Sample

SAMPLE_REQUIRED = ("sample_id", "line_name", "role")
MARKER_REQUIRED = ("marker_id", "chrom", "pos_bp")


def _read_rows(path: str | Path, required: tuple[str, ...]) -> list[dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise DataContractError(f"{path}: empty file")
        names = [n.strip().lower() for n in reader.fieldnames]
        missing = [c for c in required if c not in names]
        if missing:
            raise DataContractError(f"{path}: missing required columns {missing}")
        rows = []
        for row in reader:
            clean = {k.strip().lower(): (v or "").strip() for k, v in row.items() if k is not None}
            if any(clean.values()):
                rows.append(clean)
    return rows


def read_samples(path: str | Path) -> list[Sample]:
    rows = _read_rows(path, SAMPLE_REQUIRED)
    samples: list[Sample] = []
    seen: set[str] = set()
    for row in rows:
        sid = row["sample_id"]
        if not sid:
            raise DataContractError(f"{path}: empty sample_id")
        if sid in seen:
            raise DataContractError(f"{path}: duplicate sample_id {sid!r}")
        seen.add(sid)
        role = row["role"].lower()
        if role not in SAMPLE_ROLES:
            raise DataContractError(f"{path}: sample {sid!r} has role {row['role']!r}; expected one of {SAMPLE_ROLES}")
        generation = row.get("generation") or None  # kept verbatim; QC flags 'generation_unparsed' when unreadable
        samples.append(
            Sample(
                sample_id=sid,
                line_name=row.get("line_name") or sid,
                role=role,
                generation=generation,
                family_id=row.get("family_id") or None,
                notes=row.get("notes") or None,
            )
        )
    roles = [s.role for s in samples]
    if roles.count("recurrent_parent") != 1:
        raise DataContractError(f"{path}: exactly one recurrent_parent required, found {roles.count('recurrent_parent')}")
    if roles.count("donor_parent") != 1:
        raise DataContractError(f"{path}: exactly one donor_parent required, found {roles.count('donor_parent')}")
    if not any(r in ("progeny", "candidate") for r in roles):
        raise DataContractError(f"{path}: no progeny or candidate samples")
    return samples


def read_markers(path: str | Path) -> dict[str, Marker]:
    rows = _read_rows(path, MARKER_REQUIRED)
    out: dict[str, Marker] = {}
    for row in rows:
        mid = row["marker_id"]
        if mid in out:
            raise DataContractError(f"{path}: duplicate marker_id {mid!r}")
        cm_text = row.get("cm", "")
        cm = float(cm_text) if cm_text not in ("", "NA", "na", ".") else None
        out[mid] = Marker(marker_id=mid, chrom=normalize_chrom(row["chrom"]), pos_bp=int(float(row["pos_bp"])), cm=cm)
    if not out:
        raise DataContractError(f"{path}: no marker rows")
    return out


def apply_marker_map(gm: GenotypeMatrix, marker_map: dict[str, Marker]) -> tuple[GenotypeMatrix, list[str]]:
    """Overlay chrom/pos/cm from markers.csv; genotype-file positions win only when the map lacks the marker."""
    warnings: list[str] = []
    new_markers: list[Marker] = []
    n_missing = 0
    for m in gm.markers:
        mm = marker_map.get(m.marker_id)
        if mm is None:
            n_missing += 1
            new_markers.append(m)
            continue
        if mm.chrom != m.chrom or mm.pos_bp != m.pos_bp:
            warnings.append(f"marker {m.marker_id}: map position {mm.chrom}:{mm.pos_bp} overrides genotype file {m.chrom}:{m.pos_bp}")
        new_markers.append(mm)
    if n_missing:
        warnings.append(f"{n_missing} genotype markers absent from markers.csv keep their genotype-file positions and have no cM")
    gm2 = GenotypeMatrix(markers=new_markers, sample_ids=gm.sample_ids, alleles=gm.alleles, calls=gm.calls, coded=gm.coded)
    return gm2, warnings[:50]


def build_dataset(gm: GenotypeMatrix, samples: list[Sample]) -> Dataset:
    """Join manifest and genotypes; synthesise parents for coded input; report manifest/genotype mismatches."""
    warnings: list[str] = []
    rp = next(s for s in samples if s.role == "recurrent_parent")
    dp = next(s for s in samples if s.role == "donor_parent")
    if gm.coded:
        before = gm.n_samples
        gm = add_synthetic_parents(gm, rp.sample_id, dp.sample_id)
        if gm.n_samples != before:
            warnings.append("A/B/H-coded input: parent columns synthesised (recurrent = A, donor = B)")
    manifest_ids = {s.sample_id for s in samples}
    absent = [s.sample_id for s in samples if s.sample_id not in gm.sample_ids]
    if absent:
        raise DataContractError(f"samples in manifest but not in genotype file: {absent[:10]}{'...' if len(absent) > 10 else ''}")
    extra = [s for s in gm.sample_ids if s not in manifest_ids]
    if extra:
        warnings.append(f"{len(extra)} genotype columns are not in samples.csv and are ignored")
    return Dataset(genotypes=gm, samples=samples, warnings=warnings)
