"""samples.csv and markers.csv readers with boundary validation.

Responsibility: read the sample manifest (sample_id, line_name, role,
generation, family_id, notes) and the optional marker map (marker_id, chrom,
pos_bp via position.py, cm), validate them against contract/data-contract.md, and join the map onto
a GenotypeMatrix. Validation errors raise DataContractError; recoverable
issues are returned as warning strings.

Interface:
    read_samples(path) -> list[Sample]
    read_markers(path, scheme=SOYBEAN) -> dict[marker_id, Marker]
    apply_marker_map(gm, marker_map) -> (GenotypeMatrix, warnings)
    build_dataset(gm, samples) -> Dataset
"""

from __future__ import annotations

from pathlib import Path

from progeny_selector.constants import SAMPLE_ROLES
from progeny_selector.core.chrom import SOYBEAN, CompiledScheme, normalize_chrom
from progeny_selector.io.delimited import csv_rows, read_text
from progeny_selector.io.position import parse_position
from progeny_selector.io.wide_csv import add_synthetic_parents
from progeny_selector.model.dataset import DataContractError, Dataset, GenotypeMatrix, Marker, Sample

SAMPLE_REQUIRED = ("sample_id", "role")
MARKER_REQUIRED = ("marker_id", "chrom", "pos_bp")


def _read_rows(path: str | Path, required: tuple[str, ...]) -> list[tuple[int, dict[str, str]]]:
    text = read_text(path)
    all_rows = csv_rows(text, str(path))
    if not all_rows:
        raise DataContractError(f"{path}: empty file")
    names = [n.strip().lower() for n in all_rows[0][1]]
    missing = [c for c in required if c not in names]
    if missing:
        raise DataContractError(f"{path}: missing required columns {missing}")
    rows: list[tuple[int, dict[str, str]]] = []
    for line_num, row in all_rows[1:]:
        clean = {name: (row[i].strip() if i < len(row) else "") for i, name in enumerate(names)}
        rows.append((line_num, clean))
    return rows


def read_samples(path: str | Path) -> list[Sample]:
    rows = _read_rows(path, SAMPLE_REQUIRED)
    samples: list[Sample] = []
    seen: set[str] = set()
    for _line_no, row in rows:
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


def read_markers(path: str | Path, scheme: CompiledScheme = SOYBEAN) -> dict[str, Marker]:
    rows = _read_rows(path, MARKER_REQUIRED)
    out: dict[str, Marker] = {}
    for line_no, row in rows:
        mid = row["marker_id"]
        if not mid:
            raise DataContractError(f"{path}: line {line_no}: empty marker_id")
        if mid in out:
            raise DataContractError(f"{path}: line {line_no}: duplicate marker_id {mid!r}")
        cm_text = row.get("cm", "")
        try:
            cm = float(cm_text) if cm_text not in ("", "NA", "na", ".") else None
        except ValueError as exc:
            raise DataContractError(f"{path}: line {line_no}: invalid cm {cm_text!r}") from exc
        try:
            pos_bp = parse_position(row["pos_bp"])
        except ValueError as exc:
            raise DataContractError(f"{path}: line {line_no}: {exc}") from exc
        out[mid] = Marker(marker_id=mid, chrom=normalize_chrom(row["chrom"], scheme), pos_bp=pos_bp, cm=cm)
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
    """Join manifest and genotypes: synthesise parents for coded input and record them, drop genotype
    columns absent from the manifest, order the columns as the manifest, report mismatches."""
    warnings: list[str] = []
    rp = next(s for s in samples if s.role == "recurrent_parent")
    dp = next(s for s in samples if s.role == "donor_parent")
    synthetic: tuple[str, ...] = ()
    if gm.coded:
        before = set(gm.sample_ids)
        gm = add_synthetic_parents(gm, rp.sample_id, dp.sample_id)
        synthetic = tuple(s for s in (rp.sample_id, dp.sample_id) if s not in before)
        if synthetic:
            warnings.append("A/B/H-coded input: parent columns synthesised (recurrent = A, donor = B)")
    manifest_ids = [s.sample_id for s in samples]
    absent = [s for s in manifest_ids if s not in gm.sample_ids]
    if absent:
        raise DataContractError(f"samples in manifest but not in genotype file: {absent[:10]}{'...' if len(absent) > 10 else ''}")
    extra = [s for s in gm.sample_ids if s not in set(manifest_ids)]
    if extra:
        warnings.append(
            f"{len(extra)} genotype column(s) not in samples.csv dropped: {', '.join(extra[:5])}{', ...' if len(extra) > 5 else ''}"
        )
    gm = gm.select_samples(manifest_ids)
    return Dataset(genotypes=gm, samples=samples, warnings=warnings, synthetic_sample_ids=synthetic)
