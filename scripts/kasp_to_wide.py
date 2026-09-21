"""Convert an LGC KASP export to the wide-CSV genotype contract (marker_id,chrom,pos_bp,<samples...>).

Two input shapes are accepted (docs/adr/0019, decision 9; docs/m2-phases.md Phase 10):

- long: one row per (sample, marker) call, detected case-insensitively by the headers
  ``SubjectID``, ``SNPID``, ``Call`` (this shape is unverified against a primary source).
- grid: SNPviewer's export, one row per sample or one row per marker, detected by matching
  the header row or the first column against the marker ids in ``--markers`` in either
  orientation.

A call matches ``^[ACGT][:\\-.]?[ACGT]$`` (case-insensitive; separator ``:``, ``-``, ``.`` or
none) and becomes the nucleotide pair of its first and last character. ``Uncallable``,
``Missing``, ``?``, ``Bad``, ``Dupe``, ``NTC``, empty, and any call not matching the pattern
become ``N`` (missing), counted per reason. Rows/columns whose sample id is in ``--controls``
(default ``NTC``) are dropped. Conflicting non-identical calls for one sample x SNP are an
error naming both source rows. A SNP absent from ``--markers`` is an error listing the ids,
unless ``--drop-unplaced`` (counted). Sample columns keep first-seen order; marker rows are
sorted by chromosome then position.

Usage: python3 scripts/kasp_to_wide.py --kasp export.csv --markers markers.csv --out genotypes.csv
    [--drop-unplaced] [--controls NTC,H2O]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from progeny_selector.core.chrom import chrom_sort_key
from progeny_selector.io.delimited import csv_rows, read_text
from progeny_selector.io.manifest import read_markers
from progeny_selector.model.dataset import DataContractError

CALL_RE = re.compile(r"^[ACGT][:\-.]?[ACGT]$", re.IGNORECASE)
# Missing tokens named by decision 9 / Phase 10 (case-insensitive), plus '?' and empty.
NAMED_MISSING_TOKENS = ("Uncallable", "Missing", "Bad", "Dupe", "NTC")
LONG_HEADERS = ("subjectid", "snpid", "call")


def classify_call(text: str) -> tuple[tuple[str, str] | None, str | None]:
    """(nucleotide pair or None, missing reason or None). Never raises: unrecognised text is missing."""
    t = str(text).strip()
    if t == "":
        return None, "empty"
    if t == "?":
        return None, "?"
    upper = t.upper()
    for token in NAMED_MISSING_TOKENS:
        if upper == token.upper():
            return None, token
    if CALL_RE.match(upper):
        return (upper[0], upper[-1]), None
    return None, "unrecognized"


def _read_grid_or_long(path: Path) -> list[tuple[int, list[str]]]:
    return csv_rows(read_text(path), str(path))


def detect_shape(header: list[str]) -> str:
    lowered = [h.strip().lower() for h in header]
    if all(h in lowered for h in LONG_HEADERS):
        return "long"
    return "grid"


def _unordered(value: tuple[str, str] | None) -> tuple[str, ...] | None:
    """A call as an allele multiset: ``A:G`` and ``G:A`` are one genotype (decision 9 makes order insignificant)."""
    return None if value is None else tuple(sorted(value))


def _record(
    calls: dict[tuple[str, str], tuple[tuple[str, str] | None, int]],
    sample: str,
    marker: str,
    value: tuple[str, str] | None,
    line_no: int,
) -> None:
    key = (sample, marker)
    if key in calls:
        prev_value, prev_line = calls[key]
        if _unordered(prev_value) != _unordered(value):
            raise DataContractError(
                f"conflicting calls for sample {sample!r} x marker {marker!r}: line {prev_line} says "
                f"{prev_value or 'N'!r}, line {line_no} says {value or 'N'!r}"
            )
        return
    calls[key] = (value, line_no)


def parse_long(
    rows: list[tuple[int, list[str]]], controls: set[str]
) -> tuple[dict[tuple[str, str], tuple[tuple[str, str] | None, int]], list[str], dict[str, int]]:
    header = [h.strip().lower() for h in rows[0][1]]
    i_subject = header.index("subjectid")
    i_snp = header.index("snpid")
    i_call = header.index("call")
    calls: dict[tuple[str, str], tuple[tuple[str, str] | None, int]] = {}
    samples: list[str] = []
    seen_samples: set[str] = set()
    reasons: dict[str, int] = {}
    for line_no, row in rows[1:]:
        subject = row[i_subject].strip() if i_subject < len(row) else ""
        snp = row[i_snp].strip() if i_snp < len(row) else ""
        call_text = row[i_call] if i_call < len(row) else ""
        if subject in controls:
            continue
        if not subject or not snp:
            raise DataContractError(f"line {line_no}: empty SubjectID or SNPID")
        if subject not in seen_samples:
            seen_samples.add(subject)
            samples.append(subject)
        value, reason = classify_call(call_text)
        if reason is not None:
            reasons[reason] = reasons.get(reason, 0) + 1
        _record(calls, subject, snp, value, line_no)
    return calls, samples, reasons


def _match_count(ids: list[str], markers: dict) -> int:
    return sum(1 for i in ids if i.strip() in markers)


def parse_grid(
    rows: list[tuple[int, list[str]]], markers: dict, controls: set[str]
) -> tuple[dict[tuple[str, str], tuple[tuple[str, str] | None, int]], list[str], dict[str, int]]:
    header = [h.strip() for h in rows[0][1]]
    header_ids = header[1:]
    first_col_ids = [row[0].strip() for _line_no, row in rows[1:] if row]
    markers_as_columns = _match_count(header_ids, markers)
    markers_as_rows = _match_count(first_col_ids, markers)
    if markers_as_columns == 0 and markers_as_rows == 0:
        raise DataContractError("could not detect KASP grid shape: neither the header row nor the first column matches --markers")
    calls: dict[tuple[str, str], tuple[tuple[str, str] | None, int]] = {}
    samples: list[str] = []
    seen_samples: set[str] = set()
    reasons: dict[str, int] = {}

    if markers_as_columns >= markers_as_rows:
        # header (minus first cell) = markers; first column of each row = sample id
        marker_ids = header_ids
        for line_no, row in rows[1:]:
            sample = row[0].strip() if row else ""
            if not sample or sample in controls:
                continue
            if sample not in seen_samples:
                seen_samples.add(sample)
                samples.append(sample)
            for j, marker_id in enumerate(marker_ids, start=1):
                marker_id = marker_id.strip()
                if not marker_id:
                    continue
                cell = row[j] if j < len(row) else ""
                value, reason = classify_call(cell)
                if reason is not None:
                    reasons[reason] = reasons.get(reason, 0) + 1
                _record(calls, sample, marker_id, value, line_no)
    else:
        # first column of each row = marker id; header (minus first cell) = sample ids
        sample_ids = [s.strip() for s in header_ids]
        for s in sample_ids:
            if s and s not in controls and s not in seen_samples:
                seen_samples.add(s)
                samples.append(s)
        for line_no, row in rows[1:]:
            marker_id = row[0].strip() if row else ""
            if not marker_id:
                continue
            for j, sample in enumerate(sample_ids, start=1):
                if not sample or sample in controls:
                    continue
                cell = row[j] if j < len(row) else ""
                value, reason = classify_call(cell)
                if reason is not None:
                    reasons[reason] = reasons.get(reason, 0) + 1
                _record(calls, sample, marker_id, value, line_no)
    return calls, samples, reasons


def build_wide(
    calls: dict[tuple[str, str], tuple[tuple[str, str] | None, int]],
    samples: list[str],
    markers: dict,
    drop_unplaced: bool,
) -> tuple[list[str], dict[str, dict[str, str]], int]:
    """(ordered marker ids present in the map, {marker_id: {sample: cell text}}, count of dropped unplaced SNPs)."""
    snp_ids = sorted({marker for _sample, marker in calls})
    unplaced = [s for s in snp_ids if s not in markers]
    if unplaced and not drop_unplaced:
        raise DataContractError(f"{len(unplaced)} SNP(s) absent from --markers: {', '.join(unplaced[:20])}")
    placed = [s for s in snp_ids if s in markers]
    placed.sort(key=lambda mid: (chrom_sort_key(markers[mid].chrom), markers[mid].pos_bp, mid))
    grid: dict[str, dict[str, str]] = {mid: {} for mid in placed}
    for (sample, marker), (value, _line_no) in calls.items():
        if marker not in grid:
            continue
        grid[marker][sample] = "N" if value is None else value[0] + value[1]
    return placed, grid, len(unplaced)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--kasp", required=True, type=Path)
    parser.add_argument("--markers", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--drop-unplaced", action="store_true")
    parser.add_argument("--controls", default="NTC")
    args = parser.parse_args(argv)

    controls = {c.strip() for c in args.controls.split(",") if c.strip()}

    try:
        markers = read_markers(args.markers)
        rows = _read_grid_or_long(args.kasp)
        if not rows:
            raise DataContractError(f"{args.kasp}: empty file")
        shape = detect_shape(rows[0][1])
        if shape == "long":
            calls, samples, reasons = parse_long(rows, controls)
        else:
            calls, samples, reasons = parse_grid(rows, markers, controls)
        placed, grid, n_unplaced = build_wide(calls, samples, markers, args.drop_unplaced)
    except DataContractError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        fh.write("marker_id,chrom,pos_bp," + ",".join(samples) + "\n")
        for marker_id in placed:
            m = markers[marker_id]
            row_calls = [grid[marker_id].get(sample, "N") for sample in samples]
            fh.write(f"{marker_id},{m.chrom},{m.pos_bp}," + ",".join(row_calls) + "\n")

    print(f"shape: {shape}; {len(placed)} markers, {len(samples)} samples written to {args.out}")
    if n_unplaced:
        print(f"warning: {n_unplaced} unplaced SNP(s) dropped (--drop-unplaced)", file=sys.stderr)
    for reason, count in sorted(reasons.items()):
        print(f"missing ({reason}): {count}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
