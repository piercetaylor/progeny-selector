"""Convert Song et al. 2016 Table S1 (soybean genetic maps) to markers.csv, with Marey-map interpolation.

Source: Song et al. 2016, BMC Genomics 17:33, Table S1, CC BY 4.0
(https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12864-015-2344-0/MediaObjects/12864_2015_2344_MOESM1_ESM.xls),
exported by hand to CSV. The export has a title row, then the header on row 2 with the columns
`ss ID`, `SNP ID`, `Glyma1.01 Chromosome`, `Glyma1.01 Coordinate`, `Wm82.a2.v1 Chromosome`,
`Wm82.a2.v1 Coordinate`, `WP Linkage Group`, `WP linkage position`, `EW Linkage Group `
(trailing space) and `EW linkage position`, among others. Header names are matched exactly, then
with surrounding spaces removed. The WP and EW maps are never averaged (docs/adr/0019).

Steps, each counted in the `.log` written beside `--out`:
  1. rows read;
  2. no position on the chosen assembly: empty or unparseable coordinate, or a chromosome cell
     (whitespace stripped) that is not Gm01..Gm20 after core.chrom.normalize_chrom (scaffolds);
  3. no linkage group or linkage position on the chosen map;
  4. LG mismatch: the linkage group, a whole number written `1`, `1.0`, is not the chromosome number;
  5. non-monotone dropped: markers outside the longest non-decreasing subsequence of cM along bp,
     per chromosome, candidates sorted by (bp, cM) (core.genetic_map.longest_nondecreasing_mask);
  6. mapped: the markers kept. None kept on any chromosome stops the script.
With `--markers-in` (a markers.csv-shaped file on the same assembly), each marker is matched to the
table by `ss ID`, else by `SNP ID`; at least one must match, and a matched marker whose position
differs from the table's means the assemblies are mixed, so the script refuses. A marker with a cM
value keeps it (cM kept from --markers-in); a marker matched to a mapped row takes that cM; every
other marker gets cM by core.genetic_map.interpolate_cm over its chromosome's mapped markers
(interpolated; clamped, a subset of interpolated, when outside the mapped bp range); markers on a
chromosome with no mapped marker get an empty cm (unplaced), which disables cM mode for a dataset
that uses them.

Usage: python scripts/song2016_map.py --table data/song2016/table_s1.csv [--map WP|EW]
    [--assembly a2|a1] --out markers.csv [--markers-in markers_in.csv]
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from progeny_selector.constants import SOYBEAN_CHROMOSOMES
from progeny_selector.core.chrom import chrom_sort_key, normalize_chrom
from progeny_selector.core.genetic_map import interpolate_cm, longest_nondecreasing_mask
from progeny_selector.io.delimited import csv_rows, read_text
from progeny_selector.io.manifest import read_markers
from progeny_selector.io.position import parse_position
from progeny_selector.model.dataset import DataContractError

CITATION = "Song et al. 2016, BMC Genomics 17:33, Table S1, CC BY 4.0"
MARKER_COLUMN = "ss ID"
SNP_COLUMN = "SNP ID"
ASSEMBLY_COLUMNS = {
    "a2": ("Wm82.a2.v1 Chromosome", "Wm82.a2.v1 Coordinate"),
    "a1": ("Glyma1.01 Chromosome", "Glyma1.01 Coordinate"),
}
MAP_COLUMNS = {
    "WP": ("WP Linkage Group", "WP linkage position"),
    "EW": ("EW Linkage Group ", "EW linkage position"),
}
DECIMAL_RE = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$")
DROP_KEYS = ("no position on assembly", "no linkage position on map", "LG mismatch", "non-monotone dropped")
LOG_KEYS = (
    "rows read",
    *DROP_KEYS,
    "mapped",
    "cM kept from --markers-in",
    "interpolated",
    "clamped",
    "unplaced",
)

Row = tuple[str, str, int, float]


def parse_decimal(text: str) -> float | None:
    """A finite decimal number written plainly (no `_`, no `nan`/`inf`), else None."""
    t = text.strip()
    if not DECIMAL_RE.match(t):
        return None
    value = float(t)
    return value if np.isfinite(value) else None


def lg_to_chrom(lg: str) -> str | None:
    """Chromosome named by a linkage-group cell holding a whole number 1..20 (`1`, `1.0`), else None."""
    value = parse_decimal(lg)
    if value is None or value != int(value) or not 1 <= int(value) <= len(SOYBEAN_CHROMOSOMES):
        return None
    return SOYBEAN_CHROMOSOMES[int(value) - 1]


def _find_column(header: list[str], name: str) -> int:
    if name in header:
        return header.index(name)
    stripped = [h.strip() for h in header]
    if name.strip() in stripped:
        return stripped.index(name.strip())
    raise DataContractError(f"Table S1: column {name!r} not found in the header on row 2")


def read_table(
    path: Path, assembly: str, map_name: str
) -> tuple[list[Row], dict[str, int], dict[str, tuple[str, int]], dict[str, str], set[str]]:
    """Rows kept after steps 1-4 in table order, the counts, every placed row's position, SNP ID -> ss id, all ss ids."""
    rows = csv_rows(read_text(path), str(path))
    if len(rows) < 2:
        raise DataContractError(f"{path}: expected a title row and a header row")
    header = rows[1][1]
    i_id = _find_column(header, MARKER_COLUMN)
    i_snp = _find_column(header, SNP_COLUMN)
    i_chrom, i_pos = (_find_column(header, c) for c in ASSEMBLY_COLUMNS[assembly])
    i_lg, i_cm = (_find_column(header, c) for c in MAP_COLUMNS[map_name])
    counts = dict.fromkeys(LOG_KEYS, 0)
    kept: list[Row] = []
    seen: set[str] = set()
    placed: dict[str, tuple[str, int]] = {}
    by_snp: dict[str, str] = {}

    def cell(row: list[str], i: int) -> str:
        return row[i].strip() if i < len(row) else ""

    for line_no, row in rows[2:]:
        counts["rows read"] += 1
        raw_id = cell(row, i_id)
        try:
            marker_id = "ss" + str(int(float(raw_id)))
        except ValueError as exc:
            raise DataContractError(f"{path}: line {line_no}: invalid ss ID {raw_id!r}") from exc
        if marker_id in seen:
            raise DataContractError(f"{path}: line {line_no}: duplicate marker {marker_id}")
        seen.add(marker_id)
        snp = cell(row, i_snp)
        if snp:
            by_snp.setdefault(snp, marker_id)
        chrom = normalize_chrom(cell(row, i_chrom))
        try:
            pos_bp = parse_position(cell(row, i_pos))
        except ValueError:
            pos_bp = None
        if pos_bp is None or chrom not in SOYBEAN_CHROMOSOMES:
            counts["no position on assembly"] += 1
            continue
        placed[marker_id] = (chrom, pos_bp)
        lg_text = cell(row, i_lg)
        cm = parse_decimal(cell(row, i_cm))
        if not lg_text or cm is None:
            counts["no linkage position on map"] += 1
            continue
        if lg_to_chrom(lg_text) != chrom:
            counts["LG mismatch"] += 1
            continue
        kept.append((marker_id, chrom, pos_bp, cm))
    return kept, counts, placed, by_snp, seen


def clean_monotone(rows: list[Row], counts: dict[str, int]) -> list[Row]:
    """Per chromosome, keep the longest non-decreasing subsequence of cM along bp, candidates sorted by (bp, cM)."""
    by_chrom: dict[str, list[Row]] = defaultdict(list)
    for r in rows:
        by_chrom[r[1]].append(r)
    out: list[Row] = []
    for chrom_rows in by_chrom.values():
        chrom_rows.sort(key=lambda r: (r[2], r[3]))
        mask = longest_nondecreasing_mask(np.array([r[3] for r in chrom_rows], dtype=float))
        counts["non-monotone dropped"] += int((~mask).sum())
        kept = [r for r, k in zip(chrom_rows, mask, strict=True) if k]
        # markers sharing a bp are written at their mean cM, the anchor interpolate_cm uses
        by_bp: dict[int, list[float]] = defaultdict(list)
        for r in kept:
            by_bp[r[2]].append(r[3])
        out.extend((r[0], r[1], r[2], float(np.mean(by_bp[r[2]]))) for r in kept)
    counts["mapped"] = len(out)
    return out


def check_strict_cm(path: Path) -> None:
    """Reject a --markers-in cm cell that is not a plain decimal (read_markers alone accepts `1_000`)."""
    rows = csv_rows(read_text(path), str(path))
    if not rows:
        return
    names = [n.strip().lower() for n in rows[0][1]]
    if "cm" not in names:
        return
    i = names.index("cm")
    for line_no, row in rows[1:]:
        text = row[i].strip() if i < len(row) else ""
        if text not in ("", "NA", "na", ".") and parse_decimal(text) is None:
            raise DataContractError(f"{path}: line {line_no}: invalid cm {text!r}")


def _fmt_cm(cm: float | None) -> str:
    return "" if cm is None else repr(round(float(cm), 6))


def _write_log(args: argparse.Namespace, counts: dict[str, int], warnings: list[str]) -> Path:
    log_path = args.out.with_suffix(".log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", newline="", encoding="utf-8") as fh:
        fh.write(f"source: {CITATION}\n")
        fh.write(f"table: {args.table}\nmap: {args.map_name}\nassembly: {args.assembly}\n")
        for key in LOG_KEYS:
            fh.write(f"{key}: {counts[key]}\n")
        for w in warnings:
            fh.write(f"warning: {w}\n")
    return log_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--map", dest="map_name", choices=sorted(MAP_COLUMNS), default="WP")
    parser.add_argument("--assembly", choices=sorted(ASSEMBLY_COLUMNS), default="a2")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--markers-in", type=Path, default=None)
    args = parser.parse_args(argv)

    try:
        rows, counts, placed, by_snp, table_ids = read_table(args.table, args.assembly, args.map_name)
        mapped = clean_monotone(rows, counts)
        markers_in = {}
        if args.markers_in is not None:
            check_strict_cm(args.markers_in)
            markers_in = read_markers(args.markers_in)
    except DataContractError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if not mapped:
        _write_log(args, counts, [])
        worst = max(DROP_KEYS, key=lambda k: counts[k])
        print(
            f"error: no marker kept on any chromosome from {counts['rows read']} rows (largest drop: {worst} {counts[worst]}); "
            f"check --map and --assembly; {args.out} not written",
            file=sys.stderr,
        )
        return 2

    # markers-in id -> table ss id, by ss ID then SNP ID
    matched: dict[str, str] = {}
    for mid in markers_in:
        if mid in table_ids:
            matched[mid] = mid
        elif mid in by_snp:
            matched[mid] = by_snp[mid]
    if markers_in and not matched:
        print("error: no --markers-in marker matches Table S1 by ss ID or SNP ID; cannot check the assembly", file=sys.stderr)
        return 2
    mixed = [mid for mid, ss in matched.items() if ss in placed and (markers_in[mid].chrom, markers_in[mid].pos_bp) != placed[ss]]
    if mixed:
        print(
            f"error: {len(mixed)} markers in --markers-in have positions that differ from Table S1 on assembly {args.assembly} "
            f"(first: {', '.join(mixed[:5])}); positions must be on one assembly",
            file=sys.stderr,
        )
        return 2

    mapped_by_id = {r[0]: r for r in mapped}
    out: dict[str, tuple[str, str, int, float | None]] = {r[0]: r for r in mapped}
    warnings: list[str] = []
    by_chrom: dict[str, list[Row]] = defaultdict(list)
    for r in mapped:
        by_chrom[r[1]].append(r)
    queries: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for mid, m in markers_in.items():
        ss = matched.get(mid)
        if ss is not None and ss != mid and ss not in markers_in:
            out.pop(ss, None)  # matched by SNP ID: one row, under the --markers-in id
        if m.cm is not None:
            out[mid] = (mid, m.chrom, m.pos_bp, m.cm)
            counts["cM kept from --markers-in"] += 1
        elif mid in matched and matched[mid] in mapped_by_id:
            out[mid] = (mid, m.chrom, m.pos_bp, mapped_by_id[matched[mid]][3])
        else:
            queries[m.chrom].append((mid, m.pos_bp))
    for chrom in sorted(queries, key=chrom_sort_key):
        q = queries[chrom]
        ref = by_chrom.get(chrom, [])
        if not ref:
            counts["unplaced"] += len(q)
            warnings.append(f"chromosome {chrom}: no mapped marker; {len(q)} markers have no cM, so cM mode will be disabled")
            out.update((mid, (mid, chrom, pos, None)) for mid, pos in q)
            continue
        ref_bp = np.array([r[2] for r in ref], dtype=float)
        q_bp = np.array([pos for _, pos in q], dtype=float)
        cms = interpolate_cm(ref_bp, np.array([r[3] for r in ref], dtype=float), q_bp)
        counts["interpolated"] += len(q)
        counts["clamped"] += int(((q_bp < ref_bp.min()) | (q_bp > ref_bp.max())).sum())
        out.update((mid, (mid, chrom, pos, float(c))) for (mid, pos), c in zip(q, cms, strict=True))

    out_rows = sorted(out.values(), key=lambda r: (chrom_sort_key(r[1]), r[2], r[3] is None, r[3] or 0.0, r[0]))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        fh.write("marker_id,chrom,pos_bp,cm\n")
        for mid, chrom, pos, cm in out_rows:
            fh.write(f"{mid},{chrom},{pos},{_fmt_cm(cm)}\n")

    log_path = _write_log(args, counts, warnings)
    for w in warnings:
        print(f"warning: {w}", file=sys.stderr)
    print(f"{len(out_rows)} markers written to {args.out}; counts in {log_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
