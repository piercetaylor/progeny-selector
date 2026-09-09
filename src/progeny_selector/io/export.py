"""Writers for results, selection lists and the next-round sample manifest.

Responsibility: serialise analysis rows to CSV with a stable column order,
write the selected-ID list with free-text notes, and emit a samples.csv
skeleton for the next genotyping round (same contract as the input manifest,
so the file round-trips into this tool and the sibling isoline browser).

Interface:
    write_results_csv(rows, path) -> None
    write_selection_csv(rows, path, notes: dict[sample_id, str] | None) -> None
    write_next_round_manifest(rows, path, next_generation: str, rp: Sample, donor: Sample) -> None
"""

from __future__ import annotations

import csv
from pathlib import Path

from progeny_selector.model.dataset import Sample

LEADING_COLUMNS = (
    "rank_overall",
    "rank_in_family",
    "sample_id",
    "line_name",
    "family_id",
    "generation",
    "passes_filters",
    "exclusion_reason",
    "composite_score",
    "foreground_all_pass",
    "avoid_all_pass",
    "rpp_total",
    "rpp_carrier",
    "rpp_noncarrier",
    "expected_rpp",
    "drag_total_est",
    "drag_total_max",
    "drag_unit",
    "ibs_rp",
    "ibs_donor",
    "missing_rate",
    "het_rate",
    "expected_het",
    "qc_flags",
)


def _columns(rows: list[dict]) -> list[str]:
    seen: list[str] = [c for c in LEADING_COLUMNS if rows and c in rows[0]]
    for row in rows:
        for c in row:
            if c not in seen:
                seen.append(c)
    return seen


def _fmt(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "TRUE" if v else "FALSE"
    if isinstance(v, float):
        return f"{v:.6g}" if abs(v) < 1e6 else f"{v:.0f}"
    return v


def write_results_csv(rows: list[dict], path: str | Path) -> None:
    cols = _columns(rows)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(cols)
        for row in rows:
            writer.writerow([_fmt(row.get(c)) for c in cols])


def write_selection_csv(rows: list[dict], path: str | Path, notes: dict[str, str] | None = None) -> None:
    notes = notes or {}
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ["sample_id", "line_name", "family_id", "generation", "rank_overall", "rank_in_family", "composite_score", "rpp_total", "notes"]
        )
        for r in rows:
            writer.writerow(
                [
                    r["sample_id"],
                    r.get("line_name", ""),
                    r.get("family_id") or "",
                    r.get("generation") or "",
                    _fmt(r.get("rank_overall")),
                    _fmt(r.get("rank_in_family")),
                    _fmt(r.get("composite_score")),
                    _fmt(r.get("rpp_total")),
                    notes.get(r["sample_id"], ""),
                ]
            )


def write_next_round_manifest(
    rows: list[dict], path: str | Path, next_generation: str, rp: Sample, donor: Sample, n_per_selected: int = 1
) -> None:
    """samples.csv for the next round: both parents plus placeholder progeny ids derived from each selected line."""
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_id", "line_name", "role", "generation", "family_id", "notes"])
        writer.writerow([rp.sample_id, rp.line_name, "recurrent_parent", "", "", rp.notes or ""])
        writer.writerow([donor.sample_id, donor.line_name, "donor_parent", "", "", donor.notes or ""])
        for r in rows:
            for k in range(1, n_per_selected + 1):
                writer.writerow(
                    [
                        f"{r['sample_id']}-{next_generation}-{k:03d}",
                        r.get("line_name", ""),
                        "progeny",
                        next_generation,
                        r["sample_id"],
                        f"derived from {r['sample_id']}",
                    ]
                )
