"""Writers for results, selection lists and the next-round sample manifest.

Responsibility: serialise analysis rows to CSV with a stable column order,
write the selected-ID list with free-text notes, and emit a samples.csv
skeleton for the next genotyping round (same contract as the input manifest,
so the file round-trips into this tool and the sibling Backcross tool).

Interface:
    write_results_csv(rows, path) -> None
    write_selection_csv(rows, path, notes: dict[sample_id, str] | None) -> None
    write_next_round_manifest(rows, path, next_generation: str, rp: Sample, donor: Sample) -> None
    results_csv_text(rows) -> str
    selection_csv_text(rows, notes) -> str
    next_round_manifest_text(rows, next_generation, rp, donor) -> str

The ``*_text`` wrappers return exactly what the path writers put in the file
(``csv.writer`` line endings, CRLF), so a download encoded as UTF-8 is
byte-identical to the CLI's output. Both results.csv and selected.csv end with the
``token_profile`` column (contract 1.4.0): results rows carry it as their last key, so
``_columns`` puts it last; selected.csv appends it after ``results_schema``.

results.csv is schema 1.0.0 (docs/adr/0016): a missing value is written ``NA`` and empty
text (``exclusion_reason``, ``qc_flags``, ``notes``) stays an empty cell; with no rows the
header is exactly ``FIXED_COLUMNS``. selected.csv carries ``results_schema`` too;
next_samples.csv does not, and keeps empty cells, because it is the input contract's
samples.csv. A placeholder-row count below 1 is refused with ``DataContractError``: a manifest
holding the two parents and no progeny is a silently truncated file rather than a smaller one.
"""

from __future__ import annotations

import csv
from io import StringIO
from pathlib import Path
from typing import TextIO

from progeny_selector.constants import RESULTS_SCHEMA
from progeny_selector.model.dataset import DataContractError, Sample

__all__ = [
    "FIXED_COLUMNS",
    "LEADING_COLUMNS",
    "MANIFEST_COLUMNS",
    "NA",
    "RESULTS_SCHEMA",
    "SELECTION_COLUMNS",
    "check_per_selected",
    "next_round_manifest_text",
    "results_csv_text",
    "selection_csv_text",
    "write_next_round_manifest",
    "write_results_csv",
    "write_selection_csv",
]

# The text written for a missing value (docs/adr/0016). Empty text stays an empty cell.
NA = "NA"

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

# The fixed prefix of results.csv: the leading columns, the composition columns, then the metadata
# columns of schema 1.0.0 (docs/adr/0016), and last ``token_profile`` (contract 1.4.0), which the
# dynamic per-locus and per-chromosome columns are written before.
FIXED_COLUMNS = (
    *LEADING_COLUMNS,
    "role",
    "n_informative_called",
    "frac_a",
    "frac_h",
    "frac_b",
    "background_model",
    "background_unit",
    "rank_mode",
    "assembly",
    "results_schema",
    "token_profile",
)

SELECTION_COLUMNS = (
    "sample_id",
    "line_name",
    "family_id",
    "generation",
    "rank_overall",
    "rank_in_family",
    "composite_score",
    "rpp_total",
    "notes",
    "results_schema",
    "token_profile",
)


# samples.csv contract columns, in the order docs/data-formats.md lists them.
MANIFEST_COLUMNS = ("sample_id", "line_name", "role", "generation", "family_id", "notes")


def _columns(rows: list[dict]) -> list[str]:
    if not rows:
        return list(FIXED_COLUMNS)
    # ``token_profile`` is left to dict order so it stays last, after the dynamic columns.
    seen: list[str] = [c for c in FIXED_COLUMNS[:-1] if c in rows[0]]
    for row in rows:
        for c in row:
            if c not in seen:
                seen.append(c)
    return seen


def _fmt(v):
    if v is None:
        return NA
    if isinstance(v, bool):
        return "TRUE" if v else "FALSE"
    if isinstance(v, float):
        if v != v:
            return NA
        return f"{v:.6g}" if abs(v) < 1e6 else f"{v:.0f}"
    return v


def _write_results_csv(fh: TextIO, rows: list[dict]) -> None:
    cols = _columns(rows)
    writer = csv.writer(fh)
    writer.writerow(cols)
    for row in rows:
        writer.writerow([_fmt(row.get(c)) for c in cols])


def write_results_csv(rows: list[dict], path: str | Path) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        _write_results_csv(fh, rows)


def results_csv_text(rows: list[dict]) -> str:
    buf = StringIO(newline="")
    _write_results_csv(buf, rows)
    return buf.getvalue()


def _write_selection_csv(fh: TextIO, rows: list[dict], notes: dict[str, str] | None = None) -> None:
    notes = notes or {}
    writer = csv.writer(fh)
    writer.writerow(SELECTION_COLUMNS)
    for r in rows:
        writer.writerow(
            [
                r["sample_id"],
                r.get("line_name", ""),
                r.get("family_id") or NA,
                r.get("generation") or NA,
                _fmt(r.get("rank_overall")),
                _fmt(r.get("rank_in_family")),
                _fmt(r.get("composite_score")),
                _fmt(r.get("rpp_total")),
                notes.get(r["sample_id"], ""),
                RESULTS_SCHEMA,
                r.get("token_profile", ""),
            ]
        )


def write_selection_csv(rows: list[dict], path: str | Path, notes: dict[str, str] | None = None) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        _write_selection_csv(fh, rows, notes)


def selection_csv_text(rows: list[dict], notes: dict[str, str] | None = None) -> str:
    buf = StringIO(newline="")
    _write_selection_csv(buf, rows, notes)
    return buf.getvalue()


def check_per_selected(n_per_selected: int) -> None:
    """Every manifest writer's check on the placeholder-row count: below 1 there is nothing to write.

    Zero or less would emit the two parent rows and no progeny — a manifest that looks written and
    carries no next generation — so it is refused here, where both the CLI and the Export screen
    reach it, rather than in either caller. The CLI also rejects it at the argparse level, which is
    what gives ``--per-selected 0`` a usage exit code.
    """
    if n_per_selected < 1:
        raise DataContractError(f"placeholder rows per selected individual must be 1 or more, got {n_per_selected}")


def _write_next_round_manifest(
    fh: TextIO, rows: list[dict], next_generation: str, rp: Sample, donor: Sample, n_per_selected: int = 1
) -> None:
    check_per_selected(n_per_selected)
    writer = csv.writer(fh)
    writer.writerow(MANIFEST_COLUMNS)
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


def write_next_round_manifest(
    rows: list[dict], path: str | Path, next_generation: str, rp: Sample, donor: Sample, n_per_selected: int = 1
) -> None:
    """samples.csv for the next round: both parents plus placeholder progeny ids derived from each selected line."""
    with open(path, "w", encoding="utf-8", newline="") as fh:
        _write_next_round_manifest(fh, rows, next_generation, rp, donor, n_per_selected)


def next_round_manifest_text(rows: list[dict], next_generation: str, rp: Sample, donor: Sample, n_per_selected: int = 1) -> str:
    buf = StringIO(newline="")
    _write_next_round_manifest(buf, rows, next_generation, rp, donor, n_per_selected)
    return buf.getvalue()
