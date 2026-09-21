"""Command-line interface: validate, rank, select.

Responsibility: thin argparse layer over ``progeny_selector.io`` and
``progeny_selector.core``; exit codes 0 (ok), 1 (data-contract or criteria
error), 2 (usage). Useful for batch runs on an HPC node and for wiring
outputs into the R reader ``scripts/read_results.R``.

Interface (subcommands):
    progeny-selector validate --genotypes G --samples S [--markers M] [--criteria C] [--profile ID_OR_FILE]
    progeny-selector rank     --genotypes G --samples S --criteria C [--markers M] [--profile ID_OR_FILE] --out results.csv
    progeny-selector select   --results results.csv --top N [--overall] --out selected.csv
                              [--next-manifest next_samples.csv --next-generation BC3F1 --samples S
                               --per-selected N]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys

from progeny_selector.core.pipeline import run_analysis
from progeny_selector.core.selection import project_next_generation, select_top_n
from progeny_selector.io import load_dataset, read_criteria, read_samples
from progeny_selector.io.export import write_next_round_manifest, write_results_csv, write_selection_csv
from progeny_selector.model.criteria import CriteriaError
from progeny_selector.model.dataset import DataContractError


def _add_data_args(p: argparse.ArgumentParser, criteria_required: bool) -> None:
    p.add_argument("--genotypes", required=True, help="VCF(.gz), HapMap (.hmp.txt) or wide CSV")
    p.add_argument("--samples", required=True, help="samples.csv manifest")
    p.add_argument("--markers", help="optional markers.csv with chrom, pos_bp, cm")
    p.add_argument("--criteria", required=criteria_required, help="criteria.yaml")
    p.add_argument("--coding", default="auto", choices=("auto", "nucleotide", "abh"), help="wide CSV call coding")
    p.add_argument(
        "--profile",
        metavar="ID_OR_FILE",
        help="token profile for HapMap and wide CSV cells: a built-in id (tassel, soybase-report, dart, axiom, kasp) or a JSON file",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="progeny-selector", description="MABC progeny ranking and selection")
    sub = parser.add_subparsers(dest="command", required=True)

    v = sub.add_parser("validate", help="check files against the data contract and print QC warnings")
    _add_data_args(v, criteria_required=False)

    r = sub.add_parser("rank", help="compute statuses, scores and ranks; write results.csv")
    _add_data_args(r, criteria_required=True)
    r.add_argument("--out", required=True, help="output results.csv")

    s = sub.add_parser("select", help="pick top N per family (or overall) from results.csv")
    s.add_argument("--results", required=True)
    s.add_argument("--top", type=int, required=True)
    s.add_argument("--overall", action="store_true", help="top N overall instead of per family")
    s.add_argument("--out", required=True, help="selected.csv")
    s.add_argument("--next-manifest", help="write a samples.csv skeleton for the next round")
    s.add_argument("--next-generation", default="BC3F1")
    s.add_argument(
        "--per-selected",
        type=_positive_int,
        default=1,
        metavar="N",
        help="placeholder progeny rows per selected individual in --next-manifest (default 1)",
    )
    s.add_argument("--samples", help="current samples.csv (needed for --next-manifest to copy the parents)")
    s.add_argument("--project", choices=("backcross", "self"), default="backcross")
    return parser


def _positive_int(value: str) -> int:
    """argparse type for a count of at least 1; anything else is a usage error (exit 2).

    The rule itself lives in ``io.export.check_per_selected``, which every manifest writer runs;
    this is the CLI's earlier, exit-code-shaped report of the same rule.
    """
    try:
        n = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected an integer, got {value!r}") from None
    if n < 1:
        raise argparse.ArgumentTypeError(f"must be 1 or more, got {n}")
    return n


def _profile_ref(ref: str | None) -> str | dict | None:
    """A built-in id as given; a value containing / or \\ or ending .json is read as a profile JSON file."""
    if ref is None or not ("/" in ref or "\\" in ref or ref.lower().endswith(".json")):
        return ref
    try:
        with open(ref, encoding="utf-8") as fh:
            obj = json.load(fh)
    except (OSError, ValueError) as exc:
        raise DataContractError(f"token profile file {ref}: {exc}") from exc
    if not isinstance(obj, dict):
        raise DataContractError("token profile: must be a JSON object")
    return obj


def _load(args: argparse.Namespace):
    dataset = load_dataset(args.genotypes, args.samples, args.markers, coding=args.coding, profile=_profile_ref(args.profile))
    criteria = read_criteria(args.criteria) if args.criteria else None
    return dataset, criteria


def cmd_validate(args: argparse.Namespace) -> int:
    dataset, criteria = _load(args)
    gm = dataset.genotypes
    print(f"genotypes: {gm.n_markers} markers x {gm.n_samples} samples; coded={gm.coded}; cM map={'yes' if gm.has_cm() else 'no'}")
    print(f"samples.csv: {len(dataset.samples)} rows; progeny/candidates={len(dataset.progeny)}")
    for w in dataset.warnings:
        print(f"warning: {w}")
    if criteria is not None:
        res = run_analysis(dataset, criteria)
        print(f"informative markers: {res.classification.n_informative}")
        for w in res.warnings:
            print(f"warning: {w}")
        flagged = [q for q in res.qc if q.flags]
        print(f"QC-flagged individuals: {len(flagged)}")
        for q in flagged[:20]:
            print(f"  {q.sample_id}: {'|'.join(q.flags)}")
        print(f"duplicate pairs: {len(res.duplicates)}")
        for a, b, ibs in res.duplicates[:20]:
            print(f"  {a}, {b}: ibs {ibs:.3f}")
    return 0


def cmd_rank(args: argparse.Namespace) -> int:
    dataset, criteria = _load(args)
    res = run_analysis(dataset, criteria)
    write_results_csv(res.rows, args.out)
    n_pass = sum(1 for r in res.rows if r["passes_filters"])
    print(f"wrote {args.out}: {len(res.rows)} individuals, {n_pass} pass hard filters")
    for w in res.warnings:
        print(f"warning: {w}")
    return 0


NUMERIC_RESULT_COLUMNS = ("rank_overall", "rank_in_family", "composite_score", "rpp_total", "frac_a", "frac_h", "frac_b")
# Columns where a missing value is meaningful and ``NA`` therefore means "missing" rather than the text "NA".
# Identifier and free-text columns (sample_id, line_name, notes, qc_flags, ...) keep whatever the file says,
# so a line genuinely named NA survives ``rank`` -> ``select`` -> next_samples.csv (docs/adr/0016).
NA_RESULT_COLUMNS = frozenset({*NUMERIC_RESULT_COLUMNS, "family_id", "generation"})


def _read_results(path: str) -> list[dict]:
    rows: list[dict] = []
    with open(path, encoding="utf-8", newline="") as fh:
        for raw in csv.DictReader(fh):
            # results.csv writes NA for a missing value (docs/adr/0016); a pre-freeze file writes an empty cell.
            row: dict = {k: (None if (v == "NA" and k in NA_RESULT_COLUMNS) else v) for k, v in raw.items()}
            for k in NUMERIC_RESULT_COLUMNS:
                row[k] = float(row[k]) if row.get(k) not in (None, "") else None
            row["passes_filters"] = row.get("passes_filters") == "TRUE"
            rows.append(row)
    return rows


def cmd_select(args: argparse.Namespace) -> int:
    rows = _read_results(args.results)
    chosen = select_top_n(rows, args.top, per_family=not args.overall)
    write_selection_csv(chosen, args.out)
    print(f"wrote {args.out}: {len(chosen)} selected")
    if chosen:
        proj = project_next_generation(chosen, args.project)
        print(
            f"projection ({proj.step}): mean RPP now {proj.mean_rpp_selected:.3f}; "
            f"expected RPP next {proj.expected_rpp_next:.3f}; expected het next {proj.expected_het_next:.3f}"
        )
        print(f"  {proj.note}")
    if args.next_manifest:
        if not args.samples:
            print("error: --next-manifest needs --samples to copy the parents", file=sys.stderr)
            return 2
        samples = read_samples(args.samples)
        rp = next(s for s in samples if s.role == "recurrent_parent")
        dp = next(s for s in samples if s.role == "donor_parent")
        write_next_round_manifest(chosen, args.next_manifest, args.next_generation, rp, dp, args.per_selected)
        print(f"wrote {args.next_manifest}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "validate":
            return cmd_validate(args)
        if args.command == "rank":
            return cmd_rank(args)
        if args.command == "select":
            return cmd_select(args)
    except (DataContractError, CriteriaError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 2
