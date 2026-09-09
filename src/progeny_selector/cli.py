"""Command-line interface: validate, rank, select.

Responsibility: thin argparse layer over ``progeny_selector.io`` and
``progeny_selector.core``; exit codes 0 (ok), 1 (data-contract or criteria
error), 2 (usage). Useful for batch runs on an HPC node and for wiring
outputs into downstream Shiny dashboards.

Interface (subcommands):
    progeny-selector validate --genotypes G --samples S [--markers M] [--criteria C]
    progeny-selector rank     --genotypes G --samples S --criteria C [--markers M] --out results.csv
    progeny-selector select   --results results.csv --top N [--overall] --out selected.csv
                              [--next-manifest next_samples.csv --next-generation BC3F1 --samples S]
"""

from __future__ import annotations

import argparse
import csv
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
    s.add_argument("--samples", help="current samples.csv (needed for --next-manifest to copy the parents)")
    s.add_argument("--project", choices=("backcross", "self"), default="backcross")
    return parser


def _load(args: argparse.Namespace):
    dataset = load_dataset(args.genotypes, args.samples, args.markers, coding=args.coding)
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


def _read_results(path: str) -> list[dict]:
    rows: list[dict] = []
    with open(path, encoding="utf-8", newline="") as fh:
        for raw in csv.DictReader(fh):
            row = dict(raw)
            for k in ("rank_overall", "rank_in_family", "composite_score", "rpp_total", "frac_a", "frac_h", "frac_b"):
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
        write_next_round_manifest(chosen, args.next_manifest, args.next_generation, rp, dp)
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
