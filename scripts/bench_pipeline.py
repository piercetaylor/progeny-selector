"""Generate a synthetic dataset at a given size and measure wall-clock and peak memory
through the CPython pipeline (docs/adr/0021, docs/limits.md).

Usage:
    python scripts/bench_pipeline.py --markers 50000 --progeny 2000 --out data/bench/50k_x_2000 \
        [--seed 0] [--no-gzip] [--regenerate]

Two parts. ``generate_case`` writes a deterministic synthetic VCF, samples.csv, markers.csv
and criteria.yaml under ``out/`` (skipped when ``out/genotypes.vcf`` already exists, unless
``--regenerate``). ``measure_case`` loads and analyses the case, timing each stage with
``time.perf_counter`` and taking the CPython memory figure from ``tracemalloc.get_traced_memory()[1]``
around each stage separately (numpy reports its allocations to ``tracemalloc``), so each figure
excludes whatever was already live from the previous stage. Prints one Markdown row for
``docs/limits.md`` and writes ``out/bench.json`` with every figure plus the machine and package
versions; the commit hash is written by the main session, not this script (doers never run git).

Duplicate detection (``core.qc.duplicate_pairs``) runs only at 2,000 progeny or fewer
(``core/pipeline.py`` skips the quadratic scan above that); this script marks the row
accordingly and, run with more than 2,000 progeny, prints "duplicate scan: skipped".
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import platform
import shutil
import sys
import tracemalloc
from pathlib import Path
from time import perf_counter

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from progeny_selector.constants import SOYBEAN_CHROM_LENGTHS_BP_WM82A4
from progeny_selector.core.pipeline import run_analysis
from progeny_selector.core.qc import duplicate_pairs
from progeny_selector.io import load_dataset, read_criteria
from progeny_selector.io.export import write_results_csv

CHROMS = [f"Gm{i:02d}" for i in range(1, 21)]
RP_ID = "RP_bench"
DONOR_ID = "DONOR_bench"
_MAX_DUPLICATE_SAMPLES = 2000  # mirrors core/pipeline.py::_MAX_DUPLICATE_SAMPLES

_STATE_TOKENS = {"A": "0/0", "H": "0/1", "B": "1/1", "N": "./."}
_STATE_ORDER = ("A", "H", "B", "N")
_STATE_PROBS = (0.85, 0.12, 0.01, 0.02)


def _fixture_criteria_doc() -> dict:
    fixture = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "synthetic_bc2f1" / "criteria.yaml"
    with open(fixture, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def generate_case(out: Path, n_markers: int, n_progeny: int, seed: int = 0, gzip_too: bool = True) -> None:
    """Write genotypes.vcf[.gz], samples.csv, markers.csv and criteria.yaml under ``out``.

    Deterministic (``np.random.default_rng(seed)``). Skipped entirely when ``out/genotypes.vcf``
    already exists; the caller passes a fresh directory to force regeneration.
    """
    out.mkdir(parents=True, exist_ok=True)
    vcf_path = out / "genotypes.vcf"
    if vcf_path.exists():
        return

    rng = np.random.default_rng(seed)
    per = math.ceil(n_markers / 20)

    marker_ids: list[str] = []
    marker_chroms: list[str] = []
    marker_pos: list[int] = []
    i = 0
    for chrom in CHROMS:
        length = SOYBEAN_CHROM_LENGTHS_BP_WM82A4[chrom]
        for k in range(per):
            pos = round((k + 1) * length / (per + 1))
            marker_ids.append(f"bm{i:06d}")
            marker_chroms.append(chrom)
            marker_pos.append(pos)
            i += 1
    n_total = len(marker_ids)

    # Donor is uninformative (0/0, matching RP) at every 20th marker (global index).
    donor_uninformative = [idx % 20 == 0 for idx in range(n_total)]

    sample_ids = [f"P{j + 1:05d}" for j in range(n_progeny)]

    # Per-marker progeny state draws, iid, overridden to A (0/0) at monomorphic markers.
    codes = np.arange(len(_STATE_ORDER))
    draws = rng.choice(codes, size=(n_total, n_progeny), p=_STATE_PROBS)

    with open(vcf_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join([RP_ID, DONOR_ID, *sample_ids]) + "\n")
        for m in range(n_total):
            donor_gt = "0/0" if donor_uninformative[m] else "1/1"
            row_tokens = ["0/0", donor_gt]
            for j in range(n_progeny):
                if donor_uninformative[m]:
                    row_tokens.append("0/0")
                else:
                    row_tokens.append(_STATE_TOKENS[_STATE_ORDER[draws[m, j]]])
            fh.write(f"{marker_chroms[m]}\t{marker_pos[m]}\t{marker_ids[m]}\tA\tT\t.\t.\t.\tGT\t" + "\t".join(row_tokens) + "\n")

    if gzip_too:
        gz_path = out / "genotypes.vcf.gz"
        with open(vcf_path, "rb") as src, gzip.open(gz_path, "wb", compresslevel=6) as dst:
            shutil.copyfileobj(src, dst)

    with open(out / "samples.csv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("sample_id,line_name,role,generation,family_id,notes\n")
        fh.write(f"{RP_ID},{RP_ID},recurrent_parent,,,\n")
        fh.write(f"{DONOR_ID},{DONOR_ID},donor_parent,,,\n")
        for j, sid in enumerate(sample_ids):
            fh.write(f"{sid},{sid},progeny,BC2F1,F{j % 20:02d},\n")

    with open(out / "markers.csv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("marker_id,chrom,pos_bp,cm\n")
        for mid, chrom, pos in zip(marker_ids, marker_chroms, marker_pos, strict=True):
            cm = pos / 1e6 * 2.5
            fh.write(f"{mid},{chrom},{pos},{cm}\n")

    doc = _fixture_criteria_doc()
    gm06_len = SOYBEAN_CHROM_LENGTHS_BP_WM82A4["Gm06"]
    mid_gm06 = gm06_len // 2
    # Half-width is max(1_000_000, d), not a flat 1 Mb: markers on a chromosome are evenly
    # spaced at interval d = L / (per + 1), so every point on the chromosome lies within d/2
    # of a marker; a flat 1 Mb window is empty whenever d exceeds it (every small case,
    # including this generator's own smoke-test size), which is not what "one region on Gm06"
    # was meant to test. At the measured sizes (6,000+ markers) d is well under 1 Mb, so the
    # window is unchanged there.
    half_width = max(1_000_000, math.ceil(gm06_len / (per + 1)))
    doc["targets"] = [
        {
            "locus_id": "T1",
            "chrom": "Gm06",
            "start_bp": mid_gm06 - half_width,
            "end_bp": mid_gm06 + half_width,
            "required_state": "het",
        }
    ]
    gm10_len = SOYBEAN_CHROM_LENGTHS_BP_WM82A4["Gm10"]
    mid_gm10 = gm10_len // 2
    gm10_markers = [(mid, pos) for mid, chrom, pos in zip(marker_ids, marker_chroms, marker_pos, strict=True) if chrom == "Gm10"]
    nearest_id, _ = min(gm10_markers, key=lambda mp: abs(mp[1] - mid_gm10))
    doc["avoid"] = [{"locus_id": "AV1", "marker_id": nearest_id}]
    doc["background"]["model"] = "weighted"
    doc["filters"]["exclude_qc_flagged"] = False

    with open(out / "criteria.yaml", "w", encoding="utf-8", newline="\n") as fh:
        yaml.safe_dump(doc, fh, sort_keys=False, default_flow_style=False)


def measure_case(case: Path) -> dict:
    """Load, analyse and export the case under ``case``, timing and memory-profiling each stage."""
    vcf = case / "genotypes.vcf"
    vcf_gz = case / "genotypes.vcf.gz"
    samples = case / "samples.csv"
    markers = case / "markers.csv"
    criteria_path = case / "criteria.yaml"

    result: dict = {
        "case": case.name,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "python_version": sys.version,
        "numpy_version": np.__version__,
    }

    stage = "parse"
    try:
        tracemalloc.start()
        t0 = perf_counter()
        dataset = load_dataset(vcf, samples, markers)
        parse_s = perf_counter() - t0
        parse_peak_bytes = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()

        resident_bytes = dataset.genotypes.calls.nbytes
        n_markers = dataset.genotypes.n_markers
        n_progeny = len(dataset.progeny)

        stage = "parse (gz)"
        parse_gz_s = None
        if vcf_gz.exists():
            t0 = perf_counter()
            load_dataset(vcf_gz, samples, markers)
            parse_gz_s = perf_counter() - t0

        criteria = read_criteria(criteria_path)

        stage = "analysis"
        tracemalloc.start()
        t0 = perf_counter()
        analysis = run_analysis(dataset, criteria)
        analysis_s = perf_counter() - t0
        analysis_peak_bytes = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()

        stage = "duplicate scan"
        dedup_ran = n_progeny <= _MAX_DUPLICATE_SAMPLES
        sorted_gm = dataset.genotypes.sorted_by_position(dataset.scheme)
        progeny_ids = [s.sample_id for s in dataset.progeny]
        marker_idx = np.flatnonzero(analysis.classification.informative)
        t0 = perf_counter()
        duplicate_pairs(sorted_gm, progeny_ids, max_samples=_MAX_DUPLICATE_SAMPLES, marker_idx=marker_idx)
        dedup_s = perf_counter() - t0

        stage = "export"
        t0 = perf_counter()
        write_results_csv(analysis.rows, case / "results.csv")
        export_s = perf_counter() - t0

        n_pass = sum(1 for r in analysis.rows if r["passes_filters"])
        total_s = parse_s + analysis_s + dedup_s + export_s

        result.update(
            {
                "markers": n_markers,
                "progeny": n_progeny,
                "vcf_bytes_plain": vcf.stat().st_size,
                "vcf_bytes_gz": vcf_gz.stat().st_size if vcf_gz.exists() else None,
                "parse_s": parse_s,
                "parse_gz_s": parse_gz_s,
                "parse_peak_bytes": parse_peak_bytes,
                "resident_bytes": resident_bytes,
                "analysis_s": analysis_s,
                "analysis_peak_bytes": analysis_peak_bytes,
                "dedup_s": dedup_s,
                "dedup_ran": dedup_ran,
                "export_s": export_s,
                "total_s": total_s,
                "n_pass": n_pass,
                "outcome": "ok",
            }
        )
        if not dedup_ran:
            print("duplicate scan: skipped")
    except MemoryError:
        result["outcome"] = f"MemoryError at {stage}"
        print(f"outcome: MemoryError at {stage}")
        _write_json(case, result)
        sys.exit(5)

    _print_row(result)
    _write_json(case, result)
    return result


def _mib(n: int | None) -> str:
    if n is None:
        return "-"
    return f"{n / (1024 * 1024):.1f}"


def _mb(n: int | None) -> str:
    if n is None:
        return "-"
    return f"{n / 1_000_000:.1f}"


def _s(n: float | None) -> str:
    if n is None:
        return "-"
    return f"{n:.2f}"


def _print_row(r: dict) -> None:
    dedup_col = f"{_s(r['dedup_s'])} ({'ran' if r['dedup_ran'] else 'skipped'})"
    row = (
        f"| {r['case']} | {r['markers']} | {r['progeny']} | {_mb(r['vcf_bytes_plain'])} / {_mb(r['vcf_bytes_gz'])} "
        f"| {_s(r['parse_s'])} / {_s(r['parse_gz_s'])} | {_mib(r['parse_peak_bytes'])} | {_mib(r['resident_bytes'])} "
        f"| {_s(r['analysis_s'])} | {_mib(r['analysis_peak_bytes'])} | {dedup_col} | {_s(r['export_s'])} "
        f"| {_s(r['total_s'])} | {r['outcome']} |"
    )
    print(row)


def _write_json(case: Path, result: dict) -> None:
    import datetime

    result = dict(result)
    result["date"] = datetime.datetime.now(datetime.UTC).date().isoformat()
    with open(case / "bench.json", "w", encoding="utf-8", newline="\n") as fh:
        json.dump(result, fh, indent=2, sort_keys=True)
        fh.write("\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--markers", type=int, required=True)
    parser.add_argument("--progeny", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--no-gzip", action="store_true")
    parser.add_argument("--regenerate", action="store_true")
    args = parser.parse_args(argv)

    if args.regenerate and (args.out / "genotypes.vcf").exists():
        (args.out / "genotypes.vcf").unlink()
        gz = args.out / "genotypes.vcf.gz"
        if gz.exists():
            gz.unlink()

    generate_case(args.out, args.markers, args.progeny, seed=args.seed, gzip_too=not args.no_gzip)
    measure_case(args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
