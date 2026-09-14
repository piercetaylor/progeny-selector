"""Build a progeny-selector dataset for one backcross-derived NIL family from SoySNP50K + PATRIOT pedigrees.

Reads the SoyBase SoySNP50K VCF (streamed, never loaded into memory) and PATRIOT's
``GS Pedigrees.csv``, selects the NILs belonging to one recurrent x donor backcross
family, and writes ``genotypes.vcf``, ``samples.csv`` and ``README.md`` under
``--out`` following docs/data-formats.md ("Genotype file > VCF" and "samples.csv").

Usage: python scripts/soysnp50k_nils.py --vcf PATH --pedigrees PATH --recurrent
    Clark|Harosoy --donor NAME [--donor-pi PI] [--generation TEMPLATE] --out DIR
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path
from typing import TextIO

RECURRENT_PI = {"Clark": "PI548533", "Harosoy": "PI548573"}
PI_FORM_RE = re.compile(r"^PI\d+[A-Z]?(-\d+)?$")
GM_CHROM_RE = re.compile(r"^Gm(0[1-9]|1[0-9]|20)$")
CITATION = "Song et al. 2015, G3 5(10):1999-2006, doi:10.1534/g3.115.019000 (data licence: Open, per the SoyBase README)"


def open_maybe_gzip(path: Path) -> TextIO:
    name = str(path).lower()
    if name.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path)


def find_family(pedigrees_path: Path, recurrent: str, donor: str) -> tuple[list[tuple[str, str, int]], int]:
    """Return (rows of (line, cross_string, n) sorted by line), n_in_pedigree."""
    recurrent_pi = RECURRENT_PI[recurrent]
    cross_re = re.compile(rf"^{re.escape(recurrent)} \((\d+)\) x (.+)$")
    seen: dict[str, tuple[str, str, int]] = {}
    with open(pedigrees_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            line = (row.get("Line") or "").strip()
            female = (row.get("Female Parent") or "").strip()
            male = (row.get("Male Parent") or "").strip()
            if not line.startswith("PI"):
                continue
            if line in seen:
                continue
            for a, b in ((female, male), (male, female)):
                if a in (recurrent_pi, recurrent):
                    m = cross_re.match(b)
                    if m and m.group(2).strip() == donor:
                        seen[line] = (line, b, int(m.group(1)))
                        break
    rows = sorted(seen.values(), key=lambda r: r[0])
    return rows, len(rows)


def strip_chrom(chrom: str) -> str:
    idx = chrom.rfind(".")
    return chrom[idx + 1 :] if idx != -1 else chrom


def gt_allele_index_exceeds_alt(gt: str, n_alt: int) -> bool:
    """True when a GT subfield names an allele index greater than the number of ALT alleles."""
    for allele in re.split(r"[/|]", gt):
        if allele == ".":
            continue
        if int(allele) > n_alt:
            return True
    return False


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--pedigrees", required=True, type=Path)
    parser.add_argument("--recurrent", required=True, choices=sorted(RECURRENT_PI))
    parser.add_argument("--donor", required=True)
    parser.add_argument("--donor-pi", default=None)
    parser.add_argument("--generation", default="")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)

    donor_pi = args.donor_pi
    if donor_pi is None:
        if PI_FORM_RE.match(args.donor):
            donor_pi = args.donor
        else:
            print(
                f"error: --donor {args.donor!r} is not PI-form and --donor-pi was not given",
                file=sys.stderr,
            )
            return 2

    recurrent_pi = RECURRENT_PI[args.recurrent]

    family_rows, n_in_pedigree = find_family(args.pedigrees, args.recurrent, args.donor)
    if n_in_pedigree == 0:
        print("error: no NILs found in pedigree for this recurrent x donor cross", file=sys.stderr)
        return 2

    out_dir: Path = args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    with open_maybe_gzip(args.vcf) as fh:
        header_cols: list[str] | None = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_cols = line.rstrip("\n").split("\t")
                break
        if header_cols is None:
            print("error: no #CHROM header line found in VCF", file=sys.stderr)
            return 2

        sample_cols = header_cols[9:]
        col_index = {sid: i for i, sid in enumerate(sample_cols)}

        if recurrent_pi not in col_index:
            print(f"error: recurrent {recurrent_pi} not found in VCF header", file=sys.stderr)
            return 2
        if donor_pi not in col_index:
            print(f"error: donor {donor_pi} not found in VCF header", file=sys.stderr)
            return 2

        nils_present: list[tuple[str, str, int]] = []
        n_missing_from_vcf = 0
        for line_id, cross_str, n in family_rows:
            if line_id not in col_index:
                print(f"warning: NIL {line_id} not found in VCF header, dropping", file=sys.stderr)
                n_missing_from_vcf += 1
                continue
            nils_present.append((line_id, cross_str, n))

        if not nils_present:
            print("error: zero NILs remaining after matching against VCF header", file=sys.stderr)
            return 2

        selected_samples = [recurrent_pi, donor_pi] + [r[0] for r in nils_present]
        selected_idx = [9 + col_index[sid] for sid in selected_samples]

        n_records_written = 0
        n_dropped_non_gm = 0
        n_dropped_bad_gt = 0

        genotypes_path = out_dir / "genotypes.vcf"
        with open(genotypes_path, "w", newline="\n") as out_fh:
            out_fh.write("##fileformat=VCFv4.2\n")
            out_fh.write(f"##source=soysnp50k_nils.py --recurrent {args.recurrent} --donor {args.donor}\n")
            out_fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(selected_samples) + "\n")
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                chrom = strip_chrom(fields[0])
                if not GM_CHROM_RE.match(chrom):
                    n_dropped_non_gm += 1
                    continue
                pos, rid, ref, alt = fields[1], fields[2], fields[3], fields[4]
                n_alt = 0 if alt == "." else len(alt.split(","))
                gts = []
                bad_gt = False
                for idx in selected_idx:
                    cell = fields[idx]
                    gt = cell.split(":", 1)[0]
                    if gt_allele_index_exceeds_alt(gt, n_alt):
                        bad_gt = True
                        break
                    gts.append(gt)
                if bad_gt:
                    n_dropped_bad_gt += 1
                    continue
                out_fh.write(f"{chrom}\t{pos}\t{rid}\t{ref}\t{alt}\t.\t.\t.\tGT\t" + "\t".join(gts) + "\n")
                n_records_written += 1

    samples_path = out_dir / "samples.csv"
    family_id = f"{args.recurrent}_x_{args.donor}"
    with open(samples_path, "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["sample_id", "line_name", "role", "generation", "family_id", "notes"])
        writer.writerow([recurrent_pi, args.recurrent, "recurrent_parent", "", "", "recurrent parent"])
        writer.writerow([donor_pi, args.donor, "donor_parent", "", "", "donor parent"])
        for line_id, cross_str, n in nils_present:
            generation = args.generation.format(n=n) if args.generation else ""
            writer.writerow([line_id, line_id, "candidate", generation, family_id, cross_str])

    readme_path = out_dir / "README.md"
    readme_lines = [
        "# SoySNP50K NIL dataset",
        "",
        "## Source files",
        f"- VCF: {args.vcf}",
        f"- Pedigrees: {args.pedigrees}",
        "",
        "## Command",
        "```",
        " ".join(sys.argv),
        "```",
        "",
        "## Counts",
        f"- NILs in pedigree: {n_in_pedigree}",
        f"- NILs written: {len(nils_present)}",
        f"- NILs missing from VCF: {n_missing_from_vcf}",
        f"- Markers written: {n_records_written}",
        f"- Records dropped as non-Gm01..Gm20: {n_dropped_non_gm}",
        f"- Records dropped: GT allele index exceeds ALT count: {n_dropped_bad_gt}",
        "",
        "Positions are Wm82.a2 (gnm2); progeny-selector's chromosome-length table is "
        "Wm82.a4.v1, so drag bounds near chromosome ends are approximate.",
        "",
        "Records dropped for GT allele index exceeding the ALT count are malformed in the "
        'SoyBase source (ALT "." with a non-reference call) and cannot be repaired without '
        "the true alternate allele.",
        "",
        "NILs are advanced inbred lines.",
        "",
        f"Citation: {CITATION}",
        "",
    ]
    readme_path.write_text("\n".join(readme_lines), newline="\n")

    print(f"NILs in pedigree: {n_in_pedigree}")
    print(f"NILs written: {len(nils_present)}")
    print(f"NILs missing from VCF: {n_missing_from_vcf}")
    print(f"Markers written: {n_records_written}")
    print(f"Records dropped as non-Gm01..Gm20: {n_dropped_non_gm}")
    print(f"Records dropped: GT allele index exceeds ALT count: {n_dropped_bad_gt}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
