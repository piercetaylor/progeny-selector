"""Build a SoySNP50K/6K marker position table across Wm82 assemblies from SoyBase GFF3s.

Source: SoyBase Data Store, `https://data.soybase.org/Glycine/max/markers/<dir>/glyma.<dir>.gff3.gz`
with `<dir>` in `Wm82.gnm{1,2,4,5,6}.mrk.SoySNP50K` and `Wm82.gnm{1,2,4,5,6}.mrk.SoySNP6K`
(licence Open, per each directory's `README.*.yml`; docs/adr/0019, decision 7). GFF3 rows are
tab-separated `seqid source type start end score strand phase attributes`; `marker_id` is taken
from the `Name=` attribute (not `ID=`, which is `glyma.Wm82.gnmN.ss...`; docs/m2-phases.md decision
7's "Spec corrections"). The chromosome is the `seqid` with the `glyma.Wm82.gnmN.` prefix stripped
the way `scripts/soysnp50k_nils.py::strip_chrom` does (the text after the last `.`), then normalised
by `core.chrom.normalize_chrom`; a gnm5 seqid such as `glyma.Wm82.gnm5.Chr01` strips to `Chr01` and
normalises to `Gm01` the same way a gnm2 `...Gm01` seqid does. Scaffolds keep their stripped name.
There are no `##sequence-region` pragmas to rely on (decision 7).

The panel (SoySNP50K or SoySNP6K) is read from each `--gff3` file's path, since the SoyBase
directory names carry it and the phase does not define a separate `--panel` flag.

Usage:
    python3 scripts/soysnp_positions.py --gff3 Wm82.a2=data/soybase/glyma.Wm82.gnm2.mrk.SoySNP50K.gff3.gz \\
        [--gff3 ...] --out data/soybase/soysnp_positions.csv \\
        [--emit-markers-csv Wm82.a2 --markers-out markers.csv] [--download]

`--download` fetches every file in `SOYBASE_URLS` into `data/soybase/`, sends a browser
User-Agent, and verifies each directory's `CHECKSUM.*.md5` before use; it is never run in CI
and never exercised by the test suite.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path
from urllib.parse import urljoin

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from progeny_selector.core.chrom import chrom_sort_key, normalize_chrom
from progeny_selector.model.dataset import DataContractError

ASSEMBLIES: tuple[str, ...] = ("Wm82.a1", "Wm82.a2", "Wm82.a4", "Wm82.a5", "Wm82.a6")
GM01_20 = {f"Gm{i:02d}" for i in range(1, 21)}
PANELS = ("SoySNP50K", "SoySNP6K")
BASE_URL = "https://data.soybase.org/Glycine/max/markers/"
BROWSER_USER_AGENT = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"

# Directory name (data.soybase.org/Glycine/max/markers/<dir>/) -> gff3 filename, per decision 7.
SOYBASE_URLS: dict[str, str] = {
    f"Wm82.gnm{gnm}.mrk.{panel}": f"glyma.Wm82.gnm{gnm}.mrk.{panel}.gff3.gz" for gnm in (1, 2, 4, 5, 6) for panel in PANELS
}


def strip_chrom(chrom: str) -> str:
    """The text after the last '.', mirroring scripts/soysnp50k_nils.py::strip_chrom."""
    idx = chrom.rfind(".")
    return chrom[idx + 1 :] if idx != -1 else chrom


def _panel_from_path(path: Path) -> str:
    name = str(path)
    for panel in PANELS:
        if panel in name:
            return panel
    raise DataContractError(f"{path}: cannot tell SoySNP50K from SoySNP6K in the file name")


def _open_maybe_gzip(path: Path):
    if str(path).lower().endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path)


def _parse_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for field in text.strip().split(";"):
        field = field.strip()
        if not field or "=" not in field:
            continue
        key, value = field.split("=", 1)
        attrs[key.strip()] = value.strip()
    return attrs


def parse_gff3(path: Path) -> list[tuple[str, str, int]]:
    """Return (marker_id, chrom, pos_bp) for every feature row; comment and blank lines skipped."""
    rows: list[tuple[str, str, int]] = []
    with _open_maybe_gzip(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise DataContractError(f"{path}: line {line_no}: expected 9 tab-separated fields, found {len(fields)}")
            seqid, _source, _type, start, _end, _score, _strand, _phase, attributes = fields
            attrs = _parse_attributes(attributes)
            marker_id = attrs.get("Name")
            if not marker_id:
                raise DataContractError(f"{path}: line {line_no}: no Name= attribute")
            chrom = normalize_chrom(strip_chrom(seqid))
            try:
                pos_bp = int(start)
            except ValueError as exc:
                raise DataContractError(f"{path}: line {line_no}: invalid start {start!r}") from exc
            rows.append((marker_id, chrom, pos_bp))
    return rows


def build_table(
    files: list[tuple[str, Path]],
) -> tuple[dict[str, dict[str, object]], dict[str, dict[str, int]]]:
    """From [(assembly, path), ...], the joined per-marker table and per-assembly per-chrom max position."""
    table: dict[str, dict[str, object]] = {}
    maxima: dict[str, dict[str, int]] = defaultdict(dict)
    for assembly, path in files:
        if assembly not in ASSEMBLIES:
            raise DataContractError(f"{path}: assembly {assembly!r} must be one of {ASSEMBLIES}")
        panel = _panel_from_path(path)
        for marker_id, chrom, pos_bp in parse_gff3(path):
            blank: dict[str, object] = {"in_SoySNP50K": False, "in_SoySNP6K": False}
            blank.update({f"chrom_{a}": "" for a in ASSEMBLIES})
            blank.update({f"pos_bp_{a}": "" for a in ASSEMBLIES})
            row = table.setdefault(marker_id, blank)
            row[f"in_{panel}"] = True
            if row[f"chrom_{assembly}"] == "":
                row[f"chrom_{assembly}"] = chrom
                row[f"pos_bp_{assembly}"] = pos_bp
            if chrom in GM01_20:
                maxima[assembly][chrom] = max(maxima[assembly].get(chrom, 0), pos_bp)
    return table, maxima


def write_table(table: dict[str, dict[str, object]], out: Path) -> None:
    columns = ["marker_id", "in_SoySNP50K", "in_SoySNP6K"]
    for a in ASSEMBLIES:
        columns += [f"chrom_{a}", f"pos_bp_{a}"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        fh.write(",".join(columns) + "\n")
        for marker_id in sorted(table):
            row = table[marker_id]
            cells = [marker_id, "TRUE" if row["in_SoySNP50K"] else "FALSE", "TRUE" if row["in_SoySNP6K"] else "FALSE"]
            for a in ASSEMBLIES:
                cells += [str(row[f"chrom_{a}"]), str(row[f"pos_bp_{a}"])]
            fh.write(",".join(cells) + "\n")


def write_markers_csv(table: dict[str, dict[str, object]], assembly: str, out: Path) -> None:
    """marker_id,chrom,pos_bp for one assembly, Gm01..Gm20 only, sorted by chromosome then position."""
    rows: list[tuple[str, str, int]] = []
    for marker_id, row in table.items():
        chrom = row[f"chrom_{assembly}"]
        pos = row[f"pos_bp_{assembly}"]
        if chrom in GM01_20 and pos != "":
            rows.append((marker_id, str(chrom), int(pos)))
    rows.sort(key=lambda r: (chrom_sort_key(r[1]), r[2], r[0]))
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        fh.write("marker_id,chrom,pos_bp\n")
        for marker_id, chrom, pos in rows:
            fh.write(f"{marker_id},{chrom},{pos}\n")


def _fetch(url: str, dest: Path, opener: urllib.request.OpenerDirector) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": BROWSER_USER_AGENT})
    with opener.open(req) as resp, open(dest, "wb") as fh:
        fh.write(resp.read())


def _verify_checksum(dir_url: str, dest_dir: Path, filename: str, opener: urllib.request.OpenerDirector) -> None:
    """Fetch the directory's CHECKSUM.*.md5, find filename's line, and compare against the downloaded file."""
    import hashlib

    index_req = urllib.request.Request(dir_url, headers={"User-Agent": BROWSER_USER_AGENT})
    with opener.open(index_req) as resp:
        index_html = resp.read().decode("utf-8", errors="replace")
    match = re.search(r'href="(CHECKSUM\.[^"]*\.md5)"', index_html)
    if not match:
        print(f"warning: no CHECKSUM.*.md5 listed at {dir_url}; skipping verification for {filename}", file=sys.stderr)
        return
    checksum_url = urljoin(dir_url, match.group(1))
    checksum_req = urllib.request.Request(checksum_url, headers={"User-Agent": BROWSER_USER_AGENT})
    with opener.open(checksum_req) as resp:
        checksum_text = resp.read().decode("utf-8", errors="replace")
    expected = None
    for line in checksum_text.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[1].lstrip("*") == filename:
            expected = parts[0]
            break
    if expected is None:
        print(f"warning: {filename} not listed in {checksum_url}; skipping verification", file=sys.stderr)
        return
    actual = hashlib.md5((dest_dir / filename).read_bytes()).hexdigest()
    if actual != expected:
        raise DataContractError(f"{filename}: md5 {actual} does not match {checksum_url} ({expected})")


def download(out_dir: Path) -> list[tuple[str, Path]]:
    """Fetch every SOYBASE_URLS entry into out_dir, verified by its directory's CHECKSUM.*.md5. Never run in CI."""
    opener = urllib.request.build_opener()
    fetched: list[tuple[str, Path]] = []
    for directory, filename in SOYBASE_URLS.items():
        dir_url = BASE_URL + directory + "/"
        file_url = dir_url + filename
        dest = out_dir / filename
        _fetch(file_url, dest, opener)
        _verify_checksum(dir_url, out_dir, filename, opener)
        fetched.append((directory, dest))
    return fetched


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff3", action="append", default=[], metavar="ASSEMBLY=PATH")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--emit-markers-csv", metavar="ASSEMBLY", default=None)
    parser.add_argument("--markers-out", type=Path, default=None)
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args(argv)

    if args.download:
        fetched = download(Path("data/soybase"))
        print(f"downloaded {len(fetched)} files to data/soybase/", file=sys.stderr)

    files: list[tuple[str, Path]] = []
    for entry in args.gff3:
        if "=" not in entry:
            print(f"error: --gff3 must be ASSEMBLY=PATH, got {entry!r}", file=sys.stderr)
            return 2
        assembly, path_text = entry.split("=", 1)
        files.append((assembly, Path(path_text)))
    if not files:
        print("error: at least one --gff3 is required", file=sys.stderr)
        return 2

    try:
        table, maxima = build_table(files)
    except DataContractError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    write_table(table, args.out)
    print(f"{len(table)} markers written to {args.out}")

    if args.emit_markers_csv:
        if args.markers_out is None:
            print("error: --emit-markers-csv requires --markers-out", file=sys.stderr)
            return 2
        write_markers_csv(table, args.emit_markers_csv, args.markers_out)
        print(f"markers written to {args.markers_out}")

    for assembly in ASSEMBLIES:
        if assembly not in maxima:
            continue
        print(f"{assembly} max position per chromosome:")
        for chrom in sorted(maxima[assembly], key=chrom_sort_key):
            print(f"  {chrom}: {maxima[assembly][chrom]}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
