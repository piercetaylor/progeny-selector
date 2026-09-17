"""Tests for scripts/song2016_map.py on a hand-built Table S1 export (the real table is never read here)."""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from progeny_selector.io import load_dataset

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "song2016_map.py"


def _load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location("song2016_map", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


song2016_map = _load_module()

HEADER = [
    "ss ID",
    "SNP ID",
    "Glyma1.01 Chromosome",
    "Glyma1.01 Coordinate",
    "Wm82.a2.v1 Chromosome",
    "Wm82.a2.v1 Coordinate",
    "Glyma1.01 Gene",
    "WP Linkage Group",
    "WP linkage position",
    "EW Linkage Group ",
    "EW linkage position",
]
# ss, SNP ID, a1 chrom, a1 coord, a2 chrom, a2 coord, WP LG, WP cM, EW LG, EW cM
ROWS = [
    ("1001.0", "BARC_A", "Gm01", "1000007", "Chr01", "1000000", "1.0", "0.0", "", ""),  # mapped
    ("1002", "BARC_B", "Gm01", "2000007", "Chr01", "2000000", "1.0", "5.0", "", ""),  # mapped
    ("1003", "BARC_C", "Gm01", "3000007", "Chr01", "3000000", "1.0", "4.0", "", ""),  # non-monotone
    ("1004", "BARC_D", "Gm01", "4000007", "Chr01", "4000000", "1", "10.0", "", ""),  # mapped
    ("1005", "BARC_E", "Gm01", "5000007", "Chr01", "5000000", "2.0", "12.0", "", ""),  # LG mismatch
    ("1006", "BARC_F", "Gm02", "1000007", "\nChr02", "1000000", "2.0", "1.0", "", ""),  # mapped, leading newline
    ("1007", "BARC_G", "Gm02", "3000007", "\nChr02 ", "3000000", "2.0", "3.0", "", ""),  # mapped
    ("1008", "BARC_H", "Gm02", "5000007", "Chr02", "5000000", "2.0", "9.0", "2", "7.0"),  # mapped
    ("1009", "BARC_I", "Gm01", "9000000", "", "", "1.0", "3.0", "", ""),  # no position on a2
    ("1010", "BARC_J", "", "", "\nscaffold_298 ", "5000", "1.0", "3.0", "", ""),  # scaffold
    ("1011", "BARC_K", "Gm02", "6000007", "Chr02", "6000000", "", "", "2.0", "11.0"),  # not on the WP map
    ("1012", "BARC_L", "Gm01", "6000007", "Chr01", "6000000", "1.0", "20.0", "", ""),  # mapped
    ("1013", "BARC_M", "Gm01", "7000007", "Chr01", "7000000", "D1a", "25.0", "", ""),  # LG not a number
]


def _write_table(path: Path, rows: list[tuple[str, ...]] = ROWS) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Table S1. SNP positions and genetic maps"])
        w.writerow(HEADER)
        for ss, snp, a1c, a1p, a2c, a2p, wp_lg, wp_cm, ew_lg, ew_cm in rows:
            w.writerow([ss, snp, a1c, a1p, a2c, a2p, "", wp_lg, wp_cm, ew_lg, ew_cm])


def _read_out(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _log_counts(path: Path) -> dict[str, str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    return dict(line.split(": ", 1) for line in lines if ": " in line)


def _run(tmp_path: Path, *extra: str, markers_in: str | None = None, rows: list[tuple[str, ...]] = ROWS) -> tuple[int, Path]:
    table = tmp_path / "table_s1.csv"
    _write_table(table, rows)
    out = tmp_path / "markers.csv"
    argv = ["--table", str(table), "--out", str(out), *extra]
    if markers_in is not None:
        path = tmp_path / "in.csv"
        path.write_text(markers_in, encoding="utf-8")
        argv += ["--markers-in", str(path)]
    return song2016_map.main(argv), out


def test_table_only(tmp_path: Path):
    code, out = _run(tmp_path)
    assert code == 0
    assert [(r["marker_id"], r["chrom"], r["pos_bp"], float(r["cm"])) for r in _read_out(out)] == [
        ("ss1001", "Gm01", "1000000", 0.0),
        ("ss1002", "Gm01", "2000000", 5.0),
        ("ss1004", "Gm01", "4000000", 10.0),
        ("ss1012", "Gm01", "6000000", 20.0),
        ("ss1006", "Gm02", "1000000", 1.0),
        ("ss1007", "Gm02", "3000000", 3.0),
        ("ss1008", "Gm02", "5000000", 9.0),
    ]
    counts = _log_counts(tmp_path / "markers.log")
    assert counts["rows read"] == "13"
    assert counts["no position on assembly"] == "2"
    assert counts["no linkage position on map"] == "1"
    assert counts["LG mismatch"] == "2"
    assert counts["non-monotone dropped"] == "1"
    assert counts["mapped"] == "7"
    assert counts["interpolated"] == "0"
    assert counts["map"] == "WP"


def test_assembly_a1_uses_a1_columns(tmp_path: Path):
    code, out = _run(tmp_path, "--assembly", "a1")
    assert code == 0
    rows = {r["marker_id"]: r["pos_bp"] for r in _read_out(out)}
    assert rows["ss1001"] == "1000007"
    assert rows["ss1008"] == "5000007"
    counts = _log_counts(tmp_path / "markers.log")
    # ss1009 is placed on a1 (9000000, cM 3) and dropped as non-monotone after cM 20; ss1010 has no a1 position
    assert (counts["no position on assembly"], counts["non-monotone dropped"], counts["mapped"]) == ("1", "2", "7")


def test_ew_map_uses_trailing_space_column(tmp_path: Path):
    code, out = _run(tmp_path, "--map", "EW")
    assert code == 0
    assert [(r["marker_id"], float(r["cm"])) for r in _read_out(out)] == [("ss1008", 7.0), ("ss1011", 11.0)]
    counts = _log_counts(tmp_path / "markers.log")
    assert counts["no linkage position on map"] == "9"
    assert counts["mapped"] == "2"


def test_markers_in_interpolated_clamped_unplaced(tmp_path: Path):
    code, out = _run(
        tmp_path,
        markers_in="marker_id,chrom,pos_bp,cm\n"
        "ss1001,Gm01,1000000,\n"
        "ss1003,Gm01,3000000,\n"
        "ss1005,Gm01,5000000,\n"
        "BARC_B,Gm01,2000000,\n"
        "x_a,Gm01,500000,\n"
        "x_b,Chr02,7000000,\n"
        "x_c,Gm03,100,\n"
        "x_d,Gm01,2500000,6.5\n",
    )
    assert code == 0
    rows = _read_out(out)
    cm = {r["marker_id"]: r["cm"] for r in rows}
    assert float(cm["ss1003"]) == pytest.approx(7.5)
    assert float(cm["ss1005"]) == pytest.approx(15.0)
    assert float(cm["BARC_B"]) == pytest.approx(5.0)  # matched by SNP ID to a mapped row
    assert float(cm["x_a"]) == pytest.approx(0.0)
    assert float(cm["x_b"]) == pytest.approx(9.0)
    assert float(cm["x_d"]) == pytest.approx(6.5)
    assert cm["x_c"] == ""
    assert len(cm) == 13  # 7 mapped + 7 markers-in ids not among them - ss1002 written once as BARC_B
    assert "ss1002" not in cm
    counts = _log_counts(tmp_path / "markers.log")
    assert (counts["interpolated"], counts["clamped"], counts["unplaced"], counts["cM kept from --markers-in"]) == ("4", "2", "1", "1")
    assert "cM mode will be disabled" in (tmp_path / "markers.log").read_text(encoding="utf-8")


def test_mixed_assembly_refused(tmp_path: Path):
    code, out = _run(tmp_path, markers_in="marker_id,chrom,pos_bp\nss1002,Gm01,2000001\n")
    assert code == 2
    assert not out.exists()
    # a table marker dropped by cleaning (ss1003, non-monotone) is still checked
    code, out = _run(tmp_path, markers_in="marker_id,chrom,pos_bp\nss1003,Gm01,3000001\n")
    assert code == 2
    assert not out.exists()


def test_snp_id_named_markers_on_wrong_assembly_refused(tmp_path: Path):
    # BARC names on a1 positions while the table is read on a2
    code, out = _run(tmp_path, markers_in="marker_id,chrom,pos_bp\nBARC_A,Gm01,1000007\nBARC_B,Gm01,2000007\n")
    assert code == 2
    assert not out.exists()


def test_markers_in_without_any_match_refused(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    code, out = _run(tmp_path, markers_in="marker_id,chrom,pos_bp\nfoo_1,Gm01,1500000\n")
    assert code == 2
    assert "no --markers-in marker matches Table S1 by ss ID or SNP ID; cannot check the assembly" in capsys.readouterr().err
    assert not out.exists()


def test_markers_in_cm_strict_decimal(tmp_path: Path):
    code, out = _run(tmp_path, markers_in="marker_id,chrom,pos_bp,cm\nss1001,Gm01,1000000,1_000\n")
    assert code == 2
    assert not out.exists()


def test_equal_bp_output_cm_non_decreasing(tmp_path: Path):
    rows = [
        ("1002", "S2", "", "", "Chr01", "1000000", "1", "0.0", "", ""),
        ("1001", "S1", "", "", "Chr01", "1000000", "1", "5.0", "", ""),
        ("1003", "S3", "", "", "Chr01", "2000000", "1", "6.0", "", ""),
    ]
    code, out = _run(tmp_path, rows=rows)
    assert code == 0
    written = _read_out(out)
    assert [(r["marker_id"], float(r["cm"])) for r in written] == [("ss1001", 2.5), ("ss1002", 2.5), ("ss1003", 6.0)]
    code, out = _run(tmp_path, rows=rows, markers_in="marker_id,chrom,pos_bp\nss1001,Gm01,1000000\nq,Gm01,1500000\n")
    assert code == 0
    written = _read_out(out)
    cms = [float(r["cm"]) for r in written]
    assert cms == sorted(cms)
    assert [r["marker_id"] for r in written] == ["ss1001", "ss1002", "q", "ss1003"]
    assert cms[2] == pytest.approx(4.25)  # halfway between the mean anchor 2.5 at 1e6 and 6.0 at 2e6


def test_no_repeated_marker_id(tmp_path: Path):
    code, out = _run(
        tmp_path,
        markers_in="marker_id,chrom,pos_bp,cm\nss1001,Gm01,1000000,0.5\nss1002,Gm01,2000000,\nBARC_D,Gm01,4000000,\nBARC_F,Gm02,1000000,1.5\n",
    )
    assert code == 0
    written = _read_out(out)
    ids = [r["marker_id"] for r in written]
    assert len(ids) == len(set(ids))
    cm = {r["marker_id"]: float(r["cm"]) for r in written}
    assert cm["ss1001"] == 0.5  # own cM kept over the table's
    assert cm["ss1002"] == 5.0  # table cM
    assert cm["BARC_D"] == 10.0 and "ss1004" not in cm  # matched by SNP ID: one row under the --markers-in id
    assert cm["BARC_F"] == 1.5 and "ss1006" not in cm


def test_nothing_kept_stops(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    rows = [(r[0], r[1], r[2], r[3], r[4], r[5], "7", r[7], r[8], r[9]) for r in ROWS]
    code, out = _run(tmp_path, rows=rows)
    assert code == 2
    assert "largest drop: LG mismatch 10" in capsys.readouterr().err
    assert not out.exists()


def test_lg_to_chrom():
    assert song2016_map.lg_to_chrom("1") == "Gm01"
    assert song2016_map.lg_to_chrom("1.0") == "Gm01"
    assert song2016_map.lg_to_chrom(" 15.0 ") == "Gm15"
    assert song2016_map.lg_to_chrom("Gm01") is None
    assert song2016_map.lg_to_chrom("D1a") is None
    assert song2016_map.lg_to_chrom("1.5") is None
    assert song2016_map.lg_to_chrom("21") is None
    assert song2016_map.lg_to_chrom("") is None


def test_written_map_enables_cm(tmp_path: Path):
    code, out = _run(tmp_path, markers_in="marker_id,chrom,pos_bp\nss1001,Gm01,1000000\nss1003,Gm01,3000000\nx_b,Gm02,7000000\n")
    assert code == 0
    geno = tmp_path / "genotypes.csv"
    geno.write_text(
        "marker_id,chrom,pos_bp,RP,DONOR,P1\nss1001,Gm01,1000000,AA,GG,AG\nss1003,Gm01,3000000,CC,TT,CC\nx_b,Gm02,7000000,AA,CC,AC\n",
        encoding="utf-8",
    )
    samples = tmp_path / "samples.csv"
    samples.write_text("sample_id,role\nRP,recurrent_parent\nDONOR,donor_parent\nP1,progeny\n", encoding="utf-8")
    dataset = load_dataset(geno, samples, out)
    assert dataset.genotypes.has_cm()
