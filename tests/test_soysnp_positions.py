"""Tests for scripts/soysnp_positions.py on hand-built GFF3 snippets (synthetic; no real SoyBase data)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from progeny_selector.io.manifest import read_markers

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "soysnp_positions.py"


def _load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location("soysnp_positions", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


soysnp_positions = _load_module()

# Synthetic 6-line GFF3 snippets, one for SoySNP50K on Wm82.a2 and one for SoySNP50K on Wm82.a4.
# marker.1 appears in both; marker.only_a2 appears only in a2; marker.scaffold is on a scaffold.
A2_GFF3 = """\
##gff-version 3
glyma.Wm82.gnm2.Gm01\tsoybase\tSNP\t1000\t1000\t.\t+\t.\tID=glyma.Wm82.gnm2.ss715579458;Name=marker.1
glyma.Wm82.gnm2.Gm02\tsoybase\tSNP\t2000\t2000\t.\t+\t.\tID=glyma.Wm82.gnm2.ss715579459;Name=marker.only_a2
glyma.Wm82.gnm2.scaffold_12\tsoybase\tSNP\t50\t50\t.\t+\t.\tID=glyma.Wm82.gnm2.ss715579460;Name=marker.scaffold
"""

A4_GFF3 = """\
##gff-version 3
glyma.Wm82.gnm4.Gm01\tsoybase\tSNP\t1100\t1100\t.\t+\t.\tID=glyma.Wm82.gnm4.ss715579458;Name=marker.1
"""


def test_build_table_joins_on_marker_id_and_leaves_blanks(tmp_path):
    a2_path = tmp_path / "glyma.Wm82.gnm2.mrk.SoySNP50K.gff3"
    a4_path = tmp_path / "glyma.Wm82.gnm4.mrk.SoySNP50K.gff3"
    a2_path.write_text(A2_GFF3, encoding="utf-8")
    a4_path.write_text(A4_GFF3, encoding="utf-8")

    table, maxima = soysnp_positions.build_table([("Wm82.a2", a2_path), ("Wm82.a4", a4_path)])

    assert table["marker.1"]["chrom_Wm82.a2"] == "Gm01"
    assert table["marker.1"]["pos_bp_Wm82.a2"] == 1000
    assert table["marker.1"]["chrom_Wm82.a4"] == "Gm01"
    assert table["marker.1"]["pos_bp_Wm82.a4"] == 1100
    assert table["marker.1"]["in_SoySNP50K"] is True
    assert table["marker.1"]["in_SoySNP6K"] is False

    # only in a2: a4 cells are blank
    assert table["marker.only_a2"]["chrom_Wm82.a2"] == "Gm02"
    assert table["marker.only_a2"]["chrom_Wm82.a4"] == ""
    assert table["marker.only_a2"]["pos_bp_Wm82.a4"] == ""

    # scaffold kept with its own name, not folded into a chromosome
    assert table["marker.scaffold"]["chrom_Wm82.a2"] == "scaffold_12"

    # max-position sanity check excludes the scaffold
    assert maxima["Wm82.a2"] == {"Gm01": 1000, "Gm02": 2000}
    assert maxima["Wm82.a4"] == {"Gm01": 1100}


def test_write_table_full_pipeline(tmp_path):
    a2_path = tmp_path / "glyma.Wm82.gnm2.mrk.SoySNP50K.gff3"
    a2_path.write_text(A2_GFF3, encoding="utf-8")
    out = tmp_path / "soysnp_positions.csv"

    rc = soysnp_positions.main(["--gff3", f"Wm82.a2={a2_path}", "--out", str(out)])
    assert rc == 0
    header = out.read_text(encoding="utf-8").splitlines()[0]
    assert header == (
        "marker_id,in_SoySNP50K,in_SoySNP6K,"
        "chrom_Wm82.a1,pos_bp_Wm82.a1,chrom_Wm82.a2,pos_bp_Wm82.a2,"
        "chrom_Wm82.a4,pos_bp_Wm82.a4,chrom_Wm82.a5,pos_bp_Wm82.a5,chrom_Wm82.a6,pos_bp_Wm82.a6"
    )


def test_emit_markers_csv_excludes_scaffold_and_loads_through_read_markers(tmp_path):
    a4_path = tmp_path / "glyma.Wm82.gnm4.mrk.SoySNP50K.gff3"
    a4_gff3 = A4_GFF3 + "glyma.Wm82.gnm4.scaffold_9\tsoybase\tSNP\t77\t77\t.\t+\t.\tID=x;Name=marker.scaffold4\n"
    a4_path.write_text(a4_gff3, encoding="utf-8")
    out = tmp_path / "positions.csv"
    markers_out = tmp_path / "markers.csv"

    rc = soysnp_positions.main(
        ["--gff3", f"Wm82.a4={a4_path}", "--out", str(out), "--emit-markers-csv", "Wm82.a4", "--markers-out", str(markers_out)]
    )
    assert rc == 0

    text = markers_out.read_text(encoding="utf-8")
    assert "marker.scaffold4" not in text

    loaded = read_markers(markers_out)
    assert loaded["marker.1"].chrom == "Gm01"
    assert loaded["marker.1"].pos_bp == 1100


def test_name_attribute_required(tmp_path):
    bad = tmp_path / "bad.gff3"
    bad.write_text("glyma.Wm82.gnm2.Gm01\tsoybase\tSNP\t1\t1\t.\t+\t.\tID=only_id\n", encoding="utf-8")
    from progeny_selector.model.dataset import DataContractError

    try:
        soysnp_positions.parse_gff3(bad)
        raise AssertionError("expected DataContractError")
    except DataContractError:
        pass
