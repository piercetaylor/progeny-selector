"""Tests for scripts/kasp_to_wide.py on hand-built KASP exports (synthetic; no real KASP data)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from progeny_selector.core.classify import classify, state_labels
from progeny_selector.model.dataset import Sample

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "kasp_to_wide.py"


def _load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location("kasp_to_wide", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


kasp_to_wide = _load_module()

MARKERS_CSV = "marker_id,chrom,pos_bp\nsnp.1,Gm01,1000\nsnp.2,Gm01,2000\n"

# 10-row long-format export: two samples (RP, DONOR) x snp.1/snp.2, one NTC control row,
# one Uncallable call, one A:G het call.
LONG_KASP = """SubjectID,SNPID,Call
RP,snp.1,A:A
RP,snp.2,A:A
DONOR,snp.1,G:G
DONOR,snp.2,G:G
PROGENY1,snp.1,A:G
PROGENY1,snp.2,Uncallable
NTC,snp.1,A:A
NTC,snp.2,A:A
PROGENY2,snp.1,?
PROGENY2,snp.2,G:G
"""


def _write(tmp_path: Path, name: str, text: str) -> Path:
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def test_long_format_loads_and_classifies(tmp_path, capsys):
    markers = _write(tmp_path, "markers.csv", MARKERS_CSV)
    kasp = _write(tmp_path, "export.csv", LONG_KASP)
    out = tmp_path / "genotypes.csv"

    rc = kasp_to_wide.main(["--kasp", str(kasp), "--markers", str(markers), "--out", str(out)])
    assert rc == 0
    # Uncallable is counted under its own reason, distinct from a generic "unrecognized" bucket
    stderr = capsys.readouterr().err
    assert "missing (Uncallable): 1" in stderr
    assert "missing (?): 1" in stderr

    from progeny_selector.io import load_genotypes

    genotype_matrix = load_genotypes(str(out))
    # NTC dropped as a control; RP, DONOR, PROGENY1, PROGENY2 remain, first-seen order
    assert genotype_matrix.sample_ids == ["RP", "DONOR", "PROGENY1", "PROGENY2"]

    samples = [
        Sample(sample_id="RP", line_name="RP", role="recurrent_parent", generation=None, family_id=None, notes=None),
        Sample(sample_id="DONOR", line_name="DONOR", role="donor_parent", generation=None, family_id=None, notes=None),
        Sample(sample_id="PROGENY1", line_name="PROGENY1", role="progeny", generation=None, family_id=None, notes=None),
        Sample(sample_id="PROGENY2", line_name="PROGENY2", role="progeny", generation=None, family_id=None, notes=None),
    ]
    from progeny_selector.io.manifest import build_dataset

    dataset = build_dataset(genotype_matrix, samples)
    classification = classify(dataset.genotypes, "RP", "DONOR")
    labels = state_labels(classification.states)
    marker_index = {m.marker_id: i for i, m in enumerate(dataset.genotypes.markers)}
    sample_index = {s: i for i, s in enumerate(dataset.genotypes.sample_ids)}

    # PROGENY1 x snp.1 was A:G -> heterozygous -> H
    assert labels[marker_index["snp.1"], sample_index["PROGENY1"]] == "H"
    # PROGENY1 x snp.2 was Uncallable -> missing -> N
    assert labels[marker_index["snp.2"], sample_index["PROGENY1"]] == "N"
    # PROGENY2 x snp.1 was '?' -> missing -> N
    assert labels[marker_index["snp.1"], sample_index["PROGENY2"]] == "N"


def test_conflicting_duplicate_rows_raise(tmp_path):
    markers = _write(tmp_path, "markers.csv", MARKERS_CSV)
    conflicting = LONG_KASP + "PROGENY1,snp.1,G:G\n"  # PROGENY1 x snp.1 was A:G above, now G:G
    kasp = _write(tmp_path, "export_conflict.csv", conflicting)
    out = tmp_path / "genotypes.csv"

    rc = kasp_to_wide.main(["--kasp", str(kasp), "--markers", str(markers), "--out", str(out)])
    assert rc == 2


def test_unplaced_snp_raises_without_drop_unplaced_and_is_dropped_with_it(tmp_path):
    markers = _write(tmp_path, "markers.csv", MARKERS_CSV)
    with_unplaced = LONG_KASP.replace("PROGENY2,snp.2,G:G", "PROGENY2,snp.3,G:G")
    kasp = _write(tmp_path, "export_unplaced.csv", with_unplaced)
    out = tmp_path / "genotypes.csv"

    rc = kasp_to_wide.main(["--kasp", str(kasp), "--markers", str(markers), "--out", str(out)])
    assert rc == 2

    rc = kasp_to_wide.main(["--kasp", str(kasp), "--markers", str(markers), "--out", str(out), "--drop-unplaced"])
    assert rc == 0
    text = out.read_text(encoding="utf-8")
    assert "snp.3" not in text


# --- grid shape (SNPviewer export), required by decision 9 though not listed in the Phase 10 table row ---

GRID_MARKERS_AS_COLUMNS = """SampleID,snp.1,snp.2
RP,A:A,A:A
DONOR,G:G,G:G
PROGENY1,A:G,Missing
"""

GRID_MARKERS_AS_ROWS = """MarkerID,RP,DONOR,PROGENY1
snp.1,A:A,G:G,A:G
snp.2,A:A,G:G,Missing
"""


def test_grid_markers_as_columns(tmp_path):
    markers = _write(tmp_path, "markers.csv", MARKERS_CSV)
    kasp = _write(tmp_path, "grid_cols.csv", GRID_MARKERS_AS_COLUMNS)
    out = tmp_path / "genotypes.csv"

    rc = kasp_to_wide.main(["--kasp", str(kasp), "--markers", str(markers), "--out", str(out)])
    assert rc == 0
    text = out.read_text(encoding="utf-8")
    lines = {line.split(",")[0]: line for line in text.splitlines()}
    assert lines["snp.1"].split(",")[3:] in (["AA", "GG", "AG"], ["AA", "GG", "GA"])
    assert lines["snp.2"].split(",")[3:][-1] == "N"


def test_grid_markers_as_rows(tmp_path):
    markers = _write(tmp_path, "markers.csv", MARKERS_CSV)
    kasp = _write(tmp_path, "grid_rows.csv", GRID_MARKERS_AS_ROWS)
    out = tmp_path / "genotypes.csv"

    rc = kasp_to_wide.main(["--kasp", str(kasp), "--markers", str(markers), "--out", str(out)])
    assert rc == 0
    text = out.read_text(encoding="utf-8")
    lines = {line.split(",")[0]: line for line in text.splitlines()}
    assert lines["snp.2"].split(",")[3:][-1] == "N"
