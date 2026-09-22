"""Crop chromosome schemes (contract 1.5.0; src/progeny_selector/io/crops.py, core/chrom.py).

For each contract/crops/*.json the Python literal in BUILTIN_CROPS equals the file (the equality that
keeps the runtime copy honest, since Shinylive stages the package and not contract/), the key
derivation and the ordering match backcross src/core/chromosomes.ts, and the chosen scheme reaches
every reader, ``read_markers``, ``Dataset.crop`` and both exports.

Every threading test loads a NON-soybean crop: soybean is the default, so a soybean dataset cannot
tell a threaded scheme apart from the fallback.
"""

from __future__ import annotations

import csv
import dataclasses
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

from progeny_selector.app.screens.load import CROP_CHOICES
from progeny_selector.constants import DEFAULT_ASSEMBLY, RESULTS_SCHEMA, STATE_A, STATE_B, STATE_H
from progeny_selector.core.chrom import SOYBEAN, CropScheme, chrom_length_bp, chrom_sort_key, compile_scheme, normalize_chrom
from progeny_selector.core.classify import classify
from progeny_selector.core.pipeline import run_analysis
from progeny_selector.core.strip import chromosome_strips
from progeny_selector.io import load_dataset, read_criteria
from progeny_selector.io.crops import BUILTIN_CROPS, DEFAULT_CROP_ID, resolve_crop, validate_scheme
from progeny_selector.io.export import results_csv_text, selection_csv_text
from progeny_selector.io.hapmap import read_hapmap
from progeny_selector.io.manifest import read_markers
from progeny_selector.io.vcf import read_vcf
from progeny_selector.io.wide_csv import read_wide_csv
from progeny_selector.model.dataset import DataContractError, Dataset

ROOT = Path(__file__).resolve().parent.parent
CROPS = ROOT / "contract" / "crops"
CASES = ROOT / "contract" / "cases"
CROP_FILES = sorted(CROPS.glob("*.json"))

MAIZE = resolve_crop("maize")
RICE = resolve_crop("rice")
COTTON = resolve_crop("cotton")
OAT = resolve_crop("oat")

MINIMAL_CRITERIA = """name: crop case
targets:
  - locus_id: T1
    marker_id: r1
    required_state: either
background:
  model: count
  map_unit: bp
filters:
  max_missing_rate: 0.5
"""


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_bytes(text.encode("utf-8"))
    return path


# --- the nine files -------------------------------------------------------------------------


def test_nine_builtins() -> None:
    assert [p.stem for p in CROP_FILES] == sorted(BUILTIN_CROPS)
    assert list(BUILTIN_CROPS) == ["soybean", "maize", "rice", "sorghum", "wheat", "barley", "oat", "common-bean", "cotton"]


@pytest.mark.parametrize("path", CROP_FILES, ids=lambda p: p.stem)
def test_literal_equals_file(path: Path) -> None:
    obj = json.loads(path.read_bytes().decode("utf-8"))
    assert obj["id"] == path.stem
    assert obj == dataclasses.asdict(BUILTIN_CROPS[obj["id"]])
    assert validate_scheme(obj) == BUILTIN_CROPS[obj["id"]]


@pytest.mark.parametrize("crop_id", list(BUILTIN_CROPS), ids=list(BUILTIN_CROPS))
def test_every_canonical_name_normalises_to_itself(crop_id: str) -> None:
    scheme = resolve_crop(crop_id)
    for name in scheme.scheme.chromosomes:
        assert normalize_chrom(name, scheme) == name


# --- key derivation and ordering -------------------------------------------------------------


def test_leading_zeros_are_stripped_from_each_digit_run() -> None:
    assert normalize_chrom("01", MAIZE) == "chr1"
    assert normalize_chrom("chr_010", MAIZE) == "chr10"
    assert normalize_chrom("A01", COTTON) == "A01"
    assert normalize_chrom("chrA1", COTTON) == "A01"
    assert normalize_chrom("chr01", RICE) == "Chr1"
    # Every shipped pattern absorbs the zeros in `0*`, so the key-level strip needs a pattern that does not.
    padded = compile_scheme(
        CropScheme(
            id="padded",
            name="Padded",
            species="Genus species",
            ploidy=2,
            assembly="v1",
            chromosomes=["chr1", "chr10"],
            keys=["1", "10"],
            pattern=r"^chr([0-9]+)$",
            sources=[],
        )
    )
    assert normalize_chrom("chr001", padded) == "chr1"
    assert normalize_chrom("chr010", padded) == "chr10"


def test_keys_are_derived_case_insensitively() -> None:
    """The key is upper-cased before the lookup, so a lower-case letter group still matches."""
    assert normalize_chrom("chra1", COTTON) == "A01"
    assert normalize_chrom("d13", COTTON) == "D13"
    assert normalize_chrom("1h", resolve_crop("barley")) == "chr1H"
    assert normalize_chrom("chr7d", resolve_crop("wheat")) == "Chr7D"


def test_lg_prefix_is_soybean_only() -> None:
    assert normalize_chrom("LG7", SOYBEAN) == "Gm07"
    assert normalize_chrom("LG7", MAIZE) == "LG7"


def test_unmatched_names_are_kept_and_ordered_after_the_canonical_names() -> None:
    for name in ("ChrUn", "chrUn", "Un0", "chr00", "MT", "Pltd"):
        assert normalize_chrom(name, MAIZE) == name
    assert sorted(["scaffold_10", "chr2", "scaffold_2", "chr10"], key=lambda c: chrom_sort_key(c, MAIZE)) == [
        "chr2",
        "chr10",
        "scaffold_2",
        "scaffold_10",
    ]
    # Under soybean's key an unplaced "Un0" would sort before "chr1A" ("U" < "c"); under oat's it does not.
    assert sorted(["Un0", "chr7D", "chr1A"], key=lambda c: chrom_sort_key(c, OAT)) == ["chr1A", "chr7D", "Un0"]


def test_resolve_crop_defaults_to_soybean_and_names_the_built_ins() -> None:
    assert resolve_crop(None) is SOYBEAN
    assert resolve_crop(DEFAULT_CROP_ID) is SOYBEAN
    with pytest.raises(DataContractError, match=r'unknown crop "cowpea"; built-in crops are soybean, maize'):
        resolve_crop("cowpea")


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"id": "Maize"}, '"id" must match'),
        ({"name": " "}, '"name" must be a non-empty string'),
        ({"species": ""}, '"species" must be a non-empty string'),
        ({"ploidy": 0}, '"ploidy" must be a positive integer'),
        ({"assembly": ""}, '"assembly" must be a non-empty string'),
        ({"keys": ["1"]}, '"chromosomes" has 2 entries and "keys" 1'),
        ({"chromosomes": ["chr1", "chr1"]}, 'chromosome "chr1" appears twice'),
        ({"chromosomes": ["chr1", "chr2"], "keys": ["1", "1"]}, 'key "1" appears twice'),
        ({"keys": ["1", "02"]}, "every key must be upper-case with no leading zeros"),
        ({"keys": ["1", "a"]}, "every key must be upper-case with no leading zeros"),
        ({"pattern": "^([1-2$"}, '"pattern" does not compile'),
        ({"sources": [1]}, '"sources" must be an array of strings'),
        ({"extra": 1}, 'unknown key(s) "extra"'),
    ],
)
def test_validate_scheme_rejects(change: dict, message: str) -> None:
    valid = {
        "id": "mycrop",
        "name": "My crop",
        "species": "Genus species",
        "ploidy": 2,
        "assembly": "v1",
        "chromosomes": ["chr1", "chr2"],
        "keys": ["1", "2"],
        "pattern": "^(?:chr)?([1-2])$",
        "sources": [],
    }
    assert validate_scheme(valid).id == "mycrop"
    with pytest.raises(DataContractError, match=re.escape(f"crop scheme: {message}")):
        validate_scheme({**valid, **change})


# --- the scheme reaches every reader ----------------------------------------------------------


def test_vcf_reader_uses_the_scheme(tmp_path: Path) -> None:
    text = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tRP\tDONOR\tL1\n"
        "1\t1000\tr1\tA\tG\t.\t.\t.\tGT\t0/0\t1/1\t0/1\n"
    )
    path = _write(tmp_path, "g.vcf", text)
    assert read_vcf(path, scheme=MAIZE).markers[0].chrom == "chr1"
    assert read_vcf(path).markers[0].chrom == "Gm01"


def test_hapmap_reader_uses_the_scheme(tmp_path: Path) -> None:
    header = "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\tRP\tDONOR\tL1\n"
    path = _write(tmp_path, "g.hmp.txt", header + "r1\tA/G\t1\t1000\t+\tNA\tNA\tNA\tNA\tNA\tNA\tAA\tGG\tAG\n")
    assert read_hapmap(path, scheme=MAIZE).markers[0].chrom == "chr1"
    assert read_hapmap(path).markers[0].chrom == "Gm01"


def test_wide_csv_reader_uses_the_scheme(tmp_path: Path) -> None:
    path = _write(tmp_path, "g.csv", "marker_id,chrom,pos_bp,RP,DONOR,L1\nr1,1,1000,A,G,R\n")
    assert read_wide_csv(path, scheme=MAIZE).markers[0].chrom == "chr1"
    assert read_wide_csv(path).markers[0].chrom == "Gm01"


def test_read_markers_uses_the_scheme(tmp_path: Path) -> None:
    path = _write(tmp_path, "markers.csv", "marker_id,chrom,pos_bp\nr1,1,1000\nr2,chr12,2000\n")
    assert [m.chrom for m in read_markers(path, scheme=RICE).values()] == ["Chr1", "Chr12"]
    assert [m.chrom for m in read_markers(path).values()] == ["Gm01", "Gm12"]


def test_load_dataset_threads_the_crop_to_genotypes_and_markers(tmp_path: Path) -> None:
    genotypes = _write(tmp_path, "g.csv", "marker_id,chrom,pos_bp,RP,DONOR,L1\nr1,1,1000,A,G,R\nr2,10,2000,A,G,R\n")
    samples = _write(
        tmp_path,
        "samples.csv",
        "sample_id,line_name,role,generation,family_id,notes\nRP,RP,recurrent_parent,,,\nDONOR,D,donor_parent,,,\nL1,L1,candidate,,,\n",
    )
    markers = _write(tmp_path, "markers.csv", "marker_id,chrom,pos_bp\nr1,Chromosome_1,1500\nr2,chr10,2500\n")
    dataset = load_dataset(genotypes, samples, markers, crop="maize")
    assert dataset.crop == "maize"
    assert [m.chrom for m in dataset.genotypes.markers] == ["chr1", "chr10"]
    assert load_dataset(genotypes, samples, markers).crop == "soybean"
    with pytest.raises(DataContractError, match='unknown crop "cowpea"'):
        load_dataset(genotypes, samples, markers, crop="cowpea")


# --- the crop reaches the exports ---------------------------------------------------------------


@pytest.fixture(scope="module")
def maize_result(tmp_path_factory):
    case = CASES / "crop-maize-spellings"
    criteria_path = tmp_path_factory.mktemp("crop") / "criteria.yaml"
    criteria_path.write_text(MINIMAL_CRITERIA, encoding="utf-8")
    dataset = load_dataset(case / "genotypes.csv", case / "samples.csv", crop="maize")
    return run_analysis(dataset, read_criteria(criteria_path))


def test_results_schema_is_1_1_0(maize_result) -> None:
    assert RESULTS_SCHEMA == "1.1.0"
    assert maize_result.rows
    assert all(row["results_schema"] == "1.1.0" for row in maize_result.rows)


def test_crop_column_sits_between_results_schema_and_token_profile(maize_result) -> None:
    header = results_csv_text(maize_result.rows).split("\r\n")[0].split(",")
    assert header.index("crop") == header.index("results_schema") + 1
    assert header[-1] == "token_profile"
    rows = list(csv.DictReader(results_csv_text(maize_result.rows).splitlines()))
    assert rows
    assert all(row["crop"] == "maize" for row in rows)


def test_crop_column_in_selected_csv(maize_result) -> None:
    lines = selection_csv_text(maize_result.rows[:1]).split("\r\n")
    assert lines[0].endswith(",notes,results_schema,crop,token_profile")
    assert lines[1].endswith(",1.1.0,maize,default")


def test_empty_results_header_carries_crop_before_token_profile() -> None:
    header = results_csv_text([]).split("\r\n")[0].split(",")
    assert header[-2:] == ["crop", "token_profile"]


def test_load_screen_offers_the_nine_crops() -> None:
    assert list(CROP_CHOICES) == list(BUILTIN_CROPS)
    assert CROP_CHOICES["common-bean"] == "Common bean"


def test_cli_rank_and_select_carry_the_crop(tmp_path: Path) -> None:
    """The crop the CLI read reaches the `crop` cell of results.csv and of selected.csv."""
    case = CASES / "crop-maize-spellings"
    criteria = tmp_path / "criteria.yaml"
    criteria.write_text(MINIMAL_CRITERIA, encoding="utf-8")
    results, selected = tmp_path / "results.csv", tmp_path / "selected.csv"
    env = {"PYTHONPATH": str(ROOT / "src")}
    rank = [
        sys.executable,
        "-m",
        "progeny_selector",
        "rank",
        "--genotypes",
        str(case / "genotypes.csv"),
        "--samples",
        str(case / "samples.csv"),
        "--criteria",
        str(criteria),
        "--crop",
        "maize",
        "--out",
        str(results),
    ]
    proc = subprocess.run(rank, capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stderr
    rows = list(csv.DictReader(results.read_text(encoding="utf-8").splitlines()))
    assert rows and all(row["crop"] == "maize" for row in rows)

    sel = [sys.executable, "-m", "progeny_selector", "select", "--results", str(results), "--top", "1", "--out", str(selected)]
    proc = subprocess.run(sel, capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stderr
    chosen = list(csv.DictReader(selected.read_text(encoding="utf-8").splitlines()))
    assert chosen and all(row["crop"] == "maize" for row in chosen)


def test_cli_rejects_an_unknown_crop(tmp_path: Path) -> None:
    case = CASES / "crop-maize-spellings"
    env = {"PYTHONPATH": str(ROOT / "src")}
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "progeny_selector",
            "validate",
            "--genotypes",
            str(case / "genotypes.csv"),
            "--samples",
            str(case / "samples.csv"),
            "--crop",
            "cowpea",
        ],
        capture_output=True,
        text=True,
        env=env,
    )
    assert proc.returncode == 2
    assert "cowpea" in proc.stderr


# --- the scheme reaches core: lengths, foreground, ordering -------------------------------------


def test_chrom_length_bp_is_soybean_only() -> None:
    """The length tables are Williams 82's, so every other crop gets the caller's fallback."""
    gm01 = chrom_length_bp("Gm01", None, DEFAULT_ASSEMBLY)
    assert gm01 == chrom_length_bp("chr1", None, DEFAULT_ASSEMBLY) == chrom_length_bp("Gm01", None, DEFAULT_ASSEMBLY, SOYBEAN)
    assert gm01 is not None
    assert chrom_length_bp("chr1", None, DEFAULT_ASSEMBLY, MAIZE) is None
    assert chrom_length_bp("chr1", 42, DEFAULT_ASSEMBLY, MAIZE) == 42
    assert chrom_length_bp("Chr1", None, DEFAULT_ASSEMBLY, RICE) is None
    assert chrom_length_bp("chr1A", None, DEFAULT_ASSEMBLY, OAT) is None
    # A soybean spelling inside a non-soybean dataset takes no length either: the table is soybean's.
    assert chrom_length_bp("Gm01", None, DEFAULT_ASSEMBLY, MAIZE) is None
    assert chrom_length_bp("Gm01", 42, DEFAULT_ASSEMBLY, OAT) == 42
    # The gate is identity with the built-in scheme, not its id: a look-alike must take no length.
    lookalike = compile_scheme(
        CropScheme(
            id="soybean",
            name="Not Williams 82",
            species="Genus species",
            ploidy=2,
            assembly="v1",
            chromosomes=["Gm01"],
            keys=["1"],
            pattern=r"^Gm0*([1-9])$",
            sources=[],
        )
    )
    assert chrom_length_bp("Gm01", None, DEFAULT_ASSEMBLY, lookalike) is None


MAIZE_GENOTYPES = """marker_id,chrom,pos_bp,RP,DONOR,L1,L2
m1,chr2,1000000,A,G,R,A
m2,chr7,20000000,A,G,R,R
m3,chr7,21000000,A,G,R,R
m4,chr10,200000000,A,G,A,A
"""
MAIZE_SAMPLES = """sample_id,line_name,role,generation,family_id,notes
RP,Recurrent,recurrent_parent,,,
DONOR,Donor,donor_parent,,,
L1,Line 1,candidate,,,
L2,Line 2,candidate,,,
"""
MAIZE_CRITERIA = """name: maize region target
targets:
  - locus_id: T1
    chrom: chr7
    start_bp: 19000000
    end_bp: 22000000
    required_state: either
background:
  model: count
  map_unit: bp
filters:
  max_missing_rate: 0.5
"""


def _maize_dataset(tmp_path: Path):
    genotypes = _write(tmp_path, "g.csv", MAIZE_GENOTYPES)
    samples = _write(tmp_path, "samples.csv", MAIZE_SAMPLES)
    return load_dataset(genotypes, samples, crop="maize")


def test_a_maize_region_target_resolves_and_takes_no_soybean_length(tmp_path: Path) -> None:
    """End to end under maize: `chr7` stays `chr7` in the target, no chromosome takes a Williams 82 length."""
    criteria = _write(tmp_path, "criteria.yaml", MAIZE_CRITERIA)
    result = run_analysis(_maize_dataset(tmp_path), read_criteria(criteria))
    assert result.resolved_targets["T1"].chrom == "chr7"
    assert list(result.resolved_targets["T1"].marker_idx) != []
    assert [row["target_T1_status"] for row in result.rows if row["sample_id"] == "L1"] == ["pass"]
    # m4 sits at 200 Mb, past every Williams 82 chromosome: under soybean `chr10` would be `Gm10` and
    # the run would warn that the markers reach beyond the assembly length.
    assert not [w for w in result.warnings if "beyond" in w]
    assert [w for w in result.warnings if w.startswith("chromosome chr10: the assembly")]
    header = results_csv_text(result.rows).split("\r\n")[0].split(",")
    # Not a crop discriminator (soybean normalises "chr" names too); the oat test below is the one that is.
    assert header.index("rpp_chr2") < header.index("rpp_chr7") < header.index("rpp_chr10")
    assert all(row["crop"] == "maize" for row in result.rows)


OAT_GENOTYPES = """marker_id,chrom,pos_bp,RP,DONOR,L1,L2
m1,1A,1000000,A,G,R,A
m2,chr7D,2000000,A,G,R,R
m3,Un0,3000000,A,G,A,A
"""
OAT_CRITERIA = """name: oat marker target
targets:
  - locus_id: T1
    marker_id: m1
    required_state: either
background:
  model: count
  map_unit: bp
filters:
  max_missing_rate: 0.5
"""


def test_an_oat_dataset_orders_its_canonical_chromosomes_before_an_unplaced_contig(tmp_path: Path) -> None:
    """End to end under oat: `Un0` sorts after `chr1A` and `chr7D`; under soybean's key it would lead."""
    genotypes = _write(tmp_path, "g.csv", OAT_GENOTYPES)
    samples = _write(tmp_path, "samples.csv", MAIZE_SAMPLES)
    criteria = _write(tmp_path, "criteria.yaml", OAT_CRITERIA)
    dataset = load_dataset(genotypes, samples, crop="oat")
    assert [m.chrom for m in dataset.genotypes.markers] == ["chr1A", "chr7D", "Un0"]
    result = run_analysis(dataset, read_criteria(criteria))
    header = results_csv_text(result.rows).split("\r\n")[0].split(",")
    assert header.index("rpp_chr1A") < header.index("rpp_chr7D") < header.index("rpp_Un0")


def test_dataset_crop_is_derived_from_its_scheme(tmp_path: Path) -> None:
    dataset = _maize_dataset(tmp_path)
    assert dataset.scheme is MAIZE
    assert dataset.crop == "maize"
    assert Dataset(genotypes=dataset.genotypes, samples=dataset.samples).crop == "soybean"


def test_chromosome_strips_take_the_scheme(tmp_path: Path) -> None:
    """The Compare screen's strips order and end under the dataset's crop, not under soybean."""
    dataset = _maize_dataset(tmp_path)
    gm = dataset.genotypes.sorted_by_position(dataset.scheme)
    states = classify(gm, "RP", "DONOR").states[:, gm.sample_index("L1")]
    strips = chromosome_strips(states, gm, DEFAULT_ASSEMBLY, dataset.scheme)
    assert [s.chrom for s in strips] == ["chr2", "chr7", "chr10"]
    # chr2's markers end at 1 Mb; under soybean the strip would stretch to Gm02's 50.4 Mb.
    assert next(s.length_bp for s in strips if s.chrom == "chr2") == 1_000_000
    soybean_strips = chromosome_strips(states, gm, DEFAULT_ASSEMBLY)
    assert next(s.length_bp for s in soybean_strips if s.chrom == "chr2") > 1_000_000

    # Strip order is the scheme's too: under soybean's key the unplaced "Un0" would lead.
    oat = load_dataset(_write(tmp_path, "oat.csv", OAT_GENOTYPES), _write(tmp_path, "oat-samples.csv", MAIZE_SAMPLES), crop="oat")
    oat_gm = oat.genotypes.sorted_by_position(oat.scheme)
    oat_states = classify(oat_gm, "RP", "DONOR").states[:, oat_gm.sample_index("L1")]
    assert [s.chrom for s in chromosome_strips(oat_states, oat_gm, DEFAULT_ASSEMBLY, oat.scheme)] == ["chr1A", "chr7D", "Un0"]


# `chr1A` = A, `chr7D` = H, `Un0` = B for L1, so a row order that is not the scheme's is visible.
OAT_ALIGN_GENOTYPES = """marker_id,chrom,pos_bp,RP,DONOR,L1,L2
m1,1A,1000000,A,G,A,A
m2,chr7D,2000000,A,G,R,A
m3,Un0,3000000,A,G,G,A
"""


def test_pipeline_and_compare_sort_the_same_way_under_a_crop(tmp_path: Path) -> None:
    """run_analysis classifies its own sorted copy and the Compare screen sorts the loaded matrix again.

    Both must sort under the dataset's scheme: the screen indexes ``result.classification.states`` by
    row, so a pipeline sorted under soybean and a view sorted under oat would misalign the genotype
    strips silently. `Un0` leads under soybean's key and trails under oat's, so the two orders differ.
    """
    genotypes = _write(tmp_path, "g.csv", OAT_ALIGN_GENOTYPES)
    samples = _write(tmp_path, "samples.csv", MAIZE_SAMPLES)
    criteria = _write(tmp_path, "criteria.yaml", OAT_CRITERIA)
    dataset = load_dataset(genotypes, samples, crop="oat")
    result = run_analysis(dataset, read_criteria(criteria))

    view_gm = dataset.genotypes.sorted_by_position(dataset.scheme)  # app/screens/compare.py
    assert [m.chrom for m in view_gm.markers] == ["chr1A", "chr7D", "Un0"]
    states = result.classification.states[:, view_gm.sample_index("L1")]
    assert list(states) == [STATE_A, STATE_H, STATE_B]
    # The same rows, read through the view's own matrix: sorted under soybean the pipeline would give
    # [STATE_B, STATE_A, STATE_H] against these markers.
    assert [m.marker_id for m in view_gm.markers] == ["m1", "m2", "m3"]


def test_an_unanchored_pattern_searches_as_backcross_does() -> None:
    """`normalize_chrom` uses `re.search`, the semantics of backcross's `RegExp.exec`.

    Every shipped pattern anchors itself with `^...$`, and neither `validate_scheme` nor its
    TypeScript counterpart requires anchors, so an unanchored pattern must behave identically in
    both implementations: `re.match` would anchor it at the start and diverge.
    """
    unanchored = compile_scheme(
        CropScheme(
            id="unanchored",
            name="Unanchored",
            species="Genus species",
            ploidy=2,
            assembly="v1",
            chromosomes=["chr1", "chr2"],
            keys=["1", "2"],
            pattern=r"chr([12])",
            sources=[],
        )
    )
    assert normalize_chrom("chr1", unanchored) == "chr1"
    # `re.exec` in backcross finds "chr2" inside this name; `re.match` would not.
    assert normalize_chrom("scaffold_chr2", unanchored) == "chr2"
