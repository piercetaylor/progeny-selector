"""Token profiles (contract 1.4.0; src/progeny_selector/io/profiles.py and the profile paths of calls.py).

For each contract/profiles/*.json the Python literal in BUILTIN_PROFILES equals the file (the equality that
keeps the runtime copy honest, since Shinylive stages the package and not contract/), and the schema rules
reject what backcross tests/profiles.test.ts rejects.
"""

from __future__ import annotations

import csv
import dataclasses
import json
import subprocess
import sys
from pathlib import Path

import pytest

from progeny_selector.app.screens.load import profile_select_tag, read_profile_json
from progeny_selector.cli import _profile_ref
from progeny_selector.core.pipeline import run_analysis
from progeny_selector.io import load_dataset, load_genotypes, read_criteria
from progeny_selector.io.calls import parse_nucleotide_call
from progeny_selector.io.export import results_csv_text, selection_csv_text
from progeny_selector.io.hapmap import read_hapmap
from progeny_selector.io.profiles import BUILTIN_PROFILES, compile_profile, profile_label, resolve_profile, validate_profile
from progeny_selector.io.wide_csv import read_wide_csv
from progeny_selector.model.dataset import DataContractError

PROFILES = Path(__file__).resolve().parent.parent / "contract" / "profiles"
PROFILE_FILES = sorted(PROFILES.glob("*.json"))

VALID = {
    "id": "mylab",
    "name": "My lab",
    "description": "",
    "base": "nucleotide",
    "homozygous": {},
    "heterozygous": {"H": "*"},
    "missing": ["U"],
    "sources": [],
}


def test_five_builtins() -> None:
    assert [p.stem for p in PROFILE_FILES] == sorted(BUILTIN_PROFILES)


@pytest.mark.parametrize("path", PROFILE_FILES, ids=lambda p: p.stem)
def test_literal_equals_file(path: Path) -> None:
    obj = json.loads(path.read_bytes().decode("utf-8"))
    assert obj["id"] == path.stem
    assert obj == dataclasses.asdict(BUILTIN_PROFILES[obj["id"]])
    assert validate_profile(obj) == BUILTIN_PROFILES[obj["id"]]


def test_accepts_valid_custom() -> None:
    assert validate_profile(VALID).id == "mylab"


@pytest.mark.parametrize(
    ("change", "pattern"),
    [
        ({"id": "My Lab"}, r"id"),
        ({"base": None}, r"base"),
        ({"missing": ["U", "h"]}, r'"heterozygous" and "missing"'),
        ({"homozygous": {"Z": "*"}}, r'"\*"'),
        ({"base": "none", "homozygous": {"0": "0"}, "heterozygous": {}, "missing": ["-"]}, r'""'),
    ],
)
def test_rejections(change: dict, pattern: str) -> None:
    obj = {**VALID, **change}
    if change.get("base", "") is None:
        del obj["base"]
    with pytest.raises(DataContractError, match=rf"^token profile: .*{pattern}"):
        validate_profile(obj)


def test_resolve_and_label() -> None:
    assert resolve_profile(None) is None
    assert resolve_profile("default") is None
    assert profile_label(None) == "default"
    with pytest.raises(DataContractError, match="built-in profiles are tassel, soybase-report, dart, axiom, kasp"):
        resolve_profile("nope")
    assert profile_label(BUILTIN_PROFILES["dart"]) == "dart"
    assert profile_label(validate_profile(VALID)) == "custom:mylab"


def test_axiom_cells() -> None:
    axiom = compile_profile(BUILTIN_PROFILES["axiom"])
    assert parse_nucleotide_call("ab", profile=axiom) == ("A", "B")
    assert parse_nucleotide_call("nocall", profile=axiom) is None
    with pytest.raises(ValueError, match="unrecognised nucleotide call"):
        parse_nucleotide_call("A", profile=axiom)


# Review findings, backcross docs/adr/0019 amendment of 2026-09-16.

WIDE_HEADER = "marker_id,chrom,pos_bp,RP,DONOR,L1\n"
HAPMAP_HEADER = "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\tRP\tDONOR\tL1\n"
HAPMAP_FIXED = "+\tNA\tNA\tNA\tNA\tNA\tNA"


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_bytes(text.encode("utf-8"))
    return path


def test_claimed_tokens_do_not_decide_coding(tmp_path: Path) -> None:
    path = _write(tmp_path, "g.csv", WIDE_HEADER + "r1,Gm02,1000,A,A,H\n")
    assert read_wide_csv(path).coded is True
    # Under soybase-report H does not vote, so the file is nucleotide and H at a one-allele marker is ambiguous.
    with pytest.raises(DataContractError, match="heterozygote token but the marker shows 1 allele"):
        read_wide_csv(path, profile=BUILTIN_PROFILES["soybase-report"])


def test_profile_on_coded_file_is_an_error(tmp_path: Path) -> None:
    path = _write(tmp_path, "g.csv", WIDE_HEADER + "r1,Gm02,1000,A,B,H\nr2,Gm02,2000,A,B,B\n")
    with pytest.raises(
        DataContractError,
        match=(
            r'^token profile "soybase-report" applies to HapMap and wide CSV nucleotide calls; '
            r"the genotype file is coded A/B/H \(detected\)$"
        ),
    ):
        read_wide_csv(path, profile=BUILTIN_PROFILES["soybase-report"])
    with pytest.raises(DataContractError, match=r'"dart" applies to HapMap and wide CSV .*\(requested\)'):
        read_wide_csv(path, coding="abh", profile=BUILTIN_PROFILES["dart"])


def test_custom_file_with_builtin_id_is_custom() -> None:
    copy = validate_profile(json.loads(json.dumps(dataclasses.asdict(BUILTIN_PROFILES["dart"]))))
    assert profile_label(copy) == "custom:dart"
    assert profile_label(resolve_profile(dataclasses.asdict(BUILTIN_PROFILES["dart"]))) == "custom:dart"
    assert profile_label(resolve_profile("dart")) == "dart"


@pytest.mark.parametrize(
    ("change", "pattern"),
    [
        ({"missing": ["U", "u"]}, r'appears twice in "missing"'),
        ({"homozygous": {"aa": "A", "AA": "A"}}, r'appears twice in "homozygous"'),
        ({"heterozygous": {"AB": ["A", "A"]}}, r"identical"),
        ({"tokens": {}}, r'unknown key\(s\) "tokens"'),
    ],
)
def test_more_rejections(change: dict, pattern: str) -> None:
    with pytest.raises(DataContractError, match=rf"^token profile: .*{pattern}"):
        validate_profile({**VALID, **change})


def test_hapmap_alleles_upper_cased(tmp_path: Path) -> None:
    path = _write(tmp_path, "g.hmp.txt", HAPMAP_HEADER + f"r1\ta/g\tGm02\t1000\t{HAPMAP_FIXED}\tAA\tAA\tH\n")
    gm = read_hapmap(path, profile=BUILTIN_PROFILES["soybase-report"])
    col = gm.sample_index("L1")
    assert sorted(gm.alleles[0][int(i)] for i in gm.calls[0, col]) == ["A", "G"]


def test_hapmap_indel_allele_rejects_h(tmp_path: Path) -> None:
    path = _write(tmp_path, "g.hmp.txt", HAPMAP_HEADER + f"r1\tA/-\tGm02\t1000\t{HAPMAP_FIXED}\tAA\tAA\tH\n")
    with pytest.raises(DataContractError, match="heterozygote token but the marker shows an indel allele"):
        read_hapmap(path, profile=BUILTIN_PROFILES["soybase-report"])


def test_token_profile_in_results_rows_and_selected_csv(fixture_dir: Path) -> None:
    dataset = load_dataset(fixture_dir / "genotypes.vcf", fixture_dir / "samples.csv", fixture_dir / "markers.csv")
    result = run_analysis(dataset, read_criteria(fixture_dir / "criteria.yaml"))
    assert result.rows
    for row in result.rows:
        assert list(row)[-1] == "token_profile"
        assert row["token_profile"] == "default"
    assert results_csv_text(result.rows).split("\r\n")[0].endswith(",token_profile")
    selected = selection_csv_text(result.rows[:2]).split("\r\n")
    assert selected[0].endswith(",notes,results_schema,token_profile")
    assert all(line.endswith(",default") for line in selected[1:3])


def test_profile_select_disabled_while_custom_file_set() -> None:
    assert 'disabled="disabled"' in str(profile_select_tag("default", True))
    assert "disabled" not in str(profile_select_tag("dart", False))


CASES = Path(__file__).resolve().parent.parent / "contract" / "cases"
MINIMAL_CRITERIA = """name: profile case
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


def test_hapmap_abh_with_profile_is_an_error(tmp_path: Path) -> None:
    """A profile on a HapMap file read as A/B/H is genotypes.profile_format, as for wide CSV."""
    path = _write(tmp_path, "g.hmp.txt", HAPMAP_HEADER + f"r1\ta/g\tGm02\t1000\t{HAPMAP_FIXED}\tAA\tAA\tH\n")
    with pytest.raises(
        DataContractError, match=r"applies to HapMap and wide CSV nucleotide calls; the genotype file is coded A/B/H \(requested\)"
    ):
        load_genotypes(path, coding="abh", profile="soybase-report")
    with pytest.raises(DataContractError, match="applies to HapMap and wide CSV"):
        read_hapmap(path, coding="abh", profile=BUILTIN_PROFILES["soybase-report"])


def test_trailing_newline_rejected_in_id_and_symbol() -> None:
    """The patterns are anchored on the whole string, so a trailing newline is not a valid id or symbol."""
    with pytest.raises(DataContractError, match=r'^token profile: "id" must match'):
        validate_profile({**VALID, "id": "dart\n"})
    with pytest.raises(DataContractError, match="must map to a non-blank allele symbol"):
        validate_profile({**VALID, "homozygous": {"0": "A\n"}})
    with pytest.raises(DataContractError, match="must map to two allele symbols"):
        validate_profile({**VALID, "heterozygous": {"2": ["A\n", "T"]}})


def test_profile_ref_file_branch(tmp_path: Path) -> None:
    assert _profile_ref(None) is None
    assert _profile_ref("dart") == "dart"
    path = tmp_path / "p.json"
    path.write_text(json.dumps(VALID), encoding="utf-8")
    assert _profile_ref(str(path)) == VALID
    with pytest.raises(DataContractError, match="token profile file"):
        _profile_ref(str(tmp_path / "absent.json"))
    bad = tmp_path / "bad.json"
    bad.write_text("[1]", encoding="utf-8")
    with pytest.raises(DataContractError, match="must be a JSON object"):
        _profile_ref(str(bad))


def test_read_profile_json_for_the_load_screen(tmp_path: Path) -> None:
    path = tmp_path / "p.json"
    path.write_text(json.dumps(VALID), encoding="utf-8")
    assert read_profile_json(str(path)) == VALID
    with pytest.raises(DataContractError, match="token profile file"):
        read_profile_json(str(tmp_path / "absent.json"))
    notjson = tmp_path / "bad.json"
    notjson.write_text("nope", encoding="utf-8")
    with pytest.raises(DataContractError, match="token profile file"):
        read_profile_json(str(notjson))
    arr = tmp_path / "arr.json"
    arr.write_text("[1]", encoding="utf-8")
    with pytest.raises(DataContractError, match="must be a JSON object"):
        read_profile_json(str(arr))


def test_load_dataset_records_the_profile_from_an_id_and_from_a_dict() -> None:
    case = CASES / "profile-dart-wide"
    by_id = load_dataset(case / "genotypes.csv", case / "samples.csv", profile="dart")
    assert by_id.token_profile == "dart"
    obj = json.loads((PROFILES / "dart.json").read_bytes().decode("utf-8"))
    by_dict = load_dataset(case / "genotypes.csv", case / "samples.csv", profile=obj)
    assert by_dict.token_profile == "custom:dart"


def test_cli_rank_and_select_carry_the_profile(tmp_path: Path) -> None:
    """The profile the CLI read reaches the last cell of results.csv and of selected.csv."""
    case = CASES / "profile-dart-wide"
    criteria = tmp_path / "criteria.yaml"
    criteria.write_text(MINIMAL_CRITERIA, encoding="utf-8")
    results, selected = tmp_path / "results.csv", tmp_path / "selected.csv"
    env = {"PYTHONPATH": str(Path(__file__).resolve().parents[1] / "src")}
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
        "--profile",
        "dart",
        "--out",
        str(results),
    ]
    proc = subprocess.run(rank, capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stderr
    rows = list(csv.reader(results.read_text(encoding="utf-8").splitlines()))
    assert rows[0][-1] == "token_profile"
    assert [r[-1] for r in rows[1:]] == ["dart"] * (len(rows) - 1)

    sel = [sys.executable, "-m", "progeny_selector", "select", "--results", str(results), "--top", "1", "--out", str(selected)]
    proc = subprocess.run(sel, capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stderr
    chosen = list(csv.reader(selected.read_text(encoding="utf-8").splitlines()))
    assert chosen[0][-1] == "token_profile"
    assert len(chosen) > 1
    assert [r[-1] for r in chosen[1:]] == ["dart"] * (len(chosen) - 1)
