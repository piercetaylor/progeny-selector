"""Contract 1.11.0: every text input is strict UTF-8 (docs/adr/0031).

A byte sequence that is not well-formed UTF-8 is a DataContractError (CriteriaError for
criteria.yaml) naming the file, the physical line (LF-counted, 1-based) and the 1-based byte
position of the ill-formed sequence's lead byte within that line. The files are written with
`write_bytes`, so the bytes under test are exactly the bytes on disk.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from progeny_selector.cli import _profile_ref, main
from progeny_selector.io import load_dataset, load_genotypes, read_criteria
from progeny_selector.io.delimited import invalid_utf8_message
from progeny_selector.io.manifest import read_markers, read_samples
from progeny_selector.model.criteria import CriteriaError
from progeny_selector.model.dataset import DataContractError
from tests.test_io import HAPMAP_FIXED, HAPMAP_HEADER

VCF_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
VCF_UTF8_ROWS = [
    "##fileformat=VCFv4.2",
    f"{VCF_HEADER}\tRP\tDONOR\tL1",
    "Gm06\t1000\tg1\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\t0/1",
]
WIDE = b"marker_id,chrom,pos_bp,RP,DONOR,L1\ne1,Gm04,1000,A,G,A/G\ne2,Gm04,2000,C,T,C\n"


def vcf_bad_id(bad: bytes) -> bytes:
    """The shared case's VCF_BAD_ID: the bad bytes sit in the ID of line 4, at byte position 12."""
    head = "\n".join(VCF_UTF8_ROWS).encode("utf-8") + b"\n"
    return head + b"Gm06\t2000\tg" + bad + b"2\tC\tT\t.\tPASS\t.\tGT\t0/0\t1/1\t1/1\n"


def write_vcf(path: Path, data: bytes) -> Path:
    if path.name.endswith((".gz", ".bgz")):
        with gzip.open(path, "wb") as fh:
            fh.write(data)
    else:
        path.write_bytes(data)
    return path


# (a) the Latin-1 byte through the genotype loader, plain, gzip and bgzip ---------------------------


@pytest.mark.parametrize("name", ["g.vcf", "g.vcf.gz", "g.vcf.bgz"])
def test_vcf_latin1_byte_names_line_and_position(tmp_path: Path, name: str) -> None:
    path = write_vcf(tmp_path / name, vcf_bad_id(b"\xe9"))
    with pytest.raises(DataContractError, match=r"line 4: not valid UTF-8 \(byte 0xE9 at position 12\)") as info:
        load_genotypes(path)
    assert str(info.value).startswith(f"{path}: line 4: ")
    assert info.value.__cause__ is None and info.value.__suppress_context__


# (b) overlong and surrogate encodings in the same slot ----------------------------------------------


@pytest.mark.parametrize(("bad", "lead"), [(b"\xc0\x80", "0xC0"), (b"\xed\xa0\x80", "0xED")])
def test_vcf_overlong_and_surrogate(tmp_path: Path, bad: bytes, lead: str) -> None:
    path = write_vcf(tmp_path / "g.vcf", vcf_bad_id(bad))
    with pytest.raises(DataContractError, match=rf"line 4: not valid UTF-8 \(byte {lead} at position 12\)"):
        load_genotypes(path)


# (c) a sequence cut short by the end of the file ----------------------------------------------------


def test_vcf_truncated_sequence_at_eof(tmp_path: Path) -> None:
    rows = [*VCF_UTF8_ROWS, "Gm06\t2000\tg2\tC\tT\t.\tPASS\t.\tGT\t0/0\t1/1\t1/1"]
    path = write_vcf(tmp_path / "g.vcf", "\n".join(rows).encode("utf-8") + b"\xe2\x82")
    with pytest.raises(DataContractError, match=r"line 4: not valid UTF-8 \(byte 0xE2 at position 41\)"):
        load_genotypes(path)


# (d) the other readers --------------------------------------------------------------------------------


def test_hapmap_latin1_byte(tmp_path: Path) -> None:
    text = (
        f"{HAPMAP_HEADER}\tRP\tDONOR\tL1\n".encode()
        + f"h1\tA/G\tGm02\t100\t{HAPMAP_FIXED}\tAA\tGG\tAG\n".encode()
        + f"h\xe92\tC/T\tGm02\t200\t{HAPMAP_FIXED}\tCC\tTT\tCT\n".encode("latin-1")
    )
    path = tmp_path / "g.hmp.txt"
    path.write_bytes(text)
    with pytest.raises(DataContractError, match=r"line 3: not valid UTF-8 \(byte 0xE9 at position 2\)"):
        load_genotypes(path)


def test_wide_csv_latin1_byte(tmp_path: Path) -> None:
    path = tmp_path / "g.csv"
    path.write_bytes(b"marker_id,chrom,pos_bp,RP,DONOR,L1\ne1,Gm04,1000,A,G,A/G\ne\xe92,Gm04,2000,C,T,C\n")
    with pytest.raises(DataContractError, match=r"line 3: not valid UTF-8 \(byte 0xE9 at position 2\)"):
        load_genotypes(path)


def test_samples_latin1_byte(tmp_path: Path) -> None:
    path = tmp_path / "samples.csv"
    path.write_bytes(
        b"sample_id,line_name,role,generation,family_id,notes\n"
        b"RP,Recurrent,recurrent_parent,,,\nDONOR,Donor,donor_parent,,,\nL1,Line 1,candidate,,,s\xe9lection\n"
    )
    with pytest.raises(DataContractError, match=r"line 4: not valid UTF-8 \(byte 0xE9 at position 24\)"):
        read_samples(path)


def test_markers_latin1_byte(tmp_path: Path) -> None:
    path = tmp_path / "markers.csv"
    path.write_bytes(b"marker_id,chrom,pos_bp,cm\ne1,Gm04,1000,0.5\ne\xe92,Gm04,2000,1.5\n")
    with pytest.raises(DataContractError, match=r"line 3: not valid UTF-8 \(byte 0xE9 at position 2\)"):
        read_markers(path)


def test_criteria_latin1_byte(tmp_path: Path) -> None:
    path = tmp_path / "criteria.yaml"
    path.write_bytes(b"name: s\xe9lection\n")
    with pytest.raises(CriteriaError, match=r"not valid UTF-8"):
        read_criteria(path)


def test_bom_bytes_count_in_the_position(tmp_path: Path) -> None:
    path = tmp_path / "g.csv"
    path.write_bytes(b"\xef\xbb\xbfmarker_id,chrom,pos_bp,R\xe9\n")
    with pytest.raises(DataContractError, match=r"line 1: not valid UTF-8 \(byte 0xE9 at position 28\)"):
        load_genotypes(path)


def test_results_csv_latin1_byte_through_select(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """results.csv is outside the contract, but `select` reports a bad byte as the same error, exit 1."""
    path = tmp_path / "results.csv"
    path.write_bytes(b"sample_id,notes\nL1,s\xe9lection\n")
    assert main(["select", "--results", str(path), "--top", "1", "--out", str(tmp_path / "selected.csv")]) == 1
    err = capsys.readouterr().err
    assert f"{path}: line 2: not valid UTF-8 (byte 0xE9 at position 5)" in err


def test_custom_token_profile_latin1_byte_through_the_cli(tmp_path: Path) -> None:
    """Outside the contract (D8), but named with the same line and position as backcross's `token profile`."""
    path = tmp_path / "profile.json"
    path.write_bytes(b'{\n  "id": "caf\xe9"\n}\n')
    with pytest.raises(DataContractError) as info:
        _profile_ref(str(path))
    assert str(info.value) == f"token profile file {path}: line 2: not valid UTF-8 (byte 0xE9 at position 13)"


# (e) positives ----------------------------------------------------------------------------------------


def test_bom_prefixed_wide_csv_loads(tmp_path: Path) -> None:
    path = tmp_path / "g.csv"
    path.write_bytes(b"\xef\xbb\xbf" + WIDE)
    gm = load_genotypes(path)
    assert gm.sample_ids == ["RP", "DONOR", "L1"]
    assert [m.marker_id for m in gm.markers] == ["e1", "e2"]


def test_non_ascii_ids_round_trip(tmp_path: Path) -> None:
    geno = tmp_path / "genotypes.csv"
    geno.write_bytes("marker_id,chrom,pos_bp,RP,DONOR,Lé1\nré1,Gm04,1000,A,G,A/G\nr2,Gm04,2000,C,T,C\n".encode())
    samples = tmp_path / "samples.csv"
    samples.write_bytes("sample_id,line_name,role\nRP,,recurrent_parent\nDONOR,,donor_parent\nLé1,Línea 1,candidate\n".encode())
    ds = load_dataset(geno, samples)
    assert ds.genotypes.sample_ids == ["RP", "DONOR", "Lé1"]
    assert [s.sample_id for s in ds.samples] == ["RP", "DONOR", "Lé1"]
    assert [m.marker_id for m in ds.genotypes.markers] == ["ré1", "r2"]


def test_bgzip_character_split_across_members(tmp_path: Path) -> None:
    head = "\n".join(VCF_UTF8_ROWS).encode("utf-8") + b"\nGm06\t2000\tz2\xe2"
    tail = b"\x82\xac\tC\tT\t.\tPASS\t.\tGT\t0/0\t1/1\t1/1\n"
    path = tmp_path / "g.vcf.bgz"
    path.write_bytes(gzip.compress(head) + gzip.compress(tail))
    gm = load_genotypes(path)
    assert [m.marker_id for m in gm.markers] == ["g1", "z2€"]


# (f) the fallback message -----------------------------------------------------------------------------


def test_message_on_a_well_formed_file_names_the_file_only(tmp_path: Path) -> None:
    path = tmp_path / "g.csv"
    path.write_bytes(WIDE)
    assert invalid_utf8_message(path) == f"{path}: not valid UTF-8"


def test_message_when_the_reread_fails(tmp_path: Path) -> None:
    path = tmp_path / "g.vcf.gz"
    path.write_bytes(b"not gzip at all")
    assert invalid_utf8_message(path) == f"{path}: not valid UTF-8"
