"""Boundary parsers: VCF (plain and gzip), HapMap, wide CSV (nucleotide and coded), manifests, criteria."""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pytest

from progeny_selector.constants import HAPMAP_MISSING, WIDE_NUCLEOTIDE_MISSING
from progeny_selector.io import load_dataset, load_genotypes, read_criteria
from progeny_selector.io.calls import detect_coding, parse_nucleotide_call
from progeny_selector.io.criteria import criteria_from_dict
from progeny_selector.io.manifest import read_markers, read_samples
from progeny_selector.io.profiles import BUILTIN_PROFILES, compile_profile
from progeny_selector.io.wide_csv import read_wide_csv
from progeny_selector.model.criteria import CriteriaError
from progeny_selector.model.dataset import DataContractError

VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tRP\tDONOR\tP1
chr6\t1000\tm1\tA\tT\t.\tPASS\t.\tGT:DP\t0/0:10\t1/1:9\t0|1:8
6\t2000\t.\tC\tG,T\t.\t.\t.\tGT\t0/0\t1/1\t2/2
Gm06\t3000\tm3\tA\tT\t.\t.\t.\tGT\t0/0\t1/1\t./.
"""

HAPMAP = (
    "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\tRP\tDONOR\tP1\n"
    + "m1\tA/T\t6\t1000\t+\tNA\tNA\tNA\tNA\tNA\tNA\tAA\tTT\tW\n"
    + "m2\tC/G\t6\t2000\t+\tNA\tNA\tNA\tNA\tNA\tNA\tC\tG\tNN\n"
)

WIDE_NUC = "marker_id,chrom,pos_bp,RP,DONOR,P1\nm1,Gm06,1000,A,T,A/T\nm2,chr6,2000,C,G,GG\nm3,6,3000,A,T,NA\n"
WIDE_CODED = "marker_id,chrom,pos_bp,P1,P2\nm1,Gm06,1000,H,A\nm2,Gm06,2000,B,N\nm3,Gm06,3000,A,H\n"
SAMPLES = (
    "sample_id,line_name,role,generation,family_id,notes\n"
    "RP,RP line,recurrent_parent,,,\nDONOR,Donor,donor_parent,,,\nP1,P1,progeny,BC2F1,F1,\n"
)


def test_vcf_plain_and_gzip(tmp_path: Path):
    p = tmp_path / "g.vcf"
    p.write_text(VCF)
    gm = load_genotypes(p)
    assert [m.chrom for m in gm.markers] == ["Gm06"] * 3
    assert gm.markers[1].marker_id == "6_2000"
    assert gm.alleles[1] == ["C", "G", "T"]
    assert tuple(gm.calls[0, 2]) == (0, 1) and tuple(gm.calls[1, 2]) == (2, 2) and tuple(gm.calls[2, 2]) == (-1, -1)
    gz = tmp_path / "g.vcf.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(VCF)
    assert np.array_equal(load_genotypes(gz).calls, gm.calls)


def test_vcf_bgz_extension(tmp_path: Path):
    p = tmp_path / "g.vcf"
    p.write_text(VCF)
    gm = load_genotypes(p)
    bgz = tmp_path / "g.vcf.bgz"
    with gzip.open(bgz, "wt") as fh:
        fh.write(VCF)
    assert np.array_equal(load_genotypes(bgz).calls, gm.calls)


def test_bcf_extension_rejected(tmp_path: Path):
    p = tmp_path / "g.bcf"
    with pytest.raises(DataContractError, match="cannot infer genotype format"):
        load_genotypes(p)


def test_hapmap(tmp_path: Path):
    p = tmp_path / "g.hmp.txt"
    p.write_text(HAPMAP)
    gm = load_genotypes(p)
    assert gm.alleles[0] == ["A", "T"] and tuple(gm.calls[0, 2]) == (0, 1)
    assert tuple(gm.calls[1, 2]) == (-1, -1)


def test_hapmap_blank_line_before_header(tmp_path: Path):
    p = tmp_path / "g.hmp.txt"
    p.write_text("  \n" + HAPMAP)
    gm = load_genotypes(p)
    assert gm.alleles[0] == ["A", "T"] and tuple(gm.calls[0, 2]) == (0, 1)


def test_hapmap_hash_line_before_header_is_an_error(tmp_path: Path):
    p = tmp_path / "g.hmp.txt"
    p.write_text("# note\n" + HAPMAP)
    with pytest.raises(DataContractError, match="must start with rs#"):
        load_genotypes(p)


def test_wide_csv_tab_blank_line_before_header(tmp_path: Path):
    p = tmp_path / "g.csv"
    p.write_text("   \nmarker_id\tchrom\tpos_bp\tRP\tDONOR\nm1\tGm01\t100\tA\tT\n")
    gm = read_wide_csv(p)
    assert gm.sample_ids == ["RP", "DONOR"]
    assert [m.marker_id for m in gm.markers] == ["m1"]


def test_samples_blank_line_before_header(tmp_path: Path):
    s = tmp_path / "samples.csv"
    s.write_text("\n  \nsample_id,role\nRP,recurrent_parent\nDONOR,donor_parent\nL1,candidate\n")
    assert [x.sample_id for x in read_samples(s)] == ["RP", "DONOR", "L1"]


def test_hapmap_nbsp_line_is_not_blank(tmp_path: Path):
    p = tmp_path / "g.hmp.txt"
    p.write_bytes((chr(0xA0) + "\n" + HAPMAP).encode("utf-8"))
    with pytest.raises(DataContractError, match="must start with rs#"):
        load_genotypes(p)


def test_samples_quoted_line_break_crlf_read_as_lf(tmp_path: Path):
    s = tmp_path / "samples.csv"
    s.write_bytes(b'sample_id,role,notes\r\nRP,recurrent_parent,"first\r\nsecond"\r\nDONOR,donor_parent,\r\nL1,candidate,\r\n')
    samples = read_samples(s)
    assert [x.sample_id for x in samples] == ["RP", "DONOR", "L1"]
    assert samples[0].notes == "first\nsecond"


def test_samples_quoted_line_break_lf(tmp_path: Path):
    s = tmp_path / "samples.csv"
    s.write_bytes(b'sample_id,role,notes\nRP,recurrent_parent,"first\nsecond"\nDONOR,donor_parent,\nL1,candidate,\n')
    assert read_samples(s)[0].notes == "first\nsecond"


def test_midfield_quote_is_literal(tmp_path: Path):
    s = tmp_path / "samples.csv"
    s.write_bytes(b'sample_id,line_name,role\nRP,,recurrent_parent\nL1,6" pot,candidate\nDONOR,,donor_parent\n')
    samples = read_samples(s)
    assert [x.sample_id for x in samples] == ["RP", "L1", "DONOR"]
    assert samples[1].line_name == '6" pot'


def test_unterminated_quote_names_opening_line(tmp_path: Path):
    m = tmp_path / "markers.csv"
    m.write_bytes(b'marker_id,chrom,pos_bp,cm\r\nr1,Gm02,1000,1.5\r\nr2,Gm02,"2000,2.5\r\nr3,Gm02,3000,3.5\r\n')
    with pytest.raises(DataContractError, match="line 3: unterminated quoted field"):
        read_markers(m)


def test_markers_invalid_cm_names_line(tmp_path: Path):
    m = tmp_path / "markers.csv"
    m.write_bytes(b"marker_id,chrom,pos_bp,cm\nr1,Gm02,1000,1.5\nr2,Gm02,2000,abc\n")
    with pytest.raises(DataContractError, match=r"line 3: invalid cm 'abc'"):
        read_markers(m)


def test_wide_csv_nucleotide_and_coded(tmp_path: Path):
    p = tmp_path / "g.csv"
    p.write_text(WIDE_NUC)
    gm = read_wide_csv(p)
    assert not gm.coded and tuple(gm.calls[0, 2]) == (0, 1) and tuple(gm.calls[1, 2]) == (1, 1)
    q = tmp_path / "c.csv"
    q.write_text(WIDE_CODED)
    gm2 = read_wide_csv(q)
    assert gm2.coded and tuple(gm2.calls[0, 0]) == (0, 1) and tuple(gm2.calls[1, 0]) == (1, 1) and tuple(gm2.calls[1, 1]) == (-1, -1)
    s = tmp_path / "samples.csv"
    s.write_text("sample_id,line_name,role\nRP,RP,recurrent_parent\nDONOR,D,donor_parent\nP1,P1,progeny\nP2,P2,progeny\n")
    ds = load_dataset(q, s)
    assert ds.genotypes.n_samples == 4 and any("synthesised" in w for w in ds.warnings)
    assert ds.synthetic_sample_ids == ("RP", "DONOR") and ds.genotypes.sample_ids == ["RP", "DONOR", "P1", "P2"]


def test_manifest_rules(tmp_path: Path):
    s = tmp_path / "s.csv"
    s.write_text(SAMPLES)
    assert [x.role for x in read_samples(s)] == ["recurrent_parent", "donor_parent", "progeny"]
    s.write_text(SAMPLES.replace("donor_parent", "recurrent_parent"))
    with pytest.raises(DataContractError, match="exactly one recurrent_parent"):
        read_samples(s)
    s.write_text(SAMPLES.replace("progeny", "child"))
    with pytest.raises(DataContractError, match="role"):
        read_samples(s)
    m = tmp_path / "m.csv"
    m.write_text("marker_id,chrom,pos_bp,cm\nm1,6,1000,0.5\nm2,6,2000,\n")
    mm = read_markers(m)
    assert mm["m1"].cm == 0.5 and mm["m2"].cm is None and mm["m2"].chrom == "Gm06"
    v = tmp_path / "g.vcf"
    v.write_text(VCF)
    s.write_text(SAMPLES + "P9,P9,progeny,BC2F1,F1,\n")
    with pytest.raises(DataContractError, match="not in genotype file"):
        load_dataset(v, s)


def test_criteria_parsing_and_validation(tmp_path: Path):
    doc = {
        "targets": [{"id": "T1", "region": "Gm06:1,000-3,000", "required_state": "hom_donor"}],
        "avoid": [{"locus_id": "A1", "marker_id": "m2"}],
        "weights": {"drag": 0.3},
    }
    crit = criteria_from_dict(doc)
    assert crit.targets[0].chrom == "Gm06" and crit.targets[0].start_bp == 1000 and crit.targets[0].kind() == "region"
    assert crit.weights.drag == 0.3 and crit.weights.rpp_noncarrier == 0.5
    with pytest.raises(CriteriaError, match="unknown keys"):
        criteria_from_dict({"targets": [{"locus_id": "T1", "marker_id": "m1", "reqd": "x"}]})
    with pytest.raises(CriteriaError, match="required_state"):
        criteria_from_dict({"targets": [{"locus_id": "T1", "marker_id": "m1", "required_state": "homo"}]})
    with pytest.raises(CriteriaError, match="at least one target"):
        criteria_from_dict({"avoid": []})
    y = tmp_path / "c.yaml"
    y.write_text("targets:\n  - locus_id: T1\n    left_marker: m1\n    right_marker: m3\nflank_unit: bp\nflank_window: 1500\n")
    c2 = read_criteria(y)
    assert c2.targets[0].kind() == "flanking" and c2.flank_unit == "bp"


HAPMAP_HEADER = "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode"
HAPMAP_FIXED = "+\tNA\tNA\tNA\tNA\tNA\tNA"


def test_manifest_tab_no_line_name(tmp_path: Path):
    s = tmp_path / "s.tsv"
    s.write_bytes("\ufeffsample_id\trole\r\nRP\trecurrent_parent\r\nDONOR\tdonor_parent\r\nP1\tprogeny\r\n".encode())
    rows = read_samples(s)
    assert [x.sample_id for x in rows] == ["RP", "DONOR", "P1"] and rows[2].line_name == "P1"


def test_wide_csv_tab_quoted_bom_crlf(tmp_path: Path):
    p = tmp_path / "g.tsv"
    p.write_bytes(
        '\ufeffmarker_id\tchrom\tpos_bp\tRP\tDONOR\t"P1"\r\nm1\tGm06\t1000\tA\tT\t"A/T"\r\nm2\tGm06\t2000\tC\tG\t.|.\r\n'.encode()
    )
    gm = read_wide_csv(p)
    assert gm.sample_ids == ["RP", "DONOR", "P1"]
    assert tuple(gm.calls[0, 2]) == (0, 1) and tuple(gm.calls[1, 2]) == (-1, -1)


def test_extra_columns_dropped_manifest_order(tmp_path: Path):
    p = tmp_path / "g.csv"
    p.write_text("marker_id,chrom,pos_bp,EXTRA,P1,RP,DONOR\nm1,Gm06,1000,A,A/T,A,T\n")
    s = tmp_path / "s.csv"
    s.write_text("sample_id,role\nDONOR,donor_parent\nP1,progeny\nRP,recurrent_parent\n")
    ds = load_dataset(p, s)
    assert ds.genotypes.sample_ids == ["DONOR", "P1", "RP"]
    assert tuple(ds.genotypes.calls[0, 1]) == (0, 1)
    assert any(w.startswith("1 genotype column(s) not in samples.csv dropped: EXTRA") for w in ds.warnings)


def test_coded_synthetic_parents_tracked(tmp_path: Path):
    p = tmp_path / "c.csv"
    p.write_text("marker_id,chrom,pos_bp,P2,P1\nm1,Gm06,1000,B,H\n")
    s = tmp_path / "s.csv"
    s.write_text("sample_id,role\nRP,recurrent_parent\nP1,progeny\nDONOR,donor_parent\nP2,progeny\n")
    ds = load_dataset(p, s)
    assert ds.synthetic_sample_ids == ("RP", "DONOR")
    assert ds.genotypes.sample_ids == ["RP", "P1", "DONOR", "P2"]


def test_wide_missing_tokens_per_mode(tmp_path: Path):
    p = tmp_path / "g.csv"
    p.write_text("marker_id,chrom,pos_bp,S1,S2,S3\nm1,Gm06,100,.|.,NN,--\nm2,Gm06,200,.,,A\n")
    gm = read_wide_csv(p)
    assert (gm.calls[0] == -1).all() and tuple(gm.calls[1, 2]) == (0, 0)
    for bad in ("--", ".", "./.", ".|.", "NN"):
        p.write_text(f"marker_id,chrom,pos_bp,S1,S2\nm1,Gm06,100,A,B\nm2,Gm06,200,H,{bad}\n")
        with pytest.raises(DataContractError, match=rf"line 3: unrecognised coded call '{re.escape(bad)}'"):
            read_wide_csv(p)


def test_nucleotide_iupac_expands(tmp_path: Path):
    p = tmp_path / "g.csv"
    p.write_text("marker_id,chrom,pos_bp,S1,S2\nm1,Gm06,100,T,R\nm2,Gm06,200,y,S\n")
    gm = read_wide_csv(p)
    assert gm.alleles[0] == ["A", "G", "T"] and tuple(gm.calls[0, 1]) == (0, 1)
    assert gm.alleles[1] == ["C", "G", "T"] and tuple(gm.calls[1, 0]) == (0, 2) and tuple(gm.calls[1, 1]) == (0, 1)


def test_nucleotide_rejects_stray_cells(tmp_path: Path):
    p = tmp_path / "g.csv"
    for bad in ("?", "B", "H", "X", "XX", "0", "+", "A?", "N?", "RR"):  # T in S1 forces nucleotide detection
        p.write_text(f"marker_id,chrom,pos_bp,S1,S2\nm1,Gm06,100,T,{bad}\n")
        with pytest.raises(DataContractError, match=rf"line 2: unrecognised nucleotide call '{re.escape(bad)}'"):
            read_wide_csv(p)


def _half_missing_cells() -> list[str]:
    """Every pair of one nucleotide and one of N, -, . , in both orders and all three spellings."""
    cells: list[str] = []
    for n in "ACGT":
        for h in ("N", "-", "."):
            for a, b in ((n, h), (h, n)):
                cells += [f"{a}{b}", f"{a}/{b}", f"{a}|{b}"]
    return cells


def test_half_missing_pairs_read_as_missing(tmp_path: Path):
    """Contract 1.6.0: a pair of one nucleotide and one of N, - or . is missing, either order, both formats."""
    for cell in _half_missing_cells():
        assert parse_nucleotide_call(cell, WIDE_NUCLEOTIDE_MISSING) is None, cell
        assert parse_nucleotide_call(cell, HAPMAP_MISSING) is None, cell
        lower = cell.lower()  # cells are compared case-insensitively
        assert parse_nucleotide_call(lower, WIDE_NUCLEOTIDE_MISSING) is None, lower
        assert parse_nucleotide_call(f" {lower} ", HAPMAP_MISSING) is None, lower
    p = tmp_path / "g.csv"
    p.write_text("marker_id,chrom,pos_bp,S1,S2,S3,S4\nm1,Gm06,100,T,AN,A-,./A\n")
    gm = read_wide_csv(p)
    assert (gm.calls[0, 1:] == -1).all() and gm.alleles[0] == ["T"]


def test_half_missing_pair_needs_one_of_the_three_tokens():
    """Contract 1.6.0: only N, - and . pair with a nucleotide; another missing token does not."""
    for bad in ("AX", "XA", "A/X", "AU", "UA", "A?"):
        for missing in (WIDE_NUCLEOTIDE_MISSING, HAPMAP_MISSING):
            with pytest.raises(ValueError, match="unrecognised nucleotide call"):
                parse_nucleotide_call(bad, missing)


def test_half_missing_pairs_under_profiles():
    """Contract 1.6.0: the rule is part of the nucleotide grammar, so it applies under base nucleotide only."""
    for pid in ("tassel", "soybase-report"):
        compiled = compile_profile(BUILTIN_PROFILES[pid])
        for cell in _half_missing_cells():
            assert parse_nucleotide_call(cell, WIDE_NUCLEOTIDE_MISSING, compiled) is None, (pid, cell)
            assert parse_nucleotide_call(cell, HAPMAP_MISSING, compiled) is None, (pid, cell)
        extra = "AX" if pid == "tassel" else "AU"
        with pytest.raises(ValueError, match="unrecognised nucleotide call"):
            parse_nucleotide_call(extra, WIDE_NUCLEOTIDE_MISSING, compiled)
    for pid in ("dart", "axiom", "kasp"):
        compiled = compile_profile(BUILTIN_PROFILES[pid])
        for cell in _half_missing_cells():
            if cell in compiled.missing:  # NA, asserted below
                continue
            for text in (cell, cell.lower()):
                with pytest.raises(ValueError, match="unrecognised nucleotide call"):
                    parse_nucleotide_call(text, WIDE_NUCLEOTIDE_MISSING, compiled)
        # "unless the profile lists that exact token": NA is in all three missing lists.
        assert "NA" in compiled.missing, pid
        assert parse_nucleotide_call("NA", WIDE_NUCLEOTIDE_MISSING, compiled) is None, pid
        assert parse_nucleotide_call("na", WIDE_NUCLEOTIDE_MISSING, compiled) is None, pid


def test_half_missing_pair_decides_nucleotide_detection(tmp_path: Path):
    """Contract 1.6.0: the pair is not a missing token, so a file whose only non-A/B cell is the pair
    is detected as nucleotide, and the parse then fails on the coded letter B, not on the pair."""
    assert detect_coding(["A", "B", "N/A"]) == "nucleotide"
    p = tmp_path / "g.csv"
    p.write_text("marker_id,chrom,pos_bp,RP,DONOR,L1\nm1,Gm01,100,A,B,N/A\n")
    with pytest.raises(DataContractError, match="line 2: unrecognised nucleotide call 'B'"):
        read_wide_csv(p)


def test_hapmap_missing_tokens(tmp_path: Path):
    p = tmp_path / "g.hmp.txt"
    p.write_text(
        f"{HAPMAP_HEADER}\tS1\tS2\tS3\tS4\n"
        f"m1\tA/T\t6\t1000\t{HAPMAP_FIXED}\tNA\t./.\t.\tAT\n"
        f"m2\tA/T\t6\t2000\t{HAPMAP_FIXED}\tN\tNN\t-\t--\n"
        f"m3\tA/T\t6\t3000\t{HAPMAP_FIXED}\t\tA\tT\tW\n"
        f"m4\tA/T\t6\t4000\t{HAPMAP_FIXED}\t.|.\tX\tXX\tA\n"
    )
    gm = load_genotypes(p)
    assert (gm.calls[0, :3] == -1).all() and tuple(gm.calls[0, 3]) == (0, 1)
    assert (gm.calls[1] == -1).all()
    assert tuple(gm.calls[2, 0]) == (-1, -1) and tuple(gm.calls[2, 1]) == (0, 0) and tuple(gm.calls[2, 3]) == (0, 1)
    assert (gm.calls[3, :3] == -1).all() and tuple(gm.calls[3, 3]) == (0, 0)


def test_hapmap_rejects_stray_cells(tmp_path: Path):
    p = tmp_path / "g.hmp.txt"
    for bad in ("?", "+", "0", "B", "H", "A?"):
        p.write_text(f"{HAPMAP_HEADER}\tS1\tS2\nm1\tA/T\t6\t1000\t{HAPMAP_FIXED}\tAA\t{bad}\n")
        with pytest.raises(DataContractError, match=rf"line 2: unrecognised nucleotide call '{re.escape(bad)}'"):
            load_genotypes(p)


def test_coding_detection_scans_every_row(tmp_path: Path):
    header = "marker_id,chrom,pos_bp,S1\n"
    body = "".join(f"m{i},Gm01,{i + 1},A\n" for i in range(250))
    p = tmp_path / "late.csv"
    p.write_text(header + body + "late,Gm01,999,B\n")
    assert read_wide_csv(p).coded
    p.write_text(header + "h,Gm01,1,H\n" + body + "late,Gm01,999,T\n")
    with pytest.raises(DataContractError, match="unrecognised nucleotide call 'H'"):
        read_wide_csv(p)  # detected as nucleotide because of the late T; a single H is then an error


def test_vcf_crlf_bom(tmp_path: Path):
    p = tmp_path / "g.vcf"
    p.write_bytes(("\ufeff" + VCF.replace("\n", "\r\n")).encode())
    gm = load_genotypes(p)
    assert gm.sample_ids == ["RP", "DONOR", "P1"] and tuple(gm.calls[2, 2]) == (-1, -1)


def test_separator_shapes_match_contract(tmp_path: Path):
    wide = tmp_path / "g.csv"
    hmp = tmp_path / "g.hmp.txt"
    for bad in ("A/", "/A", "AT/", "A//T", "A|/T", "NA/"):
        wide.write_text(f"marker_id,chrom,pos_bp,S1,S2\nm1,Gm06,100,T,{bad}\n")
        with pytest.raises(DataContractError, match=rf"line 2: unrecognised nucleotide call '{re.escape(bad)}'"):
            read_wide_csv(wide)
        hmp.write_text(f"{HAPMAP_HEADER}\tS1\tS2\nm1\tA/T\t6\t1000\t{HAPMAP_FIXED}\tAA\t{bad}\n")
        with pytest.raises(DataContractError, match=rf"line 2: unrecognised nucleotide call '{re.escape(bad)}'"):
            load_genotypes(hmp)
    wide.write_text("marker_id,chrom,pos_bp,S1,S2,S3\nm1,Gm06,100,A/T,A|T,N/A\n")
    gm = read_wide_csv(wide)
    assert tuple(gm.calls[0, 0]) == (0, 1) and tuple(gm.calls[0, 1]) == (0, 1) and tuple(gm.calls[0, 2]) == (-1, -1)
    hmp.write_text(f"{HAPMAP_HEADER}\tS1\tS2\tS3\nm1\tA/T\t6\t1000\t{HAPMAP_FIXED}\tA/T\tA|T\tN/A\n")
    gm = load_genotypes(hmp)
    assert tuple(gm.calls[0, 0]) == (0, 1) and tuple(gm.calls[0, 1]) == (0, 1) and tuple(gm.calls[0, 2]) == (-1, -1)


def test_hapmap_invalid_position(tmp_path: Path):
    p = tmp_path / "g.hmp.txt"
    p.write_text(f"{HAPMAP_HEADER}\tS1\tS2\nm1\tA/T\t6\t1000.5\t{HAPMAP_FIXED}\tAA\tTT\n")
    with pytest.raises(DataContractError, match=r"line 2: invalid position '1000\.5'"):
        load_genotypes(p)


def test_hapmap_header_requires_rs_hash(tmp_path: Path):
    p = tmp_path / "g.hmp.txt"
    p.write_text(f"{HAPMAP_HEADER.replace('rs#', 'rsid', 1)}\tS1\tS2\nm1\tA/T\t6\t1000\t{HAPMAP_FIXED}\tAA\tTT\n")
    with pytest.raises(DataContractError, match="HapMap header must start with rs#"):
        load_genotypes(p)


def test_parse_position_grammar():
    """Contract 1.2.0: whole-valued text is the integer; everything else is a ValueError naming the cell."""
    from progeny_selector.io.position import parse_position

    ok = [
        ("1000", 1000),
        (" 1000 ", 1000),
        ("+1000", 1000),
        ("1000.", 1000),
        ("1000.0", 1000),
        ("1e3", 1000),
        ("1.0E3", 1000),
        ("1.9E+07", 19000000),
        ("1.5e3", 1500),
        ("0", 0),
        ("-0", 0),
        ("-0.0", 0),
        ("\u00a01000", 1000),
        ("\ufeff1000", 1000),
    ]
    for text, value in ok:
        assert parse_position(text) == value, text
    for bad in ["", " ", "100.7", ".5", "1e-3", "-5", "nan", "inf", "Infinity", "1e400", "0x10", "1_000", "1,000", "abc", "\u0661\u0660"]:
        with pytest.raises(ValueError, match="invalid position"):
            parse_position(bad)
    for text, value in [("1000", 1000), ("01000", 1000)]:
        assert parse_position(text, grammar="digits") == value, text
    for bad in ["+1000", "1000.0", "1e3", "", "-5", "1" * 5000]:
        with pytest.raises(ValueError, match="invalid position"):
            parse_position(bad, grammar="digits")
    for bad in ["\x851000", "\x1f1000"]:
        with pytest.raises(ValueError, match="invalid position"):
            parse_position(bad)


def test_positions_whole_floats_accepted_fractions_rejected(tmp_path: Path):
    wide = tmp_path / "g.csv"
    wide.write_text("marker_id,chrom,pos_bp,S1,S2\nm1,Gm06,1000.0,A,T\nm2,Gm06,2e3,A,T\nm3,Gm06,3.0E3,A,T\nm4,Gm06,+4000,A,T\n")
    assert [m.pos_bp for m in read_wide_csv(wide).markers] == [1000, 2000, 3000, 4000]
    wide.write_text("marker_id,chrom,pos_bp,S1,S2\nm1,Gm06,1000,A,T\nm2,Gm06,100.7,A,T\n")
    with pytest.raises(DataContractError, match=r"line 3: invalid position '100\.7'"):
        read_wide_csv(wide)
    hm = tmp_path / "g.hmp.txt"
    hm.write_text(f"{HAPMAP_HEADER}\tS1\tS2\nm1\tA/T\t6\t1e3\t{HAPMAP_FIXED}\tAA\tTT\nm2\tA/T\t6\t2000.0\t{HAPMAP_FIXED}\tAA\tTT\n")
    assert [m.pos_bp for m in load_genotypes(hm).markers] == [1000, 2000]
    vcf = tmp_path / "g.vcf"
    vcf.write_text(VCF.replace("\t1000\t", "\t1e3\t", 1))
    with pytest.raises(DataContractError, match=r"line 3: invalid position '1e3'"):
        load_genotypes(vcf)
    vcf.write_text(VCF.replace("chr6\t1000\tm1\t", "chr6\t01000\t.\t", 1))
    gm = load_genotypes(vcf)
    assert gm.markers[0].marker_id == "chr6_1000" and gm.markers[0].pos_bp == 1000
    m = tmp_path / "m.csv"
    m.write_text("marker_id,chrom,pos_bp,cm\nm1,6,1000.0,0.5\nm2,6,2e3,1.25\n")
    mm = read_markers(m)
    assert mm["m1"].pos_bp == 1000 and mm["m2"].pos_bp == 2000 and mm["m2"].cm == 1.25
    m.write_text("marker_id,chrom,pos_bp,cm\nm1,6,1000,0.5\nm2,6,2000.25,\n")
    with pytest.raises(DataContractError, match=r"line 3: invalid position '2000\.25'"):
        read_markers(m)


def test_wide_csv_skips_all_empty_rows_and_names_physical_lines(tmp_path: Path):
    p = tmp_path / "g.csv"
    p.write_text("marker_id,chrom,pos_bp,S1,S2\nm1,Gm06,100,A,T\n,,,,\n , , , , \nm2,Gm06,200,A,T\n")
    assert [m.marker_id for m in read_wide_csv(p).markers] == ["m1", "m2"]
    p.write_text("marker_id,chrom,pos_bp,S1,S2\nm1,Gm06,100,A,T\n,Gm06,200,A,T\n")
    with pytest.raises(DataContractError, match="line 3: empty marker_id"):
        read_wide_csv(p)
    p.write_text("marker_id,chrom,pos_bp,S1,S2\nm1,Gm06,100,A,T\n\nm2,Gm06,,A,T\n")
    with pytest.raises(DataContractError, match=r"line 4: invalid position ''"):
        read_wide_csv(p)


def test_markers_csv_empty_marker_id_names_line(tmp_path: Path):
    m = tmp_path / "m.csv"
    m.write_text("marker_id,chrom,pos_bp,cm\n,6,1000,0.5\n")
    with pytest.raises(DataContractError, match="line 2: empty marker_id"):
        read_markers(m)
