"""Boundary parsers: VCF (plain and gzip), HapMap, wide CSV (nucleotide and coded), manifests, criteria."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pytest

from progeny_selector.io import load_dataset, load_genotypes, read_criteria
from progeny_selector.io.criteria import criteria_from_dict
from progeny_selector.io.manifest import read_markers, read_samples
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
    assert gm.markers[1].marker_id == "Gm06_2000"
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
