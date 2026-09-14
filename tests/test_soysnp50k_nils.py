"""Tests for scripts/soysnp50k_nils.py."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from progeny_selector.io import load_dataset

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "soysnp50k_nils.py"


def _load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location("soysnp50k_nils", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


soysnp50k_nils = _load_module()


def _write_vcf(path: Path) -> None:
    samples = ["PI548533", "PI86024", "PI547416", "PI999999", "PIOTHER"]
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples])
    lines = [
        "##fileformat=VCFv4.2",
        header,
        "\t".join(
            [
                "glyma.Wm82.gnm2.Gm01",
                "100",
                "snp1",
                "A",
                "T",
                ".",
                ".",
                ".",
                "GT",
                "0/0:12",
                "1/1:12",
                "0/1:12",
                "0/0:12",
                "0/0:12",
            ]
        ),
        "\t".join(
            [
                "glyma.Wm82.gnm2.Gm20",
                "200",
                "snp2",
                "C",
                "G",
                ".",
                ".",
                ".",
                "GT",
                "0/0:5",
                "1/1:5",
                "1/1:5",
                "0/0:5",
                "0/0:5",
            ]
        ),
        "\t".join(
            [
                "glyma.Wm82.gnm2.scaffold_105",
                "300",
                "snp3",
                "A",
                "C",
                ".",
                ".",
                ".",
                "GT",
                "0/0:7",
                "1/1:7",
                "0/1:7",
                "0/0:7",
                "0/0:7",
            ]
        ),
        "\t".join(
            [
                "glyma.Wm82.gnm2.Gm01",
                "150",
                "snp4",
                "A",
                ".",
                ".",
                ".",
                ".",
                "GT",
                "0/0:9",
                "1/1:9",
                "0/0:9",
                "0/0:9",
                "0/0:9",
            ]
        ),
    ]
    with gzip.open(path, "wt", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")


def _write_pedigree(path: Path) -> None:
    rows = [
        ("Line", "Female Parent", "Male Parent"),
        ("PI547416", "PI548533", "Clark (5) x PI86024"),
        ("PI547417", "PI548533", "Clark (5) x PI86024"),
        ("PI547418", "PI548533", "Clark (5) x Higan"),
        ("Clark (5) x PI86024", "PI548533", "PI86024"),
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerows(rows)


def test_full_run(tmp_path: Path) -> None:
    vcf_path = tmp_path / "in.vcf.gz"
    pedigree_path = tmp_path / "GS_Pedigrees.csv"
    out_dir = tmp_path / "out"
    _write_vcf(vcf_path)
    _write_pedigree(pedigree_path)

    rc = soysnp50k_nils.main(
        [
            "--vcf",
            str(vcf_path),
            "--pedigrees",
            str(pedigree_path),
            "--recurrent",
            "Clark",
            "--donor",
            "PI86024",
            "--generation",
            "BC{n}F6",
            "--out",
            str(out_dir),
        ]
    )
    assert rc == 0

    genotypes_lines = (out_dir / "genotypes.vcf").read_text().splitlines()
    header_line = next(line_ for line_ in genotypes_lines if line_.startswith("#CHROM"))
    header_cols = header_line.split("\t")
    assert header_cols[9:] == ["PI548533", "PI86024", "PI547416"]

    record_lines = [line_ for line_ in genotypes_lines if not line_.startswith("#")]
    assert len(record_lines) == 2
    chroms = [line_.split("\t")[0] for line_ in record_lines]
    assert chroms == ["Gm01", "Gm20"]
    positions = [line_.split("\t")[1] for line_ in record_lines]
    assert "150" not in positions, "malformed record (ALT '.' with non-reference GT) must be dropped"
    for line_ in record_lines:
        gts = line_.split("\t")[9:]
        for gt in gts:
            assert ":" not in gt

    with open(out_dir / "samples.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    by_id = {r["sample_id"]: r for r in rows}
    assert by_id["PI548533"]["role"] == "recurrent_parent"
    assert by_id["PI86024"]["role"] == "donor_parent"
    assert by_id["PI547416"]["role"] == "candidate"
    assert by_id["PI547416"]["generation"] == "BC5F6"
    assert by_id["PI547416"]["family_id"] == "Clark_x_PI86024"
    assert by_id["PI547416"]["notes"] == "Clark (5) x PI86024"

    dataset = load_dataset(out_dir / "genotypes.vcf", out_dir / "samples.csv")
    chrom_names = {m.chrom for m in dataset.genotypes.markers}
    assert chrom_names == {"Gm01", "Gm20"}


def test_missing_donor_returns_2(tmp_path: Path) -> None:
    vcf_path = tmp_path / "in.vcf.gz"
    pedigree_path = tmp_path / "GS_Pedigrees.csv"
    out_dir = tmp_path / "out"
    _write_vcf(vcf_path)
    _write_pedigree(pedigree_path)

    rc = soysnp50k_nils.main(
        [
            "--vcf",
            str(vcf_path),
            "--pedigrees",
            str(pedigree_path),
            "--recurrent",
            "Clark",
            "--donor",
            "NotADonor",
            "--donor-pi",
            "PI000000",
            "--out",
            str(out_dir),
        ]
    )
    assert rc == 2


def test_named_donor_without_pi_fails(tmp_path: Path) -> None:
    vcf_path = tmp_path / "in.vcf.gz"
    pedigree_path = tmp_path / "GS_Pedigrees.csv"
    out_dir = tmp_path / "out"
    _write_vcf(vcf_path)
    _write_pedigree(pedigree_path)

    try:
        rc = soysnp50k_nils.main(
            [
                "--vcf",
                str(vcf_path),
                "--pedigrees",
                str(pedigree_path),
                "--recurrent",
                "Clark",
                "--donor",
                "Higan",
                "--out",
                str(out_dir),
            ]
        )
    except SystemExit as exc:
        rc = exc.code
    assert rc not in (0, None)
