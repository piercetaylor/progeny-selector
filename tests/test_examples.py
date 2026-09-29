"""The shipped examples: byte copies of the fixture, and the ``example`` command (docs/m4-phases.md, phase 2)."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from progeny_selector.cli import main
from progeny_selector.examples import example_files

REPO = Path(__file__).resolve().parents[1]
FIXTURES = REPO / "tests" / "fixtures"

# (example name, file, fixture directory the bytes must match)
COPY_CASES = [
    ("synthetic_bc2f1", "genotypes.vcf", "synthetic_bc2f1"),
    ("synthetic_bc2f1", "samples.csv", "synthetic_bc2f1"),
    ("synthetic_bc2f1", "markers.csv", "synthetic_bc2f1"),
    ("synthetic_bc2f1", "criteria.yaml", "synthetic_bc2f1"),
    ("synthetic_bc3f1", "genotypes.vcf", "synthetic_bc3f1"),
    ("synthetic_bc3f1", "samples.csv", "synthetic_bc3f1"),
    # shared: bc3f1 reuses bc2f1's markers.csv and criteria.yaml
    ("synthetic_bc3f1", "markers.csv", "synthetic_bc2f1"),
    ("synthetic_bc3f1", "criteria.yaml", "synthetic_bc2f1"),
]


@pytest.mark.parametrize(("name", "file", "src_name"), COPY_CASES)
def test_package_copy_is_byte_identical_to_the_fixture(name: str, file: str, src_name: str) -> None:
    assert example_files(name)[file].read_bytes() == (FIXTURES / src_name / file).read_bytes()


def test_example_command_writes_self_contained_directories(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    ex = tmp_path / "ex"
    assert main(["example", "--out", str(ex)]) == 0
    assert sorted(os.listdir(ex / "synthetic_bc3f1")) == ["criteria.yaml", "genotypes.vcf", "markers.csv", "samples.csv"]
    capsys.readouterr()
    d = ex / "synthetic_bc3f1"
    rc = main(
        [
            "rank",
            "--genotypes",
            str(d / "genotypes.vcf"),
            "--samples",
            str(d / "samples.csv"),
            "--markers",
            str(d / "markers.csv"),
            "--criteria",
            str(d / "criteria.yaml"),
            "--out",
            str(tmp_path / "results.csv"),
        ]
    )
    assert rc == 0
    assert "40 individuals" in capsys.readouterr().out


def test_example_refuses_to_overwrite_unless_forced(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    out = str(tmp_path / "ex")
    assert main(["example", "--out", out]) == 0
    capsys.readouterr()
    assert main(["example", "--out", out]) == 1
    assert "exists; pass --force" in capsys.readouterr().err
    assert main(["example", "--out", out, "--force"]) == 0


def test_refusal_under_all_writes_nothing(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    out = tmp_path / "ex"
    (out / "synthetic_bc3f1").mkdir(parents=True)
    (out / "synthetic_bc3f1" / "samples.csv").write_text("existing\n", encoding="utf-8")
    assert main(["example", "--out", str(out), "--name", "all"]) != 0
    assert "exists; pass --force" in capsys.readouterr().err
    assert not (out / "synthetic_bc2f1").exists()
    assert (out / "synthetic_bc3f1" / "samples.csv").read_text(encoding="utf-8") == "existing\n"


def test_unknown_example_name_is_a_usage_error() -> None:
    with pytest.raises(SystemExit) as exc:
        main(["example", "--name", "nope"])
    assert exc.value.code == 2


def test_example_output_is_ignored_by_gitignore_but_package_data_is_not() -> None:
    lines = [line.strip() for line in (REPO / ".gitignore").read_text(encoding="utf-8").splitlines()]
    assert "*.vcf" in lines
    assert "!src/progeny_selector/examples/**" in lines
    assert lines.index("!src/progeny_selector/examples/**") > lines.index("*.vcf")
