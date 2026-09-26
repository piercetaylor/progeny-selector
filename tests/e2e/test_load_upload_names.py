"""Load screen: uploads whose Shiny-assigned temp name keeps only the last suffix.

Shiny stores every upload at ``<tempdir>/<index><suffix>`` where ``suffix`` is
``pathlib.Path(name).suffix`` (the last one only). A genotype file with a
compound extension -- ``genotypes.vcf.gz``, ``genotypes.hmp.txt`` -- therefore
arrives on disk as ``0.gz`` or ``0.txt`` unless the Load screen restores the
original name before handing the path to ``io.load_genotypes``.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest
from playwright.sync_api import Page, expect
from shiny.playwright import controller
from shiny.pytest import create_app_fixture
from shiny.run import ShinyAppProc

from tests.e2e.helpers import FIXTURE, LOAD_STATUS

app = create_app_fixture(["../../src/progeny_selector/app/app.py"])

HAPMAP_HEADER = "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode"
HAPMAP_FIXED = "+\tNA\tNA\tNA\tNA\tNA\tNA"


def _vcf_to_hapmap(vcf_text: str) -> str:
    """The fixture VCF (biallelic REF/ALT, GT-only FORMAT) rewritten as HapMap."""
    header_fields: list[str] | None = None
    sample_ids: list[str] = []
    lines = [HAPMAP_HEADER]
    for line in vcf_text.splitlines():
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header_fields = line.lstrip("#").split("\t")
            sample_ids = header_fields[9:]
            lines[0] = HAPMAP_HEADER + "\t" + "\t".join(sample_ids)
            continue
        fields = line.split("\t")
        chrom, pos, marker_id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        gts = fields[9:]
        calls = []
        for gt in gts:
            if gt == "./.":
                calls.append("NN")
                continue
            alleles = [ref if a == "0" else alt for a in gt.split("/")]
            calls.append("".join(sorted(alleles)))
        row = [marker_id, f"{ref}/{alt}", chrom, pos, HAPMAP_FIXED, *calls]
        lines.append("\t".join(row))
    return "\n".join(lines) + "\n"


def _make_gz(tmp_path: Path) -> Path:
    vcf_text = (FIXTURE / "genotypes.vcf").read_text()
    out = tmp_path / "genotypes.vcf.gz"
    with gzip.open(out, "wt") as fh:
        fh.write(vcf_text)
    return out


def _make_hapmap(tmp_path: Path) -> Path:
    vcf_text = (FIXTURE / "genotypes.vcf").read_text()
    out = tmp_path / "genotypes.hmp.txt"
    out.write_text(_vcf_to_hapmap(vcf_text))
    return out


@pytest.mark.parametrize("make_genotypes", [_make_gz, _make_hapmap], ids=["vcf.gz", "hmp.txt"])
def test_compound_suffix_upload_names(page: Page, app: ShinyAppProc, tmp_path: Path, make_genotypes) -> None:
    genotypes_path = make_genotypes(tmp_path)
    page.goto(app.url)
    controller.InputFile(page, "load-genotypes").set(genotypes_path)
    controller.InputFile(page, "load-samples").set(FIXTURE / "samples.csv")
    controller.InputFile(page, "load-markers").set(FIXTURE / "markers.csv")
    controller.InputFile(page, "load-criteria").set(FIXTURE / "criteria.yaml")
    controller.InputActionButton(page, "load-run").click()
    expect(page.locator("#load-status")).to_contain_text(LOAD_STATUS, timeout=60_000)
