"""The two-pass `read_vcf` (docs/adr/0022): the memory bound, and every rule unchanged.

The discriminating test is `test_parse_peak_is_bounded_by_the_matrix`: the pre-phase
reader held a list of per-record `(s, 2)` arrays and then `np.stack`ed them, so its
`tracemalloc` peak was about twice the matrix; the two-pass reader allocates the matrix
once. The rest pin what must not change: line numbers, messages, blank-line handling,
FORMAT order, and that a truncated `.gz` is a contract error naming the file.
"""

from __future__ import annotations

import gzip
import tracemalloc
from pathlib import Path

import numpy as np
import pytest

from progeny_selector.core.chrom import SOYBEAN
from progeny_selector.io import vcf as vcf_module
from progeny_selector.io.vcf import read_vcf
from progeny_selector.model.dataset import DataContractError

HEADER_FIXED = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
GT_POOL = ("0/0", "0/1", "1/1", "./.")


def _write_vcf(
    path: Path,
    n_markers: int,
    n_samples: int,
    rng: np.random.Generator,
    fmt: str = "GT",
    n_depths: int = 1,
) -> None:
    """A VCF of `n_markers` biallelic records over `n_samples` samples, LF, one FORMAT."""
    sample_ids = [f"S{j:05d}" for j in range(n_samples)]
    keys = fmt.split(":")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("\t".join([*HEADER_FIXED, *sample_ids]) + "\n")
        for i in range(n_markers):
            gts = [GT_POOL[int(k)] for k in rng.integers(0, len(GT_POOL), n_samples)]
            if keys == ["GT"]:
                cells = gts
            else:
                cells = []
                for j, gt in enumerate(gts):
                    dp = str((i * n_samples + j) % n_depths)
                    cells.append(":".join(gt if k == "GT" else dp for k in keys))
            fh.write("\t".join([f"Gm{i % 20 + 1:02d}", str(1000 * (i + 1)), f"m{i:06d}", "A", "T", ".", ".", ".", fmt, *cells]) + "\n")


def test_parse_peak_is_bounded_by_the_matrix(tmp_path: Path) -> None:
    p = tmp_path / "g.vcf"
    _write_vcf(p, 2_000, 1_000, np.random.default_rng(0))
    tracemalloc.start()
    gm = read_vcf(p)
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert gm.calls.nbytes == 4_000_000
    assert peak < 1.5 * gm.calls.nbytes, f"peak={peak} matrix={gm.calls.nbytes} ratio={peak / gm.calls.nbytes:.2f}"


def test_gz_matches_plain(tmp_path: Path) -> None:
    p = tmp_path / "g.vcf"
    _write_vcf(p, 40, 6, np.random.default_rng(1))
    gz = tmp_path / "g.vcf.gz"
    with gzip.open(gz, "wb") as out:
        out.write(p.read_bytes())
    plain, zipped = read_vcf(p), read_vcf(gz)
    assert np.array_equal(plain.calls, zipped.calls)
    assert plain.markers == zipped.markers
    assert plain.sample_ids == zipped.sample_ids
    assert plain.alleles == zipped.alleles


def test_blank_lines_are_skipped_in_both_passes(tmp_path: Path) -> None:
    p = tmp_path / "g.vcf"
    _write_vcf(p, 5, 3, np.random.default_rng(2))
    lines = p.read_text(encoding="utf-8").split("\n")
    padded = tmp_path / "padded.vcf"
    with open(padded, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(lines[0] + "\n" + lines[1] + "\n")
        for line in lines[2:7]:
            fh.write("\n   \t \n" + line + "\n")
        fh.write("\n \n")
    plain, with_blanks = read_vcf(p), read_vcf(padded)
    assert with_blanks.n_markers == plain.n_markers == 5
    assert np.array_equal(with_blanks.calls, plain.calls)


def test_format_order_and_multidigit_indices(tmp_path: Path) -> None:
    p = tmp_path / "wide.vcf"
    alts = ",".join(["T", "G", "C", "AA", "AT", "AG", "AC", "TA", "TT", "TG"])  # 11 alleles with REF
    with open(p, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("\t".join([*HEADER_FIXED, "S1"]) + "\n")
        fh.write("\t".join(["Gm01", "1000", "m1", "A", alts, ".", ".", ".", "DP:GT", "7:10/3"]) + "\n")
    gm = read_vcf(p)
    assert tuple(int(x) for x in gm.calls[0, 0]) == (10, 3)

    plain = tmp_path / "plain.vcf"
    withdp = tmp_path / "withdp.vcf"
    _write_vcf(plain, 30, 10, np.random.default_rng(3))
    _write_vcf(withdp, 30, 10, np.random.default_rng(3), fmt="GT:DP", n_depths=300)
    assert np.array_equal(read_vcf(plain).calls, read_vcf(withdp).calls)


def test_error_line_numbers_are_unchanged(tmp_path: Path) -> None:
    rng = np.random.default_rng(4)
    good = tmp_path / "good.vcf"
    _write_vcf(good, 4, 3, rng)
    lines = good.read_text(encoding="utf-8").rstrip("\n").split("\n")

    short = tmp_path / "short.vcf"
    body = list(lines)
    body[4] = "\t".join(body[4].split("\t")[:-1])  # physical line 5 loses one sample column
    short.write_text("\n".join(body) + "\n", encoding="utf-8", newline="")
    with pytest.raises(DataContractError, match=r"line 5: expected \d+ columns, found \d+"):
        read_vcf(short)

    hashed = tmp_path / "hashed.vcf"
    hashed.write_text("\n".join([*lines[:4], "# a note", *lines[4:]]) + "\n", encoding="utf-8", newline="")
    with pytest.raises(DataContractError, match=r"line 5: line beginning with '#' after the #CHROM header"):
        read_vcf(hashed)

    bad_index = tmp_path / "bad_index.vcf"
    fields = lines[3].split("\t")
    fields[9] = "2/0"  # ALT has one allele, so index 2 does not exist
    bad_index.write_text("\n".join([*lines[:3], "\t".join(fields), *lines[4:]]) + "\n", encoding="utf-8", newline="")
    with pytest.raises(DataContractError, match=r"line 4: GT allele index exceeds ALT count"):
        read_vcf(bad_index)

    early = tmp_path / "early.vcf"
    early.write_text("\n".join([lines[0], lines[2], *lines[1:]]) + "\n", encoding="utf-8", newline="")
    with pytest.raises(DataContractError, match="VCF data line before #CHROM header"):
        read_vcf(early)

    header_only = tmp_path / "header_only.vcf"
    header_only.write_text("\n".join(lines[:2]) + "\n", encoding="utf-8", newline="")
    with pytest.raises(DataContractError, match="VCF contains no variant records"):
        read_vcf(header_only)

    # No #CHROM anywhere: the file is not empty, and pass 2 says so rather than calling it empty.
    no_header = tmp_path / "no_header.vcf"
    no_header.write_text("\n".join([lines[0], *lines[2:]]) + "\n", encoding="utf-8", newline="")
    with pytest.raises(DataContractError, match="VCF data line before #CHROM header"):
        read_vcf(no_header)

    # A header followed only by a '#' line: the '#'-after-header message at its line, not "no records".
    hash_only = tmp_path / "hash_only.vcf"
    hash_only.write_text("\n".join([*lines[:2], "# a note"]) + "\n", encoding="utf-8", newline="")
    with pytest.raises(DataContractError, match=r"line 3: line beginning with '#' after the #CHROM header"):
        read_vcf(hash_only)


def test_a_malformed_header_never_outranks_an_earlier_error(tmp_path: Path) -> None:
    """The #CHROM shape is judged in pass 2, so an error on an earlier line still wins.

    Pass 1 must also keep counting records after a malformed header (its `seen_header` flag is
    explicit, not the truthiness of the sample ids), or a file that should be refused would
    instead parse into a short matrix.
    """
    bad_header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
    row = "\t".join(["Gm01", "1000", "m1", "A", "T", ".", ".", ".", "GT", "0/0"])

    data_first = tmp_path / "data_first.vcf"
    data_first.write_text("\n".join(["##fileformat=VCFv4.2", row, bad_header, row]) + "\n", encoding="utf-8", newline="")
    with pytest.raises(DataContractError, match="VCF data line before #CHROM header"):
        read_vcf(data_first)

    note_first = tmp_path / "note_first.vcf"
    note_first.write_text("\n".join(["##fileformat=VCFv4.2", "#note", bad_header, row]) + "\n", encoding="utf-8", newline="")
    with pytest.raises(DataContractError, match="VCF data line before #CHROM header"):
        read_vcf(note_first)

    header_first = tmp_path / "header_first.vcf"
    header_first.write_text("\n".join(["##fileformat=VCFv4.2", bad_header, row]) + "\n", encoding="utf-8", newline="")
    with pytest.raises(DataContractError, match="VCF header must have FORMAT and at least one sample column"):
        read_vcf(header_first)


def test_fill_judges_a_missing_header_before_it_uses_the_sample_ids(tmp_path: Path) -> None:
    """With no header there are no sample ids: `_fill` must reach the structural raise anyway."""
    p = tmp_path / "no_header.vcf"
    _write_vcf(p, 3, 3, np.random.default_rng(6))
    body = [line for line in p.read_text(encoding="utf-8").rstrip("\n").split("\n") if not line.startswith("#CHROM")]
    p.write_text("\n".join(body) + "\n", encoding="utf-8", newline="")
    calls = np.empty((0, 0, 2), dtype=np.int8)
    with pytest.raises(DataContractError, match="VCF data line before #CHROM header"):
        vcf_module._fill(p, [], calls, SOYBEAN)


def test_a_file_that_grows_between_the_passes_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The guard before `calls[filled] = gts`: a contract error, not an `IndexError`."""
    p = tmp_path / "g.vcf"
    _write_vcf(p, 6, 4, np.random.default_rng(7))
    real_scan = vcf_module._scan

    def _one_short(path: Path) -> tuple[list[str], int]:
        sample_ids, n_records = real_scan(path)
        return sample_ids, n_records - 1

    monkeypatch.setattr(vcf_module, "_scan", _one_short)
    with pytest.raises(DataContractError, match="VCF changed while it was being read"):
        read_vcf(p)


def test_a_file_that_shrinks_between_the_passes_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The check after the loop: without it, `np.empty` garbage rows would ship as genotypes."""
    p = tmp_path / "g.vcf"
    _write_vcf(p, 6, 4, np.random.default_rng(8))
    real_scan = vcf_module._scan

    def _one_long(path: Path) -> tuple[list[str], int]:
        sample_ids, n_records = real_scan(path)
        return sample_ids, n_records + 1

    monkeypatch.setattr(vcf_module, "_scan", _one_long)
    with pytest.raises(DataContractError, match="VCF changed while it was being read"):
        read_vcf(p)


def test_a_missing_file_still_raises_file_not_found(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        read_vcf(tmp_path / "absent.vcf")


def test_truncated_gz_is_a_contract_error(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Pass 1 raises: `_scan` reads the whole compressed stream to count records, so the
    truncation surfaces there and `_fill` is never entered. Pinned so a later reordering
    cannot move the error silently."""
    p = tmp_path / "g.vcf"
    _write_vcf(p, 200, 20, np.random.default_rng(5))
    gz = tmp_path / "g.vcf.gz"
    with gzip.open(gz, "wb") as out:
        out.write(p.read_bytes())
    whole = gz.read_bytes()
    gz.write_bytes(whole[: len(whole) * 2 // 3])

    def _no_pass_two(*args: object, **kwargs: object) -> None:
        pytest.fail("pass 2 was entered; the truncation should surface in _scan")

    monkeypatch.setattr(vcf_module, "_fill", _no_pass_two)
    with pytest.raises(DataContractError, match="cannot read"):
        read_vcf(gz)
