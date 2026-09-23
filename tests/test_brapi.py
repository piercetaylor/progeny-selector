"""BrAPI v2.1 loader (src/progeny_selector/io/brapi.py, docs/adr/0024).

Every test here runs with no network: the transport is a stub that serves the recorded pages
under tests/fixtures/brapi/, which scripts/make_fixture.py generates from the BC2F1 dataset.
The acceptance test compares a BrAPI load against read_vcf of the same 25 variants and 8 samples.
"""

from __future__ import annotations

import copy
import json
import urllib.error
import urllib.request
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

import numpy as np
import pytest

from progeny_selector.cli import main
from progeny_selector.core.chrom import SOYBEAN
from progeny_selector.core.foreground import resolve_locus
from progeny_selector.io import brapi
from progeny_selector.io.brapi import (
    BrapiSource,
    assign_sample_ids,
    call_sets_csv_text,
    fetch_brapi_genotypes,
    fetch_call_sets,
    load_brapi_dataset,
    marker_id_of,
    parse_brapi_call,
    variant_position,
)
from progeny_selector.io.vcf import read_vcf
from progeny_selector.model.criteria import CriteriaError, TargetSpec
from progeny_selector.model.dataset import CallSetRef, DataContractError

BRAPI_DIR = Path(__file__).parent / "fixtures" / "brapi"
BASE_URL = "https://example.org/brapi/v2"
N_VARIANTS = 25
N_CALL_SETS = 8


def _source(**kwargs) -> BrapiSource:
    defaults = {
        "base_url": BASE_URL,
        "variant_set_db_id": "vs1",
        "page_size_variants": 13,
        "page_size_call_sets": 5,
    }
    defaults.update(kwargs)
    return BrapiSource(**defaults)  # type: ignore[arg-type]


def _stub(fixture_dir: Path, recorded: list[str], variants_prefix: str = "variants"):
    """A FetchJson serving the recorded pages, checking the query the loader sent."""

    def fetch_json(url: str, headers: dict[str, str]) -> dict:
        assert headers["Accept"] == "application/json"
        recorded.append(url)
        parts = urlsplit(url)
        query = {k: v[0] for k, v in parse_qs(parts.query).items()}
        assert query["variantSetDbId"] == "vs1"
        endpoint = parts.path.rsplit("/", 1)[-1]
        if endpoint == "callsets":
            name = f"callsets.p{query['page']}.json"
        elif endpoint == "variants":
            page = query["page"]
            if page != "0":
                assert "pageToken" in query
            name = f"{variants_prefix}.p{page}.json"
        elif endpoint == "allelematrix":
            assert query["dataMatrixAbbreviations"] == "GT"
            name = f"allelematrix.v{query['dimensionVariantPage']}.c{query['dimensionCallSetPage']}.json"
        else:
            raise AssertionError(f"unexpected endpoint {endpoint}")
        return json.loads((fixture_dir / name).read_text(encoding="utf-8"))

    return fetch_json


def _mutating_stub(fixture_dir: Path, mutate):
    """The stub above with one page rewritten by ``mutate(name, body)`` before it is returned."""
    recorded: list[str] = []
    inner = _stub(fixture_dir, recorded)

    def fetch_json(url: str, headers: dict[str, str]) -> dict:
        body = inner(url, headers)
        name = urlsplit(url).path.rsplit("/", 1)[-1]
        return mutate(name, copy.deepcopy(body))

    return fetch_json


def test_brapi_equals_vcf_on_the_fixture(fixture_dir: Path) -> None:
    recorded: list[str] = []
    gm, call_sets, warnings = fetch_brapi_genotypes(_source(), _stub(BRAPI_DIR, recorded))
    assert warnings == []
    assert len(call_sets) == N_CALL_SETS
    vcf = read_vcf(fixture_dir / "genotypes.vcf")
    marker_ids = [m.marker_id for m in gm.markers]
    assert len(marker_ids) == N_VARIANTS
    rows = [vcf.marker_index(mid) for mid in marker_ids]
    cols = [vcf.sample_index(sid) for sid in gm.sample_ids]
    assert [vcf.markers[i].chrom for i in rows] == [m.chrom for m in gm.markers]
    assert [vcf.markers[i].pos_bp for i in rows] == [m.pos_bp for m in gm.markers]
    assert [vcf.alleles[i] for i in rows] == gm.alleles
    assert np.array_equal(vcf.calls[np.ix_(rows, cols)], gm.calls)
    endpoints = [urlsplit(u).path.rsplit("/", 1)[-1] for u in recorded]
    assert endpoints == ["callsets"] * 2 + ["variants"] * 2 + ["allelematrix"] * 4
    matrix_pages = [
        (parse_qs(urlsplit(u).query)["dimensionVariantPage"][0], parse_qs(urlsplit(u).query)["dimensionCallSetPage"][0])
        for u in recorded
        if u.find("/allelematrix") >= 0
    ]
    assert matrix_pages == [("0", "0"), ("0", "1"), ("1", "0"), ("1", "1")]


def test_positions_are_one_based() -> None:
    assert variant_position({"referenceName": "Gm06", "start": 4999999, "end": 5000000}, SOYBEAN) == ("Gm06", 5000000)
    assert variant_position({"referenceName": "Gm06", "start": 0}, SOYBEAN) == ("Gm06", 1)
    assert variant_position({"referenceName": "Gm06", "start": 4999999.0}, SOYBEAN) == ("Gm06", 5000000)
    assert variant_position({"start": 5}, SOYBEAN) is None
    assert variant_position({"referenceName": "Gm06"}, SOYBEAN) is None
    for bad in (4999999.5, "4999999", True):
        with pytest.raises(DataContractError):
            variant_position({"referenceName": "Gm06", "start": bad}, SOYBEAN)


def test_region_resolves_the_shifted_marker() -> None:
    """The region is built from the recorded ``start``, not from the loaded ``pos_bp``.

    Taking both from the loaded marker would pass with or without the ``+ 1`` of D1.
    """
    recorded = json.loads((BRAPI_DIR / "variants.p0.json").read_text(encoding="utf-8"))["result"]["data"][7]
    chrom, start = recorded["referenceName"], recorded["start"]
    dataset = load_brapi_dataset(_source(), BRAPI_DIR / "samples.csv", fetch_json=_stub(BRAPI_DIR, []))
    gm = dataset.genotypes
    spec = TargetSpec(locus_id="R", chrom=chrom, start_bp=start + 1, end_bp=start + 1)
    resolved = resolve_locus(spec, gm)
    assert list(resolved.marker_idx) == [7]
    assert gm.markers[7].marker_id == recorded["variantNames"][0]
    # At the 0-based start itself the region is empty, which resolve_locus refuses.
    shifted = TargetSpec(locus_id="R", chrom=chrom, start_bp=start, end_bp=start)
    with pytest.raises(CriteriaError):
        resolve_locus(shifted, gm)


def test_parse_brapi_call() -> None:
    assert parse_brapi_call("1|0", "|", "/", ".") == (1, 0)
    assert parse_brapi_call("0/1", "|", "/", ".") == (0, 1)
    assert parse_brapi_call(".", "|", "/", ".") == (-1, -1)
    assert parse_brapi_call("./1", "|", "/", ".") == (-1, -1)
    assert parse_brapi_call("1", "|", "/", ".") == (1, 1)
    for bad in ("0/1/1", "a/b", "200/0"):
        with pytest.raises(DataContractError):
            parse_brapi_call(bad, "|", "/", ".")
    assert parse_brapi_call("NA", "|", "/", "NA") == (-1, -1)
    with pytest.raises(DataContractError):
        parse_brapi_call(".", "|", "/", "NA")


def test_sample_id_rule() -> None:
    raw = [{"callSetName": n, "callSetDbId": f"cs{k + 1}"} for k, n in enumerate(["A", "B", "B", ""])]
    call_sets, warnings = assign_sample_ids(raw)
    assert [c.sample_id for c in call_sets] == ["A", "cs2", "cs3", "cs4"]
    assert len(warnings) == 1
    assert "share the name 'B'" in warnings[0]
    collide = [{"callSetName": "cs2", "callSetDbId": "cs1"}, {"callSetName": "", "callSetDbId": "cs2"}]
    with pytest.raises(DataContractError, match="cs2"):
        assign_sample_ids(collide)


def test_marker_id_rule() -> None:
    assert marker_id_of({"variantNames": ["m1"], "variantDbId": "v1"}) == "m1"
    assert marker_id_of({"variantNames": ["", ".", "m2"], "variantDbId": "v1"}) == "m2"
    assert marker_id_of({"variantNames": [], "variantDbId": "v1"}) == "v1"
    assert marker_id_of({"variantDbId": "v1"}) == "v1"
    assert marker_id_of({"variantNames": "m3", "variantDbId": "v1"}) == "m3"


def test_nopos_uses_markers_csv(fixture_dir: Path) -> None:
    with_pos, _, _ = fetch_brapi_genotypes(_source(), _stub(BRAPI_DIR, []))
    marker_map_path = fixture_dir / "markers.csv"
    from progeny_selector.io.manifest import read_markers

    nopos, _, _ = fetch_brapi_genotypes(
        _source(), _stub(BRAPI_DIR, [], variants_prefix="variants-nopos"), marker_map=read_markers(marker_map_path)
    )
    assert [(m.chrom, m.pos_bp) for m in nopos.markers] == [(m.chrom, m.pos_bp) for m in with_pos.markers]
    with pytest.raises(DataContractError, match=with_pos.markers[0].marker_id):
        fetch_brapi_genotypes(_source(), _stub(BRAPI_DIR, [], variants_prefix="variants-nopos"))


def test_repeated_page_is_an_error() -> None:
    page0 = json.loads((BRAPI_DIR / "variants.p0.json").read_text(encoding="utf-8"))

    def fetch_json(url: str, headers: dict[str, str]) -> dict:
        path = urlsplit(url).path.rsplit("/", 1)[-1]
        if path == "variants":
            return copy.deepcopy(page0)
        return _stub(BRAPI_DIR, [])(url, headers)

    with pytest.raises(DataContractError, match="repeated"):
        fetch_brapi_genotypes(_source(), fetch_json)


def test_unknown_call_set_in_matrix_is_an_error() -> None:
    def mutate(name: str, body: dict) -> dict:
        if name == "allelematrix":
            body["result"]["callSetDbIds"][0] = "cs-nowhere"
        return body

    with pytest.raises(DataContractError, match="cs-nowhere"):
        fetch_brapi_genotypes(_source(), _mutating_stub(BRAPI_DIR, mutate))


def test_no_gt_matrix_is_an_error() -> None:
    def mutate(name: str, body: dict) -> dict:
        if name == "allelematrix":
            body["result"]["dataMatrices"][0]["dataMatrixAbbreviation"] = "DP"
        return body

    with pytest.raises(DataContractError, match="no GT data matrix"):
        fetch_brapi_genotypes(_source(), _mutating_stub(BRAPI_DIR, mutate))


def test_size_guard_warns(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(brapi, "BRAPI_MAX_CALLS", 10)
    _, _, warnings = fetch_brapi_genotypes(_source(), _stub(BRAPI_DIR, []))
    assert len(warnings) == 1
    assert warnings[0].startswith(f"BrAPI: {N_VARIANTS} variants x {N_CALL_SETS} call sets = {N_VARIANTS * N_CALL_SETS} calls")
    assert "docs/limits.md" in warnings[0]


def test_token_only_in_the_header_and_never_in_errors(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: list[urllib.request.Request] = []

    def fake_open(self, request, timeout=None):
        captured.append(request)
        raise urllib.error.HTTPError(request.full_url, 401, "Unauthorized", {}, None)

    # The transport builds its own opener so a cross-host redirect can be stripped of the
    # Authorization header, so the request is captured at OpenerDirector.open, not at urlopen.
    monkeypatch.setattr(urllib.request.OpenerDirector, "open", fake_open)
    source = _source(token="secret-token")
    with pytest.raises(DataContractError) as excinfo:
        fetch_call_sets(source, brapi.urllib_fetch_json)
    assert captured[0].get_header("Authorization") == "Bearer secret-token"
    message = str(excinfo.value)
    assert BASE_URL in message
    assert "401" in message
    assert "secret-token" not in message


def test_call_sets_csv_text() -> None:
    call_sets, warnings = fetch_call_sets(_source(), _stub(BRAPI_DIR, []))
    assert warnings == []
    text = call_sets_csv_text(call_sets)
    lines = text.split("\r\n")[:-1]
    assert lines[0] == "sample_id,call_set_name,call_set_db_id,sample_db_id"
    assert len(lines) == N_CALL_SETS + 1
    assert lines[1].startswith("RP_Williams82,RP_Williams82,cs0,smp0")
    assert call_sets[0] == CallSetRef("RP_Williams82", "RP_Williams82", "cs0", "smp0")


def test_cli_brapi(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, fixture_dir: Path) -> None:
    monkeypatch.setattr(brapi, "urllib_fetch_json", _stub(BRAPI_DIR, []))
    out = tmp_path / "results.csv"
    argv = [
        "rank",
        "--brapi-url",
        BASE_URL,
        "--variant-set",
        "vs1",
        "--samples",
        str(BRAPI_DIR / "samples.csv"),
        "--markers",
        str(fixture_dir / "markers.csv"),
        "--criteria",
        str(BRAPI_DIR / "criteria.yaml"),
        "--out",
        str(out),
    ]
    assert main(argv) == 0
    assert len(out.read_text(encoding="utf-8").strip().split("\n")) == 7  # header plus 6 progeny

    with pytest.raises(SystemExit) as exc:
        main([*argv[:1], "--genotypes", str(fixture_dir / "genotypes.vcf"), *argv[1:]])
    assert exc.value.code == 2
    with pytest.raises(SystemExit) as exc:
        main(["rank", "--samples", str(BRAPI_DIR / "samples.csv"), "--criteria", str(BRAPI_DIR / "criteria.yaml"), "--out", str(out)])
    assert exc.value.code == 2
    with pytest.raises(SystemExit) as exc:
        main([*argv, "--brapi-token-env", "PS_UNSET_VAR"])
    assert exc.value.code == 2

    table = tmp_path / "c.csv"
    assert main(["brapi-callsets", "--brapi-url", BASE_URL, "--variant-set", "vs1", "--out", str(table)]) == 0
    assert len(table.read_text(encoding="utf-8").strip().split("\n")) == N_CALL_SETS + 1


def test_repeated_matrix_page_is_an_error() -> None:
    """A server ignoring dimensionCallSetPage sends page c0 for every call-set page.

    The cells of that page are not all written on the second reading, because the fixture holds
    missing calls, so only the page's identity catches it.
    """

    def fetch_json(url: str, headers: dict[str, str]) -> dict:
        parts = urlsplit(url)
        if parts.path.endswith("/allelematrix"):
            query = {k: v[0] for k, v in parse_qs(parts.query).items()}
            name = f"allelematrix.v{query['dimensionVariantPage']}.c0.json"
            return json.loads((BRAPI_DIR / name).read_text(encoding="utf-8"))
        return _stub(BRAPI_DIR, [])(url, headers)

    with pytest.raises(DataContractError, match="repeated a page already read"):
        fetch_brapi_genotypes(_source(), fetch_json)


def test_uncovered_variants_are_warned_about() -> None:
    """A server reporting one variant page leaves the second page's variants entirely missing."""

    def mutate(name: str, body: dict) -> dict:
        if name == "allelematrix":
            for entry in body["result"]["pagination"]:
                if entry["dimension"] == "VARIANTS":
                    entry["totalPages"] = 1
        return body

    gm, _, warnings = fetch_brapi_genotypes(_source(), _mutating_stub(BRAPI_DIR, mutate))
    assert len(warnings) == 1
    assert warnings[0].startswith("BrAPI: 12 variant(s) absent from every /allelematrix page")
    assert bool((gm.calls[13:] == -1).all())
    assert not bool((gm.calls[:13] == -1).all())


def test_total_pages_must_be_an_integer() -> None:
    def mutate(name: str, body: dict) -> dict:
        if name == "allelematrix":
            for entry in body["result"]["pagination"]:
                if entry["dimension"] == "CALLSETS":
                    entry["totalPages"] = 2.5
        return body

    with pytest.raises(DataContractError, match="totalPages"):
        fetch_brapi_genotypes(_source(), _mutating_stub(BRAPI_DIR, mutate))

    def integral(name: str, body: dict) -> dict:
        if name == "allelematrix":
            for entry in body["result"]["pagination"]:
                entry["totalPages"] = float(entry["totalPages"])
        return body

    gm, _, warnings = fetch_brapi_genotypes(_source(), _mutating_stub(BRAPI_DIR, integral))
    assert warnings == []
    assert gm.n_markers == N_VARIANTS


def test_redirect_drops_the_authorization_header(monkeypatch: pytest.MonkeyPatch) -> None:
    """A 302 away from the configured host must not carry the bearer token."""
    handler = brapi._StripAuthOnCrossHostRedirect()
    original = urllib.request.Request(f"{BASE_URL}/callsets", headers={"Authorization": "Bearer secret-token"})

    class _Response:
        def __init__(self) -> None:
            self.headers: dict[str, str] = {}

        def read(self) -> bytes:
            return b""

    same_host = handler.redirect_request(original, _Response(), 302, "Found", {}, f"{BASE_URL}/callsets/2")
    assert same_host is not None
    assert same_host.get_header("Authorization") == "Bearer secret-token"

    other_host = handler.redirect_request(original, _Response(), 302, "Found", {}, "https://elsewhere.example/steal")
    assert other_host is not None
    assert other_host.get_header("Authorization") is None
    assert "secret-token" not in str(other_host.headers) + str(other_host.unredirected_hdrs)


def _one_cell(token: str, variant_patch: dict | None = None):
    """A stub whose page 0,0 holds ``token`` at the first cell, with the variant optionally patched."""

    def mutate(name: str, body: dict) -> dict:
        if name == "allelematrix" and body["result"]["variantDbIds"][0] == "var0":
            body["result"]["dataMatrices"][0]["dataMatrix"][0][0] = token
        if name == "variants" and variant_patch is not None and body["result"]["data"][0]["variantDbId"] == "var0":
            body["result"]["data"][0].update(variant_patch)
        return body

    return _mutating_stub(BRAPI_DIR, mutate)


def test_cell_paths_the_fixture_does_not_reach() -> None:
    # An allele index past the symbol list of a variant that declared its bases.
    with pytest.raises(DataContractError, match="allele index 3"):
        fetch_brapi_genotypes(_source(), _one_cell("3/3"))
    # A null cell is missing, not an error.
    gm, _, warnings = fetch_brapi_genotypes(_source(), _one_cell(None))  # type: ignore[arg-type]
    assert warnings == []
    assert list(gm.calls[0, 0]) == [-1, -1]
    # A cell that is not a string at all.
    with pytest.raises(DataContractError, match="not a string"):
        fetch_brapi_genotypes(_source(), _one_cell(7))  # type: ignore[arg-type]
    # A phased pair keeps its order, and a half-call is missing as a whole.
    gm, _, _ = fetch_brapi_genotypes(_source(), _one_cell("1|0"))
    assert list(gm.calls[0, 0]) == [1, 0]
    gm, _, _ = fetch_brapi_genotypes(_source(), _one_cell("./1"))
    assert list(gm.calls[0, 0]) == [-1, -1]
    # A multi-allelic token, against a variant carrying a second ALT.
    gm, _, _ = fetch_brapi_genotypes(_source(), _one_cell("1/2", {"alternateBases": ["C", "G"]}))
    assert list(gm.calls[0, 0]) == [1, 2]
    assert gm.alleles[0][1:] == ["C", "G"]


def test_matrix_shape_mismatch_is_an_error() -> None:
    def mutate(name: str, body: dict) -> dict:
        if name == "allelematrix":
            body["result"]["dataMatrices"][0]["dataMatrix"][0].pop()
        return body

    with pytest.raises(DataContractError, match="dataMatrix is not"):
        fetch_brapi_genotypes(_source(), _mutating_stub(BRAPI_DIR, mutate))


def test_markers_csv_overrides_a_server_position(tmp_path: Path) -> None:
    """markers.csv wins over the server, exactly as it wins over a genotype file."""
    server_only = load_brapi_dataset(_source(), BRAPI_DIR / "samples.csv", fetch_json=_stub(BRAPI_DIR, []))
    first = server_only.genotypes.markers[0]
    markers_csv = tmp_path / "markers.csv"
    lines = ["marker_id,chrom,pos_bp,cm"]
    for m in server_only.genotypes.markers:
        pos = m.pos_bp + 1000 if m.marker_id == first.marker_id else m.pos_bp
        lines.append(f"{m.marker_id},{m.chrom},{pos},{pos / 1e6 * 2.5}")
    markers_csv.write_text("\n".join(lines) + "\n", encoding="utf-8")
    overridden = load_brapi_dataset(_source(), BRAPI_DIR / "samples.csv", markers_csv, fetch_json=_stub(BRAPI_DIR, []))
    assert overridden.genotypes.markers[0].pos_bp == first.pos_bp + 1000
    assert any("overrides genotype file" in w for w in overridden.warnings)


def test_token_env_set_but_empty_is_a_usage_error(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr(brapi, "urllib_fetch_json", _stub(BRAPI_DIR, []))
    monkeypatch.setenv("PS_BLANK_TOKEN", "   ")
    with pytest.raises(SystemExit) as exc:
        main(
            [
                "brapi-callsets",
                "--brapi-url",
                BASE_URL,
                "--variant-set",
                "vs1",
                "--brapi-token-env",
                "PS_BLANK_TOKEN",
                "--out",
                str(tmp_path / "c.csv"),
            ]
        )
    assert exc.value.code == 2


def test_name_collision_is_explained_when_the_manifest_then_fails() -> None:
    """A colliding callSetName falls back to its callSetDbId, and samples.csv no longer matches.

    build_dataset's "absent from the genotype file" error alone would not tell the user why, so
    the collision warning is appended to it.
    """

    def mutate(name: str, body: dict) -> dict:
        if name == "callsets":
            for row in body["result"]["data"]:
                if row["callSetDbId"] == "cs3":
                    row["callSetName"] = "BC2F1-F1-001"
        return body

    with pytest.raises(DataContractError) as excinfo:
        load_brapi_dataset(_source(), BRAPI_DIR / "samples.csv", fetch_json=_mutating_stub(BRAPI_DIR, mutate))
    message = str(excinfo.value)
    assert "samples in manifest but not in genotype file" in message
    assert "2 call sets share the name 'BC2F1-F1-001'" in message
    assert "cs2, cs3" in message
