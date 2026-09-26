"""BrAPI v2.1 allele-matrix loader (docs/adr/0024; mirrors backcross src/io/brapi.ts and its ADR 0015).

Responsibility: page one variant set's ``/callsets``, ``/variants`` and ``/allelematrix`` (GT only)
from a BrAPI v2.1 Genotyping server into the same ``Dataset`` a genotype file produces, so
samples.csv, markers.csv and every downstream check apply unchanged. Cells never go through
``io/calls.py`` and positions never through ``io/position.py``: GT tokens are allele indices read by
``parse_brapi_call``, and a position is BrAPI's 0-based ``start`` plus 1. HTTP lives in
``urllib_fetch_json`` alone, and every other function takes the transport as an argument, which is
what lets the tests run with no network. An access token is sent in the Authorization header only
and never appears in a URL, message, warning or output.

Format facts: BrAPI v2.1, ``Specification/BrAPI-Genotyping/AlleleMatrix/`` and
``Specification/BrAPI-Genotyping/Variants/Schemas/Variant.yaml`` (docs/reference-repos.md), where
``start`` is 0-based and ``end`` exclusive against VCF's 1-based POS.

Interface:
    BrapiSource, FetchJson, BRAPI_DEFAULT_PAGE_SIZE_VARIANTS, BRAPI_DEFAULT_PAGE_SIZE_CALL_SETS,
    BRAPI_REQUEST_TIMEOUT_S, BRAPI_MAX_CALLS
    normalise_base_url(raw) -> str
    parse_brapi_call(token, sep_phased, sep_unphased, unknown) -> tuple[int, int]
    variant_position(variant, scheme) -> tuple[str, int] | None
    marker_id_of(variant) -> str
    assign_sample_ids(raw_call_sets) -> (list[CallSetRef], warnings)
    fetch_call_sets(source, fetch_json) -> (list[CallSetRef], warnings)
    fetch_brapi_genotypes(source, fetch_json, scheme=SOYBEAN, marker_map=None)
        -> (GenotypeMatrix, list[CallSetRef], warnings)
    urllib_fetch_json(url, headers) -> dict
    load_brapi_dataset(source, samples_path, markers_path=None, crop="soybean", fetch_json=None) -> Dataset
    call_sets_csv_text(call_sets) -> str
"""

from __future__ import annotations

import csv
import http.client
import json
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

import numpy as np

from progeny_selector.core.chrom import SOYBEAN, CompiledScheme, normalize_chrom
from progeny_selector.io.crops import DEFAULT_CROP_ID, resolve_crop
from progeny_selector.io.manifest import apply_marker_map, build_dataset, read_markers, read_samples
from progeny_selector.model.dataset import CallSetRef, DataContractError, Dataset, GenotypeMatrix, Marker

BRAPI_DEFAULT_PAGE_SIZE_VARIANTS = 1000
BRAPI_DEFAULT_PAGE_SIZE_CALL_SETS = 500
BRAPI_REQUEST_TIMEOUT_S = 60
# The exact call count of the largest case in docs/limits.md that completed on CPython
# (row 6k_x_2000): 6,000 markers x 2,002 call sets, the 2,000 progeny plus both parents,
# a 22.9 MiB resident matrix. Set to that count and not below it, so the largest measured
# case does not warn. Above it the load is warned about and still attempted, because no
# server-sourced case has been measured.
BRAPI_MAX_CALLS = 6_000 * 2_002

MISSING = -1
MAX_ALLELE_INDEX = 127  # np.int8

#: ``(url, headers) -> parsed JSON body``; raises :class:`DataContractError` on any transport failure.
FetchJson = Callable[[str, dict[str, str]], dict]

CALL_SET_COLUMNS = ("sample_id", "call_set_name", "call_set_db_id", "sample_db_id")


@dataclass(frozen=True)
class BrapiSource:
    """One variant set on one server. ``base_url`` is what ``normalise_base_url`` returned."""

    base_url: str
    variant_set_db_id: str
    token: str | None = None
    page_size_variants: int = BRAPI_DEFAULT_PAGE_SIZE_VARIANTS
    page_size_call_sets: int = BRAPI_DEFAULT_PAGE_SIZE_CALL_SETS


def normalise_base_url(raw: str) -> str:
    """Trim, strip trailing slashes and refuse anything but http(s) without userinfo (D7)."""
    base = raw.strip().rstrip("/")
    parsed = urllib.parse.urlsplit(base)
    if parsed.scheme not in ("http", "https") or not parsed.netloc:
        raise DataContractError("BrAPI: base URL must start with http:// or https://")
    if parsed.username or parsed.password:
        raise DataContractError("BrAPI: base URL must not contain a user name or password; use --brapi-token-env")
    # Endpoint paths are appended to the base, so a query or fragment would end up in the middle
    # of every request URL. Neither is echoed: a query is where a careless user puts a token.
    if parsed.query or parsed.fragment or "?" in base or "#" in base:
        raise DataContractError("BrAPI: base URL must not contain a query string (?...) or a fragment (#...)")
    return base


def _text(value: object) -> str:
    return value.strip() if isinstance(value, str) else ""


def _total_pages(value: object, where: str) -> int:
    """``totalPages`` as a page count; absent or 0 reads as one page, a non-integer is an error.

    A float is accepted only when it equals its floor, so ``2.0`` is two pages and ``2.5`` is
    refused. Truncating silently would drop a page's worth of calls and leave them missing.
    """
    if value is None:
        return 1
    if isinstance(value, bool):
        raise DataContractError(f"BrAPI: {where}: totalPages {value!r} is not an integer")
    if isinstance(value, float):
        if not value.is_integer():
            raise DataContractError(f"BrAPI: {where}: totalPages {value!r} is not an integer")
        value = int(value)
    if not isinstance(value, int):
        raise DataContractError(f"BrAPI: {where}: totalPages {value!r} is not an integer")
    return max(value, 1)


def parse_brapi_call(token: str, sep_phased: str, sep_unphased: str, unknown: str) -> tuple[int, int]:
    """One GT token as a pair of allele indices; phase is ignored and order is kept (D3).

    A token equal to ``unknown`` (or empty) is missing, a token without a separator is haploid and
    reads as the homozygote, and a token missing on either side is missing as a whole.
    """
    text = token.strip()
    if text == "" or text == unknown:
        return (MISSING, MISSING)
    parts = [text]
    for sep in (sep_phased, sep_unphased):
        if sep:
            parts = [piece for part in parts for piece in part.split(sep)]
    if len(parts) > 2:
        raise DataContractError(f"BrAPI: invalid GT {token!r}")
    indices: list[int] = []
    for part in parts:
        if part == unknown:
            indices.append(MISSING)
            continue
        if not part.isdigit():
            raise DataContractError(f"BrAPI: invalid GT {token!r}")
        index = int(part)
        if index >= MAX_ALLELE_INDEX:
            raise DataContractError(f"BrAPI: invalid GT {token!r}")
        indices.append(index)
    a = indices[0]
    b = indices[1] if len(indices) == 2 else a
    if a == MISSING or b == MISSING:
        return (MISSING, MISSING)
    return (a, b)


def _start_as_int(value: object, marker_id: str) -> int | None:
    """``start`` as an integer; ``None`` when absent, an error when present and not integral (D1)."""
    if value is None:
        return None
    if isinstance(value, bool):
        raise DataContractError(f"BrAPI: variant {marker_id} has a non-integer start")
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    raise DataContractError(f"BrAPI: variant {marker_id} has a non-integer start")


def variant_position(variant: dict, scheme: CompiledScheme) -> tuple[str, int] | None:
    """``(chrom, pos_bp)`` from ``referenceName`` and ``start + 1``; ``None`` when either is absent (D1)."""
    marker_id = marker_id_of(variant)
    start = _start_as_int(variant.get("start"), marker_id)
    reference_name = _text(variant.get("referenceName"))
    if reference_name == "" or start is None or start < 0:
        return None
    return (normalize_chrom(reference_name, scheme), start + 1)


def marker_id_of(variant: dict) -> str:
    """The first usable ``variantNames`` entry, else ``variantDbId`` (D2); ``.`` counts as absent."""
    names = variant.get("variantNames")
    candidates = [names] if isinstance(names, str) else names if isinstance(names, list) else []
    for name in candidates:
        text = _text(name)
        if text not in ("", "."):
            return text
    return _text(variant.get("variantDbId"))


def assign_sample_ids(raw_call_sets: list[dict]) -> tuple[list[CallSetRef], list[str]]:
    """``callSetName`` when non-empty and unique, else ``callSetDbId`` (D4)."""
    db_ids_by_name: dict[str, list[str]] = {}
    for row in raw_call_sets:
        name = _text(row.get("callSetName"))
        if name:
            db_ids_by_name.setdefault(name, []).append(_text(row.get("callSetDbId")))
    call_sets: list[CallSetRef] = []
    for row in raw_call_sets:
        name = _text(row.get("callSetName"))
        db_id = _text(row.get("callSetDbId"))
        unique = name != "" and len(db_ids_by_name[name]) == 1
        call_sets.append(
            CallSetRef(
                sample_id=name if unique else db_id,
                call_set_name=name,
                call_set_db_id=db_id,
                sample_db_id=_text(row.get("sampleDbId")),
            )
        )
    warnings = [
        f"BrAPI: {len(ids)} call sets share the name {name!r}; loaded under their callSetDbId: {', '.join(ids)}"
        for name, ids in db_ids_by_name.items()
        if len(ids) > 1
    ]
    seen: set[str] = set()
    for ref in call_sets:
        if ref.sample_id in seen:
            raise DataContractError(f"BrAPI: sample id {ref.sample_id!r} is not unique after falling back to callSetDbId")
        seen.add(ref.sample_id)
    return call_sets, warnings


class _Session:
    """One base URL, one variant set and one transport; every request goes through ``get``."""

    def __init__(self, source: BrapiSource, fetch_json: FetchJson) -> None:
        self.base = normalise_base_url(source.base_url)
        self.variant_set_db_id = source.variant_set_db_id.strip()
        if self.variant_set_db_id == "":
            raise DataContractError("BrAPI: variant set id is required")
        self.page_size_variants = int(source.page_size_variants)
        self.page_size_call_sets = int(source.page_size_call_sets)
        token = (source.token or "").strip()
        self.headers = {"Accept": "application/json"}
        if token:
            self.headers["Authorization"] = f"Bearer {token}"
        self._fetch_json = fetch_json

    def url(self, endpoint: str, params: list[tuple[str, str]]) -> str:
        return f"{self.base}/{endpoint}?{urllib.parse.urlencode(params)}"

    def get(self, url: str) -> dict:
        body = self._fetch_json(url, dict(self.headers))
        if not isinstance(body, dict):
            raise DataContractError(f"BrAPI: GET {url}: response is not a JSON object")
        result = body.get("result")
        if not isinstance(result, dict):
            raise DataContractError(f"BrAPI: GET {url}: response has no result object")
        return body


def _fetch_list(session: _Session, endpoint: str, page_size: int) -> list[dict]:
    """Page ``/callsets`` or ``/variants`` (D6); ``/variants`` follows ``nextPageToken``."""
    rows: list[dict] = []
    seen_ids: set[str] = set()
    id_key = "variantDbId" if endpoint == "variants" else "callSetDbId"
    total_pages = 1
    next_page_token = ""
    page = 0
    while page < total_pages:
        params: list[tuple[str, str]] = [("variantSetDbId", session.variant_set_db_id)]
        if endpoint == "variants" and page > 0 and next_page_token != "":
            params.append(("pageToken", next_page_token))
        else:
            params.append(("page", str(page)))
            if endpoint == "variants" and page > 0:
                params.append(("pageToken", str(page)))
        params.append(("pageSize", str(page_size)))
        url = session.url(endpoint, params)
        body = session.get(url)
        metadata = body.get("metadata")
        pagination = metadata.get("pagination") if isinstance(metadata, dict) else None
        pagination = pagination if isinstance(pagination, dict) else {}
        if page == 0:
            total_pages = _total_pages(pagination.get("totalPages"), f"/{endpoint} page 0")
        next_page_token = _text(pagination.get("nextPageToken"))
        data = body["result"].get("data")
        if not isinstance(data, list):
            raise DataContractError(f"BrAPI: GET {url}: result.data is not an array")
        if page > 0 and not data:
            break
        for row in data:
            if not isinstance(row, dict):
                raise DataContractError(f"BrAPI: /{endpoint} page {page} has a row that is not an object")
            db_id = _text(row.get(id_key))
            if db_id == "":
                raise DataContractError(f"BrAPI: /{endpoint} page {page} has a row with no {id_key}")
            if db_id in seen_ids:
                raise DataContractError(f"BrAPI: /{endpoint} page {page} repeated {id_key} {db_id!r}")
            seen_ids.add(db_id)
            rows.append(row)
        page += 1
    return rows


def fetch_call_sets(source: BrapiSource, fetch_json: FetchJson) -> tuple[list[CallSetRef], list[str]]:
    """Every call set of the variant set, in server order, with their sample ids (D4, D6)."""
    session = _Session(source, fetch_json)
    rows = _fetch_list(session, "callsets", session.page_size_call_sets)
    if not rows:
        raise DataContractError(f"BrAPI: variant set {session.variant_set_db_id} has no call sets")
    return assign_sample_ids(rows)


def _allele_symbols(variant: dict) -> list[str] | None:
    """``[referenceBases, *alternateBases]``, or ``None`` when the server gave no bases (D3)."""
    reference_bases = _text(variant.get("referenceBases"))
    if reference_bases == "":
        return None
    alternate = variant.get("alternateBases")
    alternate = alternate if isinstance(alternate, list) else []
    return [reference_bases, *[_text(base) or str(k + 1) for k, base in enumerate(alternate)]]


def _read_variants(
    session: _Session, scheme: CompiledScheme, marker_map: dict[str, Marker] | None
) -> tuple[list[Marker], list[list[str] | None], dict[str, int]]:
    rows = _fetch_list(session, "variants", session.page_size_variants)
    if not rows:
        raise DataContractError(f"BrAPI: variant set {session.variant_set_db_id} has no variants")
    markers: list[Marker] = []
    symbols: list[list[str] | None] = []
    row_by_variant_db_id: dict[str, int] = {}
    for row in rows:
        marker_id = marker_id_of(row)
        position = variant_position(row, scheme)
        if position is None:
            mapped = marker_map.get(marker_id) if marker_map else None
            if mapped is None:
                raise DataContractError(f"BrAPI: variant {marker_id} has no position and markers.csv does not list it")
            position = (mapped.chrom, mapped.pos_bp)
        markers.append(Marker(marker_id=marker_id, chrom=position[0], pos_bp=position[1]))
        symbols.append(_allele_symbols(row))
        row_by_variant_db_id[_text(row.get("variantDbId"))] = len(markers) - 1
    return markers, symbols, row_by_variant_db_id


def _gt_matrix(result: dict, page_label: str) -> dict:
    matrices = result.get("dataMatrices")
    matrices = matrices if isinstance(matrices, list) else []
    for matrix in matrices:
        if isinstance(matrix, dict) and matrix.get("dataMatrixAbbreviation") == "GT":
            return matrix
    raise DataContractError(f"BrAPI: no GT data matrix on /allelematrix page {page_label}")


def _string_list(result: dict, key: str, page_label: str) -> list[str]:
    value = result.get(key)
    if not isinstance(value, list) or not all(isinstance(v, str) for v in value):
        raise DataContractError(f"BrAPI: /allelematrix page {page_label}: {key} is not an array of strings")
    ids = [v.strip() for v in value]
    seen: set[str] = set()
    for db_id in ids:
        if db_id in seen:
            raise DataContractError(f"BrAPI: /allelematrix page {page_label} repeated {key[:-1]} {db_id!r}")
        seen.add(db_id)
    return ids


def _dimension_pages(result: dict, dimension: str) -> int:
    pagination = result.get("pagination")
    pagination = pagination if isinstance(pagination, list) else []
    for entry in pagination:
        if isinstance(entry, dict) and entry.get("dimension") == dimension:
            return _total_pages(entry.get("totalPages"), f"/allelematrix {dimension}")
    return 1


def _read_matrix_page(
    result: dict,
    page_label: str,
    calls: np.ndarray,
    markers: list[Marker],
    symbols: list[list[str] | None],
    row_by_variant_db_id: dict[str, int],
    col_by_call_set_db_id: dict[str, int],
    max_index: list[int],
    seen_pages: set[tuple[tuple[str, ...], tuple[str, ...]]],
    covered_rows: np.ndarray,
    covered_cols: np.ndarray,
) -> None:
    matrix = _gt_matrix(result, page_label)
    variant_db_ids = _string_list(result, "variantDbIds", page_label)
    call_set_db_ids = _string_list(result, "callSetDbIds", page_label)
    sep_phased = _text(result.get("sepPhased")) or "|"
    sep_unphased = _text(result.get("sepUnphased")) or "/"
    unknown = _text(result.get("unknownString")) or "."
    data = matrix.get("dataMatrix")
    if (
        not isinstance(data, list)
        or len(data) != len(variant_db_ids)
        or not all(isinstance(row, list) and len(row) == len(call_set_db_ids) for row in data)
    ):
        raise DataContractError(f"BrAPI: /allelematrix page {page_label}: dataMatrix is not {len(variant_db_ids)} x {len(call_set_db_ids)}")
    rows = []
    for variant_db_id in variant_db_ids:
        if variant_db_id not in row_by_variant_db_id:
            raise DataContractError(f"BrAPI: /allelematrix page {page_label} names variant {variant_db_id!r} that /variants did not list")
        rows.append(row_by_variant_db_id[variant_db_id])
    cols = []
    for call_set_db_id in call_set_db_ids:
        if call_set_db_id not in col_by_call_set_db_id:
            raise DataContractError(f"BrAPI: /allelematrix page {page_label} names call set {call_set_db_id!r} that /callsets did not list")
        cols.append(col_by_call_set_db_id[call_set_db_id])
    # A page is identified by the ids it carries, not by the cells it would write. Comparing
    # written cells against the missing sentinel cannot tell a repeat carrying missing calls from
    # a first reading, and a server that ignores a dimension page parameter then silently loads
    # one page over the whole matrix.
    signature = (tuple(variant_db_ids), tuple(call_set_db_ids))
    if signature in seen_pages:
        raise DataContractError(f"BrAPI: /allelematrix page {page_label} repeated a page already read")
    seen_pages.add(signature)
    covered_rows[rows] = True
    covered_cols[cols] = True
    for i, row_index in enumerate(rows):
        marker = markers[row_index]
        allele_symbols = symbols[row_index]
        for j, col_index in enumerate(cols):
            cell = data[i][j]
            if cell is None:
                continue
            if not isinstance(cell, str):
                raise DataContractError(
                    f"BrAPI: /allelematrix page {page_label}: variant {marker.marker_id} has a call that is not a string"
                )
            a, b = parse_brapi_call(cell, sep_phased, sep_unphased, unknown)
            if a == MISSING:
                continue
            if allele_symbols is not None and max(a, b) >= len(allele_symbols):
                raise DataContractError(
                    f"BrAPI: variant {marker.marker_id}: allele index {max(a, b)} but only {len(allele_symbols)} allele(s) listed"
                )
            calls[row_index, col_index, 0] = a
            calls[row_index, col_index, 1] = b
            max_index[row_index] = max(max_index[row_index], a, b)


def _coverage_warnings(covered: np.ndarray, names: list[str], kind: str) -> list[str]:
    """Name what no /allelematrix page carried; those rows or columns stayed entirely missing."""
    absent = [name for name, seen in zip(names, covered.tolist(), strict=True) if not seen]
    if not absent:
        return []
    listed = ", ".join(absent[:5]) + (", ..." if len(absent) > 5 else "")
    return [f"BrAPI: {len(absent)} {kind}(s) absent from every /allelematrix page were loaded with all calls missing: {listed}"]


def fetch_brapi_genotypes(
    source: BrapiSource,
    fetch_json: FetchJson,
    scheme: CompiledScheme = SOYBEAN,
    marker_map: dict[str, Marker] | None = None,
) -> tuple[GenotypeMatrix, list[CallSetRef], list[str]]:
    """The variant set as a ``GenotypeMatrix``, its call sets in server order, and the warnings."""
    session = _Session(source, fetch_json)
    call_set_rows = _fetch_list(session, "callsets", session.page_size_call_sets)
    if not call_set_rows:
        raise DataContractError(f"BrAPI: variant set {session.variant_set_db_id} has no call sets")
    call_sets, warnings = assign_sample_ids(call_set_rows)
    markers, symbols, row_by_variant_db_id = _read_variants(session, scheme, marker_map)
    n = len(markers) * len(call_sets)
    if n > BRAPI_MAX_CALLS:
        warnings.append(
            f"BrAPI: {len(markers)} variants x {len(call_sets)} call sets = {n} calls, "
            f"above the largest measured CPython case {BRAPI_MAX_CALLS} (docs/limits.md)"
        )
    calls = np.full((len(markers), len(call_sets), 2), MISSING, dtype=np.int8)
    col_by_call_set_db_id = {ref.call_set_db_id: i for i, ref in enumerate(call_sets)}
    max_index = [-1] * len(markers)

    def matrix_url(v: int, c: int) -> str:
        return session.url(
            "allelematrix",
            [
                ("variantSetDbId", session.variant_set_db_id),
                ("dataMatrixAbbreviations", "GT"),
                ("dimensionVariantPage", str(v)),
                ("dimensionVariantPageSize", str(session.page_size_variants)),
                ("dimensionCallSetPage", str(c)),
                ("dimensionCallSetPageSize", str(session.page_size_call_sets)),
            ],
        )

    seen_pages: set[tuple[tuple[str, ...], tuple[str, ...]]] = set()
    covered_rows = np.zeros(len(markers), dtype=bool)
    covered_cols = np.zeros(len(call_sets), dtype=bool)
    page0 = session.get(matrix_url(0, 0))["result"]
    variant_pages = _dimension_pages(page0, "VARIANTS")
    call_set_pages = _dimension_pages(page0, "CALLSETS")
    first: dict | None = page0
    del page0
    for v in range(variant_pages):
        for c in range(call_set_pages):
            # Page (0, 0) is handed over and its name dropped, and each page's name is dropped
            # once it is written in, so no earlier page is referenced while the next one is
            # fetched: the peak is the matrix plus one page.
            if first is not None:
                result, first = first, None
            else:
                result = session.get(matrix_url(v, c))["result"]
            _read_matrix_page(
                result,
                f"{v},{c}",
                calls,
                markers,
                symbols,
                row_by_variant_db_id,
                col_by_call_set_db_id,
                max_index,
                seen_pages,
                covered_rows,
                covered_cols,
            )
            del result
    warnings.extend(_coverage_warnings(covered_rows, [m.marker_id for m in markers], "variant"))
    warnings.extend(_coverage_warnings(covered_cols, [c.call_set_db_id for c in call_sets], "call set"))
    alleles = [s if s is not None else [str(k) for k in range(max(max_index[i] + 1, 1))] for i, s in enumerate(symbols)]
    gm = GenotypeMatrix(
        markers=markers,
        sample_ids=[ref.sample_id for ref in call_sets],
        alleles=alleles,
        calls=calls,
        coded=False,
    )
    return gm, call_sets, warnings


class _StripAuthOnCrossHostRedirect(urllib.request.HTTPRedirectHandler):
    """Drop ``Authorization`` when a redirect leaves the scheme or host the token was meant for.

    ``urllib.request`` copies ``Request.headers`` onto the redirected request, so without this a
    302 from the configured server to any other host would carry the bearer token with it.
    Browsers strip the header on a cross-origin redirect, which is why the sibling's fetch-based
    loader needs no equivalent.
    """

    def redirect_request(self, req, fp, code, msg, headers, newurl):  # type: ignore[no-untyped-def]
        new = super().redirect_request(req, fp, code, msg, headers, newurl)
        if new is None:
            return None
        before = urllib.parse.urlsplit(req.full_url)
        after = urllib.parse.urlsplit(new.full_url)
        if (before.scheme, before.netloc) != (after.scheme, after.netloc):
            # Both spellings: Request.add_header title-cases the name, callers may not.
            for name in ("Authorization", "authorization"):
                new.headers.pop(name, None)
                new.unredirected_hdrs.pop(name, None)
        return new


def urllib_fetch_json(url: str, headers: dict[str, str]) -> dict:
    """The CPython transport: one GET, ``BRAPI_REQUEST_TIMEOUT_S`` seconds, JSON body (D7)."""
    request = urllib.request.Request(url, headers=headers, method="GET")
    opener = urllib.request.build_opener(_StripAuthOnCrossHostRedirect)
    try:
        with opener.open(request, timeout=BRAPI_REQUEST_TIMEOUT_S) as response:
            raw = response.read()
    except urllib.error.HTTPError as exc:
        # The token is in the header the caller built, never in the URL or the code.
        raise DataContractError(f"BrAPI: GET {url} returned HTTP {exc.code}") from None
    except (urllib.error.URLError, TimeoutError, OSError, http.client.HTTPException) as exc:
        # HTTPException covers IncompleteRead, raised by response.read() when the connection
        # closes early; it is not an OSError and would otherwise escape as a traceback.
        reason = getattr(exc, "reason", exc)
        raise DataContractError(f"BrAPI: could not reach {url}: {reason}") from None
    try:
        body = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, ValueError):
        raise DataContractError(f"BrAPI: GET {url}: response is not a JSON object") from None
    if not isinstance(body, dict):
        raise DataContractError(f"BrAPI: GET {url}: response is not a JSON object")
    return body


def load_brapi_dataset(
    source: BrapiSource,
    samples_path: str | Path,
    markers_path: str | Path | None = None,
    crop: str | None = DEFAULT_CROP_ID,
    fetch_json: FetchJson | None = None,
) -> Dataset:
    """A BrAPI variant set as the ``Dataset`` a file produces; roles come from samples.csv."""
    scheme = resolve_crop(crop)
    marker_map = read_markers(markers_path, scheme=scheme) if markers_path is not None else None
    # Resolved here, so a test that replaces the module's transport is honoured.
    transport: FetchJson = fetch_json if fetch_json is not None else urllib_fetch_json
    gm, call_sets, warnings = fetch_brapi_genotypes(source, transport, scheme, marker_map)
    map_warnings: list[str] = []
    if marker_map is not None:
        gm, map_warnings = apply_marker_map(gm, marker_map)
    try:
        dataset = build_dataset(gm, read_samples(samples_path))
    except DataContractError as exc:
        # A sample "absent from the genotype file" is usually a call set whose name collided and
        # was therefore loaded under its callSetDbId. Say so, or the user cannot act on it.
        collisions = [w for w in warnings if "call sets share the name" in w]
        if collisions and "samples in manifest but not in genotype file" in str(exc):
            raise DataContractError(f"{exc}; note: {'; '.join(collisions)}") from exc
        raise
    dataset.warnings = warnings + map_warnings + dataset.warnings
    dataset.scheme = scheme
    dataset.token_profile = "default"
    dataset.call_sets = tuple(call_sets)
    return dataset


def call_sets_csv_text(call_sets: list[CallSetRef]) -> str:
    """The ``brapi-callsets.csv`` table, CRLF as the other writers (D5)."""
    buf = StringIO(newline="")
    writer = csv.writer(buf)
    writer.writerow(CALL_SET_COLUMNS)
    for ref in call_sets:
        writer.writerow([ref.sample_id, ref.call_set_name, ref.call_set_db_id, ref.sample_db_id])
    return buf.getvalue()


__all__ = [
    "BRAPI_DEFAULT_PAGE_SIZE_CALL_SETS",
    "BRAPI_DEFAULT_PAGE_SIZE_VARIANTS",
    "BRAPI_MAX_CALLS",
    "BRAPI_REQUEST_TIMEOUT_S",
    "BrapiSource",
    "FetchJson",
    "assign_sample_ids",
    "call_sets_csv_text",
    "fetch_brapi_genotypes",
    "fetch_call_sets",
    "load_brapi_dataset",
    "marker_id_of",
    "normalise_base_url",
    "parse_brapi_call",
    "urllib_fetch_json",
    "variant_position",
]
