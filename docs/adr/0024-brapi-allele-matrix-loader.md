# BrAPI allele-matrix loader: one variant set on the command line, call-set names as sample ids, server ids on the Dataset

Status: accepted. Date: 2026-09-23. Format: MADR 4.0.0 [web] https://adr.github.io/madr/.

## Context and Problem Statement

PLAN.md names BrAPI Genotyping as the second source of genotypes after files. A programme whose SNP
data sits in Gigwa, Germinate or Breedbase should not have to export a VCF to rank its backcross
progeny. The sibling project backcross settled the reading rules for a BrAPI v2.1 variant set in its
ADR 0015 and implements them in `src/io/brapi.ts`. The two tools are meant to read one server the
same way, so the question here is not what the rules should be. It is how they land in a Python
package whose transport, threading model and outputs differ from a browser worker's.

## Decision Drivers

One `Dataset` for every source, so the compute core, the screens and the exporters stay unaware of
where the genotypes came from. samples.csv still declares the recurrent and donor parents, and
markers.csv still overrides positions. No real genotype data in the repository and no network in any
test. An access token must not reach a URL, a message or a file. The shared input contract covers
files and says nothing about a server, so the loader has to sit outside it without bending it.

## Considered Options

1. Mirror backcross ADR 0015 in `io/brapi.py` with an injected transport, loading on the command
   line only.
2. Mirror the rules and offer the source on the Load screen as well.
3. Convert BrAPI to a temporary VCF and reuse `read_vcf`.

## Decision Outcome

Option 1. `src/progeny_selector/io/brapi.py` pages `/callsets`, `/variants` and `/allelematrix` (GT
only) for one variant set, joins each matrix cell to its variant and call set by the DbIds of that
response, and returns a `GenotypeMatrix` that `build_dataset` then joins with samples.csv exactly as
a file's would be. Every function takes the transport as an argument of type `FetchJson`; HTTP lives
in `urllib_fetch_json` alone. Option 3 would write genotypes to disk and lose the server ids, and it
would run every cell through the file vocabulary in `io/calls.py`, which is not what a GT token is.
Option 2 is deferred: Q6 was answered command line only for M3, because a Shinylive build cannot
reach a server that does not allow the page's origin, and Phase 8 of docs/m3-phases.md holds the
Load-screen work for the milestone that takes it up.

The seven rules below carry backcross ADR 0015's decision numbers.

1. **Positions** (D1). `pos_bp` is BrAPI's 0-based `start` plus 1, since `end` is exclusive and VCF
   POS is 1-based. `start` must be an `int` that is not a `bool`, or a float equal to its floor. A
   variant with no `referenceName` or no `start` takes both from markers.csv; with no entry there
   the load fails naming the marker. `end` is read for nothing. `referenceName` goes through
   `normalize_chrom` under the crop scheme the load resolved.
2. **Marker id** (D2) is the first `variantNames` entry that is neither empty nor `.`, and
   `variantDbId` otherwise. A bare string counts as one entry. Gigwa and Germinate put the marker
   name in `variantNames` and an opaque key in `variantDbId`, and markers.csv is keyed by the name.
3. **Calls** (D3) are allele indices read by `parse_brapi_call`, never by the nucleotide vocabulary
   in `io/calls.py`. Both separators split, so phase is ignored; a token equal to the response's
   `unknownString` is missing; a token missing on either side is missing as a whole; a token with no
   separator is haploid and reads as the homozygote. More than two parts, a part that is neither a
   non-negative integer nor the unknown string, and an index of 127 or more are errors, the last
   because the matrix is int8. Allele symbols are `[referenceBases, *alternateBases]` where the
   bases are given, and the index strings `"0"`, `"1"` and so on where they are not. Allele order
   within a call is kept and not sorted, which is what `read_vcf` does. The two loaders agree cell
   for cell on the fixture, and one case differs in the stored pair without differing in any
   result: `read_vcf` stores `./1` as `(-1, 1)` where `parse_brapi_call` stores `(-1, -1)`. Every
   consumer treats a negative allele on either side as missing, so no state, heterozygote count,
   IBS value or RPP figure changes.
4. **Sample id** (D4) is the trimmed `callSetName` when it is non-empty and unique among the variant
   set's call sets, and `callSetDbId` otherwise, with one warning per shared name. A collision among
   the resulting ids fails the load. Replicate call sets stay separate samples.
5. **Traceability** (D5). `Dataset.call_sets` carries a `CallSetRef` per sample and is empty for a
   file. `progeny-selector brapi-callsets` writes `sample_id,call_set_name,call_set_db_id,sample_db_id`
   before samples.csv exists, because samples.csv is built from the server's own ids. No results.csv
   or selected.csv column changes in M3, and the results schema is therefore untouched (docs/adr/0016).
6. **Paging** (D6). `/callsets` pages by number through the `totalPages` of page 0. `/variants` sends
   no token on page 0, then `pageToken=<nextPageToken>` when the previous response gave a non-empty
   one and `page=N&pageToken=N` when it did not, because test-server.brapi.org ignored `page` there on
   2026-09-15 and honoured `pageToken`. `/allelematrix` pages both dimensions, variant page outermost,
   with the dimension counts read from `result.pagination[]` of the first page. A repeated id across
   pages of one endpoint, a repeated id within a matrix page, a matrix page whose
   `(variantDbIds, callSetDbIds)` pair was already read, and an id in a matrix page that the lists
   did not carry are all errors. A variant or call set that no matrix page carried is a warning
   naming it, since its calls are all missing.
7. **Transport** (D7). `urllib_fetch_json` sends `Accept: application/json` and, when a token is
   configured, `Authorization: Bearer <token>`. The token is named on the command line only as the
   environment variable holding it (`--brapi-token-env`), so it never enters a shell history, a URL,
   an error, a warning or a written file. An `HTTPError` becomes `BrAPI: GET <url> returned HTTP
   <code>`, a `URLError` or timeout becomes `BrAPI: could not reach <url>: <reason>`, and a body that
   is not a JSON object is refused as such. Each request times out after 60 s.

**Where this differs from backcross ADR 0015.** The transport is `urllib.request` on CPython and not
`fetch` in a worker, so there is no cancel message and no CORS problem to report. The token comes
from an environment variable named on the command line and not from a password field. Server ids
live on `Dataset.call_sets` and are not added to every export, because results.csv and selected.csv
are frozen under `results_schema` and the sibling's two extra columns would break the R reader and
the round-trip tests. Loading is CPython only in M3. Two further divergences are worth naming. The sibling reads a bare `.` as a missing allele whatever
`unknownString` says (`brapi.ts:158`); this loader honours `unknownString` alone, so a server
declaring `unknownString: "NA"` and then sending `./.` is refused as an invalid GT instead of read
as missing. Refusing is the safer reading of a server contradicting itself, and no observed server
does it. The sibling also appends its call-set name collisions to the "absent from the genotype
file" error (`brapi.ts:135-142`); `load_brapi_dataset` does the same, catching that failure from
`build_dataset` and naming the collisions, so a user whose samples.csv lists a colliding name learns
the call set was loaded under its `callSetDbId`.

One rule is new here: after `/callsets` and
`/variants` and before the matrix is allocated, a variant set larger than `BRAPI_MAX_CALLS` adds a
warning. That constant is `6_000 * 2_002`, the exact call count of the largest CPython case in
docs/limits.md that completed: row `6k_x_2000`, 6,000 markers by 2,000 progeny plus both parents, a
22.9 MiB resident matrix. Writing 12,000,000 instead would make the measured case itself warn. The
load is still attempted, because no server-sourced case has been measured and the number bounds
nothing but the warning.

**Test data.** `scripts/make_fixture.py` writes `tests/fixtures/brapi/` from the BC2F1 dataset already
in memory: 25 variants on the target's chromosome and 8 call sets, paged 13 then 12 and 5 then 3, so
both list endpoints and all four allele-matrix pages are exercised. The directory stays under the
64 KB cap, which the generator asserts on every run rather than leaving the figure to a document. GT tokens are spelled differently per page (collapsed
homozygotes, one phased heterozygote, a bare `.` for a missing call on one page) while denoting the
same alleles, so `tests/test_brapi.py` tests the parser and not the fixture. That test compares the
BrAPI load against `read_vcf` of the same 25 marker ids and 8 sample ids and asserts the matrices are
equal. No test reaches the network: the transport is a stub over the recorded pages, and the one test
of `urllib_fetch_json` replaces `urllib.request.OpenerDirector.open`. It patches the opener and
not `urlopen`, because the transport builds its own opener to strip `Authorization` on a
cross-host redirect.

`tests/fixtures/brapi/criteria.yaml` is the BC2F1 criteria with the avoid locus dropped. That locus is
`syn_Gm13_10` and the variant set holds Gm06 alone, so `resolve_locus` refuses it.

**Questions Q6 to Q13 and their answers.** Q6, a browser source: no, command line only in M3. Q7,
traceability: `Dataset.call_sets` and `brapi-callsets.csv`, and no results column. Q8, the sample id:
`callSetName` when unique, else `callSetDbId`. Q9, the size guard: a warning above the call count of
the largest measured CPython case in docs/limits.md. Q11, the source fetched: plantbreeding/BrAPI on 2026-09-22 (docs/reference-repos.md).
Q13, a half-missing call: missing as a whole. VCF keeps the called side in the stored pair and this
loader does not, and no consumer can tell the difference, as D3 records.

### Consequences

Good:

- A BrAPI variant set becomes the same `Dataset` a file produces, so nothing downstream changes.
- Server ids survive the load and can be written out before samples.csv exists.
- The outputs and the shared contract are unchanged, so an existing reader keeps working.
- The tests need no server and no network.

Bad:

- The paging rule for `/variants` rests on one server's behaviour observed on one date. A server that
  ignores both `page` and `pageToken` fails with the repeated-page error and does not load.
- The repeated-page check is by page identity, so a server returning genuinely different ids for the
  same underlying page is not caught. Nothing cheaper would catch it, and content comparison cannot:
  a repeated page whose calls are all missing looks exactly like a first reading.
- Gigwa names call sets with numeric strings, so those sample ids are opaque until samples.csv maps
  them to line names.
- Server ids are on the `Dataset` and not in results.csv, so a row in an exported table leads back to
  the server only through `brapi-callsets.csv`.
- `apply_marker_map`'s override warning names "genotype file" even for a server-sourced position.

**Revisit when** the Load screen takes up the source (Phase 8), a server needs an authentication
scheme other than a bearer token, naming through `/samples` or `/germplasm` is wanted, or a variant
set is too large to page in one process.
