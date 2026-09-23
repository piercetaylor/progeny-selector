# BrAPI fixture (generated)

`scripts/make_fixture.py` writes these files; do not edit them by hand. They are recorded
BrAPI v2.1 pages for one variant set, `vs1`, holding exactly the calls the BC2F1
fixture's `genotypes.vcf` holds for 25 variants on Gm06 and 8 call sets
(both parents and the first 3 progeny of each family).

Paging is deliberately small so every path is exercised with no network. `/callsets` returns
8 call sets over two pages of 5, `/variants` returns 25
variants over two pages of 13, and `/allelematrix` therefore has four pages,
`allelematrix.v0.c0.json` through `allelematrix.v1.c1.json`, read variant page outermost.

GT tokens are spelled differently per page while denoting the same alleles, which is what makes
`tests/test_brapi.py` a test of the parser and not of the fixture. Page 0,0 collapses homozygotes
to a single index (`0`, `1`) and phases one heterozygote (`0|1`); the other pages expand
homozygotes (`0/0`, `1/1`). A missing call is `.` on page 1,1 and `./.` elsewhere.

`variants-nopos.p0.json` repeats the same 25 variants in one page with `referenceName`,
`start`, `end` and `referenceBases` null and `alternateBases` empty. It is the D1 fallback case:
positions then come from markers.csv, and a load without markers.csv fails naming the first marker.

`samples.csv` declares the roles of the 8 sample ids the server's `callSetName` values
produce. `criteria.yaml` is the BC2F1 criteria with the avoid locus dropped, because that locus is
on Gm13 and this variant set holds Gm06 alone.
