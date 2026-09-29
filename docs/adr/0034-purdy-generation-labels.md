# Single-cross Purdy generation labels set the backcross count only

Status: accepted. Date: 2026-09-27. Decided by research the maintainer delegated to Fable ("ask fable"), question Q-C of 2026-09-27. backcross is unchanged: there the generation is display only.

## Context and Problem Statement

`core/generation.py` read `generation` values in the BCnFm grammar (`BC2F1`, `BC3F2`, `F2`, `BC1`, `BC2S1`) and flagged anything else `generation_unparsed`. Many programmes, CIMMYT-derived wheat programmes among them, write a backcross line in Purdy notation instead: `RP*3/DONOR`, the recurrent parent with a dose. Those lines lost their expected RPP and were flagged as unparsed. The contract's `generation` column is "display and expected-value lookup only", so the accepted grammar is this tool's vocabulary, not the contract's, and no contract change is needed.

## Considered Options

1. Leave Purdy labels unparsed.
2. Parse them into a full BCnF1 generation.
3. Parse them into the backcross count only, with the filial generation marked unknown.

## Decision Outcome

Option 3. After the BCnFm grammar fails, `parse_generation` tries two patterns, a single cross with one `/` and one `*` and the dose on either side of `*`:

- `^\s*(?:(?P<rp>[^/*\s]+)\*(?P<n1>\d+)|(?P<n2>\d+)\*(?P<rp2>[^/*\s]+))/(?P<d>[^/*\s]+)\s*$`
- `^\s*(?P<d>[^/*\s]+)/(?:(?P<rp>[^/*\s]+)\*(?P<n1>\d+)|(?P<n2>\d+)\*(?P<rp2>[^/*\s]+))\s*$`

The dose counts every use of the recurrent parent including the initial cross, so a dose n >= 2 is BC(n-1) and n = 1 is the F1 (`A*2/B` = BC1; https://bio.libretexts.org/Bookshelves/Agriculture_and_Horticulture/Crop_Improvement_(Suza_and_Lamkey)/01:_Chapters/1.02:_Pedigree_Naming_Systems_and_Symbols; https://cropforge.github.io/iciswiki/articles/w/h/e/Wheat_Pedigree_ecd4.html). UPOV puts the dose before the asterisk; Purdy et al. 1968 itself was not verified. A plain `A/B`, a dose of 0 and anything else stay `generation_unparsed`.

Purdy carries no filial generation, so the parsed `Generation` has a new field `filial_known = False` (True for BCnFm labels). It sets `expected_rpp` only: `expected_het` is `NA`, and `het_rate_deviates` and `possible_self_or_outcross` are not evaluated for it, in `core/qc.py`. Asserting F1 would be silently wrong: a selfed Purdy line has heterozygosity near 0.25 against a BC1F1's 0.5 and would fire `het_rate_deviates`. `possible_outcross`, which rests on non-parental alleles rather than on the generation, is still evaluated. `Generation.label()` returns `BC{n}` with no filial part when `filial_known` is False, and a dose of 1 labels as `F1`, the cross itself.

### Consequences

Good: Purdy-labelled lines get an expected RPP and stop carrying `generation_unparsed`; no QC flag fires on a filial generation the label never stated. Bad: those lines get no heterozygosity or self/outcross check at all, even when the line is in fact a BCnF1; a breeder who wants the checks writes the BCnFm label. Only single crosses are read; a three-way or a double cross stays unparsed. `tests/test_generation_purdy.py` pins the grammar, the dose rule and the skipped checks.
