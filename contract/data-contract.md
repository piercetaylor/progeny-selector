# Data contract

Contract version: 1.2.1

The input files shared by backcross and progeny-selector. This directory is the canonical copy in backcross and is mirrored byte for byte into progeny-selector; README.md gives the version rules and the cases under cases/ are the machine-checked examples.

## Chromosome names

Chromosome names in this version are soybean-only: the 20 Glycine max chromosomes are normalised to `Gm01`..`Gm20` and every other name is kept as written; schemes for other crops are a later version. Accepted spellings for the 20 Glycine max chromosomes: an optional prefix `Gm`, `Chr`, `Chromosome` or `LG` (case-insensitive), an optional single separator `_`, space or `-`, then the number 1..20 with any number of leading zeros; so `Gm01`, `gm1`, `Chr07`, `chr7`, `Chromosome_07`, `LG7`, `Gm-7`, `7` and `07` all normalise to `Gm07`. As a pattern on the trimmed cell: `^(?:gm|chr|chromosome|lg)?[_\s-]?0*([1-9]|1[0-9]|20)$`, case-insensitive. All are normalised to `Gm01`..`Gm20` for display, ordering and export. Any other name (scaffolds, unplaced contigs, `ch7`) is kept unchanged and ordered after Gm20 in natural order: the name is split into digit runs and text runs, digit runs compare numerically and text runs lexically, so `scaffold_2` precedes `scaffold_10`. Positions are 1-based base pairs on whichever Williams 82 assembly the user's files use; the tool does not convert between assemblies. SoyBase lists the Wm82.a1.v1.1, Wm82.a2.v1 and Wm82.a4.v1 assemblies with Gm01–Gm20 naming, and explains the naming pattern: in `Wm82.a4.v1` the middle field is the assembly version and the last field the annotation version [web] https://www.soybase.org/resources/genome_info/. A newer near-gapless assembly, Wm82.a6, was built from the single-plant sub-line Wm82-ISU-01 and released through Phytozome v13 [web] https://www.biorxiv.org/content/10.1101/2024.04.26.591401v1 (published as Espina et al. 2024, Plant Journal 120:1221–1235, DOI 10.1111/tpj.17026, as listed in PubMed search results; the publisher page returned 403 when fetched). The SoyBase genome-information page fetched on 2026-09-04 did not list a6, so which assembly is "current" for a given program's marker positions is an open question recorded in PLAN.md.

## Genotype file

One of three formats, selected by file extension: `.vcf`; `.hmp.txt`, `.hmp` or `.hapmap`; `.csv`, `.tsv` or `.txt`; each optionally followed by `.gz` or `.bgz` for gzip or bgzip compression, which is inflated first, including multi-member bgzip streams. A genotype file must carry one of these extensions; every case under cases/ does. A repository may additionally detect the format from content (`##fileformat=VCF`, a leading `rs#` header, otherwise wide CSV) and may parse a compressed VCF as a stream instead of inflating it to text, but neither is required by this contract. Every genotype is a diploid call, two alleles per sample per marker or missing; polyploid dosage is out of scope. Every text input (genotype file, samples.csv, markers.csv) may use LF or CRLF line ends, and a UTF-8 byte-order mark before the first line is ignored. Every position column (VCF `POS`, HapMap `pos`, wide-CSV and markers.csv `pos_bp`) is a non-negative integer. VCF `POS` must be written as decimal digits, as the VCF specification declares it Integer. The other position columns also accept a whole-valued number written with a fraction part or an exponent (`1000.0`, `1e3`, `1.0E3`, `1.9E+07`), read as that integer, which is how spreadsheets, pandas and R export integer columns. As a pattern on the trimmed cell: `^\d+$` for VCF `POS`, `^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$` elsewhere, with the value finite, equal to its floor and not negative. Anything else is an error naming the line and the value: an empty cell, a non-zero fraction (`100.7`, `1e-3`), a negative number, `NaN`, `inf`, hexadecimal, or digits with `_` or `,` separators. Position 0 is accepted (VCF uses 0 and N+1 for telomeres); a negative position is an error in every column. `cm` is not a position and keeps its fraction.

### VCF 4.2 or later

Fixed columns `#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT` followed by one column per sample, as in the VCF 4.2 specification [web] https://samtools.github.io/hts-specs/VCFv4.2.pdf. Only CHROM, POS, ID, REF, ALT, FORMAT and the GT sub-field are read. GT allele indices refer to the REF,ALT list (0 = REF); `.` is missing; `/` and `|` are treated alike (phase ignored); haploid GT is read as homozygous; multiallelic ALT is supported. Records with ID `.` or empty get the id `<CHROM>_<POS>` from CHROM exactly as written in the file, before chromosome normalisation, and POS as the parsed integer, so a record `chr13 019000000 .` is `chr13_19000000`, the string bcftools writes for `%CHROM_%POS`. INFO, QUAL and FILTER are ignored; filter upstream with bcftools.

### HapMap (TASSEL style)

Eleven fixed columns `rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode`, then one column per taxon [web] https://statgen-esalq.github.io/Hapmap-and-VCF-formats-and-its-integration-with-onemap/; TASSEL's own description states that only `chrom` and `pos` must be populated, genotypes are two characters (`AA`) or single nucleotide codes, and `N` is missing [web] https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load. Cells are two nucleotides (`AA`, `AT`), a slash or bar pair (`A/T`, `A|T`), one nucleotide (`A`, homozygous) or one IUPAC heterozygote code (R, Y, S, W, K, M), which is read as its two nucleotides (R = A/G, Y = C/T, S = C/G, W = A/T, K = G/T, M = A/C). Missing cells are `N`, `NN`, `NA`, `-`, `--`, `.`, `./.`, `.|.`, `X`, `XX` and empty (TASSEL reads `-` as a deletion and `X` as unknown; both are missing here). Any other cell is an error naming the cell: `?`, `B`, `H`, `0`, `+`, any other single character, and a two-character or slash cell containing a character other than A, C, G, T, N, `-` or `.` (`A?`, `N?`). A pair of one nucleotide and one of N, `-`, `.` (`AN`, `A-`) is not defined by this version. Cells are trimmed and compared case-insensitively. Tab-delimited.

### Wide CSV

| column            | type    | rule                                                               |
| ----------------- | ------- | ------------------------------------------------------------------ |
| marker_id         | text    | unique                                                             |
| chrom             | text    | any accepted chromosome spelling                                   |
| pos_bp            | integer | 1-based position                                                   |
| `<sample_id>` ... | text    | one column per sample; header is the sample_id used in samples.csv |

Comma or tab delimited (sniffed from the header line: more tabs than commas means tab), RFC 4180 quoting. The first three columns are `marker_id`, `chrom` and `pos_bp` in that order, header names case-insensitive, and every column after them is a sample; a file that places the three elsewhere is outside the contract. backcross finds the three by name in any position and treats every other column as a sample, which accepts more than the contract requires for nucleotide files; its A/B/H detection still assumes the three come first, so a coded file with them elsewhere may be misread. A row whose every cell is empty or whitespace (`,,,`) is skipped; in any other row an empty `marker_id` or an empty `pos_bp` is an error naming the line. Two cell vocabularies:

| mode       | homozygous                                        | heterozygous                                     | missing                                                       |
| ---------- | ------------------------------------------------- | ------------------------------------------------ | ------------------------------------------------------------- |
| nucleotide | `A`, `AA`                                         | `A/T`, `A\|T`, `AT`, one IUPAC code R Y S W K M  | empty, `N`, `NN`, `NA`, `-`, `--`, `.`, `./.`, `.\|.`         |
| coded      | `A` (recurrent-parent allele), `B` (donor allele) | `H`                                              | empty, `N`, `NA`                                              |

Cells are trimmed and compared case-insensitively. Mode `auto` scans every cell of the file, not a leading window: it selects coded when every cell outside the nucleotide missing set is in {A, B, H} and at least one B or H occurs; otherwise nucleotide. In nucleotide mode an IUPAC code is read as its two nucleotides, as in HapMap, and any other cell is an error naming the cell: `?`, a single `B` or `H` (a nucleotide file cannot mix in coded letters), `X`, `XX`, `0`, `+`, any other single character, and a two-character or slash cell containing a character other than A, C, G, T, N, `-` or `.` (`A?`, `N?`); `X` and `XX` are missing in HapMap only. A pair of one nucleotide and one of N, `-`, `.` (`AN`, `A-`) is not defined by this version. In coded mode any cell outside {A, B, H} and the coded missing set is an error, so `-`, `--`, `.`, `./.`, `.\|.` and `NN` are errors there. `?` is not a missing token in either mode. In coded mode allele 0 is A and allele 1 is B at every marker, so the parents may be absent from the file; the manifest still names them.

Example (nucleotide):

```
marker_id,chrom,pos_bp,RP_Williams,DONOR_PI,NIL_01
syn_Gm13_10,Gm13,19000000,C,T,T
syn_Gm13_11,13,21000000,G,A,G/A
syn_Gm13_12,chr13,23000000,A,A,A
```

## samples.csv (sample manifest)

| column     | required | values                                                                |
| ---------- | -------- | --------------------------------------------------------------------- |
| sample_id  | yes      | must match a genotype column (except parents of a coded file); unique |
| line_name  | no       | display name; defaults to sample_id                                   |
| role       | yes      | `recurrent_parent`, `donor_parent`, `candidate`, `progeny`            |
| generation | no       | e.g. `BC5F3`; display and expected-value lookup only                  |
| family_id  | no       | grouping label; defaults to empty                                     |
| notes      | no       | free text                                                             |

Exactly one `recurrent_parent` and exactly one `donor_parent` per file; at least one `candidate` or `progeny`. Column names are case-insensitive; comma or tab delimited. Genotype columns absent from the manifest are dropped with a warning, and the loaded samples are in manifest order, not genotype-file order; manifest samples absent from the genotype file are an error, except the two parents of a coded file, which need no column. A loader may create parent columns internally for a coded file (all A, all B), but such parents are not part of its reported sample list: in a case's `expected.json`, `sampleIds` and `calls` hold only samples with a column in the genotype file. Pyramided lines with two donors are a later extension: the contract will add `donor_parent_2` and a per-marker donor assignment rather than change the meaning of existing columns.

Example:

```
sample_id,line_name,role,generation,family_id,notes
RP_Williams,Williams 82,recurrent_parent,,,
DONOR_PI,PI 000000,donor_parent,,,
NIL_01,NIL-01,candidate,BC5F3,FAM1,one donor segment Gm13
```

## markers.csv (marker map, optional)

| column    | required | values                                                          |
| --------- | -------- | --------------------------------------------------------------- |
| marker_id | yes      | matches the genotype file                                       |
| chrom     | yes      | any accepted spelling                                           |
| pos_bp    | yes      | integer                                                         |
| cm        | no       | genetic position; enables cM-weighted RPP and cM segment bounds |

Map positions override genotype-file positions when they differ (the count of overrides is reported as a warning); markers absent from the map keep their genotype-file position and have no cM.

## Class codes in exports

Where a class is written it uses the labels `missing, rp_hom, donor_hom, het, uninformative, nonparental` (numeric codes 0–5 in src/core/types.ts; the numbers are stable and are never renumbered).
