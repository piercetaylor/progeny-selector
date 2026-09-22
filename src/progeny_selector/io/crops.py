"""Crop chromosome schemes (contract/data-contract.md 1.5.0, "Chromosome names"; contract/crops/;
docs/adr/0020; mirrors backcross src/io/crops.ts).

Responsibility: the nine schemes shipped under contract/crops/, the schema check for one of
them, and the id-to-compiled-scheme lookup the readers, ``load_dataset``, the CLI and the
Load screen resolve a chosen crop through. ``soybean`` is the default and reproduces the
1.2.0 rule exactly. No user-supplied scheme is accepted in this version. Shinylive stages
the package, not the repository root, so the literals below are the runtime copy of
contract/crops/*.json; tests/test_crops.py checks them against the files.

``CropScheme``, ``CompiledScheme`` and ``compile_scheme`` are defined in ``core/chrom.py``,
which holds the soybean literal because ``core`` imports no ``io``; they are re-exported here.

Interface:
    CropScheme, CompiledScheme, compile_scheme, DEFAULT_CROP_ID
    BUILTIN_CROPS: dict[str, CropScheme]   (soybean, maize, rice, sorghum, wheat, barley, oat,
                                            common-bean, cotton, in that order)
    validate_scheme(obj) -> CropScheme     (DataContractError "crop scheme: <reason>")
    resolve_crop(crop: str | None) -> CompiledScheme   (None -> soybean; an unknown id is an error)
"""

from __future__ import annotations

import re

from progeny_selector.core.chrom import SOYBEAN, SOYBEAN_SCHEME, CompiledScheme, CropScheme, compile_scheme
from progeny_selector.model.dataset import DataContractError

DEFAULT_CROP_ID = "soybean"

_KEYS = frozenset({"id", "name", "species", "ploidy", "assembly", "chromosomes", "keys", "pattern", "sources"})
_ID_PATTERN = re.compile(r"^[a-z0-9][a-z0-9-]{0,31}$")
_LEADING_ZERO = re.compile(r"(?<![0-9])0(?=[0-9])")

__all__ = [
    "BUILTIN_CROPS",
    "DEFAULT_CROP_ID",
    "CompiledScheme",
    "CropScheme",
    "compile_scheme",
    "resolve_crop",
    "validate_scheme",
]


def _fail(reason: str) -> DataContractError:
    return DataContractError(f"crop scheme: {reason}")


def _string_list(value: object, field: str) -> list[str]:
    if not isinstance(value, list) or not all(isinstance(v, str) for v in value):
        raise _fail(f'"{field}" must be an array of strings')
    return [str(v) for v in value]


def validate_scheme(obj: object) -> CropScheme:
    """The schema rules of contract 1.5.0; raises DataContractError naming the first rule broken."""
    if not isinstance(obj, dict):
        raise _fail("must be a JSON object")
    unknown = [k for k in obj if k not in _KEYS]
    if unknown:
        raise _fail("unknown key(s) " + ", ".join(f'"{k}"' for k in unknown))
    cid = obj.get("id")
    if not isinstance(cid, str) or not _ID_PATTERN.fullmatch(cid):
        raise _fail(f'"id" must match {_ID_PATTERN.pattern}')
    name = obj.get("name")
    if not isinstance(name, str) or not name.strip():
        raise _fail('"name" must be a non-empty string')
    species = obj.get("species")
    if not isinstance(species, str) or not species.strip():
        raise _fail('"species" must be a non-empty string')
    ploidy = obj.get("ploidy")
    if not isinstance(ploidy, int) or isinstance(ploidy, bool) or ploidy < 1:
        raise _fail('"ploidy" must be a positive integer')
    assembly = obj.get("assembly")
    if not isinstance(assembly, str) or not assembly.strip():
        raise _fail('"assembly" must be a non-empty string')
    chroms = _string_list(obj.get("chromosomes"), "chromosomes")
    keys = _string_list(obj.get("keys"), "keys")
    if not chroms:
        raise _fail('"chromosomes" must be non-empty')
    if len(chroms) != len(keys):
        raise _fail(f'"chromosomes" has {len(chroms)} entries and "keys" {len(keys)}')
    for i, c in enumerate(chroms):
        if c in chroms[:i]:
            raise _fail(f'chromosome "{c}" appears twice')
    for i, k in enumerate(keys):
        if k in keys[:i]:
            raise _fail(f'key "{k}" appears twice')
    if any(k != k.upper() or _LEADING_ZERO.search(k) for k in keys):
        raise _fail("every key must be upper-case with no leading zeros")
    pattern = obj.get("pattern")
    if not isinstance(pattern, str):
        raise _fail('"pattern" must be a string')
    try:
        re.compile(pattern)
    except re.error as exc:
        raise _fail(f'"pattern" does not compile: {exc}') from exc
    return CropScheme(
        id=cid,
        name=name,
        species=species,
        ploidy=ploidy,
        assembly=assembly,
        chromosomes=chroms,
        keys=keys,
        pattern=pattern,
        sources=_string_list(obj.get("sources"), "sources"),
    )


BUILTIN_CROPS: dict[str, CropScheme] = {
    "soybean": SOYBEAN_SCHEME,
    "maize": CropScheme(
        id="maize",
        name="Maize",
        species="Zea mays",
        ploidy=2,
        assembly="Zm-B73-REFERENCE-NAM-5.0",
        chromosomes=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10"],
        keys=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"],
        pattern=r"^(?:chr|chromosome)?[_\s-]?0*([1-9]|10)$",
        sources=[
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_assembly_report.txt",
            "https://rest.ensembl.org/info/assembly/zea_mays",
        ],
    ),
    "rice": CropScheme(
        id="rice",
        name="Rice",
        species="Oryza sativa",
        ploidy=2,
        assembly="IRGSP-1.0 / MSU7",
        chromosomes=["Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10", "Chr11", "Chr12"],
        keys=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"],
        pattern=r"^(?:chr|chromosome)?[_\s-]?0*([1-9]|1[0-2])$",
        sources=[
            "https://rice.uga.edu/annotation_pseudo_current.shtml",
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_assembly_report.txt",
            "https://rest.ensembl.org/info/assembly/oryza_sativa",
        ],
    ),
    "sorghum": CropScheme(
        id="sorghum",
        name="Sorghum",
        species="Sorghum bicolor",
        ploidy=2,
        assembly="BTx623 v3.1.1 (NCBIv3)",
        chromosomes=["Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10"],
        keys=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"],
        pattern=r"^(?:chr|chromosome)?[_\s-]?0*([1-9]|10)$",
        sources=[
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/195/GCF_000003195.3_Sorghum_bicolor_NCBIv3/GCF_000003195.3_Sorghum_bicolor_NCBIv3_assembly_report.txt",
        ],
    ),
    "wheat": CropScheme(
        id="wheat",
        name="Wheat",
        species="Triticum aestivum",
        ploidy=6,
        assembly="IWGSC CS RefSeq v2.1",
        chromosomes=[
            "Chr1A",
            "Chr1B",
            "Chr1D",
            "Chr2A",
            "Chr2B",
            "Chr2D",
            "Chr3A",
            "Chr3B",
            "Chr3D",
            "Chr4A",
            "Chr4B",
            "Chr4D",
            "Chr5A",
            "Chr5B",
            "Chr5D",
            "Chr6A",
            "Chr6B",
            "Chr6D",
            "Chr7A",
            "Chr7B",
            "Chr7D",
        ],
        keys=["1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"],
        pattern=r"^(?:chr|chromosome)?[_\s-]?0*([1-7])[_-]?([ABD])$",
        sources=[
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/294/505/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_assembly_report.txt",
        ],
    ),
    "barley": CropScheme(
        id="barley",
        name="Barley",
        species="Hordeum vulgare",
        ploidy=2,
        assembly="MorexV3",
        chromosomes=["chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H"],
        keys=["1H", "2H", "3H", "4H", "5H", "6H", "7H"],
        pattern=r"^(?:chr|chromosome)?[_\s-]?0*([1-7])[_-]?(H)$",
        sources=[
            "https://rest.ensembl.org/info/assembly/hordeum_vulgare",
        ],
    ),
    "oat": CropScheme(
        id="oat",
        name="Oat",
        species="Avena sativa",
        ploidy=6,
        assembly="OT3098 v2",
        chromosomes=[
            "chr1A",
            "chr1C",
            "chr1D",
            "chr2A",
            "chr2C",
            "chr2D",
            "chr3A",
            "chr3C",
            "chr3D",
            "chr4A",
            "chr4C",
            "chr4D",
            "chr5A",
            "chr5C",
            "chr5D",
            "chr6A",
            "chr6C",
            "chr6D",
            "chr7A",
            "chr7C",
            "chr7D",
        ],
        keys=["1A", "1C", "1D", "2A", "2C", "2D", "3A", "3C", "3D", "4A", "4C", "4D", "5A", "5C", "5D", "6A", "6C", "6D", "7A", "7C", "7D"],
        pattern=r"^(?:chr|chromosome)?[_\s-]?0*([1-7])[_-]?([ACD])$",
        sources=[
            "https://rest.ensembl.org/info/assembly/avena_sativa_ot3098",
        ],
    ),
    "common-bean": CropScheme(
        id="common-bean",
        name="Common bean",
        species="Phaseolus vulgaris",
        ploidy=2,
        assembly="G19833 v2.1",
        chromosomes=["Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11"],
        keys=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"],
        pattern=r"^(?:chr|pv|chromosome)?[_\s-]?0*([1-9]|1[01])$",
        sources=[
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/499/845/GCF_000499845.2_P._vulgaris_v2.0/GCF_000499845.2_P._vulgaris_v2.0_assembly_report.txt",
            "https://rest.ensembl.org/info/assembly/phaseolus_vulgaris",
            "https://knowpulse.usask.ca/bio_data/2691095",
        ],
    ),
    "cotton": CropScheme(
        id="cotton",
        name="Upland cotton",
        species="Gossypium hirsutum",
        ploidy=4,
        assembly="TM-1 UTX v2.1",
        chromosomes=[
            "A01",
            "A02",
            "A03",
            "A04",
            "A05",
            "A06",
            "A07",
            "A08",
            "A09",
            "A10",
            "A11",
            "A12",
            "A13",
            "D01",
            "D02",
            "D03",
            "D04",
            "D05",
            "D06",
            "D07",
            "D08",
            "D09",
            "D10",
            "D11",
            "D12",
            "D13",
        ],
        keys=[
            "A1",
            "A2",
            "A3",
            "A4",
            "A5",
            "A6",
            "A7",
            "A8",
            "A9",
            "A10",
            "A11",
            "A12",
            "A13",
            "D1",
            "D2",
            "D3",
            "D4",
            "D5",
            "D6",
            "D7",
            "D8",
            "D9",
            "D10",
            "D11",
            "D12",
            "D13",
        ],
        pattern=r"^(?:chr|chromosome)?[_\s-]?([AD])[_-]?0*([1-9]|1[0-3])$",
        sources=[
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/990/345/GCF_007990345.1_Gossypium_hirsutum_v2.1/GCF_007990345.1_Gossypium_hirsutum_v2.1_assembly_report.txt",
        ],
    ),
}

_COMPILED: dict[str, CompiledScheme] = {
    "soybean": SOYBEAN,
    **{cid: compile_scheme(s) for cid, s in BUILTIN_CROPS.items() if cid != "soybean"},
}


def resolve_crop(crop: str | None) -> CompiledScheme:
    """The compiled scheme for a crop id; ``None`` is the default (soybean)."""
    compiled = _COMPILED.get(DEFAULT_CROP_ID if crop is None else crop)
    if compiled is None:
        raise DataContractError(f'unknown crop "{crop}"; built-in crops are {", ".join(BUILTIN_CROPS)}')
    return compiled
