"""Chromosome name normalisation and ordering under a crop scheme.

Responsibility: map the chromosome spellings a crop scheme accepts onto that crop's
canonical names, keep every other name unchanged, and provide a sort key that puts the
scheme's canonical names first, in scheme order, and everything else after them in
natural (numeric-aware) order. A scheme is a JSON document (contract/crops/<id>.json,
contract/data-contract.md 1.5.0, docs/adr/0020): an ordered list of canonical chromosome
names, a parallel list of keys, and one alias regular expression whose capture groups
yield the key. Key derivation: join the defined capture groups, upper-case, strip the
leading zeros of each digit run; look the key up in ``keys``. A miss keeps the name as
written and orders it after the canonical names, which is the 1.2.0 rule.

``core`` imports no ``io``, so ``soybean.json`` is repeated here as ``SOYBEAN_SCHEME``
(tests/test_crops.py checks the two against each other) and every function defaults to
it; ``io/crops.py`` holds the other eleven. Mirrors backcross src/core/chromosomes.ts.

Interface:
    CropScheme, CompiledScheme (``.id`` is the scheme's id), compile_scheme(scheme) -> CompiledScheme
    SOYBEAN_SCHEME: CropScheme, SOYBEAN: CompiledScheme
    normalize_chrom(name: str, scheme: CompiledScheme = SOYBEAN) -> str
    chrom_sort_key(name: str, scheme: CompiledScheme = SOYBEAN) -> tuple[int, tuple[tuple[int, str], ...]]
    chrom_length_bp(name: str, fallback: int | None, assembly: str, scheme: CompiledScheme = SOYBEAN) -> int | None
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from progeny_selector.constants import DEFAULT_ASSEMBLY, SOYBEAN_CHROM_LENGTHS_BP

_DIGIT_RUN = re.compile(r"([0-9]+)")
# Leading zeros of a digit run: the `0`s of `01`, never the `0`s inside `100`.
_LEADING_ZEROS = re.compile(r"(?<![0-9])0+(?=[0-9])")


@dataclass(frozen=True)
class CropScheme:
    """The fields of contract/crops/<id>.json."""

    id: str
    name: str
    species: str
    ploidy: int
    assembly: str
    chromosomes: list[str] = field(default_factory=list)
    keys: list[str] = field(default_factory=list)
    pattern: str = ""
    sources: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class CompiledScheme:
    scheme: CropScheme
    regex: re.Pattern[str]
    index_by_key: dict[str, int]
    index_by_canonical: dict[str, int]

    @property
    def id(self) -> str:
        """The scheme's id, so callers write ``scheme.id`` and not ``scheme.scheme.id``."""
        return self.scheme.id


def compile_scheme(scheme: CropScheme) -> CompiledScheme:
    """The lookup form of a scheme: its pattern compiled case-insensitively and its two indexes."""
    return CompiledScheme(
        scheme=scheme,
        regex=re.compile(scheme.pattern, re.IGNORECASE),
        index_by_key={k: i for i, k in enumerate(scheme.keys)},
        index_by_canonical={c: i for i, c in enumerate(scheme.chromosomes)},
    )


#: The literal of contract/crops/soybean.json (contract 1.5.0).
SOYBEAN_SCHEME = CropScheme(
    id="soybean",
    name="Soybean",
    species="Glycine max",
    ploidy=2,
    assembly="Williams 82 (Wm82.a2.v1 / a4.v1 / a6.v1 naming)",
    chromosomes=[f"Gm{i:02d}" for i in range(1, 21)],
    keys=[str(i) for i in range(1, 21)],
    pattern=r"^(?:gm|chr|chromosome|lg)?[_\s-]?0*([1-9]|1[0-9]|20)$",
    sources=[
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/GCF_000004515.6_Glycine_max_v4.0_assembly_report.txt",
        "https://rest.ensembl.org/info/assembly/glycine_max",
    ],
)

SOYBEAN = compile_scheme(SOYBEAN_SCHEME)


def _key_of(match: re.Match[str]) -> str:
    """The scheme key a match yields: the defined groups joined, upper-cased, unpadded."""
    joined = "".join(g for g in match.groups() if g is not None).upper()
    return _LEADING_ZEROS.sub("", joined)


def normalize_chrom(name: str, scheme: CompiledScheme = SOYBEAN) -> str:
    """The scheme's canonical name for any accepted spelling; any other name stripped but unchanged."""
    raw = str(name).strip()
    # ``search``, not ``match``: backcross's ``regex.exec`` searches, and every shipped pattern
    # anchors itself with ``^...$``, so the two agree on an unanchored pattern as well.
    match = scheme.regex.search(raw)
    if match is None:
        return raw
    index = scheme.index_by_key.get(_key_of(match))
    return raw if index is None else scheme.scheme.chromosomes[index]


def _natural_key(name: str) -> tuple[tuple[int, str], ...]:
    """Split on digit runs: even parts are text, odd parts are numbers, so aligned parts share a type."""
    parts = _DIGIT_RUN.split(name)
    return tuple((int(p), "") if i % 2 else (0, p) for i, p in enumerate(parts))


def chrom_sort_key(name: str, scheme: CompiledScheme = SOYBEAN) -> tuple[int, tuple[tuple[int, str], ...]]:
    """The scheme's canonical names first in scheme order, then any other name in natural order (scaffold_2 before scaffold_10)."""
    canonical = normalize_chrom(name, scheme)
    index = scheme.index_by_canonical.get(canonical)
    if index is not None:
        return (index, ())
    return (len(scheme.scheme.chromosomes), _natural_key(canonical))


def chrom_length_bp(
    name: str, fallback: int | None = None, assembly: str = DEFAULT_ASSEMBLY, scheme: CompiledScheme = SOYBEAN
) -> int | None:
    """Chromosome length in the named assembly, or ``fallback`` when the chromosome or the assembly is unknown.

    ``assembly`` is the criteria.yaml key (the contract leaves chromosome lengths to the tool); "none"
    and any name without a table return ``fallback``. The tables are Williams 82's, so every scheme but
    the built-in ``SOYBEAN`` returns ``fallback`` and its callers fall back to the last marker
    position on the chromosome (contract 1.5.0, docs/adr/0020): a maize ``chr1`` must never take
    ``Gm01``'s length, and neither must a look-alike scheme that only calls itself "soybean".
    """
    if scheme is not SOYBEAN:
        return fallback
    return SOYBEAN_CHROM_LENGTHS_BP.get(assembly, {}).get(normalize_chrom(name, scheme), fallback)
