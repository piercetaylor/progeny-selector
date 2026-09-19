"""Named token profiles (contract/data-contract.md 1.4.0, "Token profiles"; mirrors backcross src/io/profiles.ts).

Responsibility: the five built-in profiles, a schema check for a user-supplied
profile, and the compiled lookup form the HapMap and wide-CSV readers resolve
cells through (calls.py). A profile is keyed by platform, never by crop. Tokens
are trimmed and compared case-insensitively, so compiled keys are upper-cased.
Shinylive stages the package, not the repository root, so the literals below
are the runtime copy of contract/profiles/*.json; tests/test_profiles.py checks
them against the files.

Interface:
    TokenProfile (dataclass, the JSON fields), CompiledProfile, DEFAULT_PROFILE_ID
    BUILTIN_PROFILES: dict[str, TokenProfile]   (tassel, soybase-report, dart, axiom, kasp)
    validate_profile(obj) -> TokenProfile        (DataContractError "token profile: <reason>")
    resolve_profile(ref: str | dict | None) -> TokenProfile | None   (None or "default" -> None)
    profile_label(p) -> "default" | id | "custom:<id>"   (id only for the BUILTIN_PROFILES object itself, so a
                                                        user file is custom:<id> even with a built-in's id)
    compile_profile(p) -> CompiledProfile
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from progeny_selector.model.dataset import DataContractError

DEFAULT_PROFILE_ID = "default"

_KEYS = frozenset({"id", "name", "description", "base", "homozygous", "heterozygous", "missing", "sources"})
_ID_PATTERN = re.compile(r"^[a-z0-9][a-z0-9-]{0,31}$")
_SYMBOL_PATTERN = re.compile(r"^[^\s]+$")


@dataclass(frozen=True)
class TokenProfile:
    id: str
    name: str
    description: str
    base: str  # "nucleotide" | "none"
    homozygous: dict[str, str] = field(default_factory=dict)
    heterozygous: dict[str, list[str] | str] = field(default_factory=dict)  # token -> [symbol, symbol] or "*"
    missing: list[str] = field(default_factory=list)
    sources: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class CompiledProfile:
    missing: frozenset[str]
    homozygous: dict[str, str]
    heterozygous: dict[str, tuple[str, str] | str]
    base: str
    label: str


def _norm(token: str) -> str:
    return token.strip().upper()


def _fail(reason: str) -> DataContractError:
    return DataContractError(f"token profile: {reason}")


def validate_profile(obj: object) -> TokenProfile:
    if not isinstance(obj, dict):
        raise _fail("must be a JSON object")
    unknown = [k for k in obj if k not in _KEYS]
    if unknown:
        raise _fail("unknown key(s) " + ", ".join(f'"{k}"' for k in unknown))
    pid = obj.get("id")
    if not isinstance(pid, str) or not _ID_PATTERN.match(pid):
        raise _fail(f'"id" must match {_ID_PATTERN.pattern}')
    name = obj.get("name")
    if not isinstance(name, str) or not name.strip():
        raise _fail('"name" must be a non-empty string')
    description = obj.get("description")
    if not isinstance(description, str):
        raise _fail('"description" must be a string')
    base = obj.get("base")
    if base not in ("nucleotide", "none"):
        raise _fail('"base" must be "nucleotide" or "none"')

    homozygous = obj.get("homozygous")
    if not isinstance(homozygous, dict):
        raise _fail('"homozygous" must be an object of token to allele symbol')
    hom: dict[str, str] = {}
    for token, symbol in homozygous.items():
        if symbol == "*":
            raise _fail(f'"*" is allowed only in "heterozygous" (homozygous token "{token}")')
        if not isinstance(symbol, str) or not _SYMBOL_PATTERN.match(symbol):
            raise _fail(f'homozygous token "{token}" must map to a non-blank allele symbol')
        hom[str(token)] = symbol

    heterozygous = obj.get("heterozygous")
    if not isinstance(heterozygous, dict):
        raise _fail('"heterozygous" must be an object of token to two allele symbols or "*"')
    het: dict[str, list[str] | str] = {}
    for token, pair in heterozygous.items():
        if pair == "*":
            het[str(token)] = "*"
            continue
        if not isinstance(pair, list) or len(pair) != 2 or not all(isinstance(s, str) and _SYMBOL_PATTERN.match(s) for s in pair):
            raise _fail(f'heterozygous token "{token}" must map to two allele symbols or "*"')
        if pair[0] == pair[1]:
            raise _fail(f'heterozygous token "{token}" maps to two identical symbols "{pair[0]}"')
        het[str(token)] = [pair[0], pair[1]]

    missing = obj.get("missing")
    if not isinstance(missing, list) or not all(isinstance(t, str) for t in missing):
        raise _fail('"missing" must be an array of strings')
    sources = obj.get("sources")
    if not isinstance(sources, list) or not all(isinstance(t, str) for t in sources):
        raise _fail('"sources" must be an array of strings')

    seen: dict[str, str] = {}
    for tokens, set_name in ((list(hom), "homozygous"), (list(het), "heterozygous"), (list(missing), "missing")):
        for token in tokens:
            key = _norm(token)
            prior = seen.get(key)
            if prior == set_name:
                raise _fail(f'token "{token}" appears twice in "{set_name}"')
            if prior is not None:
                raise _fail(f'token "{token}" appears in both "{prior}" and "{set_name}"')
            seen[key] = set_name

    if base == "none":
        if not hom:
            raise _fail('with base "none", "homozygous" must be non-empty')
        if "" not in missing:
            raise _fail('with base "none", "missing" must contain ""')

    return TokenProfile(
        id=pid,
        name=name,
        description=description,
        base=str(base),
        homozygous=hom,
        heterozygous=het,
        missing=list(missing),
        sources=list(sources),
    )


BUILTIN_PROFILES: dict[str, TokenProfile] = {
    "tassel": TokenProfile(
        id="tassel",
        name="TASSEL HapMap vocabulary in a wide CSV",
        description=(
            "The HapMap vocabulary applied to a wide CSV converted from TASSEL: X and XX are missing, as in HapMap. Nothing else changes."
        ),
        base="nucleotide",
        homozygous={},
        heterozygous={},
        missing=["X", "XX"],
        sources=[
            "https://bitbucket.org/tasseladmin/tassel-5-source/raw/master/src/net/maizegenetics/dna/snp/NucleotideAlignmentConstants.java"
        ],
    ),
    "soybase-report": TokenProfile(
        id="soybase-report",
        name="SoyBase SNP allele report",
        description=(
            "Nucleotides for homozygotes, H for a heterozygote (resolved to the marker's two alleles on the row), U for a missing call."
        ),
        base="nucleotide",
        homozygous={},
        heterozygous={"H": "*"},
        missing=["U"],
        sources=["maintainer's note, 2026-09-14 (docs/adr/0014); SoyBase pages returned 403 on 2026-09-14"],
    ),
    "dart": TokenProfile(
        id="dart",
        name="DArTseq SNP report, one row per marker",
        description=(
            "0 is homozygous reference, 1 homozygous alternate, 2 heterozygous, - missing. Allele symbols are 0 and 1. "
            "Do not use for the two-row report or for dartR's recoding, where 1 is the heterozygote."
        ),
        base="none",
        homozygous={"0": "0", "1": "1"},
        heterozygous={"2": ["0", "1"]},
        missing=["", "-", "NA"],
        sources=["docs/research/cross-crop-genotype-conventions.md [3], search only"],
    ),
    "axiom": TokenProfile(
        id="axiom",
        name="Axiom Analysis Suite or GenomeStudio AA/AB/BB export",
        description=(
            "AA, AB, BB with NoCall or -- missing; also the numeric form 0, 1, 2 with -1 (NoCall) and -2 (OTV) missing. "
            "Allele symbols are A and B, the array's alleles, not the parents'."
        ),
        base="none",
        homozygous={"AA": "A", "BB": "B", "0": "A", "2": "B"},
        heterozygous={"AB": ["A", "B"], "BA": ["A", "B"], "1": ["A", "B"]},
        missing=["", "NoCall", "OTV", "-1", "-2", "--", "NA"],
        sources=["docs/research/cross-crop-genotype-conventions.md [5], [6], search only"],
    ),
    "kasp": TokenProfile(
        id="kasp",
        name="KASP SNPviewer export",
        description=(
            "X:X and Y:Y homozygous, X:Y heterozygous, ?, Uncallable and Missing missing; NTC, Dupe and Bad are missing too. "
            "Allele symbols are X and Y, the assay's alleles."
        ),
        base="none",
        homozygous={"X:X": "X", "Y:Y": "Y"},
        heterozygous={"X:Y": ["X", "Y"], "Y:X": ["X", "Y"]},
        missing=["", "?", "Uncallable", "Missing", "NTC", "Dupe", "Bad", "NA"],
        sources=[
            "https://www.biosearchtech.com/products/pcr-reagents-kits-and-instruments/pcr-instruments-and-software/genotyping-and-lims-software/snpviewer"
        ],
    ),
}


def resolve_profile(ref: str | dict | TokenProfile | None) -> TokenProfile | None:
    if ref is None or ref == DEFAULT_PROFILE_ID:
        return None
    if isinstance(ref, TokenProfile):
        return ref
    if isinstance(ref, str):
        if ref not in BUILTIN_PROFILES:
            raise DataContractError(f'unknown token profile "{ref}"; built-in profiles are {", ".join(BUILTIN_PROFILES)}')
        return BUILTIN_PROFILES[ref]
    return validate_profile(ref)


def profile_label(p: TokenProfile | None) -> str:
    if p is None:
        return DEFAULT_PROFILE_ID
    return p.id if BUILTIN_PROFILES.get(p.id) is p else f"custom:{p.id}"


def compile_profile(p: TokenProfile) -> CompiledProfile:
    het: dict[str, tuple[str, str] | str] = {}
    for token, pair in p.heterozygous.items():
        het[_norm(token)] = pair if isinstance(pair, str) else (pair[0], pair[1])
    return CompiledProfile(
        missing=frozenset(_norm(t) for t in p.missing),
        homozygous={_norm(t): s for t, s in p.homozygous.items()},
        heterozygous=het,
        base=p.base,
        label=profile_label(p),
    )
