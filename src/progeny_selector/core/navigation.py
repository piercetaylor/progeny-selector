"""Cross / family / generation tree for the Navigate screen, and the row filter it drives.

Responsibility: derive the navigation tree from the sample manifest (progeny and
candidates only; parents are not nodes), filter result rows by the selected
family and generation, and label the breadcrumb. Pure functions; no Shiny.

Interface:
    GenerationNode, FamilyNode, NavTree (frozen dataclasses)
    build_tree(dataset: Dataset) -> NavTree
    filter_rows(rows, family, generation) -> list[dict]   (None means no filter)
    crumb_labels(tree, family, generation) -> list[str]
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from progeny_selector.model.dataset import Dataset


@dataclass(frozen=True)
class GenerationNode:
    generation: str | None
    n: int


@dataclass(frozen=True)
class FamilyNode:
    family_id: str | None
    n: int
    generations: tuple[GenerationNode, ...]


@dataclass(frozen=True)
class NavTree:
    cross: str
    n: int
    families: tuple[FamilyNode, ...]


def _none_last(value: str | None) -> tuple[bool, str]:
    return (value is None, value or "")


def build_tree(dataset: Dataset) -> NavTree:
    """One cross node, its families, and each family's generations; None sorts last at each level."""
    progeny = dataset.progeny
    cross = f"{dataset.recurrent_parent.line_name} x {dataset.donor_parent.line_name}"
    fam_counts = Counter(s.family_id for s in progeny)
    gen_counts = Counter((s.family_id, s.generation) for s in progeny)
    families: list[FamilyNode] = []
    for fam in sorted(fam_counts, key=_none_last):
        gens = sorted((g for f, g in gen_counts if f == fam), key=_none_last)
        nodes = tuple(GenerationNode(g, gen_counts[(fam, g)]) for g in gens)
        families.append(FamilyNode(fam, fam_counts[fam], nodes))
    return NavTree(cross=cross, n=len(progeny), families=tuple(families))


def filter_rows(rows: list[dict], family: str | None, generation: str | None) -> list[dict]:
    """Rows whose family_id and generation match; a None selector does not filter. Order is preserved."""
    return [
        r for r in rows if (family is None or r.get("family_id") == family) and (generation is None or r.get("generation") == generation)
    ]


def crumb_labels(tree: NavTree, family: str | None, generation: str | None) -> list[str]:
    """Breadcrumb labels: the cross, then the selected family and generation when set."""
    labels = [tree.cross]
    if family is not None:
        labels.append(family)
    if generation is not None:
        labels.append(generation)
    return labels
