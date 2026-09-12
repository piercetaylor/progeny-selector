"""Cross / family / generation tree for the Navigate screen, and the row filter it drives.

Responsibility: derive the navigation tree from the sample manifest (progeny and
candidates only; parents are not nodes), filter result rows by the selected
family and generation, and label the breadcrumb. Pure functions; no Shiny.

Interface:
    GenerationNode, FamilyNode, NavTree (frozen dataclasses)
    build_tree(dataset: Dataset) -> NavTree
    UNASSIGNED = ""   selector meaning "only rows with no value" (None or empty string)
    filter_rows(rows, family, generation) -> list[dict]   (None means no filter)
    crumb_labels(tree, family, generation) -> list[str]
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from progeny_selector.model.dataset import Dataset

# Selector for the unassigned node: rows whose family_id (or generation) is None or "".
UNASSIGNED = ""


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


def _matches(value: object, selector: str | None) -> bool:
    if selector is None:
        return True
    if selector == UNASSIGNED:
        return value is None or value == ""
    return value == selector


def filter_rows(rows: list[dict], family: str | None, generation: str | None) -> list[dict]:
    """Rows whose family_id and generation match; None does not filter, UNASSIGNED keeps rows with no value. Order is preserved."""
    return [r for r in rows if _matches(r.get("family_id"), family) and _matches(r.get("generation"), generation)]


def crumb_labels(tree: NavTree, family: str | None, generation: str | None) -> list[str]:
    """Breadcrumb labels: the cross, then the family and generation when set; UNASSIGNED reads "(no family)" / "(no generation)"."""
    labels = [tree.cross]
    if family is not None:
        labels.append("(no family)" if family == UNASSIGNED else family)
    if generation is not None:
        labels.append("(no generation)" if generation == UNASSIGNED else generation)
    return labels
