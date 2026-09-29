"""Shipped example inputs: the synthetic BC2F1 and BC3F1 fixtures as package data.

Responsibility: locate the example files inside the installed package (and inside the
Shinylive staged copy, which is the same directory tree) and write them to a directory.
The files are byte copies of tests/fixtures/ written by scripts/make_fixture.py, the only
author of both; tests/test_examples.py proves the copies identical. Synthetic data: the
rankings are test cases, not breeding recommendations.

Interface:
    EXAMPLES: dict[str, tuple[str, ...]]          # name -> files shipped under examples/<name>/
    SHARED_FROM_BC2F1: tuple[str, ...]            # ("markers.csv", "criteria.yaml"): bc3f1 reuses bc2f1's
    EXAMPLE_FILES: tuple[str, ...]                # ("genotypes.vcf", "samples.csv", "markers.csv", "criteria.yaml")
    example_files(name) -> dict[str, Traversable] # the four inputs of an example, shared ones resolved from bc2f1
    example_paths(name) -> AbstractContextManager[dict[str, Path]]   # real paths (importlib.resources.as_file)
    write_example(name, out_dir: Path, *, force: bool = False) -> list[Path]
                                                  # writes out_dir/<name>/<four files>; FileExistsError naming the
                                                  # first existing target unless force
"""

from __future__ import annotations

import errno
import os
from collections.abc import Iterator
from contextlib import AbstractContextManager, ExitStack, contextmanager
from importlib.resources import as_file, files
from importlib.resources.abc import Traversable
from pathlib import Path

EXAMPLE_FILES: tuple[str, ...] = ("genotypes.vcf", "samples.csv", "markers.csv", "criteria.yaml")
SHARED_FROM_BC2F1: tuple[str, ...] = ("markers.csv", "criteria.yaml")
EXAMPLES: dict[str, tuple[str, ...]] = {
    "synthetic_bc2f1": EXAMPLE_FILES,
    "synthetic_bc3f1": ("genotypes.vcf", "samples.csv"),
}
_SHARED_SOURCE = "synthetic_bc2f1"


def _check_name(name: str) -> None:
    if name not in EXAMPLES:
        raise KeyError(f"unknown example {name!r}; valid names: {', '.join(EXAMPLES)}")


def example_files(name: str) -> dict[str, Traversable]:
    """The four inputs of example ``name``, keyed by file name; shared ones resolved from bc2f1."""
    _check_name(name)
    root = files("progeny_selector.examples")
    shipped = EXAMPLES[name]
    out: dict[str, Traversable] = {}
    for file in EXAMPLE_FILES:
        src_name = name if file in shipped else _SHARED_SOURCE
        out[file] = root / src_name / file
    return out


@contextmanager
def _example_paths(name: str) -> Iterator[dict[str, Path]]:
    resources = example_files(name)
    with ExitStack() as stack:
        yield {file: stack.enter_context(as_file(res)) for file, res in resources.items()}


def example_paths(name: str) -> AbstractContextManager[dict[str, Path]]:
    """Real file-system paths to the four inputs of example ``name``, valid inside the ``with`` block."""
    _check_name(name)
    return _example_paths(name)


def write_example(name: str, out_dir: Path, *, force: bool = False) -> list[Path]:
    """Write the four inputs of example ``name`` to ``out_dir/<name>/``; return the written paths.

    Every target is checked before anything is written: without ``force``, the first existing
    target raises ``FileExistsError`` naming it and nothing is written.
    """
    resources = example_files(name)
    dest = Path(out_dir) / name
    targets = {file: dest / file for file in resources}
    if not force:
        for target in targets.values():
            if target.exists():
                raise FileExistsError(errno.EEXIST, os.strerror(errno.EEXIST), str(target))
    dest.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for file, res in resources.items():
        targets[file].write_bytes(res.read_bytes())
        written.append(targets[file])
    return written
