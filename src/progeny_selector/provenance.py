"""Which build of the tool wrote a file: the ``tool``, ``tool_version`` and ``tool_commit`` columns.

Responsibility: name the tool, its version and the commit it was built from, for the provenance
columns of results.csv and selected.csv (docs/adr/0035). ``tool_commit`` is ``g`` plus the first
seven hex digits of HEAD, with ``-dirty`` when a tracked file differs from HEAD, and ``None``
(written ``NA``) when the commit cannot be known.

Lookup order, first hit wins, cached per process:

1. ``_build_commit.txt`` next to this module, stripped and used verbatim. Only
   ``scripts/build_shinylive.py`` writes it, into its staged copy of the package; it is gitignored
   in the source tree.
2. Under Pyodide (``sys.platform == "emscripten"``) there is no git: ``None``.
3. git, run in the package directory with a 5 s timeout, but only when the work tree's
   ``src/progeny_selector`` is this package. A wheel installed in a venv that sits inside some
   other repository must not report that repository's commit. Any OSError, timeout or non-zero
   exit gives ``None``.

Interface:
    TOOL -> "progeny-selector"
    tool_version() -> str
    tool_commit() -> str | None
    git_commit(package_dir) -> str | None
"""

from __future__ import annotations

import functools
import string
import subprocess
import sys
from pathlib import Path

__all__ = ["BUILD_COMMIT_FILE", "PACKAGE_DIR", "TOOL", "git_commit", "tool_commit", "tool_version"]

TOOL = "progeny-selector"
BUILD_COMMIT_FILE = "_build_commit.txt"
# Read at call time, so a test can point the lookup at another directory.
PACKAGE_DIR = Path(__file__).resolve().parent
GIT_TIMEOUT_S = 5


def tool_version() -> str:
    # Imported here: this module is imported while ``progeny_selector/__init__`` is still loading.
    from progeny_selector import __version__

    return __version__


def _git(package_dir: Path, *args: str) -> str | None:
    try:
        done = subprocess.run(
            ["git", *args],
            cwd=package_dir,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=GIT_TIMEOUT_S,
            check=False,
        )
    except (OSError, subprocess.SubprocessError, ValueError):
        return None
    return done.stdout if done.returncode == 0 else None


def git_commit(package_dir: Path) -> str | None:
    """``g<7 hex>`` of HEAD, ``-dirty`` when tracked files differ, for the repository whose package this is."""
    top = _git(package_dir, "rev-parse", "--show-toplevel")
    if top is None or not top.strip():
        return None
    try:
        if (Path(top.strip()) / "src" / "progeny_selector").resolve() != package_dir.resolve():
            return None
    except OSError:
        return None
    head = _git(package_dir, "rev-parse", "HEAD")
    if head is None:
        return None
    sha = head.strip().lower()
    if len(sha) < 7 or any(c not in string.hexdigits for c in sha):
        return None
    status = _git(package_dir, "status", "--porcelain", "--untracked-files=no")
    if status is None:
        return None
    return f"g{sha[:7]}" + ("-dirty" if status.strip() else "")


@functools.cache
def tool_commit() -> str | None:
    """The commit this build came from, or ``None`` when it cannot be known (written ``NA``)."""
    build_file = PACKAGE_DIR / BUILD_COMMIT_FILE
    try:
        if build_file.is_file():
            return build_file.read_text(encoding="utf-8").strip()
    except OSError:
        pass
    if sys.platform == "emscripten":
        return None
    return git_commit(PACKAGE_DIR)
