#!/usr/bin/env python3
"""Verify the mirrored data contract: MANIFEST.sha256 and, optionally, the canonical copy.

Responsibility: recompute the manifest of contract/ exactly as backcross's
scripts/make-contract.mjs writes it (every file except MANIFEST.sha256, posix paths
relative to contract/, sorted, "<sha256 hex>  <path>" per line, LF) and compare it with
the committed MANIFEST.sha256; when a sibling repository root is given, also byte-compare
every file with <sibling>/contract, as backcross's scripts/check-contract-mirror.mjs
does from the other side. Never writes: this repository holds the mirror, and the
canonical copy is regenerated in backcross.

Usage: python scripts/check_contract.py [--contract DIR] [SIBLING_REPO]
Exit codes: 0 manifest matches (and the mirror is identical, when checked);
1 manifest mismatch or mirror difference, one line per difference on stderr;
2 usage error or a path that does not exist.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

REPO_CONTRACT = Path(__file__).resolve().parents[1] / "contract"
MANIFEST = "MANIFEST.sha256"


def relative_files(root: Path) -> list[str]:
    """Every file under root except the manifest, as sorted posix paths relative to root."""
    return sorted(p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file() and p.name != MANIFEST)


def manifest_text(root: Path) -> str:
    return "".join(f"{hashlib.sha256((root / p).read_bytes()).hexdigest()}  {p}\n" for p in relative_files(root))


def _entries(text: str) -> dict[str, str]:
    """path -> digest for each manifest line."""
    out: dict[str, str] = {}
    for line in text.splitlines():
        digest, sep, path = line.partition("  ")
        if sep:
            out[path] = digest
    return out


def check_manifest(root: Path) -> list[str]:
    """Differences between the committed manifest and a recomputation; empty when identical."""
    committed = (root / MANIFEST).read_bytes().decode("utf-8")
    recomputed = manifest_text(root)
    if committed == recomputed:
        return []
    have, want = _entries(committed), _entries(recomputed)
    problems: list[str] = []
    for path in sorted(set(have) | set(want)):
        if path not in have:
            problems.append(f"not in manifest: {path}")
        elif path not in want:
            problems.append(f"listed but absent: {path}")
        elif have[path] != want[path]:
            problems.append(f"hash differs: {path}")
    if not problems:
        problems.append("MANIFEST.sha256 differs from the recomputation in order or line ends")
    return problems


def check_mirror(ours: Path, theirs: Path) -> list[str]:
    """Files missing on either side or differing in bytes; empty when the copies are identical."""
    a, b = {*relative_files(ours), MANIFEST}, {*relative_files(theirs), MANIFEST}
    differences: list[str] = []
    for path in sorted(a | b):
        if path not in b:
            differences.append(f"missing in sibling: {path}")
        elif path not in a:
            differences.append(f"only in sibling: {path}")
        elif (ours / path).read_bytes() != (theirs / path).read_bytes():
            differences.append(f"differs: {path}")
    return differences


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Verify contract/MANIFEST.sha256 and, optionally, the mirror against a sibling checkout.")
    parser.add_argument("sibling", nargs="?", help="root of the backcross checkout; its contract/ is byte-compared with ours")
    parser.add_argument("--contract", type=Path, default=REPO_CONTRACT, help="contract directory to check (default: this repository's)")
    args = parser.parse_args(argv)
    root = args.contract.resolve()
    if not (root / MANIFEST).is_file():
        print(f"check_contract: {root / MANIFEST} does not exist", file=sys.stderr)
        return 2
    theirs: Path | None = None
    if args.sibling is not None:
        theirs = Path(args.sibling).resolve() / "contract"
        if not theirs.is_dir():
            print(f"check_contract: {theirs} does not exist", file=sys.stderr)
            return 2
        if not (theirs / MANIFEST).is_file():
            print(f"check_contract: {theirs / MANIFEST} does not exist", file=sys.stderr)
            return 2
    problems = check_manifest(root)
    for line in problems:
        print(f"contract manifest: {line}", file=sys.stderr)
    if not problems:
        version_file = root / "VERSION"
        version = version_file.read_text(encoding="utf-8").strip() if version_file.is_file() else "?"
        print(f"contract {version}: manifest ok, {len(relative_files(root))} files")
    differences: list[str] = []
    if theirs is not None:
        differences = check_mirror(root, theirs)
        for line in differences:
            print(f"contract mirror: {line}", file=sys.stderr)
        if not differences:
            print(f"contract mirror: {len(relative_files(root)) + 1} files identical")
    return 1 if problems or differences else 0


if __name__ == "__main__":
    sys.exit(main())
