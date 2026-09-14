"""Stage a copy of the package next to a stub app.py, then run `shinylive export`.

`shinylive export` finds Pyodide packages by scanning the app directory's own
imports and its requirements.txt (docs/m1-phases.md resolution 5); it never sees
`progeny_selector` or numpy unless the package sits next to the app it exports.
This script copies `src/progeny_selector` into a staging directory, writes a stub
`app.py` that imports the real app object, and an explicit `requirements.txt` so
the scanner's handling of inline imports does not matter.

Usage: python scripts/build_shinylive.py [--staging build/shinylive-app] [--out site]
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

STUB_APP_PY = "from progeny_selector.app.app import app  # noqa: F401\n"
STUB_REQUIREMENTS = "pyyaml\npandas\nnumpy\n"
EXPECTED_PIPELINE_ENTRY = '"name": "progeny_selector/core/pipeline.py"'


def _dir_size_bytes(path: Path) -> int:
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def _run_shinylive_export(staging: Path, out: Path) -> None:
    """Run `shinylive export staging out`.

    `python -m shinylive` is the form used elsewhere in this repo's tooling, but
    shinylive 0.8.11 ships no `__main__.py`: `python -m shinylive` fails with
    "No module named shinylive.__main__". The documented form is the console
    script (https://github.com/posit-dev/py-shinylive), installed next to the
    interpreter that installed it, so that is the fallback used here.
    """
    module_cmd = [sys.executable, "-m", "shinylive", "export", str(staging), str(out)]
    result = subprocess.run(module_cmd, capture_output=True, text=True)
    if result.returncode == 0:
        sys.stdout.write(result.stdout)
        return
    if "No module named shinylive.__main__" not in result.stderr:
        # A real export failure: surface it rather than masking it with the fallback.
        sys.stdout.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise SystemExit(result.returncode)

    console_script = Path(sys.executable).with_name("shinylive" + (".exe" if os.name == "nt" else ""))
    if not console_script.exists():
        raise SystemExit(f"`python -m shinylive` is unavailable and the console script is missing at {console_script}")
    subprocess.run([str(console_script), "export", str(staging), str(out)], check=True)


def stage_app(staging: Path) -> None:
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    (staging / "app.py").write_text(STUB_APP_PY, newline="\n")
    (staging / "requirements.txt").write_text(STUB_REQUIREMENTS, newline="\n")
    shutil.copytree(
        REPO_ROOT / "src" / "progeny_selector",
        staging / "progeny_selector",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staging", default=str(REPO_ROOT / "build" / "shinylive-app"))
    parser.add_argument("--out", default=str(REPO_ROOT / "site"))
    args = parser.parse_args()
    staging = Path(args.staging).resolve()
    out = Path(args.out).resolve()

    stage_app(staging)

    start = time.monotonic()
    _run_shinylive_export(staging, out)
    elapsed = time.monotonic() - start

    app_json = out / "app.json"
    if not app_json.exists():
        raise SystemExit(f"export did not produce {app_json}")
    contents = app_json.read_text(encoding="utf-8")
    if EXPECTED_PIPELINE_ENTRY not in contents:
        raise SystemExit(f"export did not stage progeny_selector: expected {EXPECTED_PIPELINE_ENTRY} in {app_json}")

    size_mb = _dir_size_bytes(out) / 1_000_000
    print(f"exported {out} ({size_mb:.1f} MB) in {elapsed:.1f}s")


if __name__ == "__main__":
    main()
