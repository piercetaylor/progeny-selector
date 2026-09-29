"""Every fenced sh block command in the README and the tutorial runs, in document order."""

import re
import shlex
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
PREFIX = "progeny-selector "
FENCE = re.compile(r"^```sh\s*\n(.*?)^```", re.DOTALL | re.MULTILINE)


def _commands(text: str) -> list[str]:
    lines = [ln.strip() for block in FENCE.findall(text) for ln in block.splitlines()]
    return [
        ln[len(PREFIX) :] for ln in lines if ln.startswith(PREFIX) and not ln.startswith("progeny-selector app") and "--brapi-url" not in ln
    ]


@pytest.mark.parametrize("doc", ["README.md", "docs/tutorial.md"])
def test_fenced_sh_commands_run(doc: str, tmp_path: Path) -> None:
    commands = _commands((ROOT / doc).read_text(encoding="utf-8"))
    assert commands, f"no runnable progeny-selector commands in {doc}"
    for rest in commands:
        done = subprocess.run(
            [sys.executable, "-m", "progeny_selector", *shlex.split(rest)],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )
        assert done.returncode == 0, f"{doc}: progeny-selector {rest}\n{done.stderr}"
