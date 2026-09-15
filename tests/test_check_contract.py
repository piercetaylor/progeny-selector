"""scripts/check_contract.py: manifest recomputation and the optional byte-compare with a sibling checkout."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "check_contract.py"
CONTRACT = ROOT / "contract"


def run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True, check=False)


def copy_contract(repo: Path) -> Path:
    shutil.copytree(CONTRACT, repo / "contract")
    return repo / "contract"


def test_committed_contract_passes() -> None:
    result = run()
    assert result.returncode == 0, result.stderr
    assert result.stdout.startswith("contract ") and "manifest ok" in result.stdout


def test_manifest_mismatch_exits_1(tmp_path: Path) -> None:
    ours = copy_contract(tmp_path / "ours")
    (ours / "VERSION").write_bytes((ours / "VERSION").read_bytes() + b"x")
    (ours / "cases" / "extra.txt").write_bytes(b"x")
    result = run("--contract", str(ours))
    assert result.returncode == 1
    assert "contract manifest: hash differs: VERSION" in result.stderr
    assert "contract manifest: not in manifest: cases/extra.txt" in result.stderr


def test_identical_sibling_exits_0(tmp_path: Path) -> None:
    ours = copy_contract(tmp_path / "ours")
    copy_contract(tmp_path / "sibling")
    result = run("--contract", str(ours), str(tmp_path / "sibling"))
    assert result.returncode == 0, result.stderr
    assert "contract mirror:" in result.stdout and "files identical" in result.stdout


def test_differing_sibling_exits_1(tmp_path: Path) -> None:
    ours = copy_contract(tmp_path / "ours")
    theirs = copy_contract(tmp_path / "sibling")
    (theirs / "cases" / "vcf-basic" / "samples.csv").write_bytes(b"sample_id,role\n")
    (theirs / "extra.txt").write_bytes(b"x")
    (ours / "cases" / "ours-only.txt").write_bytes(b"x")
    result = run("--contract", str(ours), str(tmp_path / "sibling"))
    assert result.returncode == 1
    assert "contract mirror: differs: cases/vcf-basic/samples.csv" in result.stderr
    assert "contract mirror: only in sibling: extra.txt" in result.stderr
    assert "contract mirror: missing in sibling: cases/ours-only.txt" in result.stderr


def test_missing_paths_exit_2(tmp_path: Path) -> None:
    assert run("--contract", str(tmp_path / "nowhere")).returncode == 2
    ours = copy_contract(tmp_path / "ours")
    assert run("--contract", str(ours), str(tmp_path / "no-sibling")).returncode == 2


def test_sibling_without_manifest_exits_2(tmp_path: Path) -> None:
    ours = copy_contract(tmp_path / "ours")
    theirs = copy_contract(tmp_path / "sibling")
    (theirs / "MANIFEST.sha256").unlink()
    result = run("--contract", str(ours), str(tmp_path / "sibling"))
    assert result.returncode == 2
    assert "MANIFEST.sha256 does not exist" in result.stderr and "Traceback" not in result.stderr
