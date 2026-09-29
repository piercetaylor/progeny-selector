"""``tool_commit``: the build file, Pyodide, and the guarded git lookup (docs/adr/0035)."""

from __future__ import annotations

import re
import shutil
import subprocess
from collections.abc import Iterator
from pathlib import Path

import pytest

from progeny_selector import provenance

needs_git = pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")


@pytest.fixture(autouse=True)
def fresh_cache() -> Iterator[None]:
    provenance.tool_commit.cache_clear()
    yield
    provenance.tool_commit.cache_clear()


def git(repo: Path, *args: str) -> str:
    cmd = ["git", "-c", "user.name=t", "-c", "user.email=t@example.org", "-c", "commit.gpgsign=false", *args]
    return subprocess.run(cmd, cwd=repo, capture_output=True, text=True, check=True).stdout


def make_repo(root: Path) -> Path:
    """A git repository with one commit holding ``src/progeny_selector/__init__.py`` and a README."""
    package = root / "src" / "progeny_selector"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("", encoding="utf-8")
    (root / "README.md").write_text("tracked\n", encoding="utf-8")
    git(root, "init", "-q")
    git(root, "add", ".")
    git(root, "commit", "-q", "-m", "init")
    return package


def head(repo: Path) -> str:
    return git(repo, "rev-parse", "HEAD").strip()


@needs_git
def test_clean_checkout_gives_g_and_seven_hex(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    package = make_repo(tmp_path)
    monkeypatch.setattr(provenance, "PACKAGE_DIR", package)
    commit = provenance.tool_commit()
    assert commit == "g" + head(tmp_path)[:7]
    assert re.fullmatch(r"g[0-9a-f]{7}", commit)


@needs_git
def test_a_modified_tracked_file_adds_dirty(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    package = make_repo(tmp_path)
    monkeypatch.setattr(provenance, "PACKAGE_DIR", package)
    (tmp_path / "untracked.txt").write_text("new\n", encoding="utf-8")
    assert provenance.tool_commit() == "g" + head(tmp_path)[:7]  # an untracked file alone is not dirty
    provenance.tool_commit.cache_clear()
    (tmp_path / "README.md").write_text("changed\n", encoding="utf-8")
    assert provenance.tool_commit() == "g" + head(tmp_path)[:7] + "-dirty"


@needs_git
def test_a_package_inside_another_repository_gives_none(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A wheel installed in a venv inside some other repository must not report that repository's commit."""
    make_repo(tmp_path)
    installed = tmp_path / ".venv" / "Lib" / "site-packages" / "progeny_selector"
    installed.mkdir(parents=True)
    monkeypatch.setattr(provenance, "PACKAGE_DIR", installed)
    assert provenance.git_commit(tmp_path / "src" / "progeny_selector") is not None  # the repo itself resolves
    assert provenance.tool_commit() is None


def test_outside_any_repository_gives_none(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GIT_CEILING_DIRECTORIES", str(tmp_path))
    monkeypatch.setattr(provenance, "PACKAGE_DIR", tmp_path / "progeny_selector")
    (tmp_path / "progeny_selector").mkdir()
    assert provenance.tool_commit() is None


def test_git_missing_gives_none(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def no_git(*args, **kwargs):
        raise FileNotFoundError("git")

    monkeypatch.setattr(provenance.subprocess, "run", no_git)
    monkeypatch.setattr(provenance, "PACKAGE_DIR", tmp_path)
    assert provenance.tool_commit() is None


def test_a_git_timeout_gives_none(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def slow_git(cmd, **kwargs):
        assert kwargs["timeout"] == 5
        raise subprocess.TimeoutExpired(cmd, kwargs["timeout"])

    monkeypatch.setattr(provenance.subprocess, "run", slow_git)
    monkeypatch.setattr(provenance, "PACKAGE_DIR", tmp_path)
    assert provenance.tool_commit() is None


@needs_git
def test_the_build_file_takes_precedence_over_git(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    package = make_repo(tmp_path)
    (package / provenance.BUILD_COMMIT_FILE).write_text("g0badc0d-dirty\n", encoding="utf-8")
    monkeypatch.setattr(provenance, "PACKAGE_DIR", package)
    assert provenance.tool_commit() == "g0badc0d-dirty"  # stripped, verbatim; git would give HEAD and -dirty


def test_the_build_file_is_read_under_pyodide(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / provenance.BUILD_COMMIT_FILE).write_text("g1234567", encoding="utf-8")
    monkeypatch.setattr(provenance, "PACKAGE_DIR", tmp_path)
    monkeypatch.setattr(provenance.sys, "platform", "emscripten")
    assert provenance.tool_commit() == "g1234567"


@needs_git
def test_pyodide_without_a_build_file_gives_none_without_calling_git(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    package = make_repo(tmp_path)  # git would find a commit here
    monkeypatch.setattr(provenance, "PACKAGE_DIR", package)
    monkeypatch.setattr(provenance.sys, "platform", "emscripten")

    def forbidden(*args, **kwargs):
        raise AssertionError("git must not run under Pyodide")

    monkeypatch.setattr(provenance.subprocess, "run", forbidden)
    assert provenance.tool_commit() is None


def test_the_result_is_cached_per_process(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    build_file = tmp_path / provenance.BUILD_COMMIT_FILE
    build_file.write_text("g1111111", encoding="utf-8")
    monkeypatch.setattr(provenance, "PACKAGE_DIR", tmp_path)
    assert provenance.tool_commit() == "g1111111"
    build_file.write_text("g2222222", encoding="utf-8")
    assert provenance.tool_commit() == "g1111111"


def test_tool_and_version() -> None:
    from progeny_selector import __version__

    assert provenance.TOOL == "progeny-selector"
    assert provenance.tool_version() == __version__
