"""``progeny-selector app``: the missing-dependency message and the arguments handed to ``shiny.run_app``."""

from __future__ import annotations

import sys

import pytest

from progeny_selector.cli import APP_MISSING_MESSAGE, main


def test_app_without_shiny_prints_the_install_line(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    monkeypatch.setitem(sys.modules, "shiny", None)
    assert main(["app"]) == 2
    err = capsys.readouterr().err
    assert "is not installed" in err
    assert APP_MISSING_MESSAGE.splitlines()[1] in err


def test_app_passes_host_port_and_browser_flag(monkeypatch: pytest.MonkeyPatch) -> None:
    import shiny

    calls: list[tuple[tuple, dict]] = []

    def recorder(*args, **kwargs) -> None:
        calls.append((args, kwargs))

    monkeypatch.setattr(shiny, "run_app", recorder)
    assert main(["app", "--host", "0.0.0.0", "--port", "8123", "--no-browser"]) == 0
    assert calls == [(("progeny_selector.app.app:app",), {"host": "0.0.0.0", "port": 8123, "launch_browser": False})]
    calls.clear()
    assert main(["app"]) == 0
    assert calls == [(("progeny_selector.app.app:app",), {"host": "127.0.0.1", "port": 8000, "launch_browser": True})]
