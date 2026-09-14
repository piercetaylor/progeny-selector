# Shinylive export stages the package; a smoke test asserts data stays in the tab

Status: accepted. Date: 2026-09-14.

## Context and Problem Statement

`shinylive export` bundles every file under the app directory it is pointed at, and finds Pyodide packages by scanning that directory's own imports and its `requirements.txt`. Pointed at `src/progeny_selector/app` as ADR 0003 originally described, it never sees `progeny_selector` itself or numpy: only `app/app.py`'s imports are scanned, and `app/requirements.txt` lists only pyyaml. The CI step "shinylive export builds" therefore checked only that the command exits 0, not that the exported site can run the pipeline. How is the package bundled into the export, and what proves the export works?

## Considered Options

1. Publish a wheel to a URL and add it to `app/requirements.txt`; `shinylive export` ignores `http(s)://` lines, so micropip installs it at run time in the browser.
2. Stage a copy of the package next to a stub `app.py` before exporting, so the scanner sees `import numpy` and the package's own files are bundled directly.
3. Keep pointing `shinylive export` at `src/progeny_selector/app` and accept that the export never runs the real pipeline.

## Decision Outcome

Option 2, with option 1 as the documented fallback (resolution 5, `docs/m1-phases.md`). `scripts/build_shinylive.py` deletes and recreates `build/shinylive-app/`, copies `src/progeny_selector` into it, writes a stub `app.py` (`from progeny_selector.app.app import app`) and an explicit `requirements.txt` (`pyyaml`, `pandas`, `numpy`), then runs `shinylive export` on the staging directory. Whether Pyodide's `sys.path` includes the staged directory was not established by reading the exporter's source, so `tests/e2e/test_shinylive_export.py` is the check: it serves the export over `http.server`, loads the fixture through the real UI, and asserts every network request the page makes starts with the page's own origin. A request to PyPI or a CDN would mean the staged copy did not import and the wheel-URL fallback is needed instead; on this machine (Windows, shinylive 0.8.11, Chromium matching the installed Playwright) the staged copy imports without it — see the M1 verification block in PLAN.md for the run that established this.

The smoke test also settles a mechanical detail the resolution did not anticipate: the exported page is a Shiny app rendered inside a same-origin `iframe` (`/app_<id>/`), not injected into the top-level document. Locators in the test are scoped to that frame with `page.frame_locator("iframe")` rather than `page`.

### Consequences

Good: no release, no PyPI dependency, and the export works from any checkout without a deploy having happened first; the smoke test is a real end-to-end check rather than an exit-code check, and it fails loudly (a non-origin request) if a future dependency change breaks the staged import. Bad: the staging directory duplicates `src/progeny_selector` under `build/` for every export (gitignored, but real disk and copy time); the wheel-URL fallback, if ever needed, ties the exported site to a GitHub Pages deploy that has already happened, which is why it stays a fallback and not the default. Neutral: `app/requirements.txt` is no longer what `shinylive export` reads when built through this script; it stays for `shiny run`, where it is unused, and for anyone exporting the unstaged app directory directly.
