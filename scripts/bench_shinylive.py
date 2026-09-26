"""Measure wall-clock and memory for a bench case running in Shinylive/Pyodide, in a real
browser (docs/adr/0021, docs/limits.md, Q4).

Usage:
    python scripts/bench_shinylive.py --site site --case data/bench/6k_x_2000 \
        [--port 8010] [--timeout-min 30]

Serves ``site/`` with Cross-Origin-Opener-Policy: same-origin and Cross-Origin-Embedder-Policy:
require-corp so ``performance.measureUserAgentSpecificMemory()`` (Q4) is available; asserts
``crossOriginIsolated`` in the page and exits 3 rather than printing a fabricated memory figure
when it is not available. Uploads the case's genotypes.vcf.gz, samples.csv, markers.csv and
criteria.yaml, clicks #load-run, and waits for #load-status to report the expected markers/progeny
count within ``--timeout-min`` minutes.
"""

from __future__ import annotations

import argparse
import http.server
import json
import re
import sys
import threading
import time
from pathlib import Path

from playwright.sync_api import FrameLocator, expect, sync_playwright


class _CoiHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self) -> None:
        self.send_header("Cross-Origin-Opener-Policy", "same-origin")
        self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
        super().end_headers()

    def log_message(self, format: str, *args: object) -> None:
        pass


def _serve(site: Path, host: str, port: int) -> http.server.ThreadingHTTPServer:
    handler = lambda *a, **kw: _CoiHandler(*a, directory=str(site), **kw)  # noqa: E731
    server = http.server.ThreadingHTTPServer((host, port), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    return server


def _wait_upload_complete(frame: FrameLocator, field_id: str) -> None:
    # Copied from tests/e2e/test_shinylive_export.py::_wait_upload_complete: scripts do not
    # import tests. Shiny animates each input_file's progress bar to width: 100% once the
    # upload is registered; clicking Run before that races the upload under Pyodide.
    progress = (
        frame.locator(f"#{field_id}")
        .locator("xpath=ancestor::div[contains(@class,'shiny-input-container')][1]//div[contains(@class,'progress-bar')]")
        .first
    )
    expect(progress).to_have_attribute("style", re.compile(r"width:\s*100%"), timeout=30_000)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site", type=Path, default=Path("site"))
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--port", type=int, default=8010)
    parser.add_argument("--timeout-min", type=float, default=30)
    args = parser.parse_args(argv)

    bench_json = args.case / "bench.json"
    with open(bench_json, encoding="utf-8") as fh:
        case_info = json.load(fh)
    n_markers, n_progeny = case_info["markers"], case_info["progeny"]

    genotypes = args.case / "genotypes.vcf.gz"
    if not genotypes.exists():
        genotypes = args.case / "genotypes.vcf"
    samples = args.case / "samples.csv"
    markers = args.case / "markers.csv"
    criteria = args.case / "criteria.yaml"
    gz_bytes = genotypes.stat().st_size

    host = "127.0.0.1"
    timeout_ms = int(args.timeout_min * 60_000)
    server = _serve(args.site, host, args.port)
    origin = f"http://{host}:{args.port}/"

    try:
        with sync_playwright() as pw:
            # channel="chromium" runs the full browser in its new headless mode. Plain
            # headless=True launches chromium-headless-shell (Playwright 1.49 and later),
            # whose measureUserAgentSpecificMemory() throws SecurityError even when the
            # page is crossOriginIsolated.
            browser = pw.chromium.launch(headless=True, channel="chromium")
            page = browser.new_page()
            try:
                t_goto = time.monotonic()
                page.goto(origin, timeout=timeout_ms)
                frame = page.frame_locator("iframe").first
                frame.locator("#load-run").wait_for(state="attached", timeout=timeout_ms)
                t_load_run = time.monotonic() - t_goto

                if page.evaluate("crossOriginIsolated") is not True:
                    print("not cross-origin isolated; memory cannot be measured (Q4)")
                    return 3

                mem_before = page.evaluate("performance.measureUserAgentSpecificMemory().then(r => r.bytes)")

                t_up = time.monotonic()
                frame.locator("#load-genotypes").set_input_files(str(genotypes))
                frame.locator("#load-samples").set_input_files(str(samples))
                frame.locator("#load-markers").set_input_files(str(markers))
                frame.locator("#load-criteria").set_input_files(str(criteria))
                for field_id in ("load-genotypes", "load-samples", "load-markers", "load-criteria"):
                    _wait_upload_complete(frame, field_id)
                t_uploads = time.monotonic() - t_up

                t_click = time.monotonic()
                frame.locator("#load-run").click()
                # Wait for either outcome, so a load error ends the run when it appears
                # and not when the timeout expires.
                status = frame.locator("#load-status")
                expected = f"{n_markers} markers, {n_progeny} progeny"
                expect(status).to_contain_text(re.compile(re.escape(expected) + "|error:"), timeout=timeout_ms)
                t_status = time.monotonic() - t_click
                if expected not in status.inner_text():
                    print(f"outcome: load failed after {t_status:.1f} s: {status.inner_text()}")
                    return 4

                mem_after = page.evaluate("performance.measureUserAgentSpecificMemory().then(r => r.bytes)")

                result = {
                    "case": args.case.name,
                    "gz_bytes": gz_bytes,
                    "t_load_run_s": t_load_run,
                    "t_uploads_s": t_uploads,
                    "t_status_s": t_status,
                    "mem_before_bytes": mem_before,
                    "mem_after_bytes": mem_after,
                    "mem_delta_bytes": mem_after - mem_before,
                    "outcome": "ok",
                }
                print(
                    f"| {args.case.name} | {gz_bytes / 1_000_000:.1f} | {t_load_run:.2f} | {t_uploads:.2f} "
                    f"| {t_status:.2f} | {mem_before / (1024 * 1024):.1f} / {mem_after / (1024 * 1024):.1f} "
                    f"/ {(mem_after - mem_before) / (1024 * 1024):.1f} | ok |"
                )
                with open(args.case / "bench_browser.json", "w", encoding="utf-8", newline="\n") as fh:
                    json.dump(result, fh, indent=2, sort_keys=True)
                    fh.write("\n")
                return 0
            except TimeoutError:
                print(f"outcome: did not complete within {args.timeout_min:.0f} min")
                return 4
            except Exception as exc:
                print(f"outcome: {exc}")
                return 4
            finally:
                browser.close()
    finally:
        server.shutdown()


if __name__ == "__main__":
    sys.exit(main())
