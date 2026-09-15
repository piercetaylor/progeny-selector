"""Text opening and delimiter sniffing shared by every reader.

Responsibility: open plain or gzip/bgzip text as UTF-8 with a leading
byte-order mark removed (contract 1.1.0 accepts CRLF and a BOM on every
text input), and choose the delimiter of a CSV/TSV from its header line
by the contract's rule: more tabs than commas means tab, otherwise comma.

Interface:
    open_text(path, newline=None) -> TextIO
    read_text(path) -> str            (newline="" so the csv module sees CRLF itself)
    sniff_delimiter(text) -> str      ("\\t" or ",")
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import TextIO


def open_text(path: str | Path, newline: str | None = None) -> TextIO:
    if str(path).endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8-sig", newline=newline)
    return open(path, encoding="utf-8-sig", newline=newline)


def read_text(path: str | Path) -> str:
    with open_text(path, newline="") as fh:
        return fh.read()


def sniff_delimiter(text: str) -> str:
    first = text.split("\n", 1)[0]
    return "\t" if first.count("\t") > first.count(",") else ","
