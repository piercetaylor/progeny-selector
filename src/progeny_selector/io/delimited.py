"""Text opening, blank-line rules, delimiter sniffing and CSV rows shared by every reader.

Responsibility: open plain or gzip/bgzip text as UTF-8 with a leading
byte-order mark removed (contract 1.1.0 accepts CRLF and a BOM on every
text input), and choose the delimiter of a CSV/TSV from its header line
by the contract's rule: more tabs than commas means tab, otherwise comma.
Contract 1.3.0: a blank line or cell is empty or holds only spaces and tabs
(no other character, so a no-break space is not blank); the delimiter is
sniffed from the first physical line that is not blank, and text that is all
blank sniffs a comma. `csv_rows` reads a delimited file the way the contract
states it: a `"` opens a quoted field only at the start of a field and is
otherwise literal (the csv module's own reading), a line break inside a quoted
field is read as a single LF, a quoted field still open at the end of the
file raises DataContractError naming the physical line where it opened, and
rows whose every cell is blank are skipped.

Interface:
    open_text(path, newline=None) -> TextIO
    read_text(path) -> str            (newline="" so the csv module sees CRLF itself)
    is_blank(text) -> bool            (only spaces and tabs, after a line end is removed)
    sniff_delimiter(text) -> str      ("\\t" or ",")
    csv_rows(text, label) -> list[(line_num, row)]   (line_num: the row's last physical line)
"""

from __future__ import annotations

import csv
import gzip
import io
from pathlib import Path
from typing import TextIO

from progeny_selector.model.dataset import DataContractError


def open_text(path: str | Path, newline: str | None = None) -> TextIO:
    if str(path).endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8-sig", newline=newline)
    return open(path, encoding="utf-8-sig", newline=newline)


def read_text(path: str | Path) -> str:
    with open_text(path, newline="") as fh:
        return fh.read()


def is_blank(text: str) -> bool:
    """True when `text`, with a trailing line end removed, is empty or only spaces and tabs."""
    return text.rstrip("\r\n").strip(" \t") == ""


def sniff_delimiter(text: str) -> str:
    first = next((line for line in text.split("\n") if not is_blank(line)), "")
    return "\t" if first.count("\t") > first.count(",") else ","


def _scan_quotes(text: str, delimiter: str) -> tuple[int | None, bool]:
    """One linear pass with the contract's quote automaton.

    Returns (the physical line where a quoted field still open at the end of the
    text opened, or None; whether any quoted field holds a line break). A line
    with no `"` outside an open quoted field is skipped without a character scan.
    """
    in_quotes = False
    field_started = False
    open_line = 0
    multiline = False
    line_no = 0
    start = 0
    n = len(text)
    while start < n:
        end = text.find("\n", start)
        if end == -1:
            end = n
        line_no += 1
        if in_quotes:
            multiline = True
        elif text.find('"', start, end) == -1:
            start = end + 1
            continue
        else:
            field_started = False
        i = start
        while i < end:
            ch = text[i]
            if in_quotes:
                if ch == '"':
                    if i + 1 < end and text[i + 1] == '"':
                        i += 1
                    else:
                        in_quotes = False
            elif ch == delimiter:
                field_started = False
            elif ch == '"' and not field_started:
                in_quotes = True
                field_started = True
                open_line = line_no
            elif ch != "\r" or i + 1 != end:
                field_started = True
            i += 1
        start = end + 1
    return (open_line if in_quotes else None), multiline


def csv_rows(text: str, label: str) -> list[tuple[int, list[str]]]:
    delimiter = sniff_delimiter(text)
    open_line, multiline = _scan_quotes(text, delimiter)
    if open_line is not None:
        raise DataContractError(f"{label}: line {open_line}: unterminated quoted field (a quote opened here is never closed)")
    reader = csv.reader(io.StringIO(text, newline=""), delimiter=delimiter)
    rows: list[tuple[int, list[str]]] = []
    for row in reader:
        if all(is_blank(cell) for cell in row):
            continue
        if multiline:
            row = [cell.replace("\r\n", "\n") for cell in row]
        rows.append((reader.line_num, row))
    return rows
