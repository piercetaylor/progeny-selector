"""Position text shared by the VCF, HapMap, wide-CSV and markers.csv readers (contract/data-contract.md 1.2.0, "Genotype file";
mirrors isoline-browser src/io/position.ts).

Responsibility: read one position cell as a non-negative whole number. The cell, trimmed as JavaScript String.prototype.trim does,
must be ASCII decimal digits with an optional sign, fraction and exponent ("1000", "+1000", "1000.0", "1e3", "1.0E3",
"1.9E+07"), and the value must be finite, equal to its floor and not negative. ValueError naming the cell for anything else:
empty, a non-zero fraction ("100.7",
"1e-3"), a negative number, nan, inf, hexadecimal, "_" or "," separators. The grammar is checked before float() because
float() alone accepts "1_000", "nan", "inf" and non-ASCII digits, and int(float("inf")) raises OverflowError. Plain digits
are converted with int() so they never round-trip through float. Callers wrap the ValueError in DataContractError with the
line number. VCF POS uses the "digits" grammar: ASCII decimal digits only, as the VCF specification declares POS an Integer.

Interface:
    parse_position(text: str, grammar: str = "number") -> int   (grammar "number" or "digits")
"""

from __future__ import annotations

import math
import re

_POSITION = re.compile(r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?", re.ASCII)
_DIGITS = re.compile(r"[+-]?\d+", re.ASCII)
_DIGITS_ONLY = re.compile(r"\d+", re.ASCII)
# The characters JavaScript String.prototype.trim removes; str.strip() also removes \x1c-\x1f and \x85.
_JS_WHITESPACE = (
    " \t\n\v\f\r\u00a0\u1680\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u2028\u2029\u202f\u205f\u3000\ufeff"
)


def parse_position(text: str, grammar: str = "number") -> int:
    if grammar not in ("number", "digits"):
        raise ValueError(f"unknown position grammar {grammar!r}")
    cell = str(text).strip(_JS_WHITESPACE)
    if grammar == "digits":
        digits_message = f"invalid position {cell!r} (expected decimal digits; VCF POS is an Integer)"
        if _DIGITS_ONLY.fullmatch(cell) is None:
            raise ValueError(digits_message)
        try:
            return int(cell)
        except ValueError as exc:  # more digits than sys.get_int_max_str_digits()
            raise ValueError(digits_message) from exc
    message = f"invalid position {cell!r} (expected a whole number such as 1000, 1000.0 or 1e3)"
    if _POSITION.fullmatch(cell) is None:
        raise ValueError(message)
    if _DIGITS.fullmatch(cell):
        try:
            value = int(cell)
        except ValueError as exc:  # more digits than sys.get_int_max_str_digits()
            raise ValueError(message) from exc
    else:
        as_float = float(cell)
        if not math.isfinite(as_float) or as_float != math.floor(as_float):
            raise ValueError(message)
        value = int(as_float)
    if value < 0:
        raise ValueError(message)
    return value
