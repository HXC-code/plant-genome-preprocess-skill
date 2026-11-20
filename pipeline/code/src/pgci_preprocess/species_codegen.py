"""Derive species abbreviations from Latin names when the table has no row.

Rule (as agreed for this project):

1. 3 letters: genus initial + two letters taken from the specific epithet.
2. The default choice is the first two epithet letters (``Cercis chinensis``
   starts with ``Cch``); when that code is already taken, try the other
   two-letter combinations inside the epithet in a fixed order.
3. When every 3-letter combination is taken, widen to 4 letters
   (genus initial + three epithet letters), then 5, and so on.
4. Only when the epithet is too short to widen further do we append a
   numeric suffix, so code generation never fails.

Codes are compared case-insensitively; the first letter is upper case and
the rest lower case, matching the existing ``species_abbreviations.tsv``.
"""

from __future__ import annotations

import re
from itertools import combinations

LETTER_RE = re.compile(r"[A-Za-z]+")


def split_latin_name(latin_name: str) -> tuple[str, str]:
    """Return ``(genus, epithet)`` from a Latin binomial with optional authority."""

    tokens = [token for token in LETTER_RE.findall(str(latin_name or ""))]
    if len(tokens) < 2:
        raise ValueError(f"not a Latin binomial: {latin_name!r}")
    genus, epithet = tokens[0], tokens[1]
    return genus, epithet.lower()


def _candidates(genus: str, epithet: str):
    """Yield candidate codes in priority order, 3 letters first."""

    initial = genus[0].upper()
    for width in range(2, len(epithet) + 1):
        for indexes in combinations(range(len(epithet)), width):
            yield initial + "".join(epithet[index] for index in indexes)
    # Epithet exhausted (e.g. "Pinus taeda" after 4 letters): keep widening
    # with the epithet prefix plus a numeric suffix so the caller always gets
    # a usable, unique code.
    stem = initial + epithet
    for suffix in range(2, 100):
        yield f"{stem}{suffix}"


def derive_code(latin_name: str, occupied: set[str] | None = None) -> str:
    """Return the first free code for ``latin_name`` given ``occupied`` codes."""

    taken = {code.casefold() for code in (occupied or set())}
    genus, epithet = split_latin_name(latin_name)
    for candidate in _candidates(genus, epithet):
        if candidate.casefold() not in taken:
            return candidate
    raise ValueError(f"cannot derive a species code for {latin_name!r}")


def derive_codes(latin_names: list[str], occupied: set[str] | None = None) -> list[str]:
    """Assign codes in list order, each one avoiding earlier assignments."""

    taken = {code.casefold() for code in (occupied or set())}
    result: list[str] = []
    for name in latin_names:
        code = derive_code(name, taken)
        taken.add(code.casefold())
        result.append(code)
    return result
