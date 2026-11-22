"""Content-integrity checks for WGDI comparison outputs.

The downstream runner historically resumed on "file exists and is non-empty",
which treats a truncated PNG or an error page written in place of a table as a
finished result. These checks detect the common corruption modes without
false-failing legitimate header-only / zero-block results:

* PNG products must carry the PNG magic bytes (definitive).
* ``.ks`` must have a tab-separated header containing a ``ks_`` column.
* ``.icl`` must contain alignment blocks (or be comment-only).
* ``.bi.csv`` / ``.bi.c.csv`` must have the comma-separated ``id,...`` header
  (header-only is allowed: tandem removal can yield zero blocks).

They are deliberately structural, not biological: a passing check means the
file is a well-formed WGDI product, not that the comparison is meaningful.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


def _check_png(path: Path) -> str:
    with path.open("rb") as handle:
        head = handle.read(8)
    if head != PNG_MAGIC:
        return f"bad_png_magic:{path.name}"
    return ""


def _check_ks(path: Path) -> str:
    lines = [
        line
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.strip()
    ]
    if not lines:
        return f"empty_ks:{path.name}"
    header = lines[0].split("\t")
    if len(header) < 4 or not any("ks_" in column for column in header):
        return f"bad_ks_header:{path.name}"
    return ""


def _check_icl(path: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        return f"empty_icl:{path.name}"
    if "Alignment" in text or all(line.startswith("#") for line in lines):
        return ""
    return f"bad_icl_content:{path.name}"


def _check_bi_csv(path: Path) -> str:
    lines = [
        line
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.strip()
    ]
    if not lines:
        return f"empty_bi:{path.name}"
    header = lines[0]
    if "," not in header or not header.startswith("id"):
        return f"bad_bi_header:{path.name}"
    return ""


def check_comparison_outputs(
    folder: Path, prefix: str, ks_enabled: bool = False
) -> list[str]:
    """Return integrity problems for one comparison's products.

    Default (SOP step 4 core): dotplot + collinearity. With ``ks_enabled``
    the six Ks-stage products are required as well.
    """

    folder = Path(folder)
    specs: list[tuple[str, Callable[[Path], str]]] = [
        (f"{prefix}.png", _check_png),
        (f"{prefix}.icl", _check_icl),
    ]
    if ks_enabled:
        specs += [
            (f"{prefix}.bk.png", _check_png),
            (f"{prefix}.kp.png", _check_png),
            (f"{prefix}.pf.png", _check_png),
            (f"{prefix}.ks", _check_ks),
            (f"{prefix}.bi.csv", _check_bi_csv),
            (f"{prefix}.bi.c.csv", _check_bi_csv),
        ]
    problems: list[str] = []
    for name, checker in specs:
        path = folder / name
        if not path.is_file() or path.stat().st_size == 0:
            problems.append(f"missing_or_empty:{name}")
            continue
        try:
            problem = checker(path)
        except OSError as exc:
            problem = f"unreadable:{name}:{exc}"
        if problem:
            problems.append(problem)
    return problems