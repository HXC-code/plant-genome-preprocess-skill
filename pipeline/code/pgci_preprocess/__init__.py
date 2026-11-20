"""Source-tree launcher for ``python -m pgci_preprocess``.

The canonical implementation lives in ``src/pgci_preprocess``.  This package only
extends ``__path__`` so the module command works from the project root without an
editable install; it contains no business logic.

The launcher used to stop after setting ``__path__``.  Because the project root
can precede ``src`` on ``sys.path`` (pytest inserts it during collection), the
submodule imports worked while package-level names such as ``__version__`` and
the documented public API stayed invisible.  The canonical ``__init__`` is now
executed into this namespace so both import paths expose exactly the same names.
"""

from pathlib import Path

_SRC = Path(__file__).resolve().parent.parent / "src" / "pgci_preprocess"
__path__ = [str(_SRC)]

_init = _SRC / "__init__.py"
if not _init.is_file():
    raise ImportError(f"canonical package initialiser is missing: {_init}")

exec(compile(_init.read_text(encoding="utf-8"), str(_init), "exec"), globals())
