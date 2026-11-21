"""Single source of truth for the package version.

Every version string that reaches a QC record, a conversion manifest or the CLI
must derive from this module.  ``pyproject.toml`` repeats the same value because
setuptools needs a static field; ``tests/test_versioning.py`` fails the build if
the two ever drift apart.

The historical drift (package 2.3.0 while ``PIPELINE_VERSION`` said 2.1.0 and
``CONVERTER_VERSION`` said 2.2.0) made provenance in ``*.qc.json`` ambiguous, so
implementation versions are no longer tracked independently.
"""

from __future__ import annotations

__version__ = "3.5.3"




