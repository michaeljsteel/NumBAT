"""Compatibility alias for the legacy modes module."""

import sys

from numbat_core._deprecation import warn_legacy_import
from numbat_core.solvers import modes as _implementation

warn_legacy_import("modes", "numbat.solvers")

sys.modules[__name__] = _implementation
