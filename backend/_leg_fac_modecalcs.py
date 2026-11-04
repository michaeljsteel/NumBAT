"""Compatibility alias for the legacy mode calculation module."""

import sys

from numbat_core._deprecation import warn_legacy_import
from numbat_core.solvers import modecalcs as _implementation

warn_legacy_import("modecalcs", "numbat.solvers")

sys.modules[__name__] = _implementation
