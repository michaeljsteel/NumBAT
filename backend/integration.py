"""Compatibility alias for the legacy integration module."""

import sys

from numbat_core._deprecation import warn_legacy_import
from numbat_core.solvers import integration as _implementation

warn_legacy_import("integration", "numbat.solvers")

sys.modules[__name__] = _implementation
