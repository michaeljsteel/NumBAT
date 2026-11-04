"""Compatibility alias for the legacy SBS gain module."""

import sys

from numbat_core._deprecation import warn_legacy_import
from numbat_core.solvers import sbsgain as _implementation

warn_legacy_import("integration", "numbat.solvers.sbsgain")

sys.modules[__name__] = _implementation
