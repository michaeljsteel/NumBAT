"""Compatibility alias for the legacy structure module."""

import sys

from numbat_core._deprecation import warn_legacy_import
from numbat_core.geometry import structure as _implementation

warn_legacy_import("structure", "numbat.geometry")

sys.modules[__name__] = _implementation
