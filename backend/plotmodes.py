"""Compatibility alias for legacy mode plotting imports."""

import sys

from numbat_core._deprecation import warn_legacy_import
from numbat_core.plotting import modes as _implementation

warn_legacy_import("plotmodes", "numbat.plotting")

sys.modules[__name__] = _implementation
