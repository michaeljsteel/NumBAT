"""Compatibility alias for legacy raw-mode plotting imports."""

import sys

from numbat_core._deprecation import warn_legacy_import
from numbat_core.plotting import raw_modes as _implementation

warn_legacy_import("plotmoderaw", "numbat.plotting")

sys.modules[__name__] = _implementation
