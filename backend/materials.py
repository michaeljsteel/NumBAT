"""Compatibility alias for the legacy top-level materials module."""

import sys

from numbat_core._deprecation import warn_legacy_import
from numbat_core.materials import materials as _implementation

warn_legacy_import("materials", "numbat.materials")

sys.modules[__name__] = _implementation
