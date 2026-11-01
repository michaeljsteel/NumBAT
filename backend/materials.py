"""Compatibility alias for the legacy top-level materials module."""

import sys

from numbat_core.materials import materials as _implementation

sys.modules[__name__] = _implementation
