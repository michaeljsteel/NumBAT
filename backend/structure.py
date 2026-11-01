"""Compatibility alias for the legacy structure module."""

import sys

from numbat_core.geometry import structure as _implementation

sys.modules[__name__] = _implementation
