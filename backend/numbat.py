"""Compatibility alias for the legacy application module."""

import sys

from numbat_core.api import application as _implementation

sys.modules[__name__] = _implementation
