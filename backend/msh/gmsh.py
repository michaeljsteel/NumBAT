"""Compatibility alias for legacy mesh-rendering imports."""

import sys

from numbat_core.plotting import mesh as _implementation

sys.modules[__name__] = _implementation
