"""Compatibility alias for legacy mesh-plotting imports."""

import sys

from numbat_core.plotting import mesh as _implementation

sys.modules[__name__] = _implementation
