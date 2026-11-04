"""Compatibility alias for legacy plot-preference imports."""

import sys

from numbat_core.plotting import preferences as _implementation

sys.modules[__name__] = _implementation
