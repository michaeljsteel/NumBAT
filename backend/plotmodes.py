"""Compatibility alias for legacy mode plotting imports."""

import sys

from numbat_core.plotting import modes as _implementation

sys.modules[__name__] = _implementation
