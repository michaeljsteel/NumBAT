"""Compatibility alias for legacy raw-mode plotting imports."""

import sys

from numbat_core.plotting import raw_modes as _implementation

sys.modules[__name__] = _implementation
