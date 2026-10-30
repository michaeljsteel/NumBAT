"""Compatibility alias for legacy gain plotting imports."""

import sys

from numbat_core.plotting import gain as _implementation

sys.modules[__name__] = _implementation
