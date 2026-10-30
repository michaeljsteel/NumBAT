"""Compatibility alias for legacy axis-helper imports."""

import sys

from numbat_core.plotting import axes as _implementation

sys.modules[__name__] = _implementation
