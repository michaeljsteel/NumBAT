"""Compatibility alias for legacy material-plotting imports."""

import sys

from numbat_core.plotting import materials as _implementation

sys.modules[__name__] = _implementation
