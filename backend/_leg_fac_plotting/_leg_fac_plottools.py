"""Compatibility alias for legacy plotting utility imports."""

import sys

from numbat_core.plotting import tools as _implementation

sys.modules[__name__] = _implementation
