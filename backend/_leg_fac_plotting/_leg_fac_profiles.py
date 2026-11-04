"""Compatibility alias for legacy profile-plotting imports."""

import sys

from numbat_core.plotting import profiles as _implementation

sys.modules[__name__] = _implementation
