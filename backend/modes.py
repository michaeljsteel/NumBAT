"""Compatibility alias for the legacy modes module."""

import sys

from numbat_core.solvers import modes as _implementation

sys.modules[__name__] = _implementation
