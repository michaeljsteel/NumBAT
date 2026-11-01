"""Compatibility alias for the legacy mode calculation module."""

import sys

from numbat_core.solvers import modecalcs as _implementation

sys.modules[__name__] = _implementation
