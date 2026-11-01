"""Compatibility alias for the legacy integration module."""

import sys

from numbat_core.solvers import integration as _implementation

sys.modules[__name__] = _implementation
