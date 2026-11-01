"""Compatibility alias for the legacy simulation result module."""

import sys

from numbat_core.solvers import simresult as _implementation

sys.modules[__name__] = _implementation
