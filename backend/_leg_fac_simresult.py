"""Compatibility alias for the legacy simulation result module."""

import sys

from numbat_core._deprecation import warn_legacy_import
from numbat_core.solvers import simresult as _implementation

warn_legacy_import("simresult", "numbat.solvers")

sys.modules[__name__] = _implementation
