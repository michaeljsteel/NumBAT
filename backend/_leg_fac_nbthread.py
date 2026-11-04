"""Compatibility alias for legacy worker thread and process imports."""

import sys

from numbat_core._deprecation import warn_legacy_import
from numbat_core.utils import threading as _implementation

warn_legacy_import("nbthread", "numbat.parallel")

sys.modules[__name__] = _implementation
