"""Compatibility alias for legacy worker thread and process imports."""

import sys

from numbat_core.utils import threading as _implementation

sys.modules[__name__] = _implementation
