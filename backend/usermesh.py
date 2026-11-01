"""Compatibility alias for legacy user geometry imports."""

import sys

from numbat_core.geometry import usermesh as _implementation

sys.modules[__name__] = _implementation
