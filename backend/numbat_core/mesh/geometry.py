"""Compatibility alias for user geometry support previously in the mesh package."""

import sys

from numbat_core.geometry import usermesh as _implementation

sys.modules[__name__] = _implementation
