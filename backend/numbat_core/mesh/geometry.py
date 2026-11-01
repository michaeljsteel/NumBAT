"""Compatibility alias for user geometry support previously in the mesh package."""

import sys

from numbat_core._deprecation import warn_legacy_import
from numbat_core.geometry import usermesh as _implementation

warn_legacy_import("numbat_core.mesh.geometry", "numbat.geometry")

sys.modules[__name__] = _implementation
