"""Compatibility alias for legacy mesh-rendering imports."""

import sys

from numbat_core._deprecation import warn_legacy_import
from numbat_core.plotting import mesh as _implementation

warn_legacy_import("msh.gmsh", "numbat.mesh")

sys.modules[__name__] = _implementation
