"""Compatibility facade for the legacy FEM mesh module."""

from numbat_core._deprecation import warn_legacy_import
from numbat_core.mesh.fem import *

warn_legacy_import("femmesh", "numbat.mesh")
