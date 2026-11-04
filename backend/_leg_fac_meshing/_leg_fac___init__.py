"""Compatibility package for legacy meshing imports."""

from numbat_core._deprecation import warn_legacy_import

warn_legacy_import("meshing", "numbat.mesh")