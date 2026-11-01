"""Compatibility facade for the legacy MAIL mesh module."""

from numbat_core._deprecation import warn_legacy_import
from numbat_core.mesh.mail import *

warn_legacy_import("nbgmsh", "numbat.mesh")
