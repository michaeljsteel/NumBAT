"""Legacy facade. Implementation lives in numbat_core.materials.materialprops."""
from numbat_core._deprecation import warn_legacy_import
from numbat_core.materials.materialprops import *

warn_legacy_import("materialprops", "numbat.materials")
