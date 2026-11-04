"""Legacy facade. Implementation lives in numbat_core.materials.bulkprops."""
from numbat_core._deprecation import warn_legacy_import
from numbat_core.materials.bulkprops import *

warn_legacy_import("bulkprops", "numbat.materials")
