"""Legacy facade. Implementation lives in numbat_core.materials.crystalsyms."""
from numbat_core._deprecation import warn_legacy_import
from numbat_core.materials.crystalsyms import *

warn_legacy_import("crystalsyms", "numbat.materials")
