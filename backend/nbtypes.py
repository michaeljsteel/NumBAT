"""Legacy facade. Implementation lives in numbat_core.nbtypes."""
from numbat_core._deprecation import warn_legacy_import
from numbat_core.nbtypes import *

warn_legacy_import("nbtypes", "numbat.types")
