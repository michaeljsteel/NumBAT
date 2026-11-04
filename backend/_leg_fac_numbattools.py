"""Legacy facade. Implementation lives in numbat_core.utils.numeric."""
from numbat_core._deprecation import warn_legacy_import
from numbat_core.utils.numeric import *

warn_legacy_import("numbattools", "numbat.utils")
