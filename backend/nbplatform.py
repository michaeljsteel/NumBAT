"""Legacy facade. Implementation lives in numbat_core.nbplatform."""
from numbat_core._deprecation import warn_legacy_import
from numbat_core.nbplatform import *

warn_legacy_import("nbplatform", "numbat")
