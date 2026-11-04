"""Legacy facade. Implementation lives in numbat_core.nbversion."""
from numbat_core._deprecation import warn_legacy_import
from numbat_core.nbversion import *

warn_legacy_import("nbversion", "numbat")
