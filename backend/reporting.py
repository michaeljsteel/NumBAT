"""Legacy facade. Implementation lives in numbat_core.reporting."""
from numbat_core._deprecation import warn_legacy_import
from numbat_core.reporting import *

warn_legacy_import("reporting", "numbat")
