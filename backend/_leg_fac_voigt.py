"""Legacy facade. Implementation lives in numbat_core.materials.elastic_tensors."""
from numbat_core._deprecation import warn_legacy_import
from numbat_core.materials.elastic_tensors import *

warn_legacy_import("voigt", "numbat.materials")
