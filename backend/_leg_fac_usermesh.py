"""Compatibility alias for legacy user geometry imports."""

import sys

from numbat_core._deprecation import warn_legacy_import
from numbat_core.geometry import usermesh as _implementation

warn_legacy_import("usermesh", "numbat.geometry")

sys.modules[__name__] = _implementation
