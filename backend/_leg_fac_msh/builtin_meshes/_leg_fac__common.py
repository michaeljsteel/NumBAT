"""Compatibility alias for shared built-in geometry helpers."""

import sys

from numbat_core.mesh.builtin import _common as _internal_common

sys.modules[__name__] = _internal_common
