"""Compatibility alias for the legacy template-loading module."""

import sys

from numbat_core.mesh import templates as _templates

sys.modules[__name__] = _templates
