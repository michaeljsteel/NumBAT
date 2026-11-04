"""Compatibility alias for legacy IPython plotting imports."""

import sys

from numbat_core.plotting import ipython as _implementation

sys.modules[__name__] = _implementation
