"""Public material API: ``from numbat.materials import Material``."""

import sys

from numbat_core.materials import materials as _implementation

sys.modules[__name__] = _implementation
