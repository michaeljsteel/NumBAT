"""Material definitions: ``from numbat.materials import Material``."""

import sys

import materials as _implementation

sys.modules[__name__] = _implementation
