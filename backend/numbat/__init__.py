"""Public NumBAT package.

This package is a thin facade over the implementation in ``numbat_core``.
It replaces the former flat ``numbat.py`` module: every name that module
exposed (``NumBATApp``, ``load_simulation``, ...) is still available as
``numbat.<name>``. Preferred submodules are ``numbat.materials``,
``numbat.mesh``, ``numbat.geometry`` and ``numbat.solvers``.
"""

from numbat_core.api import application as _application
from numbat_core.api.application import *  # noqa: F401,F403
from numbat_core.api.application import (
    NumBATApp,
    NumBATPlotPrefs,
    git_hash,
    load_simulation,
    register_waveguide_directory,
    version,
)


def __getattr__(name):
    # Anything else the legacy module exposed (including private helpers and
    # live module state) is looked up on the implementation module.
    try:
        return getattr(_application, name)
    except AttributeError:
        raise AttributeError(f"module 'numbat' has no attribute {name!r}") from None


def __dir__():
    return sorted(set(globals()) | set(dir(_application)))
