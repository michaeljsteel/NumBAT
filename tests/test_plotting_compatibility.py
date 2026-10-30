"""Compatibility checks for the internal plotting package and legacy facades."""

import importlib
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1] / "backend"))


def test_legacy_plotting_modules_alias_internal_implementations():
    module_pairs = [
        ("plotmodes", "numbat_core.plotting.modes"),
        ("plotmoderaw", "numbat_core.plotting.raw_modes"),
        ("plotgain", "numbat_core.plotting.gain"),
        ("plotting.plotprefs", "numbat_core.plotting.preferences"),
        ("plotting.plottools", "numbat_core.plotting.tools"),
        ("plotting.tidyaxes", "numbat_core.plotting.axes"),
        ("plotting.materials", "numbat_core.plotting.materials"),
        ("plotting.profiles", "numbat_core.plotting.profiles"),
        ("plotting.gmsh", "numbat_core.plotting.mesh"),
        ("msh.gmsh", "numbat_core.plotting.mesh"),
    ]

    for legacy_path, internal_path in module_pairs:
        assert importlib.import_module(legacy_path) is importlib.import_module(internal_path)
