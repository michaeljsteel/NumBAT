"""Compatibility checks for the supported plotting package."""

import importlib
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1] / "backend"))


def test_public_plotting_modules_alias_internal_implementations():
    public = importlib.import_module("numbat.plotting")
    internal_modes = importlib.import_module("numbat.plotting.modes")

    assert public.Decorator is internal_modes.Decorator
