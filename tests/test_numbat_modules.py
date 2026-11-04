"""Checks for supported internal NumBAT imports."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

BACKEND_DIR = Path(__file__).resolve().parents[1] / "backend"

NUMBAT_CORE_MODULES = [
    "numbat_core.nbtypes",
    "numbat_core.reporting",
    "numbat_core.nbversion",
    "numbat_core.nbplatform",
    "numbat_core.utils.threading",
    "numbat_core.utils.numeric",
    "numbat_core.materials.materials",
    "numbat_core.materials.crystalsyms",
    "numbat_core.materials.materialprops",
    "numbat_core.materials.bulkprops",
    "numbat_core.materials.elastic_tensors",
    "numbat_core.geometry.structure",
    "numbat_core.geometry.usermesh",
    "numbat_core.solvers.modes",
    "numbat_core.solvers.modecalcs",
    "numbat_core.solvers.sbsgain",
    "numbat_core.solvers.simresult",
]

def _run(code: str) -> subprocess.CompletedProcess:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(BACKEND_DIR)
    return subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        env=env,
        cwd=BACKEND_DIR,
        timeout=120,
        check=False,
    )


@pytest.mark.parametrize("module", NUMBAT_CORE_MODULES)
def test_numbat_core_module_imports_standalone(module):
    res = _run(f"import importlib; importlib.import_module({module!r})")
    assert res.returncode == 0, res.stderr


def test_public_import_does_not_warn():
    code = """
import warnings
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    import numbat.types

assert not [
    warning for warning in caught
    if issubclass(warning.category, FutureWarning)
], caught
"""
    res = _run(code)
    assert res.returncode == 0, res.stderr
