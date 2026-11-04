"""Checks for supported internal NumBAT imports."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

BACKEND_DIR = Path(__file__).resolve().parents[1] / "backend"

NUMBAT_MODULES = [
    "numbat.nbtypes",
    "numbat.reporting",
    "numbat.nbversion",
    "numbat.nbplatform",
    "numbat.utils.threading",
    "numbat.utils.numeric",
    "numbat.materials.materials",
    "numbat.materials.crystalsyms",
    "numbat.materials.materialprops",
    "numbat.materials.bulkprops",
    "numbat.materials.elastic_tensors",
    "numbat.geometry.structure",
    "numbat.geometry.usermesh",
    "numbat.solvers.modes",
    "numbat.solvers.modecalcs",
    "numbat.solvers.sbsgain",
    "numbat.solvers.simresult",
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


@pytest.mark.parametrize("module", NUMBAT_MODULES)
def test_numbat_module_imports_standalone(module):
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
