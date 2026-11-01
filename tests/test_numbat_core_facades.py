"""Checks for supported imports and fatal legacy-import errors."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

BACKEND_DIR = Path(__file__).resolve().parents[1] / "backend"

FACADE_PAIRS = [
    ("nbtypes", "numbat_core.nbtypes"),
    ("reporting", "numbat_core.reporting"),
    ("nbversion", "numbat_core.nbversion"),
    ("nbplatform", "numbat_core.nbplatform"),
    ("nbthread", "numbat_core.utils.threading"),
    ("numbattools", "numbat_core.utils.numeric"),
    ("materials", "numbat_core.materials.materials"),
    ("crystalsyms", "numbat_core.materials.crystalsyms"),
    ("materialprops", "numbat_core.materials.materialprops"),
    ("bulkprops", "numbat_core.materials.bulkprops"),
    ("voigt", "numbat_core.materials.elastic_tensors"),
    ("structure", "numbat_core.geometry.structure"),
    ("usermesh", "numbat_core.geometry.usermesh"),
    ("modes", "numbat_core.solvers.modes"),
    ("modecalcs", "numbat_core.solvers.modecalcs"),
    ("integration", "numbat_core.solvers.sbsgain"),
    ("simresult", "numbat_core.solvers.simresult"),
    ("numbat_core.mesh.geometry", "numbat_core.geometry.usermesh"),
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
    )


@pytest.mark.parametrize("old, _new", FACADE_PAIRS)
def test_legacy_imports_are_fatal(old, _new):
    code = f"""
import importlib
try:
    importlib.import_module({old!r})
except ImportError as exc:
    assert "no longer supported" in str(exc), exc
else:
    raise AssertionError("legacy import unexpectedly succeeded")
"""
    res = _run(code)
    assert res.returncode == 0, res.stderr


@pytest.mark.parametrize("_old, new", FACADE_PAIRS)
def test_numbat_core_module_imports_standalone(_old, new):
    res = _run(f"import importlib; importlib.import_module({new!r})")
    assert res.returncode == 0, res.stderr


@pytest.mark.parametrize(
    "legacy_import, replacement",
    [
        ("import nbtypes", "numbat.types"),
        ("import materials", "numbat.materials"),
        ("import integration", "numbat.solvers.sbsgain"),
        ("import plotting.profiles", "numbat.plotting"),
        ("import meshing.templates", "numbat.mesh"),
        ("import msh.builtin_meshes.basic", "numbat.mesh"),
    ],
)
def test_legacy_import_reports_replacement(legacy_import, replacement):
    code = f"""
try:
    {legacy_import}
except ImportError as exc:
    assert {replacement!r} in str(exc), exc
else:
    raise AssertionError("legacy import unexpectedly succeeded")
"""
    res = _run(code)
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
