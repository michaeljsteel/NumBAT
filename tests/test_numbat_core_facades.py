"""Compatibility checks for moved implementations in ``numbat_core``.

Each flat module at the backend root is a facade over its ``numbat_core``
home. Every public name in the new module must resolve through the facade to
the same object, so existing ``import nbtypes`` style code keeps working.
Runs in a fresh interpreter with PYTHONPATH=backend.
"""

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
    ("numbat", "numbat_core.api.application"),
    ("numbat_core.mesh.geometry", "numbat_core.geometry.usermesh"),
]

SIMULATION_FACADE_PAIRS = [
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


@pytest.mark.parametrize("old, new", FACADE_PAIRS)
def test_facade_exposes_new_module_public_names_by_identity(old, new):
    code = f"""
import importlib
old = importlib.import_module({old!r})
new = importlib.import_module({new!r})
public = [n for n in dir(new) if not n.startswith('_')]
missing = [n for n in public if not hasattr(old, n)]
assert not missing, f'facade missing names: {{missing}}'
mismatch = [n for n in public if getattr(old, n) is not getattr(new, n)]
assert not mismatch, f'facade names not identical: {{mismatch}}'
print('ok', len(public))
"""
    res = _run(code)
    assert res.returncode == 0, res.stderr
    assert res.stdout.startswith("ok")


@pytest.mark.parametrize("_old, new", FACADE_PAIRS)
def test_numbat_core_module_imports_standalone(_old, new):
    res = _run(f"import importlib; importlib.import_module({new!r})")
    assert res.returncode == 0, res.stderr


@pytest.mark.parametrize("old, new", SIMULATION_FACADE_PAIRS)
def test_simulation_legacy_paths_are_module_aliases(old, new):
    code = f"""
import importlib
old = importlib.import_module({old!r})
new = importlib.import_module({new!r})
assert old is new
"""
    res = _run(code)
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
def test_legacy_import_warns_at_user_import(legacy_import, replacement):
    code = f"""
import warnings
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    {legacy_import}

deprecation_warnings = [
    warning for warning in caught
    if issubclass(warning.category, FutureWarning)
]
assert len(deprecation_warnings) == 1, deprecation_warnings
warning = deprecation_warnings[0]
assert warning.filename == "<string>", warning.filename
assert warning.lineno == 5, warning.lineno
assert {replacement!r} in str(warning.message), warning.message
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
