"""Tests for the public ``numbat`` package facade (Phase 6).

``backend/numbat/`` replaces the legacy flat ``numbat.py``. It must expose
every name the old module did, and the preferred submodules must resolve to
the same objects as the ``numbat_core`` implementation. Runs in a fresh
interpreter with PYTHONPATH=backend.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

BACKEND_DIR = Path(__file__).resolve().parents[1] / "backend"


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


def test_numbat_is_a_package():
    res = _run("import numbat; assert hasattr(numbat, '__path__')")
    assert res.returncode == 0, res.stderr


def test_every_legacy_application_name_is_available():
    code = """
import numbat
import numbat_core.api.application as app
names = [n for n in dir(app) if not (n.startswith('__') and n.endswith('__'))]
missing = [n for n in names if not hasattr(numbat, n)]
assert not missing, missing
mismatch = [n for n in names if getattr(numbat, n) is not getattr(app, n)]
assert not mismatch, mismatch
for n in ('NumBATApp', 'load_simulation', 'register_waveguide_directory',
          'NumBATPlotPrefs', 'version', 'git_hash'):
    assert n in vars(numbat), n
"""
    res = _run(code)
    assert res.returncode == 0, res.stderr


@pytest.mark.parametrize(
    "stmt",
    [
        "from numbat.materials import Material; import materials; assert Material is materials.Material",
        "from numbat.mesh import FemMesh, UserGeometryBase; import femmesh, usermesh; assert FemMesh is femmesh.FemMesh and UserGeometryBase is usermesh.UserGeometryBase",
        "from numbat.geometry import Structure; import structure; assert Structure is structure.Structure",
        "from numbat.solvers import EMSimulation, ACSimulation, Simulation, Mode, EMSimResult; import modecalcs; assert EMSimulation is modecalcs.EMSimulation",
        "import numbat.materials, numbat.mesh, numbat.geometry, numbat.solvers",
    ],
)
def test_public_submodules(stmt):
    res = _run(stmt)
    assert res.returncode == 0, res.stderr


def test_legacy_flat_imports_still_work_after_numbat_package():
    res = _run("import numbat, materials, structure, modecalcs, integration, usermesh, femmesh")
    assert res.returncode == 0, res.stderr
