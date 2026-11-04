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
        check=False,
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
        "from numbat.materials import Material; assert Material",
        "from numbat.mesh import FemMesh, UserGeometryBase; assert FemMesh and UserGeometryBase",
        "from numbat.geometry import Structure; assert Structure",
        "from numbat.solvers import EMSimulation, ACSimulation, Simulation, Mode, EMSimResult; assert EMSimulation and ACSimulation and Simulation and Mode and EMSimResult",
        "import numbat.materials, numbat.mesh, numbat.geometry, numbat.solvers",
        "from numbat.types import SI_GHz, twopi, PointGroup",
        "from numbat.mesh import MailData; assert MailData",
        "from numbat.parallel import launch_worker_processes_and_wait, launch_worker_threads_and_wait",
        "from numbat.plotting import Decorator; assert Decorator",
        "from numbat.plotting.ipython import img_single, img_pair, img_triple",
        "import numbat.utils",
    ],
)
def test_public_submodules(stmt):
    res = _run(stmt)
    assert res.returncode == 0, res.stderr
