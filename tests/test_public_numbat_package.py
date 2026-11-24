"""Tests for the public ``numbat`` package facade (Phase 6).

``backend/numbat/`` replaces the legacy flat ``numbat.py``. It must expose
every name the old module did, and the preferred submodules must resolve to
the same objects as the ``numbat`` implementation. Runs in a fresh
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
import numbat.api.application as app
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


def test_reload_plot_preferences_replaces_global_preferences(tmp_path):
    first_prefs = tmp_path / "first.toml"
    first_prefs.write_text(
        '[all_plots]\nem_colormap_signed = "viridis"\nvector_arrow_scale = 4\n'
        '[single_plots]\ntitle_fs = 20\n'
    )
    second_prefs = tmp_path / "second.toml"
    second_prefs.write_text('[all_plots]\nac_colormap_signed = "plasma"\n')
    code = f"""
import numbat
prefs = numbat.reload_plot_preferences({str(first_prefs)!r})
assert numbat.NumBATPlotPrefs() is prefs
assert prefs.cmap_em_field_signed == 'viridis'
assert prefs.vecplot_arrow_scale == 4
assert prefs.d_single['title_fs'] == 20
prefs = numbat.reload_plot_preferences({str(second_prefs)!r})
assert numbat.NumBATPlotPrefs() is prefs
assert prefs.cmap_ac_field_signed == 'plasma'
assert prefs.cmap_em_field_signed == 'seismic'
assert prefs.vecplot_arrow_scale == 10
assert prefs.d_single['title_fs'] == 16
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
