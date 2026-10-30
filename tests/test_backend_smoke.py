"""Smoke tests for the backend paths that the restructuring must preserve.

These exercise real structure construction, meshing, EM mode solving, built-in
template registration and geometry pickling, so they are slower than unit tests.
"""

import pickle
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1] / "backend"))

import materials
import numbat
import simresult
from meshing import templates
from usermesh import UserGeometryBase

RECT_DOMAIN_NM = 1500.0
BUILTIN_SHAPES = [
    "circular",
    "rectangular",
    "slot",
    "rib",
    "onion",
    "pedestal",
    "triangular",
    "rounded_bowtie",
    "rounded_wedge",
]


@pytest.fixture
def nbapp(tmp_path, monkeypatch):
    # NumBATApp is a singleton. Reset it so each test writes to its own tmp_path.
    monkeypatch.setattr(numbat._SingletonMeta, "_instances", {})
    prefs_file = tmp_path / "numbat.toml"
    prefs_file.write_text("")
    return numbat.NumBATApp(
        "phase0_smoke",
        outdir=str(tmp_path / "out"),
        user_settings_file=str(prefs_file),
    )


def _make_rectangular(nbapp):
    return nbapp.make_structure(
        "rectangular",
        RECT_DOMAIN_NM,
        RECT_DOMAIN_NM,
        300.0,
        280.0,
        material_bkg=materials.make_material("Vacuum"),
        material_a=materials.make_material("Si_2016_Smith"),
        lc_bkg=0.2,
        lc_refine_1=2.5,
        lc_refine_2=2.5,
    )


def test_builtin_waveguide_construction(nbapp):
    wg = _make_rectangular(nbapp)

    assert wg.inc_shape == "rectangular"
    assert wg.wg_geom._shape_name == "rectangular"
    assert type(wg.wg_geom).__module__.startswith("msh.builtin_meshes")
    assert np.real(wg.get_material("a").refindex_n) > 3.0


def test_mesh_generation_writes_nonempty_mesh_file(nbapp):
    wg = _make_rectangular(nbapp)

    mesh_file = Path(wg.mesh_mail_fname)
    assert mesh_file.is_file()
    assert mesh_file.stat().st_size > 0


def test_em_simulation_object_is_constructed(nbapp):
    wg = _make_rectangular(nbapp)
    n_eff_guess = np.real(wg.get_material("a").refindex_n) - 0.1

    sim = wg.calc_EM_modes(2, 1550.0, n_eff_guess)

    assert isinstance(sim, simresult.EMSimResult)
    neff_0 = np.real(sim.neff(0))
    neff_1 = np.real(sim.neff(1))
    assert np.isfinite(neff_0) and np.isfinite(neff_1)
    assert neff_0 > neff_1


@pytest.mark.parametrize("shape", BUILTIN_SHAPES)
def test_builtin_template_is_registered(nbapp, shape):
    registered_shapes = set()
    for entry in templates.g_waveguide_templates:
        inc_shape = entry["inc_shape"]
        registered_shapes.update([inc_shape] if isinstance(inc_shape, str) else inc_shape)
        # Built-in entries carry a resolved class. Legacy user entries carry 'wg_impl' instead.
        if "wg_template_cls" in entry:
            assert issubclass(entry["wg_template_cls"], UserGeometryBase)

    assert shape in registered_shapes


def test_builtin_geometry_pickle_round_trip(nbapp):
    geometry = _make_rectangular(nbapp).wg_geom

    rebuilt = pickle.loads(pickle.dumps(geometry))

    assert type(rebuilt) is type(geometry)
    assert rebuilt._shape_name == geometry._shape_name
    assert rebuilt._d_params["domain_x"] == geometry._d_params["domain_x"]
    assert rebuilt.__dict__.keys() == geometry.__dict__.keys()
