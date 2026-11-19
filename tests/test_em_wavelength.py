"""Compatibility tests for EM wavelength units."""

import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1] / "backend"))

from numbat.geometry.structure import Structure
from numbat.solvers import modecalcs


@pytest.fixture
def waveguide(monkeypatch):
    monkeypatch.setattr(modecalcs.FemMesh, "build_from_gmsh_mail", lambda self, wg: None)
    return SimpleNamespace(domain_x=1500.0)


@pytest.mark.parametrize(
    ("kwargs", "expected_um"),
    [({}, 1.55), ({"wl_um": 1.31}, 1.31)],
)
def test_em_simulation_micrometre_wavelength(waveguide, kwargs, expected_um):
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        sim = modecalcs.EMSimulation(waveguide, **kwargs)

    assert sim.wl_um == pytest.approx(expected_um)
    assert sim.lambda_m == pytest.approx(expected_um * 1e-6)
    assert sim.k_0 == pytest.approx(2 * np.pi / sim.lambda_m)


def test_em_simulation_nanometre_wavelength_takes_precedence(waveguide):
    with pytest.warns(DeprecationWarning, match="wl_nm.*wl_um"):
        sim = modecalcs.EMSimulation(waveguide, wl_nm=1550, wl_um=1.31)

    assert sim.wl_um == pytest.approx(1.55)
    assert sim.lambda_m == pytest.approx(1.55e-6)


@pytest.mark.parametrize(
    ("args", "kwargs", "expected_um", "deprecated"),
    [
        ((20, 1550, 1.5), {}, 1.55, True),
        ((20,), {"wl_nm": 1310, "n_eff": 1.5}, 1.31, True),
        ((20,), {"wl_um": 1.31, "n_eff": 1.5}, 1.31, False),
        ((20,), {"wl_nm": 1550, "wl_um": 1.31, "n_eff": 1.5}, 1.55, True),
        ((20,), {"n_eff": 1.5}, 1.55, False),
    ],
)
def test_structure_em_wavelength_paths(
    waveguide, monkeypatch, args, kwargs, expected_um, deprecated
):
    monkeypatch.setattr(modecalcs.EMSimulation, "calc_modes", lambda self, short: None)
    monkeypatch.setattr(modecalcs.EMSimulation, "get_sim_result", lambda self: self)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", DeprecationWarning)
        sim = Structure.calc_EM_modes(waveguide, *args, **kwargs)

    deprecations = [item for item in caught if issubclass(item.category, DeprecationWarning)]
    assert len(deprecations) == int(deprecated)
    if deprecated:
        assert "wl_um" in str(deprecations[0].message)
    assert sim.wl_um == pytest.approx(expected_um)
    assert sim.lambda_m == pytest.approx(expected_um * 1e-6)


def test_structure_requires_effective_index(waveguide):
    with pytest.raises(TypeError, match="n_eff"):
        Structure.calc_EM_modes(waveguide, 20, wl_um=1.55)
