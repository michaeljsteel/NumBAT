"""Sampling in translated mesh coordinates."""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1] / "backend"))

from numbat.io import mode_fields, strain_fields
from numbat.nbtypes import FieldType
from numbat.solvers import modes
from numbat.solvers.mode_sampling import ModeInterpolator
from numbat.solvers.strain_sampling import strain_cut_points
from numbat.strain_selection import StrainComponent, complex_storage_keys


class Mesh:
    def __init__(self):
        self.calls = []

    def get_xy_limits(self):
        return (2e-6, 4e-6, -3e-6, -2e-6)

    def make_interpolator_for_grid(self, x, y, nx, ny):
        self.calls.append((x.copy(), y.copy(), nx, ny))
        return lambda field: field


def test_2d_sampling_uses_same_physical_points_as_display_grid(monkeypatch):
    mesh = Mesh()
    sampler = ModeInterpolator(SimpleNamespace(fem_mesh=mesh))
    sampler.define_plot_grid_2D(10)

    x, y, _, _ = mesh.calls[-1]
    np.testing.assert_allclose(x, sampler.xy_raw["m_x"].flatten("F"))
    np.testing.assert_allclose(y, sampler.xy_raw["m_y"].flatten("F"))
    np.testing.assert_allclose(x * 1e6, sampler.xy_out["m_x"].flatten("F"))
    np.testing.assert_allclose(y * 1e6, sampler.xy_out["m_y"].flatten("F"))

    mode = modes.ModeAC.__new__(modes.ModeAC)
    mode.mode_interpolator = sampler
    mode.sim_result = sampler.sim_result
    mode.mode_num = 0
    sampled = {}
    monkeypatch.setattr(modes.strain_sampling, "calculate_mode_strain", lambda *args: None)

    def interpolate(points):
        sampled["points"] = points.copy()
        return np.zeros((len(points), 6), dtype=complex)

    monkeypatch.setattr(
        modes.strain_sampling, "make_native_strain_interpolator",
        lambda samples: interpolate,
    )
    mode._sample_strain_2d(10)
    np.testing.assert_allclose(sampled["points"], np.column_stack((x, y)))


@pytest.mark.parametrize(
    "cut,val1,val2,expected_points,expected_axis",
    [
        ("x", 3.0, None, [[3, -3], [3, -2.5], [3, -2]], [-3, -2.5, -2]),
        ("y", -2.5, None, [[2, -2.5], [3, -2.5], [4, -2.5]], [2, 3, 4]),
        ("line", (2, -3), (4, -2), [[2, -3], [3, -2.5], [4, -2]],
         [0, np.sqrt(5) / 2, np.sqrt(5)]),
    ],
)
def test_cut_points_and_export_axes_match_translated_mesh(
    monkeypatch, tmp_path, cut, val1, val2, expected_points, expected_axis
):
    mesh = Mesh()
    result = SimpleNamespace(fem_mesh=mesh)
    sampler = ModeInterpolator(result)
    sampler.define_plot_grid_1D(cut, val1, val2, 3)
    x, y, _, _ = mesh.calls[-1]
    strain_axis, strain_points = strain_cut_points(result, cut, val1, val2, 3)

    np.testing.assert_allclose(np.column_stack((x, y)) * 1e6, expected_points)
    np.testing.assert_allclose(strain_points * 1e6, expected_points)
    np.testing.assert_allclose(sampler.v_x_axis, expected_axis)
    np.testing.assert_allclose(strain_axis, expected_axis)

    monkeypatch.setattr(mode_fields.runtime, "outdir_fields_path", lambda prefix: tmp_path)
    monkeypatch.setattr(strain_fields, "outdir_fields_path", lambda prefix: tmp_path)
    field_keys = ("Fxr", "Fxi", "Fyr", "Fyi", "Fzr", "Fzi", "Fa")
    mode_fields.write_mode_fields_1d(
        prefix="", field_type=FieldType.EM_E, mode_index=0, cut=cut,
        x_axis=sampler.v_x_axis, fields={key: np.zeros(3) for key in field_keys},
    )
    keys = [key for component in StrainComponent for key in complex_storage_keys(component)]
    strain_fields.write_strain_fields_1d(
        prefix="", mode_index=0, cut=cut,
        x_axis=strain_axis, fields={key: np.zeros(3) for key in keys},
    )
    np.testing.assert_allclose(
        np.loadtxt(tmp_path / f"EM_mode_00_{cut}cut.txt")[:, 0], expected_axis, atol=1e-8
    )
    np.testing.assert_allclose(
        np.loadtxt(tmp_path / f"AC_strain_mode_00_{cut}cut.txt")[:, 0], expected_axis, atol=1e-8
    )


def test_r0_offset_changes_reporting_without_moving_raw_centre():
    mode = modes.ModeAC.__new__(modes.ModeAC)
    mode.r0 = np.array([3.0, -2.0])
    mode.set_r0_offset(1.0, 0.5)
    np.testing.assert_allclose(mode.center_of_mass(), [2.0, -2.5])
    np.testing.assert_allclose(mode.r0, [3.0, -2.0])
