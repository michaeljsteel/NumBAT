"""Analytical regression tests for native acoustic P2 strain sampling."""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from matplotlib.figure import Figure

sys.path.append(str(Path(__file__).resolve().parents[1] / "backend"))

import numbat.solvers.modes as modes_module
from numbat.plotting import strain_rendering
from numbat.plotting.modes import PlotParams2D
from numbat.plotting.strain_native import plot_native_strain_component
from numbat.solvers.modes import ModeAC
from numbat.solvers.simresult import ACSimResult
from numbat.solvers.strain_sampling import (
    P2_REFERENCE_NODES,
    NativeStrainSamples,
    calculate_element_strain,
    calculate_mode_strain,
    make_native_strain_interpolator,
    sample_native_strain,
)
from numbat.strain_selection import (
    StrainComponent,
    StrainSelection,
    complex_storage_keys,
    selection_from_storage_key,
    storage_key,
)


def _interface_strain_samples():
    coordinates = np.array((
        (
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 1.0),
            (0.5, 0.0),
            (0.5, 0.5),
            (0.0, 0.5),
        ),
        (
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (1.0, 0.5),
            (0.5, 1.0),
            (0.5, 0.5),
        ),
    ))
    values = np.empty((2, 6, 6), dtype=np.complex128)
    values[0, :, :] = 1 + 2j
    values[1, :, :] = 3 + 4j
    return NativeStrainSamples(
        reference_points=P2_REFERENCE_NODES,
        coordinates=coordinates,
        values=values,
        element_material_indices=np.array((1, 2)),
    )


def test_strain_selection_primary_storage_matches_acoustic_phase_convention():
    expected_primary_keys = {
        StrainComponent.XX: "Sxxr",
        StrainComponent.YY: "Syyr",
        StrainComponent.ZZ: "Szzr",
        StrainComponent.YZ: "Syzi",
        StrainComponent.XZ: "Sxzi",
        StrainComponent.XY: "Sxyr",
    }

    for component, expected_key in expected_primary_keys.items():
        selection = StrainSelection.primary(component)
        assert storage_key(selection) == expected_key
        assert selection_from_storage_key(expected_key) == selection
        assert complex_storage_keys(component) == (
            f"S{component.value}r",
            f"S{component.value}i",
        )


def test_native_p2_strain_matches_an_affine_displacement_on_curvilinear_element():
    coordinates = np.array((
        (
            (0.0, 0.0),
            (2.0, 0.0),
            (0.0, 1.5),
            (1.1, -0.1),
            (1.1, 0.8),
            (-0.1, 0.8),
        ),
    ))
    q_ac = 3.0
    gradients = np.array((
        (2.0 + 1.0j, -1.0 + 2.0j),
        (4.0 - 3.0j, 5.0 + 1.0j),
        (-2.0 + 4.0j, 3.0 - 2.0j),
    ))
    offsets = np.array((1.0 - 1.0j, -2.0 + 3.0j, 4.0 + 2.0j))
    displacements = offsets + np.einsum("eac,kc->eak", coordinates, gradients)

    strain = calculate_element_strain(displacements, coordinates, q_ac)
    points = np.einsum(
        "pa,eac->epc",
        np.array((
            (1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (0.0, 1.0, 0.0, 0.0, 0.0, 0.0),
            (0.0, 0.0, 1.0, 0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 1.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 1.0),
        )),
        coordinates,
    )
    expected = np.empty((1, len(P2_REFERENCE_NODES), 6), dtype=np.complex128)
    expected[:, :, 0] = gradients[0, 0]
    expected[:, :, 1] = gradients[1, 1]
    expected[:, :, 2] = 1j * q_ac * (
        offsets[2] + np.einsum("epc,c->ep", points, gradients[2])
    )
    expected[:, :, 3] = 0.5 * (
        gradients[2, 1]
        + 1j * q_ac * (offsets[1] + np.einsum("epc,c->ep", points, gradients[1]))
    )
    expected[:, :, 4] = 0.5 * (
        gradients[2, 0]
        + 1j * q_ac * (offsets[0] + np.einsum("epc,c->ep", points, gradients[0]))
    )
    expected[:, :, 5] = 0.5 * (gradients[0, 1] + gradients[1, 0])
    np.testing.assert_allclose(strain, expected)


def test_mode_strain_rejects_non_acoustic_simulations():
    sim_result = SimpleNamespace(is_AC=lambda: False)

    with pytest.raises(ValueError, match="Acoustic strain"):
        calculate_mode_strain(sim_result, 0)


def test_mode_strain_uses_acoustic_fem_storage_and_element_coordinates():
    coordinates = np.array((
        (0.0, 0.0),
        (1.0, 0.0),
        (0.0, 1.0),
        (0.5, 0.0),
        (0.5, 0.5),
        (0.0, 0.5),
    ))
    displacements = np.array(
        [[(2 * x + 3 * y, 5 * x + 7 * y, 11 * x + 13 * y) for x, y in coordinates]],
        dtype=np.complex128,
    )
    fem_evecs = np.zeros((3, 6, 1, 1), dtype=np.complex128)
    fem_evecs[:, :, 0, 0] = displacements[0].T
    fem_mesh = SimpleNamespace(
        m_elnd_to_mshpt=np.arange(1, 7, dtype=np.int64)[:, np.newaxis],
        v_mshpt_xy=coordinates.T,
    )
    sim_result = SimpleNamespace(
        is_AC=lambda: True,
        fem_evecs=fem_evecs,
        fem_mesh=fem_mesh,
        q_AC=0.0,
    )

    samples = calculate_mode_strain(sim_result, mode_index=0)

    expected = calculate_element_strain(displacements, coordinates[np.newaxis], 0.0)
    np.testing.assert_allclose(samples.values, expected)
    np.testing.assert_allclose(
        samples.coordinates, coordinates[np.newaxis]
    )


def test_native_strain_sampling_preserves_interface_jump_and_selects_first_side():
    samples = _interface_strain_samples()

    sampled = sample_native_strain(
        samples,
        np.array(((0.25, 0.75), (0.1, 0.1), (2.0, 2.0))),
    )

    np.testing.assert_allclose(sampled[0], np.full(6, 1 + 2j))
    np.testing.assert_allclose(sampled[1], np.full(6, 1 + 2j))
    assert np.all(np.isnan(sampled[2]))
    assert samples.subtriangles().shape == (8, 3)
    assert samples.element_material_indices.tolist() == [1, 2]


def test_native_strain_renderer_preserves_local_display_values():
    samples = _interface_strain_samples()
    fig = Figure(figsize=(4, 3))
    ax = fig.add_subplot()

    _, _, colorbar = plot_native_strain_component(
        samples, StrainSelection.primary(StrainComponent.YZ), ax=ax
    )

    assert len(ax.collections) == 1
    assert colorbar.mappable.get_array().shape == (12,)


def test_strain_renderer_applies_standard_acoustic_trim_limits():
    figure = Figure()
    axis = figure.subplots()
    plot_params = PlotParams2D()
    plot_params.update({
        "trim": (0.1, 0.25, 0.2, 0),
        "xlim_min": 0.9,
        "xlim_max": 0.9,
        "ylim_min": 0.9,
        "ylim_max": 0.9,
        "ticks": False,
        "title": False,
        "frame": False,
    })

    with pytest.deprecated_call(match="xlim_min.*deprecated"):
        strain_rendering._plot_component(
            axis,
            {"v_x": np.array((0.0, 1.0)), "v_y": np.array((0.0, 1.0))},
            {"Sxxr": np.ones((2, 2))},
            StrainSelection.primary(StrainComponent.XX),
            plot_params,
        )

    assert axis.get_xlim() == pytest.approx((0.1, 0.8))
    assert axis.get_ylim() == pytest.approx((0.25, 1.0))


def test_fast_native_interpolator_matches_reference_sampler_away_from_interface():
    samples = _interface_strain_samples()
    points = np.array(((0.1, 0.1), (0.9, 0.9), (2.0, 2.0)))

    np.testing.assert_allclose(
        make_native_strain_interpolator(samples)(points),
        sample_native_strain(samples, points),
        equal_nan=True,
    )


def test_mode_ac_strain_plot_samples_native_values_and_routes_components(monkeypatch):
    class ModeInterpolator:
        def __init__(self):
            self.xy_raw = {
                "m_x": np.array(((0.0, 1.0), (0.0, 1.0))),
                "m_y": np.array(((0.0, 0.0), (1.0, 1.0))),
            }
            self.xy_out = {
                "v_x": np.array((0.0, 1.0)),
                "v_y": np.array((0.0, 1.0)),
            }

        def define_plot_grid_2D(self, n_pts):
            self.n_pts = n_pts

    mode_interpolator = ModeInterpolator()
    mode = ModeAC.__new__(ModeAC)
    mode.mode_num = 0
    mode.mode_interpolator = mode_interpolator
    mode.sim_result = SimpleNamespace(get_xyshift=lambda: (0.0, 0.0))
    captured = {}

    monkeypatch.setattr(
        modes_module.strain_sampling,
        "calculate_mode_strain",
        lambda *_: _interface_strain_samples(),
    )
    monkeypatch.setattr(
        modes_module.strain_sampling,
        "make_native_strain_interpolator",
        lambda _: lambda points: np.tile(
            np.arange(6, dtype=np.complex128) + 1j, (len(points), 1)
        ),
    )
    monkeypatch.setattr(
        modes_module.strain_rendering,
        "render_strain_2d",
        lambda **kwargs: captured.update(kwargs),
    )

    mode.plot_strain(comps="xx,yz", n_pts=17, trim=(0.1, 0.2, 0.3, 0.4))

    assert mode_interpolator.n_pts == 17
    assert captured["components"] == (StrainComponent.XX, StrainComponent.YZ)
    assert captured["fields"]["Sxxr"].shape == (2, 2)
    np.testing.assert_allclose(captured["fields"]["Sxzi"], np.ones((2, 2)))
    assert captured["plot_params"]["trim"] == (0.1, 0.2, 0.3, 0.4)


def test_ac_simresult_strain_plot_routes_shared_plot_parameters():
    calls = []

    class Mode:
        def plot_strain(self, **kwargs):
            calls.append(kwargs)

    sim_result = SimpleNamespace(n_modes=2, get_mode=lambda _: Mode())

    ACSimResult.plot_strain_modes(
        sim_result,
        mode_indices=(1,),
        n_points=29,
        comps=(StrainComponent.XY,),
        prefix="strain-",
        suffix="test",
        trim=(0.1, 0.2, 0.3, 0.4),
        colorbar=False,
    )

    assert len(calls) == 1
    assert calls[0]["comps"] == (StrainComponent.XY,)
    assert calls[0]["n_pts"] == 29
    assert calls[0]["plot_params"]["prefix"] == "strain-"
    assert calls[0]["plot_params"]["suffix"] == "test"
    assert calls[0]["plot_params"]["trim"] == (0.1, 0.2, 0.3, 0.4)
    assert not calls[0]["plot_params"]["colorbar"]
