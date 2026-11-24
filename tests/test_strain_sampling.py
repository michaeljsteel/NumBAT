"""Analytical regression tests for native acoustic P2 strain sampling."""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from matplotlib.figure import Figure

sys.path.append(str(Path(__file__).resolve().parents[1] / "backend"))

import numbat.solvers.modes as modes_module
from numbat.io import strain_fields
from numbat.plotting import strain_rendering
from numbat.plotting.modes import PlotParams2D
from numbat.plotting.strain_native import plot_native_strain_component
from numbat.solvers.modes import Mode, ModeAC, ModeEM
from numbat.solvers.simresult import ACSimResult
from numbat.solvers.strain_sampling import (
    P2_REFERENCE_NODES,
    NativeStrainSamples,
    calculate_element_strain,
    calculate_mode_strain,
    make_native_strain_interpolator,
    sample_mode_strain_cut,
    sample_native_strain,
    strain_cut_points,
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


def test_strain_api_is_available_only_on_acoustic_modes():
    assert not hasattr(Mode, "plot_strain")
    assert not hasattr(ModeEM, "plot_strain")
    assert hasattr(ModeAC, "plot_strain")
    assert not hasattr(modes_module.ModeInterpolator, "plot_strain_mode_i")


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
        element_node_connectivity=np.arange(1, 7, dtype=np.int64)[:, np.newaxis],
        node_coordinates=coordinates.T,
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


def test_strain_renderer_uses_megastrain_scale_and_compact_colorbar():
    figure = Figure()
    axis = figure.subplots()
    plot_params = PlotParams2D()
    plot_params.update({"ticks": False, "title": False, "frame": False})

    image = strain_rendering._plot_component(
        axis,
        {"v_x": np.array((0.0, 1.0)), "v_y": np.array((0.0, 1.0))},
        {"Sxxr": np.full((2, 2), 2e6)},
        StrainSelection.primary(StrainComponent.XX),
        plot_params,
    )
    colorbar = strain_rendering._add_colorbar(figure, axis, image)
    strain_rendering._add_mode_label(figure, 7)
    figure.canvas.draw()

    np.testing.assert_allclose(image.get_array(), np.full((2, 2), 2.0))
    assert colorbar.ax.get_ylabel() == ""
    assert colorbar.ax.get_position().height == pytest.approx(
        0.8 * axis.get_position().height
    )
    assert [(text.get_position(), text.get_text()) for text in figure.texts] == [
        ((0.01, 0.99), "Mode 7\nStrain [$\\mathrm{M}\\epsilon$]")
    ]


def test_fast_native_interpolator_matches_reference_sampler_away_from_interface():
    samples = _interface_strain_samples()
    points = np.array(((0.1, 0.1), (0.9, 0.9), (2.0, 2.0)))

    np.testing.assert_allclose(
        make_native_strain_interpolator(samples)(points),
        sample_native_strain(samples, points),
        equal_nan=True,
    )


def test_strain_cut_points_match_standard_x_y_and_line_axes(monkeypatch):
    monkeypatch.setattr("numbat.solvers.strain_sampling.SI_um", 1.0)
    sim_result = SimpleNamespace(
        fem_mesh=SimpleNamespace(get_xy_limits=lambda: (0.0, 2.0, 0.0, 4.0)),
    )

    x_axis, x_points = strain_cut_points(sim_result, "x", 1.0, None, 3)
    y_axis, y_points = strain_cut_points(sim_result, "y", 2.0, None, 3)
    line_axis, line_points = strain_cut_points(
        sim_result, "line", (0.0, 0.0), (3.0, 4.0), 2
    )

    np.testing.assert_allclose(x_axis, (0.0, 2.0, 4.0))
    np.testing.assert_allclose(x_points, ((1.0, 0.0), (1.0, 2.0), (1.0, 4.0)))
    np.testing.assert_allclose(y_axis, (0.0, 1.0, 2.0))
    np.testing.assert_allclose(y_points, ((0.0, 2.0), (1.0, 2.0), (2.0, 2.0)))
    np.testing.assert_allclose(line_axis, (0.0, 5.0))
    np.testing.assert_allclose(line_points, ((0.0, 0.0), (3.0, 4.0)))


def test_native_strain_cut_preserves_first_interface_side(monkeypatch):
    monkeypatch.setattr("numbat.solvers.strain_sampling.SI_um", 1.0)
    samples = _interface_strain_samples()
    sim_result = SimpleNamespace(
        fem_mesh=SimpleNamespace(get_xy_limits=lambda: (0.0, 1.0, 0.0, 1.0)),
    )
    monkeypatch.setattr(
        "numbat.solvers.strain_sampling.calculate_mode_strain",
        lambda *_: samples,
    )

    cut_axis, values = sample_mode_strain_cut(
        sim_result, 0, "x", 0.5, n_points=3
    )

    np.testing.assert_allclose(cut_axis, (0.0, 0.5, 1.0))
    np.testing.assert_allclose(
        values,
        sample_native_strain(samples, np.column_stack((np.full(3, 0.5), cut_axis))),
    )


@pytest.mark.parametrize(
    ("s_cut", "val1", "val2", "expected_points"),
    [
        ("x", 0.25, None, ((0.25, 0.0), (0.25, 0.5), (0.25, 1.0))),
        ("y", 0.75, None, ((0.0, 0.75), (0.5, 0.75), (1.0, 0.75))),
        ("line", (0.0, 0.0), (1.0, 1.0), ((0.0, 0.0), (0.5, 0.5), (1.0, 1.0))),
    ],
)
def test_native_strain_cuts_match_affine_values(
    monkeypatch, s_cut, val1, val2, expected_points
):
    monkeypatch.setattr("numbat.solvers.strain_sampling.SI_um", 1.0)
    samples = _interface_strain_samples()
    coordinates = samples.coordinates
    for component_index in range(6):
        values = (
            coordinates[..., 0]
            + 2 * coordinates[..., 1]
            + component_index
            + 1j * (coordinates[..., 0] - coordinates[..., 1])
        )
        samples.values[..., component_index] = values
    sim_result = SimpleNamespace(
        fem_mesh=SimpleNamespace(get_xy_limits=lambda: (0.0, 1.0, 0.0, 1.0)),
    )
    monkeypatch.setattr(
        "numbat.solvers.strain_sampling.calculate_mode_strain",
        lambda *_: samples,
    )

    _, values = sample_mode_strain_cut(
        sim_result, 0, s_cut, val1, val2, n_points=3
    )

    points = np.asarray(expected_points)
    expected = np.column_stack([
        points[:, 0] + 2 * points[:, 1] + component_index
        + 1j * (points[:, 0] - points[:, 1])
        for component_index in range(6)
    ])
    np.testing.assert_allclose(values, expected)


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
    mode.sim_result = SimpleNamespace()
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


def test_ac_simresult_strain_plot_routes_shared_plot_parameters(monkeypatch):
    calls = []
    progress_calls = []

    class Mode:
        def plot_strain(self, **kwargs):
            calls.append(kwargs)

    sim_result = SimpleNamespace(n_modes=2, get_mode=lambda _: Mode())
    monkeypatch.setattr(
        "numbat.solvers.simresult.progressBar",
        lambda values, **kwargs: (
            progress_calls.append((tuple(values), kwargs)) or values
        ),
    )

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
    assert progress_calls == [((1,), {"prefix": "  Progress:", "length": 20})]


def test_ac_simresult_strain_plot_accepts_all_mode_convenience_request():
    calls = []

    class Mode:
        def plot_strain(self, **kwargs):
            calls.append(kwargs)

    sim_result = SimpleNamespace(n_modes=2, get_mode=lambda _: Mode())

    ACSimResult.plot_strain_modes(sim_result, mode_indices="All", n_points=29)

    assert len(calls) == 2
    assert all(call["n_pts"] == 29 for call in calls)


def test_strain_2d_export_writes_all_physical_parts(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "numbat.io.strain_fields.outdir_fields_path",
        lambda prefix: tmp_path,
    )
    fields = {}
    for index, component in enumerate(StrainComponent):
        real_key, imaginary_key = complex_storage_keys(component)
        fields[real_key] = np.full((2, 3), index + 1.0)
        fields[imaginary_key] = np.full((2, 3), -(index + 1.0))

    strain_fields.write_strain_fields_2d(
        prefix="ignored",
        mode_index=2,
        coordinates={"v_x": np.array((0.0, 1.0)), "v_y": np.array((0.0, 2.0, 4.0))},
        fields=fields,
    )

    expected_files = {
        f"AC_strain_mode_02_S{component.value}_{part}.txt"
        for component in StrainComponent
        for part in ("re", "im")
    }
    assert {path.name for path in tmp_path.iterdir()} == expected_files
    for component in StrainComponent:
        real_key, imaginary_key = complex_storage_keys(component)
        for suffix, key in (("re", real_key), ("im", imaginary_key)):
            output = tmp_path / f"AC_strain_mode_02_S{component.value}_{suffix}.txt"
            assert f"Mode 02, {key}" in output.read_text()
            np.testing.assert_allclose(np.loadtxt(output), fields[key])


def test_mode_ac_strain_export_uses_shared_sampled_grid(monkeypatch):
    coordinates = {"v_x": np.array((0.0, 1.0)), "v_y": np.array((0.0, 1.0))}
    fields = {"Sxxr": np.ones((2, 2))}
    mode = ModeAC.__new__(ModeAC)
    mode.mode_num = 3
    mode._sample_strain_2d = lambda n_pts: (coordinates, fields)
    captured = {}
    monkeypatch.setattr(
        modes_module.strain_fields,
        "write_strain_fields_2d",
        lambda **kwargs: captured.update(kwargs),
    )

    mode.write_strain(prefix="strain-", n_pts=17)

    assert captured["prefix"] == "strain-"
    assert captured["mode_index"] == 3
    assert captured["coordinates"] is coordinates
    assert captured["fields"] is fields


def test_ac_simresult_strain_export_routes_selected_modes():
    calls = []

    class Mode:
        def write_strain(self, **kwargs):
            calls.append(kwargs)

    sim_result = SimpleNamespace(n_modes=2, get_mode=lambda _: Mode())

    ACSimResult.write_strain_modes(
        sim_result,
        mode_indices=(1,),
        n_points=29,
        prefix="strain-",
    )

    assert calls == [{"prefix": "strain-", "n_pts": 29}]


def test_strain_1d_export_uses_documented_twelve_column_order(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "numbat.io.strain_fields.outdir_fields_path",
        lambda prefix: tmp_path,
    )
    fields = {}
    for index, component in enumerate(StrainComponent):
        real_key, imaginary_key = complex_storage_keys(component)
        fields[real_key] = np.full(2, index + 1.0)
        fields[imaginary_key] = np.full(2, -(index + 1.0))

    strain_fields.write_strain_fields_1d(
        prefix="ignored",
        mode_index=2,
        cut="line",
        x_axis=np.array((0.0, 1.0)),
        fields=fields,
    )

    output = tmp_path / "AC_strain_mode_02_linecut.txt"
    expected_keys = [
        key
        for component in StrainComponent
        for key in complex_storage_keys(component)
    ]
    assert output.read_text().splitlines()[0].endswith("x " + " ".join(expected_keys))
    data = np.loadtxt(output)
    assert data.shape == (2, 13)
    np.testing.assert_allclose(data[:, 1], fields["Sxxr"])
    np.testing.assert_allclose(data[:, -1], fields["Sxyi"])


def test_mode_ac_strain_1d_plot_routes_native_cut_fields(monkeypatch):
    mode = ModeAC.__new__(ModeAC)
    mode.mode_num = 3
    x_axis = np.array((0.0, 1.0))
    fields = {"Sxyr": np.ones(2)}
    mode._sample_strain_1d = lambda *args: (x_axis, fields)
    captured = {}
    monkeypatch.setattr(
        modes_module.strain_rendering,
        "render_strain_1d",
        lambda **kwargs: captured.update(kwargs),
    )

    mode.plot_strain_1D("y", 0.25, comps="xy", n_pts=17)

    assert captured["x_axis"] is x_axis
    assert captured["fields"] is fields
    assert captured["components"] == (StrainComponent.XY,)
    assert captured["x_label"] == "$x$ (μm)"
    assert captured["cut"] == "y"


def test_strain_1d_renderer_uses_megastrain_scale(monkeypatch):
    class Axis:
        def __init__(self):
            self.curves = []

        def plot(self, x_values, y_values, **kwargs):
            self.curves.append((x_values, y_values, kwargs))

        def set_xlabel(self, label):
            self.x_label = label

        def set_ylabel(self, label):
            self.y_label = label

        def legend(self):
            self.has_legend = True

    axis = Axis()
    monkeypatch.setattr(
        strain_rendering.plt,
        "subplots",
        lambda: (object(), axis),
    )
    monkeypatch.setattr(
        strain_rendering,
        "_strain_cut_filename",
        lambda *args: "strain-cut.png",
    )
    monkeypatch.setattr(
        strain_rendering.plottools,
        "save_and_close_figure",
        lambda *args: None,
    )

    strain_rendering.render_strain_1d(
        x_axis=np.array((0.0, 1.0)),
        fields={"Sxxr": np.array((1e6, 2e6))},
        components=(StrainComponent.XX,),
        x_label="$x$ (μm)",
        cut="y",
        mode_index=0,
        plot_params={"prefix": "", "suffix": ""},
    )

    np.testing.assert_allclose(axis.curves[0][1], (1.0, 2.0))
    assert axis.y_label == r"Strain [$\mathrm{M}\epsilon$]"
    assert axis.has_legend


def test_ac_simresult_strain_1d_facades_route_selected_modes(monkeypatch):
    plot_calls = []
    write_calls = []
    progress_calls = []

    class Mode:
        def plot_strain_1D(self, *args, **kwargs):
            plot_calls.append((args, kwargs))

        def write_strain_1D(self, *args, **kwargs):
            write_calls.append((args, kwargs))

    sim_result = SimpleNamespace(n_modes=2, get_mode=lambda _: Mode())
    monkeypatch.setattr(
        "numbat.solvers.simresult.progressBar",
        lambda values, **kwargs: (
            progress_calls.append((tuple(values), kwargs)) or values
        ),
    )

    ACSimResult.plot_strain_modes_1D(
        sim_result,
        "line",
        (0.0, 0.0),
        (1.0, 1.0),
        mode_indices=(1,),
        n_points=29,
        comps=(StrainComponent.XY,),
        prefix="strain-",
        suffix="test",
    )
    ACSimResult.write_strain_modes_1D(
        sim_result,
        "x",
        0.25,
        mode_indices=(1,),
        n_points=31,
        prefix="strain-",
    )

    assert len(plot_calls) == 1
    assert plot_calls[0][0] == ("line", (0.0, 0.0), (1.0, 1.0))
    assert plot_calls[0][1]["comps"] == (StrainComponent.XY,)
    assert plot_calls[0][1]["n_pts"] == 29
    assert plot_calls[0][1]["plot_params"]["prefix"] == "strain-"
    assert plot_calls[0][1]["plot_params"]["suffix"] == "test"
    assert progress_calls == [((1,), {"prefix": "  Progress:", "length": 20})]
    assert write_calls == [
        (("x", 0.25, None), {"n_pts": 31, "prefix": "strain-"})
    ]
