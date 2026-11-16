"""Regression tests for field selection and plot representation helpers."""

import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from matplotlib.figure import Figure

pytestmark = pytest.mark.filterwarnings("ignore:FieldTag is deprecated")

sys.path.append(str(Path(__file__).resolve().parents[1] / "backend"))

from numbat.field_selection import (
    complex_storage_keys,
    field_symbol,
    resolve_field_type,
    selection_from_storage_key,
    storage_key,
)
from numbat.io import mode_fields
from numbat.nbtypes import (
    ComplexPart,
    FieldCode,
    FieldComponent,
    FieldSelection,
    FieldTag,
    FieldType,
    SI_vacuum_impedance_Z0,
    parse_field_components,
)
from numbat.plotting.mode_rendering import (
    components_for_1d_plot,
    normalize_plot_components,
    render_mode_1d,
    render_mode_2d,
    selections_for_1d_plot,
    split_plot_components,
)
from numbat.plotting.modes import (
    PlotParams2D,
    add_contour_plot,
    configure_overview_subplots,
    field_display_name,
    field_type_to_intensity_code,
    field_type_to_vector_code,
    line_color,
    line_style,
    offset_mode_data_panel,
    overview_figure_size,
    resolve_plot_trim,
    select_overview_layout,
    tex_label,
    write_mode_data,
)


def test_plot_trim_tuple_overrides_legacy_limits_and_validates_shape():
    plot_params = PlotParams2D()
    plot_params.update({
        "xlim_min": 0.1,
        "xlim_max": 0.2,
        "ylim_min": 0.3,
        "ylim_max": 0.4,
    })

    with pytest.deprecated_call(match="xlim_min.*deprecated"):
        assert resolve_plot_trim(plot_params, -0.05) == (0.1, 0.3, 0.2, 0.4)

    plot_params["trim"] = (0.5, 0.6, 0.7, 0.8)
    with pytest.deprecated_call(match="xlim_min.*deprecated"):
        assert resolve_plot_trim(plot_params, -0.05) == (0.5, 0.6, 0.7, 0.8)

    plot_params["trim"] = (0.1, 0.2, 0.3)
    with pytest.deprecated_call(match="xlim_min.*deprecated"):
        with pytest.raises(ValueError, match="four-item"):
            resolve_plot_trim(plot_params, -0.05)


def test_plot_trim_tuple_does_not_warn_without_legacy_limits():
    plot_params = PlotParams2D()
    plot_params["trim"] = (0.1, 0.2, 0.3, 0.4)

    with warnings.catch_warnings(record=True) as caught_warnings:
        warnings.simplefilter("always")
        assert resolve_plot_trim(plot_params, -0.05) == (0.1, 0.2, 0.3, 0.4)

    assert caught_warnings == []
from numbat.solvers.mode_sampling import ModeInterpolator as SamplingModeInterpolator
from numbat.solvers.modes import Mode as BaseMode
from numbat.solvers.modes import ModeInterpolator
from numbat.solvers.simresult import SimResult


@pytest.mark.parametrize("field_type", list(FieldType))
def test_legacy_field_type_inputs_remain_supported(field_type):
    assert FieldType.from_str(field_type.value) is field_type
    with pytest.deprecated_call(match="FieldCode is deprecated"):
        assert FieldCode(field_type.value).as_field_type() is field_type
    with pytest.deprecated_call(match="FieldTag is deprecated"):
        assert FieldTag.make_from_field(field_type.value).as_field_type() is field_type


def test_field_tag_legacy_construction_warns_and_remains_functional():
    with pytest.deprecated_call(match="FieldTag is deprecated"):
        tag = FieldTag("Ex")

    assert tag.as_selection() == FieldSelection.primary(FieldType.EM_E, FieldComponent.X)


@pytest.mark.parametrize("component", ("x", "y", "z", "t", "a"))
def test_legacy_plot_component_inputs_remain_supported(component):
    assert normalize_plot_components((component,)) == (
        FieldComponent.parse(component),
    )
    assert FieldTag.make_from_field_and_component(FieldType.EM_E, component)


def test_plot_component_strings_and_enums_normalize_identically():
    compact = normalize_plot_components("xyat")
    typed = normalize_plot_components(
        (
            FieldComponent.X,
            FieldComponent.Y,
            FieldComponent.MAGNITUDE,
            FieldComponent.TRANSVERSE,
        )
    )

    assert compact == typed
    assert [
        storage_key(selection)
        for selection in selections_for_1d_plot(FieldType.EM_E, compact, include_secondary=False)
    ] == [
        storage_key(selection)
        for selection in selections_for_1d_plot(FieldType.EM_E, typed, include_secondary=False)
    ]


def test_field_component_spec_parses_compact_and_mixed_iterable_forms():
    assert parse_field_components("exyat") == (
        FieldComponent.EVERYTHING,
        FieldComponent.X,
        FieldComponent.Y,
        FieldComponent.MAGNITUDE,
        FieldComponent.TRANSVERSE,
    )
    assert parse_field_components(("x", FieldComponent.Y, "a")) == (
        FieldComponent.X,
        FieldComponent.Y,
        FieldComponent.MAGNITUDE,
    )


def test_plot_component_split_defaults_to_the_overview_only():
    assert split_plot_components(()) == (True, ())


def test_plot_component_split_separates_overview_from_individual_components():
    assert split_plot_components("exyt") == (
        True,
        (FieldComponent.X, FieldComponent.Y, FieldComponent.TRANSVERSE),
    )
    assert split_plot_components("xyt") == (
        False,
        (FieldComponent.X, FieldComponent.Y, FieldComponent.TRANSVERSE),
    )


def test_auto_overview_layout_matches_six_panel_cell_aspect_to_waveguide():
    params = PlotParams2D()
    wide_waveguide = {
        "v_x": np.array([0.0, 5.0]),
        "v_y": np.array([0.0, 1.0]),
    }
    tall_waveguide = {
        "v_x": np.array([0.0, 0.2]),
        "v_y": np.array([0.0, 1.0]),
    }

    assert (
        select_overview_layout(
            FieldType.EM_E, wide_waveguide, params, hide_secondary_comps=True
        )
        == "wide"
    )
    assert (
        select_overview_layout(
            FieldType.EM_E, tall_waveguide, params, hide_secondary_comps=True
        )
        == "tall"
    )


def test_overview_layout_honors_legacy_and_explicit_overrides():
    params = PlotParams2D()
    d_xy = {"v_x": np.array([0.0, 5.0]), "v_y": np.array([0.0, 1.0])}

    params["overview_layout"] = "tall"
    assert (
        select_overview_layout(
            FieldType.EM_E, d_xy, params, hide_secondary_comps=True
        )
        == "tall"
    )

    params["wide_mode"] = True
    assert (
        select_overview_layout(
            FieldType.EM_E, d_xy, params, hide_secondary_comps=True
        )
        == "wide"
    )

    params["wide_mode"] = None
    params["overview_layout"] = "auto"
    assert (
        select_overview_layout(
            FieldType.EM_E, d_xy, params, hide_secondary_comps=False
        )
        == "tall"
    )


def test_auto_overview_layout_requires_a_waveguide_aspect_above_two_for_wide():
    params = PlotParams2D()
    square_to_wide = {
        "v_x": np.array([0.0, 2.0]),
        "v_y": np.array([0.0, 1.0]),
    }
    clearly_wide = {
        "v_x": np.array([0.0, 2.1]),
        "v_y": np.array([0.0, 1.0]),
    }

    assert (
        select_overview_layout(
            FieldType.EM_E, square_to_wide, params, hide_secondary_comps=True
        )
        == "tall"
    )
    assert (
        select_overview_layout(
            FieldType.EM_E, clearly_wide, params, hide_secondary_comps=True
        )
        == "wide"
    )


def test_overview_figure_size_reserves_colorbar_space_with_bounded_area():
    params = PlotParams2D()
    d_xy = {"v_x": np.array([0.0, 0.6]), "v_y": np.array([0.0, 1.0])}

    with_colorbar = overview_figure_size(
        FieldType.EM_E, d_xy, params, hide_secondary_comps=True, layout="tall"
    )
    params["colorbar"] = False
    without_colorbar = overview_figure_size(
        FieldType.EM_E, d_xy, params, hide_secondary_comps=True, layout="tall"
    )

    base_width, base_height = params["decorator"]["figsize"]
    assert with_colorbar[0] > without_colorbar[0]
    assert with_colorbar[0] * with_colorbar[1] == pytest.approx(
        base_width * base_height
    )
    assert with_colorbar[0] / with_colorbar[1] <= (base_width / base_height) * 1.5


def test_overview_figure_size_rejects_negative_colorbar_allowance():
    params = PlotParams2D()
    params["colorbar_allowance"] = -0.1
    d_xy = {"v_x": np.array([0.0, 1.0]), "v_y": np.array([0.0, 1.0])}

    with pytest.raises(ValueError, match="colorbar_allowance"):
        overview_figure_size(
            FieldType.EM_E, d_xy, params, hide_secondary_comps=True, layout="tall"
        )


def test_tall_square_overviews_keep_the_original_figure_area():
    params = PlotParams2D()
    square_waveguide = {"v_x": np.array([0.0, 1.25]), "v_y": np.array([0.0, 1.0])}

    square_size = overview_figure_size(
        FieldType.EM_E,
        square_waveguide,
        params,
        hide_secondary_comps=True,
        layout="tall",
    )

    base_width, base_height = params["decorator"]["figsize"]
    assert square_size[0] * square_size[1] == pytest.approx(base_width * base_height)


def test_tall_square_overviews_enlarge_panels_by_reducing_whitespace():
    fig = Figure(figsize=(10, 8))
    fig.subplots(2, 3)
    original = fig.subplotpars
    original_axes_height = (original.top - original.bottom) / (2 + 0.2)

    configure_overview_subplots(
        fig, "tall", waveguide_aspect=1.0, rows=2, columns=3,
        hspace=0.2, wspace=1.2,
    )

    adjusted = fig.subplotpars
    adjusted_axes_height = (adjusted.top - adjusted.bottom) / (2 + adjusted.hspace)
    assert adjusted.wspace == pytest.approx(0.75)
    assert adjusted_axes_height == pytest.approx(original_axes_height * 1.2)


def test_mode_data_panel_is_offset_northwest_by_its_physical_diagonal():
    fig = Figure(figsize=(10, 8))
    ax = fig.add_subplot()
    ax.set_position((0.3, 0.3, 0.3, 0.3))
    original = ax.get_position()
    offset_mode_data_panel(fig, ax)
    moved = ax.get_position()

    panel_diagonal = np.hypot(original.width * 10, original.height * 8)
    diagonal_component = 0.1 * panel_diagonal / np.sqrt(2)

    assert moved.x0 == pytest.approx(original.x0 - diagonal_component / 10)
    assert moved.y0 == pytest.approx(original.y0 + diagonal_component / 8)
    assert moved.width == pytest.approx(original.width)
    assert moved.height == pytest.approx(original.height)


def test_wide_mode_data_uses_a_single_sequential_column():
    class Decorator:
        def get_property(self, _):
            return 10

    class Mode:
        def analyse_mode(self):
            pass

        def field_fracs(self):
            return (0.1, 0.2, 0.3, 0.4)

        def center_of_mass(self):
            return (0.0, 0.0)

        def second_moment_widths(self):
            return (1.0, 1.0, 0.0)

        def get_mode_data(self):
            return {}

    class SimResult:
        omega_EM = 1.0

        def get_all_modes(self):
            return (Mode(),)

        def is_EM(self):
            return True

        def kz_EM(self, _):
            return 1.0

        def neff(self, _):
            return 1.0

        def ngroup_EM_available(self):
            return False

        def symmetry_classification(self, _):
            return ()

    fig = Figure(figsize=(10, 8))
    ax = fig.add_subplot()
    write_mode_data(
        ax,
        {"decorator": Decorator()},
        SimResult(),
        mode_index=0,
        field_type=FieldType.EM_E,
        layout="wide",
    )

    text_positions = [text.get_position() for text in ax.texts]
    assert all(x == pytest.approx(0.05) for x, _ in text_positions[1:])
    assert [y for _, y in text_positions] == pytest.approx(
        [0.89, 0.79, 0.69, 0.59, 0.49, 0.39, 0.29, 0.19]
    )


def test_2d_plot_component_selection_controls_overview_and_individual_outputs(
    monkeypatch,
):
    class Decorator:
        def __init__(self):
            self.calls = []

        def set_for_multi(self):
            self.calls.append("multi")

        def set_for_single(self):
            self.calls.append("single")

    calls = []
    monkeypatch.setattr(
        "numbat.plotting.modes.plot_all_components",
        lambda *args: calls.append(("overview", args)),
    )
    monkeypatch.setattr(
        "numbat.plotting.modes.plot_one_component",
        lambda *args: calls.append(("component", args[5])),
    )
    decorator = Decorator()
    plot_params = {"decorator": decorator}

    render_mode_2d(
        field_type=FieldType.EM_E,
        coordinates=None,
        fields={},
        include_overview=True,
        components=(FieldComponent.X, FieldComponent.Y),
        simulation_result=None,
        mode_index=0,
        plot_params=plot_params,
    )

    assert calls[0][0] == "overview"
    assert [call[1].component for call in calls[1:]] == [
        FieldComponent.X,
        FieldComponent.Y,
    ]
    assert decorator.calls == ["multi", "single"]

    calls.clear()
    decorator.calls.clear()
    render_mode_2d(
        field_type=FieldType.EM_E,
        coordinates=None,
        fields={},
        include_overview=False,
        components=(FieldComponent.X,),
        simulation_result=None,
        mode_index=0,
        plot_params=plot_params,
    )

    assert [call[0] for call in calls] == ["component"]
    assert decorator.calls == ["single"]


def test_overview_component_is_rejected_for_1d_and_supplied_axes():
    with pytest.raises(SystemExit):
        components_for_1d_plot("e")

    with pytest.raises(SystemExit):
        render_mode_2d(
            field_type=FieldType.EM_E,
            coordinates=None,
            fields={},
            include_overview=True,
            components=(FieldComponent.X,),
            simulation_result=None,
            mode_index=0,
            plot_params={},
            axis=object(),
        )


@pytest.mark.parametrize("spec", (FieldComponent.X, 42))
def test_field_component_spec_rejects_non_iterable_values(spec):
    with pytest.raises(ValueError, match="compact string or an iterable"):
        parse_field_components(spec)


def test_plot_component_intensity_is_supported():
    assert normalize_plot_components("i") == (FieldComponent.INTENSITY,)


def test_1d_renderer_writes_requested_component_curves(monkeypatch):
    class Axis:
        def __init__(self):
            self.curves = []

        def plot(self, x_values, y_values, **kwargs):
            self.curves.append((x_values, y_values, kwargs))

        def set_xscale(self, scale):
            self.xscale = scale

        def set_yscale(self, scale):
            self.yscale = scale

        def set_xlabel(self, label):
            self.xlabel = label

        def set_ylabel(self, label):
            self.ylabel = label

        def legend(self):
            self.has_legend = True

    axis = Axis()
    saved_filenames = []
    monkeypatch.setattr(
        "numbat.plotting.mode_rendering.plt.subplots",
        lambda: (object(), axis),
    )
    monkeypatch.setattr(
        "numbat.plotting.mode_rendering.mode_presentation.modeplot_filename_1D",
        lambda *args, **kwargs: "mode-cut.png",
    )
    monkeypatch.setattr(
        "numbat.plotting.mode_rendering.plottools.save_and_close_figure",
        lambda figure, filename: saved_filenames.append(filename),
    )

    render_mode_1d(
        field_type=FieldType.EM_E,
        fields={"Fxr": np.array([1.0, 2.0]), "Fi": np.array([4.0, 9.0])},
        x_axis=np.array([0.0, 1.0]),
        components=(FieldComponent.X, FieldComponent.INTENSITY),
        x_label="$y$ (μm)",
        cut="x",
        mode_index=0,
        plot_params={
            "logx": False,
            "logy": False,
            "suppress_imimre": True,
            "prefix": "",
            "suffix": "",
        },
    )

    assert [curve[1].tolist() for curve in axis.curves] == [[1.0, 2.0], [4.0, 9.0]]
    assert axis.xlabel == "$y$ (μm)"
    assert axis.ylabel == "fields"
    assert axis.has_legend
    assert saved_filenames == ["mode-cut.png"]


def test_plot_modes_normalizes_compact_and_typed_component_requests(monkeypatch):
    class Mode:
        def plot_mode(self, comps, field_type, plot_params):
            calls.append((comps, field_type, plot_params["field_type"]))

    calls = []
    mode = Mode()
    sim_result = SimpleNamespace(
        is_AC=lambda: False,
        n_modes=1,
        get_mode=lambda index: mode,
        get_mode_interpolator=lambda: SimpleNamespace(
            define_plot_grid_2D=lambda **kwargs: None
        ),
    )
    monkeypatch.setattr(
        "numbat.solvers.simresult.progressBar",
        lambda values, **kwargs: values,
    )

    SimResult.plot_modes(
        sim_result,
        mode_indices=[0],
        field_type=FieldType.EM_E,
        comps="exi",
    )
    SimResult.plot_modes(
        sim_result,
        mode_indices=[0],
        field_type=FieldType.EM_E,
        comps=(
            FieldComponent.EVERYTHING,
            FieldComponent.X,
            FieldComponent.INTENSITY,
        ),
    )

    expected_components = (
        FieldComponent.EVERYTHING,
        FieldComponent.X,
        FieldComponent.INTENSITY,
    )
    assert [(comps, field_type) for comps, field_type, _ in calls] == [
        (expected_components, FieldType.EM_E),
        (expected_components, FieldType.EM_E),
    ]
    assert all(plot_field_type is FieldType.EM_E for _, _, plot_field_type in calls)


def test_plot_modes_forwards_canonical_trim_tuple(monkeypatch):
    calls = []

    class Mode:
        def plot_mode(self, comps, field_type, plot_params):
            calls.append(plot_params["trim"])

    sim_result = SimpleNamespace(
        is_AC=lambda: False,
        n_modes=1,
        get_mode=lambda index: Mode(),
        get_mode_interpolator=lambda: SimpleNamespace(
            define_plot_grid_2D=lambda **kwargs: None
        ),
    )
    monkeypatch.setattr(
        "numbat.solvers.simresult.progressBar",
        lambda values, **kwargs: values,
    )

    SimResult.plot_modes(
        sim_result,
        mode_indices=[0],
        trim=(0.1, 0.2, 0.3, 0.4),
    )

    assert calls == [(0.1, 0.2, 0.3, 0.4)]


def test_plot_modes_1d_normalizes_compact_and_typed_component_requests(monkeypatch):
    class Mode:
        def plot_mode_1D(self, scut, val1, val2, comps, field_type, plot_params):
            calls.append((scut, val1, val2, comps, field_type, plot_params["field_type"]))

    calls = []
    mode = Mode()
    sim_result = SimpleNamespace(
        is_AC=lambda: False,
        n_modes=1,
        get_mode=lambda index: mode,
        get_mode_interpolator=lambda: SimpleNamespace(
            define_plot_grid_2D=lambda **kwargs: None
        ),
    )
    monkeypatch.setattr(
        "numbat.solvers.simresult.progressBar",
        lambda values, **kwargs: values,
    )

    SimResult.plot_modes_1D(
        sim_result,
        "x",
        0,
        mode_indices=[0],
        field_type=FieldType.EM_E,
        comps="ai",
    )
    SimResult.plot_modes_1D(
        sim_result,
        "x",
        0,
        mode_indices=[0],
        field_type=FieldType.EM_E,
        comps=(FieldComponent.MAGNITUDE, FieldComponent.INTENSITY),
    )

    expected_components = (FieldComponent.MAGNITUDE, FieldComponent.INTENSITY)
    assert [(scut, val1, val2, comps, field_type) for scut, val1, val2, comps, field_type, _ in calls] == [
        ("x", 0, None, expected_components, FieldType.EM_E),
        ("x", 0, None, expected_components, FieldType.EM_E),
    ]
    assert all(plot_field_type is FieldType.EM_E for *_, plot_field_type in calls)


@pytest.mark.parametrize(
    ("field_type", "component", "complex_part"),
    [
        ("EM_E", "x", "primary"),
        (FieldType.EM_H, FieldComponent.Y, ComplexPart.SECONDARY),
        ("AC", "z", ComplexPart.PRIMARY),
    ],
)
def test_field_selection_normalizes_supported_input_forms(field_type, component, complex_part):
    selection = FieldSelection(field_type, component, complex_part)

    assert isinstance(selection.field_type, FieldType)
    assert isinstance(selection.component, FieldComponent)
    assert isinstance(selection.complex_part, ComplexPart)


@pytest.mark.parametrize(
    ("parser", "value", "message"),
    [
        (FieldType.parse, "EM_D", "field_type must be one of"),
        (FieldComponent.parse, "q", "field component must be one of"),
        (ComplexPart.parse, "imaginary", "complex part must be one of"),
    ],
)
def test_field_selection_parsers_reject_unknown_values(parser, value, message):
    with pytest.raises(ValueError, match=message):
        parser(value)


@pytest.mark.parametrize(
    "selection",
    [
        FieldSelection.primary(FieldType.EM_E, FieldComponent.X),
        FieldSelection.secondary(FieldType.EM_H, FieldComponent.Z),
        FieldSelection.magnitude(FieldType.AC),
        FieldSelection.intensity(FieldType.EM_E),
        FieldSelection.transverse(FieldType.EM_H),
    ],
)
def test_field_selection_factory_methods(selection):
    assert isinstance(selection, FieldSelection)


def test_field_selection_properties_and_part_transformation():
    selection = FieldSelection.primary(FieldType.EM_E, FieldComponent.X)
    secondary = selection.with_part(ComplexPart.SECONDARY)

    assert selection.is_cartesian
    assert not selection.is_derived
    assert not selection.is_vector
    assert secondary.complex_part == ComplexPart.SECONDARY
    assert secondary != selection
    assert FieldSelection.transverse(FieldType.EM_E).is_vector
    assert FieldSelection.magnitude(FieldType.EM_E).is_derived


@pytest.mark.parametrize(
    ("component", "complex_part", "message"),
    [
        (FieldComponent.X, None, "Cartesian field components require"),
        (FieldComponent.MAGNITUDE, ComplexPart.PRIMARY, "Derived field quantities must not"),
        (FieldComponent.TRANSVERSE, ComplexPart.SECONDARY, "Derived field quantities must not"),
    ],
)
def test_field_selection_rejects_invalid_component_part_combinations(
    component, complex_part, message
):
    with pytest.raises(ValueError, match=message):
        FieldSelection(FieldType.EM_E, component, complex_part)


def test_field_selection_is_immutable():
    selection = FieldSelection.primary(FieldType.EM_E, FieldComponent.X)

    with pytest.raises(AttributeError):
        selection.component = FieldComponent.Y


@pytest.mark.parametrize(
    ("field_type", "symbol"),
    [(FieldType.EM_E, "E"), (FieldType.EM_H, "H"), (FieldType.AC, "u")],
)
def test_field_symbol_uses_the_expected_display_symbol(field_type, symbol):
    assert field_symbol(field_type) == symbol


@pytest.mark.parametrize(
    ("selection", "expected_key"),
    [
        (FieldSelection.primary(FieldType.EM_E, FieldComponent.X), "Fxr"),
        (FieldSelection.secondary(FieldType.EM_E, FieldComponent.X), "Fxi"),
        (FieldSelection.primary(FieldType.EM_H, FieldComponent.Y), "Fyr"),
        (FieldSelection.secondary(FieldType.EM_H, FieldComponent.Y), "Fyi"),
        (FieldSelection.primary(FieldType.AC, FieldComponent.Z), "Fzi"),
        (FieldSelection.secondary(FieldType.AC, FieldComponent.Z), "Fzr"),
        (FieldSelection.transverse(FieldType.EM_E), "Ft"),
        (FieldSelection.magnitude(FieldType.EM_H), "Fa"),
        (FieldSelection.intensity(FieldType.AC), "Fi"),
    ],
)
def test_field_selection_storage_keys_round_trip(selection, expected_key):
    assert storage_key(selection) == expected_key
    assert selection_from_storage_key(selection.field_type, expected_key) == selection


def test_complex_storage_keys_preserve_real_imaginary_order():
    assert complex_storage_keys(FieldType.EM_E, FieldComponent.X) == ("Fxr", "Fxi")
    assert complex_storage_keys(FieldType.EM_E, FieldComponent.Y) == ("Fyr", "Fyi")
    assert complex_storage_keys(FieldType.EM_E, FieldComponent.Z) == ("Fzr", "Fzi")


def test_interpolated_z_fields_preserve_primary_imaginary_convention():
    class SimulationResult:
        def get_mode_fields_on_fem_mesh(self, mode_index, field_type):
            return (
                np.array([1 + 2j]),
                np.array([3 + 4j]),
                np.array([5 + 6j]),
                np.array([7.0]),
            )

    interpolator = object.__new__(ModeInterpolator)
    interpolator.sim_result = SimulationResult()
    interpolator.interper_f_1D = lambda field: field

    fields = interpolator.interpolate_mode_i(0, FieldType.EM_E, dims=1)

    assert fields["Fzi"].tolist() == [6.0]
    assert fields["Fzr"].tolist() == [5.0]
    assert fields["Fi"].tolist() == [49.0]


def test_mode_interpolator_remains_import_compatible():
    assert ModeInterpolator is SamplingModeInterpolator


def test_mode_sampler_reuses_the_cached_2d_interpolation_grid():
    class FemMesh:
        def __init__(self):
            self.interpolator_calls = []

        def get_xy_limits(self):
            return (0.0, 2.0, 0.0, 1.0)

        def make_interpolator_for_grid(self, x_values, y_values, nx, ny):
            self.interpolator_calls.append((x_values, y_values, nx, ny))
            return lambda field: field

    fem_mesh = FemMesh()
    simulation_result = SimpleNamespace(
        fem_mesh=fem_mesh,
        get_xyshift=lambda: (0.0, 0.0),
    )
    sampler = SamplingModeInterpolator(simulation_result)

    sampler.define_plot_grid_2D(n_pts=10)
    sampler.define_plot_grid_2D(n_pts=10)

    assert len(fem_mesh.interpolator_calls) == 1
    _, _, nx, ny = fem_mesh.interpolator_calls[0]
    assert (nx, ny) == (14, 7)


def test_raw_fem_plotting_does_not_pass_an_unused_mode_sampler(monkeypatch):
    calls = []
    monkeypatch.setattr(
        "numbat.plotting.raw_modes.do_raw_fem_mode_plot",
        lambda *args: calls.append(args),
    )
    fem_mesh = object()
    fem_eigenvectors = object()
    mode = SimpleNamespace(
        sim_result=SimpleNamespace(
            fem_evecs_for_ft=lambda field_type: fem_eigenvectors,
            fem_mesh=fem_mesh,
        ),
        field_type=FieldType.EM_E,
        mode_num=2,
    )

    BaseMode.plot_mode_raw_fem(mode, ("x", "y"))

    assert calls == [(("x", "y"), fem_mesh, fem_eigenvectors, 2)]


def test_2d_field_export_preserves_component_filenames_and_headers(
    monkeypatch, tmp_path
):
    monkeypatch.setattr(
        "numbat.io.mode_fields.runtime.outdir_fields_path",
        lambda prefix: tmp_path,
    )
    fields = {
        "Fxr": np.array([[1.0]]),
        "Fxi": np.array([[2.0]]),
        "Fyr": np.array([[3.0]]),
        "Fyi": np.array([[4.0]]),
        "Fzr": np.array([[5.0]]),
        "Fzi": np.array([[6.0]]),
    }

    mode_fields.write_mode_fields_2d(
        prefix="ignored",
        field_type=FieldType.EM_E,
        mode_index=2,
        coordinates={"v_x": np.array([0.0, 1.0]), "v_y": np.array([0.0, 1.0])},
        fields=fields,
    )

    expected_files = {
        "EM_mode_02_Ex_re.txt",
        "EM_mode_02_Ex_im.txt",
        "EM_mode_02_Ey_re.txt",
        "EM_mode_02_Ey_im.txt",
        "EM_mode_02_Ez_re.txt",
        "EM_mode_02_Ez_im.txt",
    }
    assert {path.name for path in tmp_path.iterdir()} == expected_files
    assert "Mode 02, Fzi_re" in (tmp_path / "EM_mode_02_Ez_re.txt").read_text()
    assert "Mode 02, Fzi_im" in (tmp_path / "EM_mode_02_Ez_im.txt").read_text()


def test_1d_field_export_preserves_the_magnitude_only_derived_schema(
    monkeypatch, tmp_path
):
    monkeypatch.setattr(
        "numbat.io.mode_fields.runtime.outdir_fields_path",
        lambda prefix: tmp_path,
    )
    fields = {
        "Fxr": np.array([1.0, 2.0]),
        "Fxi": np.array([3.0, 4.0]),
        "Fyr": np.array([5.0, 6.0]),
        "Fyi": np.array([7.0, 8.0]),
        "Fzr": np.array([9.0, 10.0]),
        "Fzi": np.array([11.0, 12.0]),
        "Fa": np.array([13.0, 14.0]),
        "Fi": np.array([169.0, 196.0]),
    }

    mode_fields.write_mode_fields_1d(
        prefix="ignored",
        field_type=FieldType.EM_E,
        mode_index=2,
        cut="x",
        x_axis=np.array([0.0, 1.0]),
        fields=fields,
    )

    output = tmp_path / "EM_mode_02_xcut.txt"
    header = output.read_text().splitlines()[0]
    assert header.endswith("Fxr Fxi Fyr Fyi Fzr Fzi Fa")
    assert "Fi" not in header
    assert np.loadtxt(output).shape == (2, 8)


def test_interpolated_magnetic_intensity_uses_scaled_field_magnitude():
    class SimulationResult:
        def get_mode_fields_on_fem_mesh(self, mode_index, field_type):
            return (
                np.array([1 + 2j]),
                np.array([3 + 4j]),
                np.array([5 + 6j]),
                np.array([2.0]),
            )

    interpolator = object.__new__(ModeInterpolator)
    interpolator.sim_result = SimulationResult()
    interpolator.interper_f_1D = lambda field: field

    fields = interpolator.interpolate_mode_i(0, FieldType.EM_H, dims=1)

    assert fields["Fa"].tolist() == [2.0 * SI_vacuum_impedance_Z0]
    assert fields["Fi"].tolist() == [(2.0 * SI_vacuum_impedance_Z0) ** 2]


def test_intensity_contour_uses_interpolated_intensity_without_resquaring(
    monkeypatch,
):
    class Axis:
        def imshow(self, field, **kwargs):
            self.field = field
            return object()

    class Decorator:
        def get_cmap_limits(self, component):
            return None

    class PlotPreferences:
        def cmap_field_signed(self, selection):
            return "signed"

        def cmap_field_unsigned(self, selection):
            return "unsigned"

    monkeypatch.setattr(
        "numbat.plotting.modes.runtime.plot_prefs", lambda: PlotPreferences()
    )
    axis = Axis()
    intensity = np.array([[4.0, 9.0]])
    add_contour_plot(
        None,
        axis,
        {
            "v_x": np.array([0.0, 1.0]),
            "v_y": np.array([0.0]),
            "m_x": np.array([[0.0, 1.0]]),
            "m_y": np.array([[0.0, 0.0]]),
        },
        intensity,
        FieldSelection.intensity(FieldType.EM_E),
        {
            "logamp": False,
            "colorbar": False,
            "contours": False,
            "contour_lst": None,
        },
        Decorator(),
    )

    assert np.array_equal(axis.field, intensity.T)


@pytest.mark.parametrize(
    ("function", "argument", "message"),
    [
        (storage_key, FieldType.EM_E, "selection must be a FieldSelection"),
        (
            lambda value: selection_from_storage_key(FieldType.EM_E, value),
            "Fq",
            "Unknown field-data storage key",
        ),
        (
            lambda value: complex_storage_keys(FieldType.EM_E, value),
            FieldComponent.MAGNITUDE,
            "is not a Cartesian field component",
        ),
    ],
)
def test_field_storage_conversion_rejects_invalid_input(function, argument, message):
    with pytest.raises((TypeError, ValueError), match=message):
        function(argument)


@pytest.mark.parametrize(
    ("requested", "simulation_is_acoustic", "expected"),
    [
        (FieldType.EM_E, False, FieldType.EM_E),
        (FieldType.EM_H, False, FieldType.EM_H),
        (FieldType.AC, True, FieldType.AC),
        (FieldType.EM_E, True, FieldType.AC),
        ("EM_H", True, FieldType.AC),
    ],
)
def test_resolve_field_type_uses_the_simulation_domain(
    requested, simulation_is_acoustic, expected
):
    assert (
        resolve_field_type(requested, simulation_is_acoustic=simulation_is_acoustic)
        is expected
    )


def test_field_code_delegates_to_the_resolver():
    with pytest.deprecated_call(match="FieldCode is deprecated"):
        field_code = FieldCode(FieldType.EM_H, force_AC=True)

    assert field_code.is_AC()
    assert not field_code.is_EM()
    assert field_code.as_field_type() is FieldType.AC


@pytest.mark.parametrize(
    ("selection", "expected_label", "expected_color", "expected_style"),
    [
        (
            FieldSelection.primary(FieldType.EM_E, FieldComponent.X),
            "Re[ $E_x$ ]",
            "blue",
            "dashed",
        ),
        (
            FieldSelection.secondary(FieldType.EM_H, FieldComponent.Z),
            "Re[ $H_z$ ]",
            "brown",
            "dotted",
        ),
        (FieldSelection.magnitude(FieldType.AC), "$|\\vec u|$", "red", "solid"),
        (FieldSelection.intensity(FieldType.EM_E), "$|\\vec E|^2$", "red", "solid"),
        (FieldSelection.transverse(FieldType.EM_H), "$\\vec H_t$", None, "solid"),
    ],
)
def test_field_selection_plot_presentation(
    selection, expected_label, expected_color, expected_style
):
    assert tex_label(selection) == expected_label
    assert line_color(selection) == expected_color
    assert line_style(selection, show_magnitude=True) == expected_style


@pytest.mark.parametrize(
    ("field_type", "expected_name"),
    [
        (FieldType.EM_E, "electric"),
        (FieldType.EM_H, "magnetic"),
        (FieldType.AC, "displacement"),
    ],
)
def test_field_display_names(field_type, expected_name):
    assert field_display_name(field_type) == expected_name


@pytest.mark.parametrize(
    ("field_type", "expected_code"),
    [
        (FieldType.EM_E, "Fi"),
        (FieldType.EM_H, "Fi"),
        (FieldType.AC, "Fi"),
    ],
)
def test_field_type_to_intensity_code_uses_supported_intensity_key(
    field_type, expected_code
):
    selection = field_type_to_intensity_code(field_type)

    assert selection == FieldSelection.intensity(field_type)
    assert storage_key(selection) == expected_code


@pytest.mark.parametrize("field_type", list(FieldType))
def test_field_type_to_vector_code_uses_supported_vector_key(field_type):
    selection = field_type_to_vector_code(field_type)

    assert selection == FieldSelection.transverse(field_type)
    assert storage_key(selection) == "Ft"


def test_field_tag_normalizes_string_field_type_for_storage_key_conversion():
    tag = FieldTag.make_from_field_and_Fcode("EM_E", "Fxr")

    assert tag.field_component() == "Ex"
    assert tag.component_as_F() == "Fxr"


def test_field_tag_secondary_component_and_style_are_consistent():
    primary = FieldTag("Ex")
    secondary = FieldTag("Ex")
    secondary.set_to_minor()

    assert not primary.is_minor()
    assert not primary.is_secondary_component()
    assert secondary.is_minor()
    assert secondary.is_secondary_component()
    assert primary.linestyle(("x", "a")) == "dashed"
    assert secondary.linestyle(("x", "a")) == "dotted"


@pytest.mark.parametrize("component", ("a", "i", "t"))
def test_field_tag_rejects_minor_part_for_derived_quantities(component):
    tag = FieldTag.make_from_field_and_component(FieldType.EM_E, component)

    with pytest.raises(ValueError, match="has no minor"):
        tag.set_to_minor()


def test_1d_plot_selections_keep_primary_components_and_skip_derived_secondaries():
    selections = selections_for_1d_plot(
        FieldType.EM_E, ("x", "y", "a", "i", "t"), include_secondary=True
    )

    assert [storage_key(selection) for selection in selections] == [
        "Fxr",
        "Fyr",
        "Fa",
        "Fi",
        "Ft",
        "Fxi",
        "Fyi",
    ]
