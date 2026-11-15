"""Regression tests for field selection and plot representation helpers."""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

pytestmark = pytest.mark.filterwarnings("ignore:FieldTag is deprecated")

sys.path.append(str(Path(__file__).resolve().parents[1] / "backend"))

from numbat.field_selection import (
    complex_storage_keys,
    field_symbol,
    resolve_field_type,
    selection_from_storage_key,
    storage_key,
)
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
from numbat.plotting.mode_rendering import render_mode_1d, render_mode_2d
from numbat.plotting.modes import (
    add_contour_plot,
    field_display_name,
    field_type_to_intensity_code,
    field_type_to_vector_code,
    line_color,
    line_style,
    tex_label,
)
from numbat.solvers.modes import (
    Mode,
    ModeInterpolator,
    check_allowed_plot_components,
    normalize_plot_components,
    selections_for_1d_plot,
    split_plot_components,
)
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
    assert check_allowed_plot_components((component,)) == (
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
        Mode.plot_mode_1D(object(), "x", 0, comps="e")

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
