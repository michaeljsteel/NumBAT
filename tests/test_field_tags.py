"""Regression tests for field selection and plot representation helpers."""

import sys
from pathlib import Path

import pytest

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
)
from numbat.plotting.modes import (
    field_display_name,
    field_type_to_intensity_code,
    field_type_to_vector_code,
    line_color,
    line_style,
    tex_label,
)
from numbat.solvers.modes import check_allowed_plot_components, selections_for_1d_plot


@pytest.mark.parametrize("field_type", list(FieldType))
def test_legacy_field_type_inputs_remain_supported(field_type):
    assert FieldType.from_str(field_type.value) is field_type
    with pytest.deprecated_call(match="FieldCode is deprecated"):
        assert FieldCode(field_type.value).as_field_type() is field_type
    assert FieldTag.make_from_field(field_type.value).as_field_type() is field_type


@pytest.mark.parametrize("component", ("x", "y", "z", "t", "a"))
def test_legacy_plot_component_inputs_remain_supported(component):
    check_allowed_plot_components((component,))
    assert FieldTag.make_from_field_and_component(FieldType.EM_E, component)


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
        FieldType.EM_E, ("x", "y", "a", "t"), include_secondary=True
    )

    assert [storage_key(selection) for selection in selections] == [
        "Fxr",
        "Fyr",
        "Fa",
        "Ft",
        "Fxi",
        "Fyi",
    ]
