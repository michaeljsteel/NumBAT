"""Conversions between typed field selections and generic field-data keys."""

from numbat.nbtypes import (
    ComplexPart,
    FieldComponent,
    FieldSelection,
    FieldType,
)

_FIELD_SYMBOLS = {
    FieldType.EM_E: "E",
    FieldType.EM_H: "H",
    FieldType.AC: "u",
}

_STORAGE_KEYS = {
    (FieldComponent.X, ComplexPart.PRIMARY): "Fxr",
    (FieldComponent.X, ComplexPart.SECONDARY): "Fxi",
    (FieldComponent.Y, ComplexPart.PRIMARY): "Fyr",
    (FieldComponent.Y, ComplexPart.SECONDARY): "Fyi",
    (FieldComponent.Z, ComplexPart.PRIMARY): "Fzi",
    (FieldComponent.Z, ComplexPart.SECONDARY): "Fzr",
    (FieldComponent.TRANSVERSE, None): "Ft",
    (FieldComponent.MAGNITUDE, None): "Fa",
    (FieldComponent.INTENSITY, None): "Fi",
}

_SELECTION_PARTS = {key: value for value, key in _STORAGE_KEYS.items()}

_REAL_IMAGINARY_PARTS = {
    FieldComponent.X: (ComplexPart.PRIMARY, ComplexPart.SECONDARY),
    FieldComponent.Y: (ComplexPart.PRIMARY, ComplexPart.SECONDARY),
    FieldComponent.Z: (ComplexPart.SECONDARY, ComplexPart.PRIMARY),
}


def field_symbol(field_type):
    """Return the display symbol for a field family.

    Examples:
        >>> field_symbol("EM_E")
        'E'
        >>> field_symbol(FieldType.AC)
        'u'
    """
    return _FIELD_SYMBOLS[FieldType.parse(field_type)]


def resolve_field_type(requested, *, simulation_is_acoustic):
    """Resolve the field family available for a simulation domain.

    Acoustic simulations always expose displacement fields, irrespective of an
    electric or magnetic field request.

    Examples:
        >>> resolve_field_type("EM_H", simulation_is_acoustic=False)
        <FieldType.EM_H: 'EM_H'>
        >>> resolve_field_type(FieldType.EM_E, simulation_is_acoustic=True)
        <FieldType.AC: 'AC'>
    """
    field_type = FieldType.parse(requested)
    return FieldType.AC if simulation_is_acoustic else field_type


def storage_key(selection):
    """Return the field-data storage key for a typed selection.

    Examples:
        >>> storage_key(FieldSelection.primary("EM_E", FieldComponent.X))
        'Fxr'
        >>> storage_key(FieldSelection.magnitude(FieldType.AC))
        'Fa'
    """
    if not isinstance(selection, FieldSelection):
        raise TypeError("selection must be a FieldSelection.")

    return _STORAGE_KEYS[(selection.component, selection.complex_part)]


def selection_from_storage_key(field_type, key):
    """Return the typed selection represented by a generic field-data key.

    Examples:
        >>> selection = selection_from_storage_key("EM_H", "Fyi")
        >>> (selection.field_type.value, selection.component.value, selection.complex_part.value)
        ('EM_H', 'y', 'secondary')
    """
    try:
        component, complex_part = _SELECTION_PARTS[key]
    except KeyError as exc:
        allowed = ", ".join(_SELECTION_PARTS)
        raise ValueError(f"Unknown field-data storage key {key!r}; expected one of: {allowed}.") from exc

    return FieldSelection(field_type, component, complex_part)


def complex_storage_keys(field_type, component):
    """Return the real and imaginary storage keys for a Cartesian component.

    The result is ordered as ``(real_key, imaginary_key)``. This order differs
    from the display-primary convention for ``z``, whose primary part is
    imaginary.

    Examples:
        >>> complex_storage_keys("EM_E", FieldComponent.X)
        ('Fxr', 'Fxi')
        >>> complex_storage_keys(FieldType.AC, "z")
        ('Fzr', 'Fzi')
    """
    component = FieldComponent.parse(component)
    try:
        real_part, imaginary_part = _REAL_IMAGINARY_PARTS[component]
    except KeyError as exc:
        raise ValueError(f"{component.value!r} is not a Cartesian field component.") from exc

    return (
        storage_key(FieldSelection(field_type, component, real_part)),
        storage_key(FieldSelection(field_type, component, imaginary_part)),
    )
