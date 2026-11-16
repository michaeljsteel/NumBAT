"""Typed selection and storage helpers for acoustic strain components."""

from dataclasses import dataclass
from enum import Enum

from numbat.nbtypes import ComplexPart


class StrainComponent(Enum):
    """The six independent components of a symmetric strain tensor."""

    XX = "xx"
    YY = "yy"
    ZZ = "zz"
    YZ = "yz"
    XZ = "xz"
    XY = "xy"

    @classmethod
    def parse(cls, value):
        if isinstance(value, cls):
            return value

        try:
            return cls(value)
        except ValueError as exc:
            allowed = ", ".join(member.value for member in cls)
            raise ValueError(f"strain component must be one of: {allowed}.") from exc


def parse_strain_components(spec):
    """Normalize strain components from typed values or comma-separated names."""
    if isinstance(spec, str):
        spec = spec.replace(",", " ").split()
    try:
        return tuple(StrainComponent.parse(component) for component in spec)
    except TypeError as exc:
        raise ValueError(
            "strain components must be an iterable of StrainComponent values "
            "or component names."
        ) from exc


_DISPLAY_PRIMARY_PART = {
    StrainComponent.XX: "r",
    StrainComponent.YY: "r",
    StrainComponent.ZZ: "r",
    StrainComponent.YZ: "i",
    StrainComponent.XZ: "i",
    StrainComponent.XY: "r",
}


def component_index(component):
    """Return the Voigt-order array index for one strain component."""
    return tuple(StrainComponent).index(StrainComponent.parse(component))


@dataclass(frozen=True)
class StrainSelection:
    """One displayed real or imaginary strain component."""

    component: StrainComponent
    complex_part: ComplexPart

    def __post_init__(self):
        object.__setattr__(self, "component", StrainComponent.parse(self.component))
        object.__setattr__(self, "complex_part", ComplexPart.parse(self.complex_part))

    @classmethod
    def primary(cls, component):
        return cls(component, ComplexPart.PRIMARY)

    @classmethod
    def secondary(cls, component):
        return cls(component, ComplexPart.SECONDARY)


def display_primary_part(component):
    """Return the physical real/imaginary part used for default display."""
    return _DISPLAY_PRIMARY_PART[StrainComponent.parse(component)]


def complex_storage_keys(component):
    """Return physical real and imaginary storage keys for one component."""
    component = StrainComponent.parse(component)
    prefix = f"S{component.value}"
    return f"{prefix}r", f"{prefix}i"


def storage_key(selection):
    """Return the storage key for a display-aware strain selection."""
    component = StrainComponent.parse(selection.component)
    primary_part = display_primary_part(component)
    if selection.complex_part == ComplexPart.PRIMARY:
        part = primary_part
    else:
        part = "i" if primary_part == "r" else "r"
    return f"S{component.value}{part}"


def selection_from_storage_key(key):
    """Return the display-aware selection represented by a physical storage key."""
    for component in StrainComponent:
        real_key, imaginary_key = complex_storage_keys(component)
        if key == real_key:
            part = (
                ComplexPart.PRIMARY
                if display_primary_part(component) == "r"
                else ComplexPart.SECONDARY
            )
            return StrainSelection(component, part)
        if key == imaginary_key:
            part = (
                ComplexPart.PRIMARY
                if display_primary_part(component) == "i"
                else ComplexPart.SECONDARY
            )
            return StrainSelection(component, part)

    allowed = ", ".join(
        key for component in StrainComponent for key in complex_storage_keys(component)
    )
    raise ValueError(f"Unknown strain storage key {key!r}; expected one of: {allowed}.")


def tex_label(selection):
    """Return the TeX label for one display-aware strain component."""
    part = "Re" if storage_key(selection).endswith("r") else "Im"
    return f"{part}[ $S_{{{selection.component.value}}}$ ]"
