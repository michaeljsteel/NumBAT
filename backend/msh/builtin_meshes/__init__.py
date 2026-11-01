"""Compatibility facade for built-in NumBAT geometry templates."""

from numbat_core._deprecation import warn_legacy_import
from numbat_core.mesh.builtin import (
    CircOnion,
    CircOnion1,
    CircOnion2,
    CircOnion3,
    Circular,
    Onion,
    Onion1,
    Onion2,
    Onion3,
    Pedestal,
    Rectangular,
    Rib,
    RibCoated,
    RibDoubleCoated,
    RibMkII,
    RoundedBowtie,
    RoundedWedge,
    RoundedWedge3,
    RoundedWedge4,
    RoundedWedge6,
    Slot,
    SlotCoated,
    TrapezoidalRib,
    Triangular,
    TwoIncl,
    TwoInclVert,
)

warn_legacy_import("msh.builtin_meshes", "numbat.mesh")
