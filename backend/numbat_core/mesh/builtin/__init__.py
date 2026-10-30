"""Built-in NumBAT geometry templates, split by family from the former builtin_meshes module."""

from .basic import Circular as Circular, Rectangular as Rectangular, Triangular as Triangular
from .two_incl import TwoIncl as TwoIncl, TwoInclVert as TwoInclVert
from .onion import (
    Onion as Onion,
    Onion1 as Onion1,
    Onion2 as Onion2,
    Onion3 as Onion3,
    CircOnion as CircOnion,
    CircOnion1 as CircOnion1,
    CircOnion2 as CircOnion2,
    CircOnion3 as CircOnion3,
)
from .rib import (
    Pedestal as Pedestal,
    TrapezoidalRib as TrapezoidalRib,
    Rib as Rib,
    RibCoated as RibCoated,
    RibDoubleCoated as RibDoubleCoated,
    RibMkII as RibMkII,
)
from .slot import Slot as Slot, SlotCoated as SlotCoated
from .rounded import (
    RoundedBowtie as RoundedBowtie,
    RoundedWedge as RoundedWedge,
    RoundedWedge3 as RoundedWedge3,
    RoundedWedge4 as RoundedWedge4,
    RoundedWedge6 as RoundedWedge6,
)
