
# types.py is a subroutine of NumBAT that contains various property datatypes.

# Copyright (C) 2017-2025  Michael Steel, Bjorn Sturmberg, Kokou Dossou.

# NumBAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


import math
from enum import Enum, IntEnum, auto
from numpy import array
import scipy.constants as const


#### Natural constants ########################################################

# These are all SI 2019 values

twopi = math.pi * 2.0
SI_Planck_h = const.Planck            # Planck's constant in Js (exact)
SI_speed_c = const.speed_of_light                 # Speed of light in vacuum in m/s (exact)
SI_charge_e = const.e                          # Charge of an electron in C (exact)
SI_permittivity_eps0 = const.epsilon_0
SI_permeability_mu0 = const.mu_0
SI_vacuum_impedance_Z0 = math.sqrt(SI_permeability_mu0/SI_permittivity_eps0)

###############################################################################

SI_ns = 1.0e-9
SI_ps = 1.0e-12
SI_fs = 1.0e-15

SI_THz = 1.0e12
SI_GHz = 1.0e9
SI_GPa = 1.0e9
SI_MHz = 1.0e6

SI_nm = 1.0e-9
SI_um = 1.0e-6
SI_km = 1.0e3
SI_kmps = 1.0e3

SI_to_gmpercc = 0.001   # 1kg/m^3 = 0.001 gm/cm^3

###########################################


unit_x = array([1.0, 0.0, 0.0])
unit_y = array([0.0, 1.0, 0.0])
unit_z = array([0.0, 0.0, 1.0])

class NormalisationConstants:
    def __init__(self):
        # adjustable
        self.t0 = SI_ns
        self.x0 = SI_um
        self.rho0 = 1000  # 1000 kg/m^3
        self.eps0 = SI_permittivity_eps0  # 8.854187817e-12  F/m

        # derived
        self.f0 = 1 / self.t0  # GHz
        self.v0 = self.x0 / self.t0  # km/s

        self.T0 = self.rho0 * self.v0**2  # GPa
        self.c0 = self.T0  # GPa

        self.E0 = math.sqrt(self.c0 / self.eps0)  # V/um
        self.V0 = self.E0 * self.x0    # V
        self.D0 = self.eps0 * self.E0  # C/m^2
        self.e0 = math.sqrt(self.c0 * self.eps0)  # C/m^2

        # self.p0 = nbtypes.SI_GPa



class SimType(Enum):
    EM = 'EM'
    AC = 'AC'


#User facing
class FieldType(Enum):
    EM_E = 'EM_E'
    EM_H = 'EM_H'
    AC = 'AC'

    @staticmethod
    def from_str(lab):
        if lab == 'EM_E':
            return FieldType.EM_E
        elif lab == 'EM_H':
            return FieldType.EM_H
        elif lab == 'AC':
            return FieldType.AC
        else:
            raise ValueError("The value of field_type must be either 'AC', 'EM_E' or 'EM_H'.")

# Internal
class FieldCode:
    def __init__(self, ft, force_AC=False):
        if not isinstance(ft, FieldType):
            ft = FieldType(ft)
        self._ft = FieldType.AC if force_AC else ft

    def is_EM(self):
        return self._ft in (FieldType.EM_E, FieldType.EM_H)

    def is_AC(self):
        return self._ft == FieldType.AC

    def is_EM_E(self):
        return self._ft == FieldType.EM_E

    def is_EM_H(self):
        return self._ft == FieldType.EM_H

    def as_str(self):
        return self._ft.name

    def mode_type_as_str(self):
        return "acoustic" if self.is_AC() else "optical"

    def as_field_type(self):
        return self._ft

class PointGroup(IntEnum):
    Unknown = 1
    C2V = auto()
    C3V = auto()
    C6V = auto()


class SymRep(IntEnum):
    # E = 1
    A = 0
    B1 = auto()
    B2 = auto()
    B3 = auto()
    Unknown = auto()


class QAcMethod(IntEnum):
    NotSet = 1
    Fixed = auto()
    Intrinsic = auto()


class CrystalGroup(IntEnum):
    Isotropic = 1
    Unknown = auto()
    Trigonal = auto()
    Cubic = auto()
    Hexagonal = auto()
    GeneralAnisotropic = auto()

class PiezoSupport(IntEnum):
    Unsupported = 1
    Disabled = auto()
    Enabled = auto()




class FieldTag:
    '''Class for transferring between user readable and code name versions of field components.'''

    @staticmethod  # make plain module method
    def make_from_field(ft, force_AC=False):
        '''Create a FieldTag from just the field type. Defaults to x component.'''

        if force_AC:
            ft=FieldType.AC
        return FieldTag.make_from_field_and_component(ft, 'x')

    @staticmethod  # make plain module method
    def make_from_field_and_component(ft, cc):
        '''Create a FieldTag from the field type and component suffix x, y, z, a, i, t.'''

        if not isinstance(ft, FieldType):
            ft = FieldType(ft)

        if ft == FieldType.EM_E:
            uc = 'E'+cc
        elif ft == FieldType.EM_H:
            uc = 'H'+cc
        else:
            uc = 'u'+cc
        return FieldTag(uc)


    @staticmethod
    def make_from_field_and_Fcode(ft, fc):
        '''Create a FieldTag from the real/imag-aware but field-agnostic component Fxr, Fxi, etc and the field type.

        The fc can override the preference for major real/imag part.'''
        if ft == FieldType.EM_E:
            uc = {'Fxr': 'Ex', 'Fxi': 'Ex', 'Fyr': 'Ey', 'Fyi': 'Ey',
                  'Fzr': 'Ez', 'Fzi': 'Ez', 'Fa': 'Ea', 'Fi': 'Ei', 'Ft': 'Et'}[fc]
        elif ft == FieldType.EM_H:
            uc = {'Fxr': 'Hx', 'Fxi': 'Hx', 'Fyr': 'Hy', 'Fyi': 'Hy',
                  'Fzr': 'Hz', 'Fzi': 'Hz', 'Fa': 'Ha', 'Fi': 'Hi', 'Ft': 'Ht'}[fc]
        else:
            uc = {'Fxr': 'ux', 'Fxi': 'ux', 'Fyr': 'uy', 'Fyi': 'uy',
                  'Fzr': 'uz', 'Fzi': 'uz', 'Fa': 'ua', 'Fi': 'ui', 'Ft': 'ut'}[fc]

        return FieldTag(uc, fc)

    @staticmethod
    def reim_major(fi):
        '''Returns the tag of the major part (real or imag) of the field component fi.'''
        try:
            return {'Ex': 'Exr', 'Ey': 'Eyr', 'Ez': 'Ezi', 'Ea': 'Ea', 'Ei': 'Ei', 'Et': 'Et',
                    'Hx': 'Hxr', 'Hy': 'Hyr', 'Hz': 'Hzi', 'Ha': 'Ha', 'Hi': 'Hi', 'Ht': 'Ht',
                    'ux': 'uxr', 'uy': 'uyr', 'uz': 'uzi', 'ua': 'ua', 'ui': 'ui', 'ut': 'ut'}[fi]
        except KeyError:
            return fi

    @staticmethod
    def reim_minor(fi):
        '''Returns the tag of the minor part (real or imag) of the field component fi.'''
        try:
            return {'Ex': 'Exi', 'Ey': 'Eyi', 'Ez': 'Ezi', 'Ea': None, 'Ei': None, 'Et': None,
                    'Hx': 'Hxi', 'Hy': 'Hyi', 'Hz': 'Hzi', 'Ha': None, 'Hi': None, 'Ht': None,
                    'ux': 'uxi', 'uy': 'uyi', 'uz': 'uzr', 'ua': None, 'ui': None, 'ut': None}[fi]
        except KeyError:
            return None


    def __init__(self, uc, fc=''):
        '''Make FieldTag knowing the actual field component Ex, Ey, Ez, Et, Ea, Ei, ux, uy etc.

           By default, the tag is set to the dominant real/imag part of the field component.
           This can be overridden or changed later using set_to_major() or set_to_minor().

           Args:
              uc: str  - User-friendly field component code.
              fc: str  - Field-agnostic symbol indicating rea/imag/a/t part of whichever field is active.

           A user code (uc) is one of Ex, Ey, Ez, Et, Ea, Ei, Hx, Hy, Hz, Ht, Ha, Hi, ux, uy, uz, ut, ua, ui.

           A field code (fc) is one of Fxr, Fxi, Fyr, Fyi, Fzr, Fzi, Ft, Fa, Fi.
        '''

        self._user_code = uc
        self._F = uc[0]    # E, H, or u

        self._field_type = {
            'E': FieldType.EM_E,
            'H': FieldType.EM_H,
            'u': FieldType.AC
            }[uc[0]]

        self._Fi = uc[:2]  # Ex, Ey, Ez, Ea, Ei, Et, Hx, Hy, etc
        self._xyz = uc[1]  # x, y, z, a, t

        # default real/imag given the _F value
        if fc:
            self._f_code = fc
        else:
            self._f_code = self.major_component_as_F()

    def __str__(self):
        return f"FieldTag: {self._user_code}, {self._f_code}"

    def __repr__(self):
        return f"FieldTag({self._user_code}, {self._f_code})"

    def field_component(self):
        return self._user_code

    def is_AC(self):
        return self._F == 'u'

    def is_EM_E(self):
        return self._F == 'E'

    def is_EM_H(self):
        return self._F == 'H'

    def is_EM(self):
        return self.is_EM_E() or self.is_EM_H()


    def is_abs(self): return self._f_code == 'Fa'

    def is_intens(self): return self._f_code == 'Fi'


    def is_x(self):
        return self._xyz == 'x'

    def is_y(self):
        return self._xyz == 'y'

    def is_z(self):
        return self._xyz == 'z'

    def is_signed_field(self): return self._f_code not in ('Ft', 'Fa', 'Fi')

    def is_transverse(self): return self._user_code in ('Et', 'Ht', 'ut')

    def is_minor(self): return self._f_code in ('Fxr', 'Fyr', 'Fzi')


    def major_component(self):
        """Returns the major component of the field in E/H/u notation."""
        return FieldTag.reim_major(self._user_code)

    def minor_component(self):
        """Returns the minor component of the field in E/H/u notation."""
        return FieldTag.reim_minor(self._user_code)

    def major_component_as_F(self):
        """Returns the major component of the field in F notation."""
        majco = self.major_component()
        majcoF = 'F' + majco[1:]
        return majcoF

    def minor_component_as_F(self):
        """Returns the minor component of the field in F notation."""
        minco = self.minor_component()
        mincoF = 'F' + minco[1:]
        return mincoF

    def component_as_F(self):
        return self._f_code

    def get_tex_plot_label(self):
        lab = {'Fx': r'Re[ $F_x$ ]', 'Fy': r'Re[ $F_y$ ]', 'Fz': r'Im[ $F_z$ ]',
               'Fxr': r'Re[ $F_x$ ]', 'Fyr': r'Re[ $F_y$ ]', 'Fzi': r'Im[ $F_z$ ]',
               'Fxi': r'Im[ $F_x$ ]', 'Fyi': r'Im[ $F_y$ ]', 'Fzr': r'Re[ $F_z$ ]',
               'Fa': r'$|\vec F|$', 'Fi': r'$|\vec F|^2$',
               'Ft': r'$\vec F_t$'}[self._f_code]
        return lab.replace('F', self._F)

    def set_to_major(self):
        self._f_code = self.major_component_as_F()

    def set_to_minor(self):
        self._f_code = self.minor_component_as_F()

    def domain_type_as_str(self):
        return "acoustic" if self.is_AC() else "optical"


    def field_type_as_str(self):
        return {FieldType.EM_E: 'electric',
                FieldType.EM_H: 'magnetic',
                FieldType.AC: 'displacement'
                }[self._field_type]

    def field_type_label(self):
        if self.is_AC():
            return 'AC'
        else:
            return 'EM'

    def linestyle(self, all_comps):
        """Field amplitudes are dashed (major) or dotted (minor) if they are plotted alongside an absolute value"""

        if self.is_abs() or self.is_intens():
            return 'solid'
        mixed = 'a' in all_comps or 'a' in all_comps

        if not mixed:
            return 'solid'

        return 'dashed' if self.is_minor() else 'dotted'

    def linecolor(self):
        """Field amplitudes have standard colours."""

        if self.is_abs() or self.is_intens():
            return 'red'
        if self.is_x():
            return 'blue'
        if self.is_y():
            return 'green'
        if self.is_z():
            return 'brown'

        return None

    def clone_as(self, cc):
        """Create new FieldTag with same field but component cc of ['x', 'y', 'z', 'a', 't']"""
        tag = FieldTag.make_from_field_and_component(self._field_type, cc)
        return tag

    def as_field_type(self):
        return self._field_type
