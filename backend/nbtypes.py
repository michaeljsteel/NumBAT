
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


#### Natural constants ########################################################

# These are all SI 2019 values

twopi = math.pi * 2.0
SI_Planck_h = 6.626_070_15e-34            # Planck's constant in Js (exact)
SI_speed_c = 299792458                     # Speed of light in vacuum in m/s (exact)
SI_charge_e = 1.602_176_634e-19            # Charge of an electron in C (exact)
SI_permittivity_eps0 = 8.854_187_8188e-12
SI_permeability_mu0 = 1.256_637_061_27e-6
SI_vacuum_impedance_Z0 = math.sqrt(SI_permeability_mu0/SI_permittivity_eps0)

###############################################################################

SI_THz = 1.0e12
SI_GHz = 1.0e9
SI_MHz = 1.0e6

SI_nm = 1.0e-9
SI_um = 1.0e-6
SI_to_gmpercc = 0.001   # 1kg/m^3 = 0.001 gm/cm^3

###########################################


unit_x = array([1.0, 0.0, 0.0])
unit_y = array([0.0, 1.0, 0.0])
unit_z = array([0.0, 0.0, 1.0])


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
    GeneralAnisotropic = auto()





class FieldTag:
    '''Class for transferring between user readable and code name versions of field components.'''

    @staticmethod  # make plain module method
    def make_from_field_and_component(ft, cc):
        '''Create a FieldTag from the field type and component suffix x, y, z, abs, t.'''
        if ft == FieldType.EM_E:
            uc = 'E'+cc
        elif ft == FieldType.EM_H:
            uc = 'H'+cc
        else:
            uc = 'u'+cc
        return FieldTag(uc)


    @staticmethod
    def make_comp_from_field_and_Fcode(ft, fc):
        '''Create a FieldTag from the real/imag-aware but field-agnostic component Fxr, Fxi, etc and the field type.

        The fc can override the preference for major real/imag part.'''
        if ft == FieldType.EM_E:
            uc = {'Fxr': 'Ex', 'Fxi': 'Ex', 'Fyr': 'Ey', 'Fyi': 'Ey',
                  'Fzr': 'Ez', 'Fzi': 'Ez', 'Fabs': 'Eabs', 'Ft': 'Et'}[fc]
        elif ft == FieldType.EM_H:
            uc = {'Fxr': 'Hx', 'Fxi': 'Hx', 'Fyr': 'Hy', 'Fyi': 'Hy',
                  'Fzr': 'Hz', 'Fzi': 'Hz', 'Fabs': 'Habs', 'Ft': 'Ht'}[fc]
        else:
            uc = {'Fxr': 'ux', 'Fxi': 'ux', 'Fyr': 'uy', 'Fyi': 'uy',
                  'Fzr': 'uz', 'Fzi': 'uz', 'Fabs': 'uabs', 'Ft': 'ut'}[fc]

        return FieldTag(uc, fc)

    def __init__(self, uc, fc=''):
        '''Make FieldTag knowing the actual field component Ex, Ey, Ez, Et, Eabs, ux, uy etc.

           By default, the tag is set to the dominant real/imag part of the field component.
           This can be overridden or changed later using set_to_major() or set_to_minor().

           Args:
              uc: str  - User-friendly field component code.
              fc: str  - Field-agnostic symbol indicating rea/imag/abs/t part of whichever field is active.

           A user code (uc) is one of Ex, Ey, Ez, Et, Eabs, Hx, Hy, Hz, Ht, Habs, ux, uy, uz, ut, uabs.

           A field code (fc) is one of Fxr, Fxi, Fyr, Fyi, Fzr, Fzi, Fabs, Ft.
        '''

        self._user_code = uc
        self._F = uc[0]    # E, H, or u
        self._Fi = uc[:2]  # Ex, Ey, Ez, Ea, Et, Hx, Hy, etc
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

    #def is_user_code(self, uc):
    #     return uc in ('Ex', 'Ey', 'Ez', 'Eabs', 'Et', 'Hx', 'Hy', 'Hz', 'Habs', 'Ht', 'ux', 'uy', 'uz', 'uabs', 'ut')

    #def is_field_code(self, fc):
    #    return fc in ('Fxr', 'Fxi', 'Fyr', 'Fyi', 'Fzr', 'Fzi', 'Fabs', 'Ft')

    def is_AC(self): return self._F == 'u'

    def is_abs(self): return self._f_code == 'Fabs'

    def is_x(self): return self._xyz == 'x'
    def is_y(self): return self._xyz == 'y'
    def is_z(self): return self._xyz == 'z'

    def is_signed_field(self): return self._f_code not in ('Ft', 'Fabs')

    def is_transverse(self): return self._user_code in ('Et', 'Ht', 'ut')

    def is_minor(self): return self._f_code in ('Fxr', 'Fyr', 'Fzi')

    @staticmethod
    def reim_major(fi):
        '''Returns the tag of the major part (real or imag) of the field component fi.'''
        try:
            return {'Ex': 'Exr', 'Ey': 'Eyr', 'Ez': 'Ezi', 'Ea': 'Ea',
                    'Hx': 'Hxr', 'Hy': 'Hyr', 'Hz': 'Hzi', 'Ha': 'Ha',
                    'ux': 'uxr', 'uy': 'uyr', 'uz': 'uzi', 'ua': 'ua'}[fi]
        except KeyError:
            return fi

    @staticmethod
    def reim_minor(fi):
        '''Returns the tag of the minor part (real or imag) of the field component fi.'''
        try:
            return {'Ex': 'Exi', 'Ey': 'Eyi', 'Ez': 'Ezi', 'Ea': 'Ea',
                    'Hx': 'Hxi', 'Hy': 'Hyi', 'Hz': 'Hzi', 'Ha': None,
                    'ux': 'uxi', 'uy': 'uyi', 'uz': 'uzr', 'ua': None}[fi]
        except KeyError:
            return None

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
        lab = {'Fx': r'Re($F_x$)', 'Fy': r'Re($F_y$)', 'Fz': r'Im($F_z$)', 'Fxr': r'Re($F_x$)',
               'Fyr': r'Re($F_y$)', 'Fzi': r'Im($F_z$)', 'Fxi': r'Im($F_x$)', 'Fyi': r'Im($F_y$)', 'Fzr': r'Re($F_z$)',
               'Fabs': r'$|\vec F|^2$', 'Ft': r'$\vec F_t$'}[self._f_code]  # adjusted so that Fabs gives |F|^2
        return lab.replace('F', self._F)

    def set_to_major(self):
        self._f_code = self.major_component_as_F()

    def set_to_minor(self):
        self._f_code = self.minor_component_as_F()


    def field_type_label(self):
        if self.is_AC:
            return 'AC'
        else:
            return 'EM'

    def linestyle(self, all_comps):
        """Field amplitudes are dashed (major) or dotted (minor) if they are plotted alongside an absolute value"""

        if self.is_abs():
            return 'solid'
        mixed = 'a' in all_comps or 'abs' in all_comps

        if not mixed:
            return 'solid'

        return 'dashed' if self.is_minor() else 'dotted'

    def linecolor(self):
        """Field amplitudes have standard colours."""

        if self.is_abs():
            return 'red'
        if self.is_x():
            return 'blue'
        if self.is_y():
            return 'green'
        if self.is_z():
            return 'brown'

        return None
