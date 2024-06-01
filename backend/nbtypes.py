
# types.py is a subroutine of NumBAT that contains various property datatypes.

# Copyright (C) 2017 Bjorn Sturmberg, Kokou Dossou.

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

from enum import Enum, IntEnum, auto
import math
import numpy as np

#### Natural constants ########################################################

# These are all SI 2019 values

twopi = np.pi * 2.0
Plancks_h = 6.626_070_15e-34            # Planck's constant in Js (exact)
speed_c = 299792458                     # Speed of light in vacuum in m/s (exact)
charge_F = 1.602_176_634e-19            # Charge of an electron in C (exact)
permittivity_Fps0 = 8.854_187_8188e-12
permeability_mu0 = 1.256_637_061_27e-6
vacuum_impedance_Z0 = math.sqrt(permeability_mu0/permittivity_Fps0)
            
###############################################################################

SI_THz = 1.0e12
SI_GHz = 1.0e9
SI_MHz = 1.0e6

SI_nm = 1.0e-9
SI_um = 1.0e-6
SI_to_gmpercc = 0.001
###########################################


unit_x = np.array([1.0, 0.0, 0.0])
unit_y = np.array([0.0, 1.0, 0.0])
unit_z = np.array([0.0, 0.0, 1.0])


class SimType(Enum):
    EM = 'EM'
    AC = 'AC'


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


class component_t(object):
    @staticmethod
    def Ecomp(c): return component_t('E'+c)

    @staticmethod
    def Hcomp(c): return component_t('H'+c)

    @staticmethod
    def ucomp(c): return component_t('u'+c)

    @staticmethod
    def make_comp(emac, fc):
        if emac == FieldType.EM_E:
            uc = {'Fxr': 'Ex', 'Fxi': 'Ex', 'Fyr': 'Ey', 'Fyi': 'Ey',
                  'Fzr': 'Ez', 'Fzi': 'Ez', 'Fabs': 'Eabs', 'Ft': 'Et'}[fc]
        elif emac == FieldType.EM_H:
            uc = {'Fxr': 'Hx', 'Fxi': 'Hx', 'Fyr': 'Hy', 'Fyi': 'Hy',
                  'Fzr': 'Hz', 'Fzi': 'Hz', 'Fabs': 'Habs', 'Ft': 'Ht'}[fc]
        else:
            uc = {'Fxr': 'ux', 'Fxi': 'ux', 'Fyr': 'uy', 'Fyi': 'uy',
                  'Fzr': 'uz', 'Fzi': 'uz', 'Fabs': 'uabs', 'Ft': 'ut'}[fc]
        cc = component_t(uc)
        cc._f_code = fc # we override the _f_code in __init__ because we may not be asking for the dominant re/im part
        return cc

    def __init__(self, uc):  
        # uc is an actual field component name:  Ex, Ey, Ez, Et, Eabs, ux, uy etc. 
        # _f_code is the field-agnostic form starting with F and in the dominant real/imag part
        self._user_code = uc
        self._F = uc[0]  # E, H, or u
        self._Fi = uc[:2]  # Ex, Ey, Ez, Ea, Et, Hx, Hy, etc
        self._xyz = uc[1]  # x, y, z, a, t
        self._Fimaj = self.reim_major(self._Fi)
        self._Fimin = self.reim_minor(self._Fi)

        # default real/imag given the _F value
        self._f_code = {'Ex': 'Fxr', 'Ey': 'Fyr', 'Ez': 'Fzi', 'Eabs': 'Fabs', 'Et': 'Ft',
                        'Hx': 'Fxr', 'Hy': 'Fyr', 'Hz': 'Fzi', 'Habs': 'Fabs', 'Ht': 'Ft',
                        'ux': 'Fxr', 'uy': 'Fyr', 'uz': 'Fzi', 'uabs': 'Fabs', 'ut': 'Ft',
                        }[self._user_code]

    def is_AC(self): return self._U == 'u'
    def get_label(self):
        c = self._F
        lab = {'Fx': r'Re($F_x$)', 'Fy': r'Re($F_y$)', 'Fz': r'Im($F_z$)', 'Fxr': r'Re($F_x$)',
               'Fyr': r'Re($F_y$)', 'Fzi': r'Im($F_z$)', 'Fxi': r'Im($F_x$)', 'Fyi': r'Im($F_y$)', 'Fzr': r'Re($F_z$)',
               'Fabs': r'$|\vec F|^2$', 'Ft': r'$\vec F_t$'}[self._f_code]  # adjusted so that Fabs gives |F|^2
        return lab.replace('F', c)

    def is_abs(self): return self._f_code == 'Fabs'
    def is_signed_field(self): return self._f_code not in ('Ft', 'Fabs')
    def is_transverse(self): return self._user_code in ('Ft', 'Ht', 'ut')
    def is_dominant(self): return self._f_code in ('Fxr', 'Fyr', 'Fzi')

    def reim_major(self, fi):
        try:
            return {'Ex': 'Exr', 'Ey': 'Eyr', 'Ez': 'Ezi', 'Ea': 'Ea',
                    'Hx': 'Hxr', 'Hy': 'Hyr', 'Hz': 'Hzi', 'Ha': 'Ha',
                    'ux': 'uxr', 'uy': 'uyr', 'uz': 'uzi', 'ua': 'ua'}[fi]
        except KeyError:
            return fi

    def reim_minor(self, fi):
        try:
            return {'Ex': 'Exi', 'Ey': 'Eyi', 'Ez': 'Ezr', 'Ea': None,
                    'Hx': 'Hxi', 'Hy': 'Hyi', 'Hz': 'Hyr', 'Ha': None,
                    'ux': 'uxi', 'uy': 'uyi', 'uz': 'uyr', 'ua': None}[fi]
        except KeyError:
            return None
