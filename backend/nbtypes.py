
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

#TODO: change name of this class. FieldCompHandle ?
class component_t(object):
    '''Class for transferring between user readable and code name versions of field components.'''

    #@staticmethod
    #def Ecomp(c): return component_t('E'+c)
#
#    @staticmethod
#    def Hcomp(c): return component_t('H'+c)
#
#    @staticmethod
#    def ucomp(c): return component_t('u'+c)

    @staticmethod
    def make_comp_noreim(emac, fc):
        '''Create a component_t from the field-agnostic component Fx, Fy, Fz, Fabs, etc and the field type.
        Set the _f_code to the default.  '''
        if emac == FieldType.EM_E:
            uc=f'E{fc[1]}'
        elif emac == FieldType.EM_H:
            uc=f'H{fc[1]}'
        else:
            uc=f'u{fc[1]}'
        cc = component_t(uc)
        return cc

    @staticmethod
    def make_comp_from_component(emac, cc):
        '''Create a component_t from the field type and component suffix x, y, z, abs, t.'''
        if emac == FieldType.EM_E:
            uc = 'E'+cc
        elif emac == FieldType.EM_H:
            uc = 'H'+cc
        else:
            uc = 'u'+cc
        t_cc = component_t(uc)
        print('making ', cc, t_cc._f_code)
        return t_cc


    @staticmethod
    def make_comp_from_Fcode(emac, fc):
        '''Create a component_t from the real/imag-aware but field-agnostic component Fxr, Fxi, etc and the field type.'''
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
        #TODO: change constructor so it has option to not choose the domnant
        return cc

    def __init__(self, uc):
        '''Make component_t knowing the actual field component Ex, Ey, Ez, Et, Eabs, ux, uy etc.

           Sets  _f_code to the field-agnostic form starting with F and with the dominant real/imag part for that component.
        '''

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

    def is_AC(self): return self._F == 'u'
    def emac(self):
        if self.is_AC:
            return 'AC'
        else:
            return 'EM'

    def get_label(self):
        c = self._F
        lab = {'Fx': r'Re($F_x$)', 'Fy': r'Re($F_y$)', 'Fz': r'Im($F_z$)', 'Fxr': r'Re($F_x$)',
               'Fyr': r'Re($F_y$)', 'Fzi': r'Im($F_z$)', 'Fxi': r'Im($F_x$)', 'Fyi': r'Im($F_y$)', 'Fzr': r'Re($F_z$)',
               'Fabs': r'$|\vec F|^2$', 'Ft': r'$\vec F_t$'}[self._f_code]  # adjusted so that Fabs gives |F|^2
        return lab.replace('F', c)

    def is_abs(self): return self._f_code == 'Fabs'
    def is_signed_field(self): return self._f_code not in ('Ft', 'Fabs')
    def is_transverse(self): return self._user_code in ('Et', 'Ht', 'ut')
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
            return {'Ex': 'Exi', 'Ey': 'Eyi', 'Ez': 'Ezi', 'Ea': 'Ea',
                    'Hx': 'Hxi', 'Hy': 'Hyi', 'Hz': 'Hzi', 'Ha': None,
                    'ux': 'uxi', 'uy': 'uyi', 'uz': 'uzr', 'ua': None}[fi]
        except KeyError:
            return None
