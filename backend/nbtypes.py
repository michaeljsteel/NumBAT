
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

#### Natural constants ########################################################
ASTM15_tot_I   = 900.084            # Integral ASTM 1.5 solar irradiance W/m**2
Plancks_h      = 6.62607015e-34     # Planck's constant in Js (exact)
speed_c        = 299792458          # Speed of light in vacuum in m/s (exact)
charge_e       = 1.602176634e-19    # Charge of an electron in C (exact)
###############################################################################

class SimType(Enum):
  EM = 'EM'
  AC = 'AC'

class FieldType(Enum):
  EM_E = 'EM_E'
  EM_H = 'EM_H'
  AC = 'AC'

  @staticmethod
  def from_str(lab):
    if lab == 'EM_E': return FieldType.EM_E
    elif lab == 'EM_H': return FieldType.EM_H
    elif lab == 'AC': return FieldType.AC
    else:
      raise NotImplementedError

class PointGroup(IntEnum):
  Unknown = 1
  C2V = auto()
  C3V = auto()
  C6V = auto()

class SymRep(IntEnum):
#E = 1
  A = 0
  B1 = auto()
  B2 = auto()
  B3 = auto()
  Unknown = auto()


class QAcMethod(IntEnum):
  NotSet=1
  Fixed=auto()
  Intrinsic=auto()

class CrystalGroup(IntEnum):
  Isotropic=1
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
    if emac==FieldType.EM_E: 
      uc = { 'Fxr':'Ex', 'Fxi':'Ex', 'Fyr':'Ey', 'Fyi':'Ey', 'Fzr':'Ez', 'Fzi':'Ez', 'Fabs':'Eabs', 'Ft':'Et'}[fc]
    elif emac==FieldType.EM_H: 
      uc = { 'Fxr':'Hx', 'Fxi':'Hx', 'Fyr':'Hy', 'Fyi':'Hy', 'Fzr':'Hz', 'Fzi':'Hz', 'Fabs':'Habs', 'Ft':'Ht'}[fc]
    else:
      uc = { 'Fxr':'ux', 'Fxi':'ux', 'Fyr':'uy', 'Fyi':'uy', 'Fzr':'uz', 'Fzi':'uz', 'Fabs':'uabs', 'Ft':'ut'}[fc]
    cc= component_t(uc)
    cc._f_code=fc
    return cc

  def __init__(self, uc): #must be an Ex, Ey style, not an EM_AC,Fxr style
    self._user_code=uc
    self._E=uc[0] # E, H, or u
    self._Ei=uc[:2] # Ex, Ey, Ez, Ea, Et, Hx, Hy, etc
    self._xyz=uc[1] # x, y, z, a, t
    self._Eimaj=self.reim_major(self._Ei)
    self._Eimin=self.reim_minor(self._Ei)

    #default real/imag given the _E value                          
    self._f_code= {'Ex':'Fxr', 'Ey':'Fyr', 'Ez':'Fzi', 'Eabs':'Fabs', 'Et':'Ft',
                'Hx':'Fxr', 'Hy':'Fyr', 'Hz':'Fzi', 'Habs':'Fabs', 'Ht':'Ft',
                'ux':'Fxr', 'uy':'Fyr', 'uz':'Fzi', 'uabs':'Fabs','ut':'Ft',
                   }[self._user_code]


  def get_label(self):
    c=self._E
    lab= { 'Fx':r'Re($E_x$)', 'Fy':r'Re($E_y$)', 'Fz':r'Im($E_z$)', 'Fxr':r'Re($E_x$)',
              'Fyr':r'Re($E_y$)', 'Fzi':r'Im($E_z$)', 'Fxi':r'Im($E_x$)', 'Fyi':r'Im($E_y$)', 'Fzr':r'Re($E_z$)',
              'Fabs':r'$|\vec E|^2$', 'Ft':r'$\vec E_t$'}[self._f_code]  #adjusted so that Fabs gives |F|^2
    return lab.replace('E', c)

  def is_abs(self): return self._f_code == 'Fabs'
  def is_signed_field(self): return self._f_code not in ('Ft', 'Fabs')
  def is_transverse(self): return self._user_code.endswith('t')
  def is_dominant(self): return self._f_code in ('Fxr', 'Fyr', 'Fzi')


  def reim_major(self, fi): 
    try:
      return { 'Ex':'Exr', 'Ey':'Eyr', 'Ez':'Ezi', 'Ea':'Ea',
             'Hx':'Hxr', 'Hy':'Hyr', 'Hz':'Hzi', 'Ha':'Ha',
             'ux':'uxr', 'uy':'uyr', 'uz':'uzi', 'ua':'ua' }[fi]
    except KeyError:
      return fi

  def reim_minor(self, fi): 
    try:
      return { 'Ex':'Exi', 'Ey':'Eyi', 'Ez':'Ezr', 'Ea':None,
             'Hx':'Hxi', 'Hy':'Hyi', 'Hz':'Hyr', 'Ha':None,
             'ux':'uxi', 'uy':'uyi', 'uz':'uyr', 'ua':None }[fi]
    except KeyError:
      return None


