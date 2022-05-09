# mode_calcs.py is a subroutine of NumBAT that contains methods to
# calculate the EM and Acoustic modes of a structure.

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


import numpy as np
import sys
import os
import copy
sys.path.append("../backend/")
from math import *

import plotting
import integration
from nbtypes import *
from fortran import NumBAT


VacCSpeed=299792458

class Mode(object):
  '''This is a base class for both EM and AC modes.''' 
  def __init__(self, simmo, m):
    self.mode_num=m
    self.owning_sim=simmo
    self.fracs=[]  # fx, fy, ft, fz
    self.r0=None  # centre of mass
    self.w2=None  # second moment width
    self.r0_offset=(0.0, 0.0)
    self.extra_data={}
    self.analysed = False

  def add_mode_data(self, d): 
    '''Adds a dictionary of user-defined information about a mode.
  
       :param dict d: Dict of (str, data) tuples of user-defined information about a mode.
  '''
    self.extra_data.update(d)

  def get_mode_data(self): 
    '''Return dictionary of user-defined information about the mode.
       
       :return: Dictionary of user-defined information about the mode.
       :rtype: dict(str, obj)
       '''
    return self.extra_data

  def field_fracs(self): 
    '''Returns tuple (*fx*, *fy*, *fz*, *ft*) of "fraction" of mode contained in *x*, *y*, *z* or *t* (sum of transverse *x+y*) components.

       Note that *fraction* is defined through a simple overlap integral. It does not necessarily represent the fraction of energy density in the component.

       :return: Tuple of mode fractions
       :rtype: tuple(float, float, float, float)
    '''
    if not self.analysed: print('mode has not being analysed')
    return self.fracs

  def __str__(self): 
    '''String representation of the mode.'''
    s='Abstract mode class'
    return s

  def is_poln_ex(self):
    '''Returns true if mode is predominantly x-polarised (ie if *fx*>0.7).
       
       :rtype: bool
       '''
    polthresh=.7
    return self.fracs[0]>polthresh

  def is_poln_ey(self):
    '''Returns true if mode is predominantly y-polarised (ie if *fy*>0.7).
       
       :rtype: bool
       '''
    polthresh=.7
    return self.fracs[1]>polthresh

  def is_poln_indeterminate(self):
    '''Returns true if transverse polarisation is neither predominantly *x* or *y* oriented.

       :rtype: bool
       '''
    return not (self.is_poln_ex() or self.is_poln_ey())

  def is_EM(self):
    '''Returns true if the mode is an electromagnetic mode.

       :rtype: bool
       '''
    return self.owning_simmo.is_EM()

  def is_AC(self):
    '''Returns true if the mode is an acoustic mode.

       :rtype: bool
       '''
    return self.owning_simmo.is_AC()

  def center_of_mass(self): 
    '''Returns the centre of mass of the mode relative to the specified origin.

       :rtype: float
    '''
    return self.r0-self.r0_offset

  def second_moment_widths(self): 
    '''Returns the second moment widths :math:`(w_x, w_y, \sqrt{w_x^2+w_y^2})` of the mode relative to the specified origin.

       :rtype: (float, float, float)
       '''
    return self.w2

  def center_of_mass_x(self): 
    '''Returns the $x$ component moment of the centre of mass  of the mode.

       :rtype: float
       '''
    return self.r0[0]-self.r0_offset[0]

  def center_of_mass_y(self): 
    '''Returns the $y$ component moment of the centre of mass  of the mode.

       :rtype: float
       '''
    return self.r0[1]-self.r0_offset[1]

  def wx(self): 
    '''Returns the $x$ component moment of the second moment width.

       :rtype: float
       '''
    return self.w2[0]

  def wy(self): 
    '''Returns the $y$ component moment of the second moment width.

       :rtype: float
       '''
    return self.w2[1]

  def w0(self): 
    '''Returns the combined second moment width :math:`\sqrt{w_x^2+w_y^2}`.

       :rtype: float
       '''
    return self.w2[2]

  def set_r0_offset(self, x0, y0):
    '''Sets the transverse position in the grid that is to be regarded as the origin for calculations of center-of-mass.
  
       This can be useful in aligning the FEM coordinate grid with a physically sensible place in the waveguide.

       :param float x0: *x* position of nominal origin.
       :param float y0: *y* position of nominal origin.
  '''
    self.r0_offset=(x0, y0)

  def analyse_mode(self, v_x, v_y, m_Refx, m_Refy, m_Refz, m_Imfx, m_Imfy, m_Imfz, m_Absf):
    '''Perform a series of measurements on the mode *f* to determine polarisation fractions, second moment widths etc.
       
       :param array v_x: Vector of x points.
       :param array v_y: Vector of y points.
       :param array m_Refx: Matrix of real part of fx.
       :param array m_Refy: Matrix of real part of fy.
       :param array m_Refz: Matrix of real part of fz.
       :param array m_Imfx: Matrix of imaginary part of fx.
       :param array m_Imfy: Matrix of imaginary part of fy.
       :param array m_Imfz: Matrix of imaginary part of fz.
       ''' 
    self.analysed = True
    s_fx=np.sum(np.sum(m_Refx*m_Refx+m_Imfx*m_Imfx))
    s_fy=np.sum(np.sum(m_Refy*m_Refy+m_Imfy*m_Imfy)) 
    s_fz=np.sum(np.sum(m_Refz*m_Refz+m_Imfz*m_Imfz)) 
    s_f=s_fx+s_fy+s_fz
    f_x=s_fx/s_f
    f_y=s_fy/s_f
    f_t=f_x+f_y
    f_z=s_fz/s_f
    self.fracs=[f_x, f_y, f_t, f_z]
    [m_x, m_y]=np.meshgrid(v_x, v_y, indexing='ij') # TODO: IS THIS RIGHT?
# print('vs', v_x[0], v_x[1], v_x[-1], v_y[0], v_y[1], v_y[-1])
#    print('sh', len(v_x), len(v_y), m_Refx.shape, m_x.shape)
#    print('m_xa', m_x[0,:5])
#    print('m_xb', m_x[:5,0])
#    print('m_xb', m_x[-5:,0])
#    print('m_ya', m_y[0,:5])
#    print('m_yb', m_y[0,-5:])
#    print('m_yc', m_y[:5,0])


    m_mod= m_Refx*m_Refx+m_Imfx*m_Imfx+m_Refy*m_Refy+m_Imfy*m_Imfy+m_Refz*m_Refz+m_Imfz*m_Imfz
    m_xmod= m_x * m_mod  # could do this by broadcasting without meshgrid?
    m_ymod= m_y * m_mod
    x0=np.sum(np.sum(m_xmod))/s_f
    y0=np.sum(np.sum(m_ymod))/s_f
    m_x2mod= np.power((m_x-x0),2) * m_mod
    m_y2mod= np.power((m_y-y0),2) * m_mod
    w2x=sqrt(np.sum(np.sum(m_x2mod))/s_f)
    w2y=sqrt(np.sum(np.sum(m_y2mod))/s_f)
    w2=sqrt(w2x*w2x+w2y*w2y)
#print ('sums', s_f, np.sum(np.sum(m_mod)), w2x, w2y, w2)
    self.r0=np.array([x0, y0])
    self.w2=np.array([w2x, w2y, w2])


class ModeEM(Mode):
  '''Class representing a single electromagnetic (EM) mode.'''
  def __init__(self, simmo, m):
    super().__init__(simmo, m)

  def __str__(self): 
    s='EM mode # {0}'.format(self.mode_num)
    return s

  def analyse_mode(self, v_x, v_y, m_Refx, m_Refy, m_Refz, m_Imfx, m_Imfy, m_Imfz, m_Absf):
    super().analyse_mode(v_x, v_y, m_Refx, m_Refy, m_Refz, m_Imfx, m_Imfy, m_Imfz, m_Absf)


class ModeAC(Mode):
  '''Class representing a single acoustic (AC) mode.'''
  def __init__(self, simmo, m):
    super().__init__(simmo, m)

    self.gain={}  # { (EM_p_i, EM_s_j): gain}
    self.gain_PE={}
    self.gain_MB={}

  def __str__(self): 
    s='AC mode # {0}'.format(self.mode_num)
    return s

  def analyse_mode(self, v_x, v_y, m_Refx, m_Refy, m_Refz, m_Imfx, m_Imfy, m_Imfz, m_Absf):
    super().analyse_mode(v_x, v_y, m_Refx, m_Refy, m_Refz, m_Imfx, m_Imfy, m_Imfz, m_Absf)


class Simmo(object):
    '''Class for calculating the electromagnetic and/or acoustic modes of a ``Struct`` object.
    '''

    def __init__(self, structure, num_modes=20, wl_nm=1, n_eff=None, shift_Hz=None, 
                 k_AC=None, EM_sim=None, Stokes=False, 
                 calc_EM_mode_energy=False, calc_AC_mode_power=False, debug=False):
        '''Sets up the problem for the mode calculation at a given optical wavelength `wl_nm` or acoustic wavenumber `k_AC`.

           For electromagnetic problems, the tool solves for the effective index or wavenumber at a given wavelength.
           For acoustic problems, the tool solves for the acoustic frequency at a given wavenumber.

             :param Simmo structure: The waveguide structure to be solved.
             :param int num_modes: The number of modes to be found.
             :param float wl_nm: For electromagnetic problems, the vacuum wavelength in nanometers.
             :param float n_eff: For electromagnetic problems, an estimated effective index to begin the eigenvalue search.
             :param float shift_Hz: For acoustic problems, an estimated frequency offset to begin the eigenvalue search.
             :param float k_AC: For acoustic problems, the acoustic wavevector of the mode.
             :param float EM_sim: For acoustic problems, the results of a previously solved EM problem to speed calculations. 
             :param bool calc_EM_mode_energy: For electromagnetic problems, whether to calculate the optical mode energy.
             :param bool calc_AC_mode_power: For acoustic problems, whether to calculate the acoustic mode power.
          '''

        self.structure = structure
        self.wl_m = wl_nm*1e-9
        self.n_eff = n_eff
        self.shift_Hz = shift_Hz

        self.k_AC = k_AC
        self.Omega_AC = None
        self.EM_sim = EM_sim

        self.num_modes = num_modes
        self.Stokes = Stokes
        self.mode_pol = None
        self.k_0 = 2 * np.pi / self.wl_m
        # just off normal incidence to avoid degeneracies
        self.k_pll = np.array([1e-16, 1e-16])
        speed_c = 299792458
        self.omega_EM = 2*np.pi*speed_c/self.wl_m # Angular freq in units of rad/s
        self.calc_EM_mode_energy = calc_EM_mode_energy
        self.calc_AC_mode_power = calc_AC_mode_power

        self.EM_mode_energy = None
        self.EM_mode_power = None

        self.AC_mode_energy = None
        self.AC_mode_power = None

        self.debug = debug
        self.EM_AC = 'EM'
        self.sym_reps = None
        self.point_group = PointGroup.Unknown
        self.Q_method = QAcMethod.NotSet
        self.ac_alpha_t = None   # temporal acoustic loss [1/s]
        self.ac_linewidth = None   # acoustic linewidth [Hz]
        self.ac_Qmech = None   # acoustic mechanical Q [dimless]

        self.mode_set=[]

    def is_EM(self): 
      '''Returns true if the solver is setup for an electromagnetic problem.'''
      return self.EM_AC == 'EM'

    def is_AC(self): 
      '''Returns true if the solver is setup for an acoustic problem.'''
      return self.EM_AC != 'EM'
    
    def get_modes(self):
      '''Returns an array of class `Mode` containing the solved electromagnetic or acoustic modes.
         
         :rtype: numarray(Mode)
         '''
      if not len(self.mode_set):
        for m in range(self.num_modes):
          if self.is_EM():
            mode=ModeEM(self, m)
          else:
            mode=ModeAC(self, m)
          self.mode_set.append(mode)
        
      return self.mode_set

    def set_r0_offset(self, rx, ry): # this is clumsy and only works if called after modes have been calced.
      for m in self.get_modes(): m.set_r0_offset(rx, ry)

    def symmetry_classification(self, m):
      '''If the point group of the structure has been specified, returns the symmetry class of the given mode.
         
         :param int m: Index of the mode of interest.
         :rtype: PointGroup
         '''
      if self.point_group == PointGroup.Unknown: return ''
      return '{0}:{1}'.format(self.point_group.name, self.sym_reps[m].name)

    def neff(self, m): 
      ''' Returns the effective index of EM mode `m`.
  
      :param int m: Index of the mode of interest.
      :rtype: float
      ''' 
      assert(self.is_EM())
      return np.real(self.Eig_values[m]*self.wl_m/(2*np.pi))

    def ngroup_EM_available(self): 
      '''Returns true if a measure of the electromagnetic group index is available.'''
      return not(self.EM_mode_energy is None or self.EM_mode_power is None)

    def vgroup_AC_available(self): 
      '''Returns true if a measure of the acoustic group velocity is available.'''
      return not(self.AC_mode_energy is None or self.AC_mode_power is None)

    def ngroup_EM(self, m):
      '''Returns the group index of electromagnetic mode `m`, if available, otherwise returns zero with a warning message.
         
         :param int m: Index of the mode of interest.
         :return: Group index of the mode.
         :rtype: float
         '''
      if not self.ngroup_EM_available():
        print('''EM group index requires calculation of mode energy and mode power when calculating EM modes. 
               Set calc_EM_mode_energy=True and calc_AC_mode_power=True in call to Simmo''')
        return 0
      vg= np.real(self.EM_mode_power[m]/self.EM_mode_energy[m])
      ng=VacCSpeed/vg
      return ng

    def ngroup_EM_all(self):
      '''Returns a numpy array of the group index of all electromagnetic modes, if available, 
         otherwise returns a zero numarray with a warning message.
         
         :return: numpy array of  index of the mode.
         :rtype: array(float)
         '''
      if not self.ngroup_EM_available():
        print('''EM group index requires calculation of mode energy and mode power when calculating EM modes. 
               Set calc_EM_mode_energy=True in call to calc_EM_modes''')
        return np.zeros(len(self.Eig_values), dtype=float)
      vg= np.real(self.EM_mode_power/self.EM_mode_energy)
      ng=VacCSpeed/vg
      return ng

    def kz_EM(self, m): 
      '''Returns the wavevector in 1/m of electromagnetic mode `m`.
         
         :param int m: Index of the mode of interest.
         :return: Wavevector k in 1/m.
         :rtype: float
         '''
      assert(self.is_EM())
      return self.Eig_values[m]

    def nu_AC(self, m): 
      '''Returns the frequency in Hz of acoustic mode `m`.
         
         :param int m: Index of the mode of interest.
         :return: Frequency :math:`\\nu` in Hz
         :rtype: float
         '''
      assert(self.is_AC())
      return self.Eig_values[m]

    def Omega_AC(self, m): 
      '''Returns the frequency in 1/s of acoustic mode `m`.
         
         :param int m: Index of the mode of interest.
         :return: Angular requency :math:`\\Omega` in Hz
         :rtype: float
         '''
      return self.Eig_values[m]*2*np.pi

    def neff_all(self): 
      ''' Return an array of the effective index of all electromagnetic modes.

         :return: numpy array of effective indices
         :rtype: array(float)
          '''
      assert(self.is_EM())
      return np.real(self.Eig_values*self.wl_m/(2*np.pi))

    def kz_EM_all(self): 
      ''' Return an array of the wavevector in 1/m of all electromagnetic modes.

         :return: numpy array of wavevectors in 1/m
         :rtype: array(float)
          '''
      assert(self.is_EM())
      return self.Eig_values

    def nu_AC_all(self): 
      ''' Return an array of the frequency in Hz of all acoustic modes.

         :return: numpy array of frequencies in Hz
         :rtype: array(float)
         '''
      assert(self.is_AC())
      return self.Eig_values

    def Omega_AC_all(self, m): 
      ''' Return an array of the angular frequency in 1/s of all acoustic modes.

         :return: numpy array of angular frequencies in 1/s
         :rtype: array(float)
         '''
      assert(self.is_AC())
      return self.Eig_values*2*np.pi

    def vp_AC(self, m): 
      """ Return phase velocity of all AC modes in m/s"""
      ''' Return the phase velocity in m/s of acoustic mode `m`.

         :return: Phase velocity of acoustic mode `m` in m/s
         :rtype: float
         '''
      assert(self.is_AC())
      return np.pi*2*np.real(self.Eig_values[m])/self.k_AC

    def vp_AC_all(self): 
      ''' Return an array of the phase velocity in m/s of all acoustic modes.

         :return: numpy array of elastic phase velocities in m/s
         :rtype: array(float)
         '''
      assert(self.is_AC())
      return np.pi*2*np.real(self.Eig_values)/self.k_AC

    def vg_AC(self, m): 
      """ Return group velocity of AC mode m in m/s"""
      if self.AC_mode_energy is None or self.AC_mode_power is None:
        print('''AC group velocity requires calculation of mode energy and mode power when calculating AC modes. 
               Set calc_AC_mode_power=True in call to calc_AC_modes''')
        return 0
      vg= np.real(self.AC_mode_power[m]/self.AC_mode_energy[m])
      return vg

    def vg_AC_all(self): 
      assert(self.is_AC())
      """ Return group velocity of all AC modes in m/s"""
      if self.AC_mode_energy is None or self.AC_mode_power is None:
        print('''AC group velocity requires calculation of mode energy and mode power when calculating AC modes. 
               Set calc_AC_mode_power=True in call to calc_AC_modes''')
        return np.zeros(len(self.Eig_values), dtype=float)
      vg= np.real(self.AC_mode_power/self.AC_mode_energy)
      return vg



    def alpha_t_AC(self, m): 
      assert(self.is_AC())
      return self.ac_alpha_t[m]

    def alpha_t_AC_all(self): 
      assert(self.is_AC())
      return self.ac_alpha_t

    def alpha_s_AC(self, m): # spatial loss [1/m]  #TODO: surely this should be the group velocity for scaling between spatial and temporal decays #Which is a problem because it requires knowing vg
      assert(self.is_AC())
      print('alpha_s:')
      print(self.alpha_t_AC(m))
      print(self.vg_AC(m))
      return self.alpha_t_AC(m)/self.vg_AC(m)

    def alpha_s_AC_all(self): # spatial loss [1/m]
      assert(self.is_AC())
      return self.alpha_t_AC_all/self.vg_AC_all

    def Qmech_AC(self, m): 
      assert(self.is_AC())
      return self.ac_Qmech[m]

    def Qmech_AC_all(self):
      assert(self.is_AC())
      return self.ac_Qmech

    def linewidth_AC(self, m):
      assert(self.is_AC())
      return self.ac_linewidth[m]

    def linewidth_AC_all(self):
      assert(self.is_AC())
      return self.ac_linewidth

    def analyse_symmetries(self, ptgrp):
      self.point_group=ptgrp
      symlist = integration.symmetries(self)
      self.sym_reps = []
      if ptgrp == PointGroup.C2V:
        for m, sym in enumerate(symlist):
          if sym == (1,1,1):     self.sym_reps.append(SymRep.A)
          elif sym == (-1,1,-1): self.sym_reps.append(SymRep.B1)
          elif sym == (1,-1,-1): self.sym_reps.append(SymRep.B2)
          elif sym == (-1,-1,1): self.sym_reps.append(SymRep.B3)
          else: 
            print('Warning: Unknown symmetry pattern', sym)
            self.sym_reps.append(SymRep.Unknown)


      else:
        print("unknown symmetry properties in mode_calcs")

    def calc_acoustic_losses(self, fixed_Q=None): # TODO: make sure this is not done more than once on the same Simmo
      alpha = None
      if fixed_Q is None: 
        self.Q_method=QAcMethod.Intrinsic

        # Calc alpha (loss) Eq. 45
        print("Acoustic loss calc")
        nnodes=6 # TODO: is this right?
        try:
            if self.EM_sim.structure.inc_shape in self.EM_sim.structure.linear_element_shapes:
                alpha = NumBAT.ac_alpha_int_v2(self.num_modes,
                    self.n_msh_el, self.n_msh_pts, nnodes,
                    self.table_nod, self.type_el, self.x_arr,
                    self.structure.nb_typ_el_AC, self.structure.eta_tensor,
                    self.k_AC, self.Omega_AC, self.sol1,
                    # sim_AC.AC_mode_power) # appropriate for alpha in [1/m]
                    self.AC_mode_energy) # appropriate for alpha in [1/s]
            else:
                if self.EM_sim.structure.inc_shape not in self.EM_sim.structure.curvilinear_element_shapes:
                    print("Warning: ac_alpha_int - not sure if mesh contains curvi-linear elements", 
                        "\n using slow quadrature integration by default.\n\n")
                Fortran_debug=0
                overlap=np.zeros(self.num_modes, dtype=complex)  # not sure why this is needed by ac_alpha_int
                alpha = NumBAT.ac_alpha_int(self.num_modes,
                    self.n_msh_el, self.n_msh_pts, nnodes,
                    self.table_nod, self.type_el, self.x_arr,
                    self.structure.nb_typ_el_AC, self.structure.eta_tensor,
                    self.k_AC, self.Omega_AC, self.sol1,
                    # sim_AC.AC_mode_power, Fortran_debug) # appropriate for alpha in [1/m]
                    self.AC_mode_energy, Fortran_debug) # appropriate for alpha in [1/s]
        except KeyboardInterrupt:
            print("\n\n Routine ac_alpha_int interrupted by keyboard.\n\n")
        self.ac_alpha_t = np.real(alpha)
        # Q_factors = 0.5*(k_AC/alpha)*np.ones(num_modes_AC) # appropriate for alpha in [1/m]
        self.ac_Qmech = 0.5*(np.real(self.Omega_AC)/self.ac_alpha_t)*np.ones(self.num_modes) # appropriate for alpha in [1/s]
      else:
        self.Q_method=QAcMethod.Fixed
        # factor of a 1/2 because alpha is for power!
        # alpha [1/m] = Omega_AC/(2*vg*fixed_Q) = k_AC/fixed_Q
        # alpha [1/s] = vg * alpha [1/m]
        # alpha [1/s] = Omega_AC/(2*fixed_Q)
        # alpha = 0.5*(k_AC/fixed_Q)*np.ones(num_modes_AC) # appropriate for alpha in [1/m]
        self.ac_Qmech = fixed_Q*np.ones(self.num_modes)
        self.ac_alpha_t = 0.5*(np.real(self.Omega_AC)/fixed_Q)*np.ones(self.num_modes) # appropriate for alpha in [1/s]

      self.ac_linewidth = self.ac_alpha_t/np.pi # SBS linewidth of each resonance in [Hz]   #TODO: not sure about the 1/pi.  
                                                                                            #If linewdith should be amplitude rate in Hz, wouldn't it be
                                                                                            # alpha/(2 * 2pi)  since alpha is a power decay rate

    def calc_EM_modes(self):
        """ Run a Fortran FEM calculation to find the optical modes.

        Returns a ``Simmo`` object that has these key values:

        Eig_values: a 1d array of Eigenvalues (propagation constants) in [1/m]

        sol1: the associated Eigenvectors, ie. the fields, stored as [field comp, node nu on element, Eig value, el nu]

        EM_mode_power: the power in the optical modes. Note this power is negative for modes travelling in the negative
                       z-direction, eg the Stokes wave in backward SBS.
        """
        self.EM_AC = 'EM'

        self.d_in_m = self.structure.unitcell_x*1e-9
        n_list = []
        n_list_tmp = np.array([self.structure.material_bkg.n, 
                               self.structure.material_a.n, self.structure.material_b.n, self.structure.material_c.n,
                               self.structure.material_d.n, self.structure.material_e.n, self.structure.material_f.n,
                               self.structure.material_g.n, self.structure.material_h.n, self.structure.material_i.n,
                               self.structure.material_j.n, self.structure.material_k.n, self.structure.material_l.n,
                               self.structure.material_m.n, self.structure.material_n.n, self.structure.material_o.n,
                               self.structure.material_p.n, self.structure.material_q.n, self.structure.material_r.n])
        self.el_conv_table_n = {}
        i = 1; j = 1
        for n in n_list_tmp:
            if n != 0:
                n_list.append(n)
                self.el_conv_table_n[i] = j
                j += 1
            i += 1
        self.n_list = np.array(n_list)
        n_list = None

        if self.structure.loss is False:
            self.n_list = self.n_list.real

        if self.num_modes < 20:
            self.num_modes = 20
            print("Warning: ARPACK needs >= 20 modes so set num_modes=20.")

        # Parameters that control how FEM routine runs
        self.E_H_field = 1  # Selected formulation (1=E-Field, 2=H-Field)
        i_cond = 2  # Boundary conditions (0=Dirichlet,1=Neumann,2=unitcell_x)
        itermax = 30  # Maximum number of iterations for convergence
        EM_FEM_debug = self.debug  # Fortran routines will display & save add. info

        # Calculate where to center the Eigenmode solver around.
        # (Shift and invert FEM method)
        shift = self.n_eff**2 * self.k_0**2

        if EM_FEM_debug == 1:
            if not os.path.exists("Normed"):
                os.mkdir("Normed")
            if not os.path.exists("Matrices"):
                os.mkdir("Matrices")
            if not os.path.exists("Output"):
                os.mkdir("Output")

        with open(self.structure.mesh_file) as f:
            self.n_msh_pts, self.n_msh_el = [int(i) for i in f.readline().split()]

        # Size of Fortran's complex superarray (scales with mesh)
        int_max, cmplx_max, real_max = NumBAT.array_size(self.n_msh_el, self.num_modes)
        if EM_FEM_debug == 1:
          print("Mesh calculated: %d nodes."%self.n_msh_el)

        try:
            resm = NumBAT.calc_em_modes(
                self.wl_m, self.num_modes,
                EM_FEM_debug, self.structure.mesh_file, self.n_msh_pts,
                self.n_msh_el, self.structure.nb_typ_el, self.n_list,
                self.k_pll, self.d_in_m, shift, self.E_H_field, i_cond, itermax,
                self.structure.plotting_fields, self.structure.plot_real,
                self.structure.plot_imag, self.structure.plot_abs,
                cmplx_max, real_max, int_max)

            self.Eig_values, self.sol1, self.mode_pol, self.table_nod, \
            self.type_el, self.type_nod, self.x_arr, self.ls_material = resm

        except KeyboardInterrupt:
            print("\n\n FEM routine calc_EM_modes",\
            "interrupted by keyboard.\n\n")

        # if not self.structure.plot_field_conc:
        #     self.mode_pol = None

        # if self.structure.plotting_fields != 1:
        #     self.sol1 = None
        #     self.n_list = None
        #     self.E_H_field = None
        #     self.table_nod = None
        #     self.type_el = None
        #     self.x_arr = None
        #     self.n_msh_pts = None
        #     self.n_msh_el = None

        if self.structure.plt_mesh:
            print("Suppressed inefficient matplotlib plotting of mesh...")
            #plotting.plot_msh(self.x_arr, prefix_str=self.structure.mesh_file, suffix_str='_EM')


### Calc unnormalised power in each EM mode Kokou equiv. of Eq. 8.
        try:
            print("Calculating EM mode powers...")
            nnodes = 6
            if self.structure.inc_shape in self.structure.linear_element_shapes:
            ## Integration using analytically evaluated basis function integrals. Fast.
                self.EM_mode_power = NumBAT.em_mode_energy_int_v2_ez(
                    self.k_0, self.num_modes, self.n_msh_el, self.n_msh_pts,
                    nnodes, self.table_nod,
                    self.x_arr, self.Eig_values, self.sol1)
            else:
                if self.structure.inc_shape not in self.structure.curvilinear_element_shapes:
                    print("Warning: em_mode_energy_int - not sure if mesh contains curvi-linear elements", 
                        "\n using slow quadrature integration by default.\n\n")
            # Integration by quadrature. Slowest.
                self.EM_mode_power = NumBAT.em_mode_energy_int_ez(
                    self.k_0, self.num_modes, self.n_msh_el, self.n_msh_pts,
                    nnodes, self.table_nod,
                    self.x_arr, self.Eig_values, self.sol1)
            # Bring Kokou's def into line with CW formulation.
            self.EM_mode_power = 2.0*self.EM_mode_power

        except KeyboardInterrupt:
            print("\n\n FEM routine EM_mode_energy_int interrupted by keyboard.\n\n")


### Calc energy (not power) in each EM mode - PRA Eq. 6.
        if self.calc_EM_mode_energy is True:
            print("Calculating EM mode energies...")
            try:
                nnodes = 6
                # import time
                # start = time.time()
                if self.structure.inc_shape in self.structure.linear_element_shapes:
                # # Semi-analytic integration. Fastest!
                # else:
                #     if self.structure.inc_shape not in self.structure.curvilinear_element_shapes:
                #         print("Warning: em_mode_e_energy_int - not sure if mesh contains curvi-linear elements", 
                #             "\n using slow quadrature integration by default.\n\n")
                # # Integration by quadrature. Slowest.
                    self.EM_mode_energy = NumBAT.em_mode_e_energy_int(
                        self.num_modes, self.n_msh_el, self.n_msh_pts, nnodes,
                        self.table_nod, self.type_el, self.structure.nb_typ_el, self.n_list,
                        self.x_arr, self.sol1)
                else:
                  print("\n\n FEM routine em_mode_e_energy_int needs work for this structure .\n\n")
                  self.EM_mode_energy=np.zeros(self.num_modes, dtype=float)
                  
            except KeyboardInterrupt:
                print("\n\n FEM routine em_mode_e_energy_int interrupted by keyboard.\n\n")

        # This group velocity calc is not accurate in the presence of dispersion!
        # self.group_velocity_EM = self.EM_mode_power/self.EM_mode_power_energy

        # If considering a the backwards propagating Stokes field.
        if self.Stokes == True:
            self.Eig_values = -1*self.Eig_values
            self.sol1 = np.conj(self.sol1)

        ### Not necessary because EM FEM mesh always normalised in area to unity.
        # print area
        # x_tmp = []
        # y_tmp = []
        # for i in np.arange(self.n_msh_pts):
        #     x_tmp.append(self.x_arr[0,i])
        #     y_tmp.append(self.x_arr[1,i])
        # x_min = np.min(x_tmp); x_max=np.max(x_tmp)
        # y_min = np.min(y_tmp); y_max=np.max(y_tmp)
        # area = abs((x_max-x_min)*(y_max-y_min))
        # print area
        # self.EM_mode_power = self.EM_mode_power*area


#def calc_EM_mode_energy(self):  # these require extraction of numerical props from Fortran. Clean that up first.
#      assert(self.is_EM())
#      if not self.EM_mode_energy is None: return  # nothing to do

#    def calc_EM_mode_power(self):
#      assert(self.is_EM())
#      if not self.EM_mode_power is None: return  # nothing to do

#    def calc_EM_mode_power_and_energy(self):
#      assert(self.is_EM())
#      self.calc_EM_mode_power()
#      self.calc_EM_mode_energy()


    def calc_AC_modes(self):
        """ Run a Fortran FEM calculation to find the acoustic modes.

        Returns a ``Simmo`` object that has these key values:

        Eig_values: a 1d array of Eigenvalues (frequencies) in [1/s]

        sol1: the associated Eigenvectors, ie. the fields, stored as
               [field comp, node nu on element, Eig value, el nu]

        AC_mode_energy: the elastic power in the acoutic modes.
        """
        self.EM_AC = 'AC'

        self.d_in_m = self.structure.inc_a_x*1e-9

        el_conv_table = {}
        i = 1; j = 1
        for matter in self.structure.acoustic_props_tmp:
            if matter.s != None:
                el_conv_table[i] = j
                j += 1
            i += 1
        final_dict = {}
        for entry in el_conv_table:
            # print entry, self.EM_sim.el_conv_table_n[entry], el_conv_table[entry]
            final_dict[self.EM_sim.el_conv_table_n[entry]] = el_conv_table[entry]
        # print final_dict
        self.typ_el_AC = final_dict

        if self.num_modes < 20:
            self.num_modes = 20
            print("Warning: ARPACK needs >= 20 modes so set num_modes=20.")

        # Parameters that control how FEM routine runs
        i_cond = 1  # Boundary conditions (0=Dirichlet,1=Neumann,2=unitcell_x)
        itermax = 30  # Maximum number of iterations for convergence
        AC_FEM_debug = 0  # Fortran routines will display & save add. info
        ARPACK_tol = 1e-10  # ARPACK accuracy (0.0 for machine precision)

        # Calculate where to center the Eigenmode solver around.
        # (Shift and invert FEM method)
        if self.shift_Hz is None:
            # For AC problem shift is a frequency; [shift] = s^-1.
            v_list = []
            for el in range(self.structure.nb_typ_el_AC):
                # Using acoustic velocity of longitudinal mode pg 215 Auld vol 1.
                v_list.append(np.sqrt(self.structure.c_tensor[0,0][el]/self.structure.rho[el]))
                # # Using acoustic velocity of shear mode pg 215 Auld vol 1.
                # v_list.append(np.sqrt(self.structure.c_tensor[3,3][el]/self.structure.rho[el]))
            AC_velocity = np.real(v_list).min()
            shift = np.real(AC_velocity*self.k_AC/(2.*np.pi))
            # print "shift", shift
            shift = 0.9*shift
            # print "shift", shift
        else:
            shift = self.shift_Hz

        # Take existing msh from EM FEM and manipulate mesh to exclude vacuum areas.
        if self.EM_sim:
            suplied_geo_flag = 1
            n_msh_el = self.EM_sim.n_msh_el
            n_msh_pts = self.EM_sim.n_msh_pts
            type_el = self.EM_sim.type_el
            type_nod = self.EM_sim.type_nod
            table_nod = self.EM_sim.table_nod
            x_arr = self.EM_sim.x_arr
            n_el_kept = 0
            n_msh_pts_AC = 0
            type_el_AC = []
            table_nod_AC_tmp = np.zeros(np.shape(table_nod))
            el_convert_tbl = {}
            el_convert_tbl_inv = {}
            node_convert_tbl = {}
            if self.structure.plt_mesh:
                plotting.plot_msh(x_arr, prefix_str=self.structure.mesh_file, suffix_str='_AC-orig')

            for el in range(n_msh_el):
                # print type_el[el]
                if type_el[el] in self.typ_el_AC:
                    # print "in", type_el[el]
                    type_el_AC.append(self.typ_el_AC[type_el[el]])
                    el_convert_tbl[n_el_kept] = el
                    el_convert_tbl_inv[el] = n_el_kept
                    for i in range(6):
                        # Leaves node numbering untouched
                        table_nod_AC_tmp[i][n_el_kept] = table_nod[i][el]
                    n_el_kept += 1
            n_msh_el_AC = n_el_kept
            # Find unique nodes
            node_lst_tmp = []
            for el in range(n_msh_el_AC):
                for i in range(6):
                    node_lst_tmp.append(table_nod_AC_tmp[i][el])
            unique_nodes = list(set(node_lst_tmp))
            n_msh_pts_AC = len(unique_nodes)
            unique_nodes = [int(j) for j in unique_nodes]
            # Mapping unique nodes to start from zero
            for i in range(n_msh_pts_AC):
                node_convert_tbl[unique_nodes[i]] = i
            # Creating finalised table_nod.
            table_nod_AC = []
            for i in range(6):
                el_tbl = []
                for el in range(n_msh_el_AC):
                    # Note table_nod needs to be adjust back to fortran indexing
                    el_tbl.append(node_convert_tbl[table_nod_AC_tmp[i][el]]+1)
                table_nod_AC.append(el_tbl)
            # Find the coordinates of chosen nodes.
            x_arr_AC = np.zeros((2,n_msh_pts_AC))
            for node in unique_nodes:
                # Note x_arr needs to be adjust back to fortran indexing
                x_arr_AC[0,node_convert_tbl[node]] = (x_arr[0,node-1])
                x_arr_AC[1,node_convert_tbl[node]] = (x_arr[1,node-1])

            self.el_convert_tbl = el_convert_tbl
            self.el_convert_tbl_inv = el_convert_tbl_inv
            self.node_convert_tbl = node_convert_tbl


            ### AC FEM uses Neumann B.C.s so type_nod is totally irrelevant!
            # # Find nodes on boundaries of materials
            # node_array = -1*np.ones(n_msh_pts)
            # interface_nodes = []
            # for el in range(n_msh_el):
            #     for i in range(6):
            #         node = table_nod[i][el]
            #         # Check if first time seen this node
            #         if node_array[node - 1] == -1: # adjust to python indexing
            #             node_array[node - 1] = type_el[el]
            #         else:
            #             if node_array[node - 1] != type_el[el]:
            #                 interface_nodes.append(node)
            # interface_nodes = list(set(interface_nodes))
            type_nod_AC = np.zeros(n_msh_pts_AC)

            # import matplotlib
            # matplotlib.use('pdf')
            # import matplotlib.pyplot as plt
            # plt.clf()
            # plt.figure(figsize=(13,13))
            # ax = plt.subplot(1,1,1)
            # for node in unique_nodes:
            #     if node in interface_nodes:
            #         type_nod_AC[node_convert_tbl[node]] = i_cond
            #         plt.plot(x_arr_AC[0,node_convert_tbl[node]], x_arr_AC[1,node_convert_tbl[node]], 'ok')
            # ax.set_aspect('equal')
            # plt.savefig('boundary.pdf', bbox_inches='tight')
            self.n_msh_pts = n_msh_pts_AC
            self.n_msh_el = n_msh_el_AC
        # Default, indicates to use geometry subroutine in FEM routine.
        else:
            suplied_geo_flag = 0
            with open("../backend/fortran/msh/"+self.structure.mesh_file) as f:
                self.n_msh_pts, self.n_msh_el = [int(i) for i in f.readline().split()]
            table_nod_AC = np.zeros((6, self.n_msh_el))
            type_el_AC = np.zeros(self.n_msh_el)
            x_arr_AC = np.zeros((2,self.n_msh_pts))
            type_nod_AC = np.zeros(self.n_msh_pts)

        if AC_FEM_debug == 1:
            print('shift', shift)
            if not os.path.exists("Normed"):
                os.mkdir("Normed")
            if not os.path.exists("Output"):
                os.mkdir("Output")
            if not os.path.exists("Matrices"):
                os.mkdir("Matrices")

        # Size of Fortran's complex superarray (scales with mesh)
        int_max, cmplx_max, real_max = NumBAT.array_size(self.n_msh_el, self.num_modes)

        try:
            resm = NumBAT.calc_ac_modes(
                self.k_AC, self.num_modes,
                AC_FEM_debug, self.structure.mesh_file, self.n_msh_pts,
                self.n_msh_el, self.structure.nb_typ_el_AC,
                self.structure.c_tensor, self.structure.rho,
                self.d_in_m, shift, i_cond, itermax, ARPACK_tol,
                self.structure.plotting_fields,
                cmplx_max, real_max, int_max, suplied_geo_flag, type_nod_AC, 
                self.structure.symmetry_flag, table_nod_AC, type_el_AC, x_arr_AC)
            table_nod_out, type_el_out, x_arr_out, \
            self.Eig_values, self.sol1, self.mode_pol = resm

            # FEM Eigenvalue is frequency, rather than angular frequency Omega
            self.Omega_AC = self.Eig_values*2*np.pi

        except KeyboardInterrupt:
            print("\n\n FEM routine calc_AC_modes",\
            "interrupted by keyboard.\n\n")

        # Retrieve the material properties of each mesh point.
        self.ls_material = NumBAT.array_material_ac(self.n_msh_pts, self.n_msh_el,
             self.structure.nb_typ_el_AC, type_el_AC,
             self.structure.rho, self.structure.c_tensor, 
             self.structure.p_tensor, self.structure.eta_tensor)

        if self.structure.plt_mesh:
            plotting.plot_msh(x_arr_AC, prefix_str=self.structure.mesh_file, suffix_str='_AC-in')
            plotting.plot_msh(x_arr_out, prefix_str=self.structure.mesh_file, suffix_str='_AC-out')

        # if self.EM_sim is None:
        #     table_nod_out = None
        #     type_el_out = None
        #     x_arr_out = None
        #     self.table_nod = table_nod_AC
        #     self.type_el = type_el_AC
        #     self.x_arr = x_arr_AC
        # else:
        self.table_nod = table_nod_out
        self.type_el = type_el_out
        self.x_arr = x_arr_out

### Calc unnormalised power in each AC mode - PRA Eq. 18.
        if self.calc_AC_mode_power is True:
            try:
                nnodes = 6
                if self.structure.inc_shape in self.structure.linear_element_shapes:
                # Semi-analytic integration following KD 9/9/16 notes. Fastest!
                    self.AC_mode_power = NumBAT.ac_mode_power_int_v4(
                        self.num_modes, self.n_msh_el, self.n_msh_pts,
                        nnodes, self.table_nod, self.type_el, self.x_arr,
                        self.structure.nb_typ_el_AC, self.structure.c_tensor,
                        self.k_AC, self.Omega_AC, self.sol1)
                else:
                    if self.structure.inc_shape not in self.structure.curvilinear_element_shapes:
                        print("Warning: ac_mode_power_int - not sure if mesh contains curvi-linear elements", 
                            "\n using slow quadrature integration by default.\n\n")
                # Integration by quadrature. Slowest.
                    self.AC_mode_power = NumBAT.ac_mode_power_int(
                        self.num_modes, self.n_msh_el, self.n_msh_pts,
                        nnodes, self.table_nod, self.type_el, self.x_arr,
                        self.structure.nb_typ_el_AC, self.structure.c_tensor_z,
                        self.k_AC, self.Omega_AC, self.sol1, AC_FEM_debug)
            except KeyboardInterrupt:
                print("\n\n FEM routine AC_mode_energy_int interrupted by keyboard.\n\n")


### Calc unnormalised elastic energy in each AC mode - PRA Eq. 16.
        try:
            nnodes = 6
            if self.structure.inc_shape in self.structure.linear_element_shapes:
            # Semi-analytic integration. Fastest!
                self.AC_mode_energy= NumBAT.ac_mode_elastic_energy_int_v4(
                    self.num_modes, self.n_msh_el, self.n_msh_pts,
                    nnodes, self.table_nod, self.type_el, self.x_arr,
                    self.structure.nb_typ_el_AC, self.structure.rho,
                    self.Omega_AC, self.sol1)
            else:
                if self.structure.inc_shape not in self.structure.curvilinear_element_shapes:
                    print("Warning: ac_mode_elastic_energy_int - not sure if mesh contains curvi-linear elements", 
                        "\n using slow quadrature integration by default.\n\n")
            # Integration by quadrature. Slowest.
                self.AC_mode_energy= NumBAT.ac_mode_elastic_energy_int(
                    self.num_modes, self.n_msh_el, self.n_msh_pts,
                    nnodes, self.table_nod, self.type_el, self.x_arr,
                    self.structure.nb_typ_el_AC, self.structure.rho,
                    self.Omega_AC, self.sol1, AC_FEM_debug)
        except KeyboardInterrupt:
            print("\n\n FEM routine AC_mode_elastic_energy_int interrupted by keyboard.\n\n")


def bkwd_Stokes_modes(EM_sim):
    """ Defines the backward travelling Stokes waves as the conjugate
        of the forward travelling pump waves.

    Returns a ``Simmo`` object that has these key values:

    Eig_values: a 1d array of Eigenvalues (propagation constants) in [1/m]

    sol1: the associated Eigenvectors, ie. the fields, stored as
           [field comp, node nu on element, Eig value, el nu]

    EM_mode_power: the power in the Stokes modes. Note this power is negative because the modes 
                   are travelling in the negative z-direction.
    """
    Stokes_modes = copy.deepcopy(EM_sim)
    Stokes_modes.sol1 = np.conj(Stokes_modes.sol1)
    Stokes_modes.Eig_values = -1.0*Stokes_modes.Eig_values
    Stokes_modes.EM_mode_power = -1.0*Stokes_modes.EM_mode_power
    return Stokes_modes



def fwd_Stokes_modes(EM_sim):
    """ Defines the forward travelling Stokes waves as a copy
        of the forward travelling pump waves.

    Returns a ``Simmo`` object that has these key values:

    """
    Stokes_modes = copy.deepcopy(EM_sim)
    return Stokes_modes
