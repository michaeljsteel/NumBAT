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

import plotting
import integration
from nbtypes import *
from fortran import NumBAT


VacCSpeed=299792458

class Simmo(object):
    """ Calculates the modes of a ``Struct`` object at a wavelength of wl_nm.
    """
    def __init__(self, structure, num_modes=20, wl_nm=1, n_eff=None, shift_Hz=None, 
                 k_AC=None, EM_sim=None, Stokes=False, 
                 calc_EM_mode_energy=False, calc_AC_mode_power=False, debug=False):
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
        self.omega_EM = 2*np.pi*speed_c/self.wl_m # Angular freq in units of Hz
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

    def is_EM(self): return self.EM_AC == 'EM'
    def is_AC(self): return self.EM_AC != 'EM'

    def symmetry_classification(self, m):
      if self.point_group == PointGroup.Unknown: return ''
      return '{0}:{1}'.format(self.point_group.name, self.sym_reps[m].name)

    def neff(self, m): 
      """ Return effective index of EM mode m""" 
      assert(self.is_EM())

      return np.real(self.Eig_values[m]*self.wl_m/(2*np.pi))

    def ngroup_EM_available(self): return not(self.EM_mode_energy is None or self.EM_mode_power is None)
    def vgroup_AC_available(self): return not(self.AC_mode_energy is None or self.AC_mode_power is None)

    def ngroup_EM(self, m):
      if self.EM_mode_energy is None or self.EM_mode_power is None:
        print('''EM group index requires calculation of mode energy and mode power when calculating EM modes. 
               Set calc_EM_mode_energy=True and calc_AC_mode_power=True in call to Simmo''')
        return 0
      vg= np.real(self.EM_mode_power[m]/self.EM_mode_energy[m])
      ng=VacCSpeed/vg
      return ng

    def ngroup_EM_all(self):
      if self.EM_mode_energy is None or self.EM_mode_power is None:
        print('''EM group index requires calculation of mode energy and mode power when calculating EM modes. 
               Set calc_EM_mode_energy=True in call to calc_EM_modes''')
        return np.zeros(len(self.Eig_values), dtype=float)
      vg= np.real(self.EM_mode_power/self.EM_mode_energy)
      ng=VacCSpeed/vg
      return ng

    def kz_EM(self, m): 
      """ Return wavevector of EM mode m in 1/m""" 
      assert(self.is_EM())
      return self.Eig_values[m]

    def nu_AC(self, m): 
      """ Return frequency of AC mode m in Hz"""
      assert(self.is_AC())
      return self.Eig_values[m]

    def om_AC(self, m): 
      """ Return angular frequency of AC mode m in 1/s"""
      return self.Eig_values[m]

    def neff_all(self): 
      """ Return effective index of all EM modes """
      assert(self.is_EM())
      return np.real(self.Eig_values*self.wl_m/(2*np.pi))

    def kz_EM_all(self): 
      """ Return wavevector of all EM modes in 1/m"""
      assert(self.is_EM())
      return self.Eig_values

    def nu_AC_all(self): 
      """ Return frequency of all AC modes in Hz"""
      assert(self.is_AC())
      return self.Eig_values

    def om_AC_all(self, m): 
      """ Return angular frequency of all AC modes in 1/s"""
      assert(self.is_AC())
      return self.Eig_values*2*np.pi

    def vp_AC(self, m): 
      """ Return phase velocity of all AC modes in m/s"""
      assert(self.is_AC())
      return np.pi*2*np.real(self.Eig_values[m])/self.k_AC

    def vp_AC_all(self): 
      """ Return phase velocity of all AC modes in m/s"""
      assert(self.is_AC())
      return np.pi*2*np.real(self.Eig_values)/self.k_AC

    def vg_AC(self, m): 
      """ Return group velocity of AC mode m in m/s"""
      if self.AC_mode_energy is None or self.AC_mode_power is None:
        print('''AC group velocity requires calculation of mode energy and mode power when calculating AC modes. 
               Set calc_AC_mode_power=True in call to calc_AC_modes''')
        return np.zeros(len(self.Eig_values), dtype=float)
      vg= np.real(self.AC_mode_power[m]/self.AC_mode_energy[m])
      return vg

    def vg_AC_all(self): 
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

    def alpha_s_AC(self, m): # spatial loss [1/m]
      assert(self.is_AC())
      return self.alpha_t_AC(m)/self.vp_AC(m)

    def alpha_s_AC_all(self): # spatial loss [1/m]
      assert(self.is_AC())
      return self.alpha_t_AC_all/self.vp_AC_all

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

      self.ac_linewidth = self.ac_alpha_t/np.pi # SBS linewidth of each resonance in [Hz]

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
