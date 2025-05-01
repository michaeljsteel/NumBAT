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



import copy
import numpy as np

from numbattools import process_fortran_return


from nbtypes import (
    QAcMethod,
    SI_nm,
    twopi,
)

from femmesh import FemMesh

from fortran import nb_fortran

from simresult import EMSimResult, ACSimResult

class Simulation:
    """Class for calculating the electromagnetic and/or acoustic modes of a ``Struct`` object."""

    def __init__(self, structure, num_modes, debug):
        """Sets up the problem for the mode calculation at a given optical wavelength `wl_nm` or acoustic wavenumber `q_AC`.

        For electromagnetic problems, the tool solves for the effective index or wavenumber at a given wavelength.
        For acoustic problems, the tool solves for the acoustic frequency at a given wavenumber.

          :param Simulation structure: The waveguide structure to be solved.
          :param int num_modes: The number of modes to be found.
          :param float wl_nm: For electromagnetic problems, the vacuum wavelength in nanometers.
          :param float n_eff: For electromagnetic problems, an estimated effective index to begin the eigenvalue search.
          :param float shift_Hz: For acoustic problems, an estimated frequency offset to begin the eigenvalue search.
          :param float q_AC: For acoustic problems, the acoustic wavevector of the mode.
          :param float simres_EM: For acoustic problems, the results of a previously solved EM problem to speed calculations.
          :param bool calc_EM_mode_energy: For electromagnetic problems, whether to calculate the optical mode energy.
          :param bool calc_AC_mode_power: For acoustic problems, whether to calculate the acoustic mode power.
        """

        self.structure = structure
        self.mode_plot_helper = None
        self.n_modes = num_modes

        self.d_in_m = (
            self.structure.domain_x * SI_nm
        )  # Scales Gmsh structure to correct size

        # seems unused
        self.mode_pol = None
        self.debug = debug

        # just off normal incidence to avoid degeneracies
        self.k_perp = np.array([1e-12, 1e-12])

        # Some of these values depend on whether EM or AC simulation because domain vacuum trimming

        # Keep linter quiet
        self.fem_mesh = None
        self.simres_EM = None
        self.sim_result = None

    def get_sim_result(self):
        return self.sim_result

    @staticmethod
    def load_simulation(prefix):  # Need to load type and instantiate correct class
        npzfile = np.load(prefix + ".npz", allow_pickle=True)
        return npzfile["simulation"].tolist()

    def clean_for_pickle(self):
        if self.mode_plot_helper is not None:
            self.mode_plot_helper.cleanup()

    def save_simulation(self, prefix):
        self.clean_for_pickle()

        # Acoustic sims can contain EM sims which must also be clean for saving
        if self.simres_EM is not None:  # TODO: to clean_for_asve9)
            self.simres_EM.clean_for_pickle()

        np.savez(prefix, simulation=self)

    def is_EM(self): # override as needed
        """Returns true if the solver is setup for an electromagnetic problem."""
        return False

    def is_AC(self): # override as needed
        """Returns true if the solver is setup for an acoustic problem."""
        return False


class EMSimulation(Simulation):
    def __init__(
        self,
        structure,
        num_modes=20,
        wl_nm=1550,
        n_eff_target=None,
        Stokes=False,
        #calc_EM_mode_energy=False,
        calc_EM_mode_energy=True,
        debug=False,
        **args  # TODO: here for some tests on concurrency, remove later
    ):

        super().__init__(structure, num_modes, debug)

        # independent frequency variable
        self.lambda_m = wl_nm * SI_nm
        self.k_0 = 2 * np.pi / self.lambda_m

        # simulation control parameters
        self.n_eff_target = n_eff_target
        self.Stokes = Stokes

        # additional measures
        self.calc_EM_mode_energy = calc_EM_mode_energy
        self.EM_mode_energy = None
        self.EM_mode_power = None

        self.simres_EM = None  # kludge to simplify save code in Simulation. Fix

        self.fem_mesh = FemMesh()
        self.fem_mesh.build_from_gmsh_mail(self.structure)

        #self.fem_mesh.report_properties(structure)

    def is_EM(self):
        """Returns true if the solver is setup for an electromagnetic problem."""
        return True

    def make_result(self):
        self.sim_result = EMSimResult(self)


    def do_main_eigensolve(self, shortrun):


        tstruc = self.structure
        fm = self.fem_mesh
        opt_props = tstruc.optical_props

        if self.n_modes < 20:
            self.n_modes = 20
            print("Warning: ARPACK needs >= 20 modes so set num_modes=20.")


        E_H_field = 1  # Selected formulation (1=E-Field, 2=H-Field)
        # bnd_cdn_i = 2       # Boundary conditions (0=Dirichlet,1=Neumann,2=domain_x)
        bnd_cdn_i = 0  # Boundary conditions (0=Dirichlet,1=Neumann,2=domain_x)  TODO: this has been set to periodic for some time?!

        itermax = 30  # Maximum number of iterations for convergence

        # Calculate where to center the Eigenmode solver around.
        # (Shift and invert FEM method)
        shift_ksqr = self.n_eff_target**2 * self.k_0**2

        EM_FEM_debug = 0

        print(
            " Boundary conditions:",
            str({0: "Dirichlet", 1: "Neumann", 2: "Periodic"}[bnd_cdn_i]),
        )

        resm = nb_fortran.calc_em_modes(
            self.n_modes,
            self.lambda_m,
            self.d_in_m,
            self.k_perp,
            shift_ksqr,     # passing single complex numbers with intel is fraught
             E_H_field,
            bnd_cdn_i,
            itermax,
            EM_FEM_debug,
            fm.mesh_mail_fname,
            fm.n_msh_pts,
            fm.n_msh_el,
            #opt_props.n_mats_em,   # f2py figures this out from v_refindex, probably because not using f2py(depends)
            opt_props.v_refindexn,
            shortrun
        )

        # self.node_physindex: GMsh physical line or surface number (a small nonneg int). Maps to fortran type_nod
        # self.type_el: material index of each element into list self.v_refindexn (unit-based)

        (
            self.eigs_kz,
            self.fem_evecs,
            self.mode_pol,
            elnd_to_mshpt,
            type_el,
            node_physindex,
            v_nd_xy,
            self.ls_material,
        ) = process_fortran_return(resm, "solving for electromagnetic modes")


        # TODO: ls_material is just refractive index of each element (13 reps for some reason)
        #       clean up and give to FemMesh

        self.fem_mesh.store_fortran_em_mesh_properties(type_el, node_physindex, elnd_to_mshpt, v_nd_xy)

    def calc_field_powers(self):
        tstruc = self.structure
        fm = self.fem_mesh
        #opt_props = tstruc.optical_props

        # Calc unnormalised power in each EM mode Kokou equiv. of Eq. 8.
        print("  Calculating EM mode powers...")

        # These quantities are of order 10^-16 W:
        #   Efield ~ .25/m, P= (\epsilon c /2) n |E|^2 A, with A~ 1 um^2
        #   If we moved to microns as unit of length, with E~0.25 V/um, c=3e14 and A~1 um^2,
        #    we would have powers of order 1 W
        if tstruc.using_linear_elements():
            print('using linear elements')
            # Integration using analytically evaluated basis function integrals. Fast.
            #self.EM_mode_power = nb_fortran.em_mode_energy_int_v2_ez(
            resm = nb_fortran.em_mode_power_sz_analytic(
                self.k_0,
                self.n_modes,
                fm.n_msh_el,
                fm.n_msh_pts,
                fm.elnd_to_mshpt,
                fm.v_nd_xy,
                self.eigs_kz,
                self.fem_evecs,
            )
            (self.EM_mode_power, ) = process_fortran_return(resm, "finding linear element EM mode power")
        else:
            if not tstruc.using_curvilinear_elements():
                print(
                    "Warning: em_mode_energy_int - not sure if mesh contains curvi-linear elements",
                    "\n using slow quadrature integration by default.\n\n",
                )
            # Integration by quadrature. Slowest.
            print('using curvilinear elements')
            resm = nb_fortran.em_mode_power_sz_quadrature(
                self.k_0,
                self.n_modes,
                fm.n_msh_el,
                fm.n_msh_pts,
                fm.elnd_to_mshpt,
                fm.v_nd_xy,
                self.eigs_kz,
                self.fem_evecs,
            )
            (self.EM_mode_power, ) = process_fortran_return(resm, "finding curvilinear element EM mode power")

        # Bring Kokou's def into line with CW formulation.
        self.EM_mode_power = 2.0 * self.EM_mode_power

        #print('EM mode powers', self.EM_mode_power)


    def calc_field_energies(self):
        # Calc linear energy density (not power) in each EM mode - PRA Eq. 6.

        tstruc = self.structure
        fm = self.fem_mesh
        opt_props = tstruc.optical_props


        # These quantities are of order 10^-24 J/m:
        #   Efield ~ .25/m, W= \epsilon n^2 |E|^2 A, with A~ 1e-12 m^2
        #   If we moved to microns as unit of length, with E~0.25 V/um, and A~1 um^2,
        #    we would have densities of order 1e-12 J/m


        if self.calc_EM_mode_energy:
            print("Calculating EM mode energies...")

            if True or tstruc.using_linear_elements():

                # # Integration by quadrature. Slowest.
                resm = nb_fortran.em_mode_act_energy_quadrature(
                    self.n_modes,
                    fm.n_msh_el,
                    fm.n_msh_pts,
                    fm.elnd_to_mshpt,
                    fm.v_nd_xy,
                    opt_props.n_mats_em,
                    fm.v_el_2_mat_idx,
                    opt_props.v_refindexn,
                    self.fem_evecs,
                )
                (self.EM_mode_energy, ) = process_fortran_return(resm,
                     "finding curvilinear element EM mode energy")
            else:
                print(
                    "\n\n FEM routine em_mode_e_energy_int is not implemented for this structure. Can't find group index. \n\n"
                )
                self.EM_mode_energy = np.zeros(self.n_modes, dtype=float)



    def convert_to_Stokes(self):
        # If considering a the backwards propagating Stokes field.
        if self.Stokes:
            self.eigs_kz = -1 * self.eigs_kz
            self.fem_evecs = np.conj(self.fem_evecs)


    def calc_modes(self, shortrun):
        """Run a Fortran FEM calculation to find the optical modes.

        Returns a ``Simulation`` object that has these key values:

        eigs_kz: a 1D array of Eigenvalues (propagation constants) in [1/m]

        fem_evecs: the associated Eigenvectors, ie. the fields, stored as [field comp, node nu on element, Eig value, el nu]

        EM_mode_power: the power in the optical modes. Note this power is negative for modes travelling in the negative
                       z-direction, eg the Stokes wave in backward SBS.
        """

        print("\n\nCalculating EM modes:")



        self.do_main_eigensolve(shortrun)

        if (shortrun): return   #TODO REMOVE ME SHORTRUN

        self.calc_field_powers()

        self.calc_field_energies()

        # This group velocity calc is not accurate in the presence of dispersion!
        # self.group_velocity_EM = self.EM_mode_power/self.EM_mode_power_energy

        if self.Stokes:
            self.convert_to_Stokes()


        self.make_result()





class ACSimulation(Simulation):
    def __init__(
        self,
        structure,
        num_modes=20,
        shift_Hz=None,
        q_AC=0,
        simres_EM=None,
        calc_AC_mode_power=False,
        debug=False,
    ):

        super().__init__(structure, num_modes, debug)

        self.shift_Hz = shift_Hz

        self.q_AC = q_AC
        self.Omega_AC = None

        self.simres_EM = simres_EM

        # Move to result

        self.calc_AC_mode_power = calc_AC_mode_power
        self.AC_mode_energy = None
        self.AC_mode_power = None

        self.ac_alpha_t = None  # temporal acoustic loss [1/s]
        self.ac_linewidth = None  # acoustic linewidth [Hz]
        self.ac_Qmech = None  # acoustic mechanical Q [dimless]
        self.Q_method = QAcMethod.NotSet

        self.ls_material = None

        self.fem_mesh = FemMesh()
        self.fem_mesh.ac_build_from_em(structure, self.simres_EM.fem_mesh)


    def is_AC(self):
        """Returns true if the solver is setup for an acoustic problem."""
        return True



    def do_main_eigensolve(self, bcs):

        if self.n_modes < 20:
            self.n_modes = 20
            print("Warning: ARPACK needs >= 20 modes so set num_modes=20.")

        # Parameters that control how FEM routine runs
        # Boundary conditions (0=Dirichlet,1=Neumann,2=domain_x)
        bnd_cdn_i = 0
        if bcs == "Open":
            print("Attempting open elastic boundary conditions.")
            # icond = 1  # TODO: DO THIS ACTUILLY WORK?

        itermax = 30  # Maximum number of iterations for convergence
        AC_FEM_debug = 0  # Fortran routines will display & save add. info
        ARPACK_tol = 1e-10  # ARPACK accuracy (0.0 for machine precision)

        shift_nu = self.choose_eigensolver_frequency_shift()

        fm = self.fem_mesh
        tstruc = self.structure
        elastic_props = self.structure.elastic_props

        #show_mem_est = False

        # TODO: rmove _AC suffixes from fm.fields_AC
        resm = nb_fortran.calc_ac_modes(
            self.n_modes,
            self.q_AC,
            self.d_in_m,
            shift_nu,  # scalar params
            bnd_cdn_i,
            itermax,
            ARPACK_tol,
            AC_FEM_debug,
            #show_mem_est,
            tstruc.symmetry_flag,
            #elastic_props.n_mats_ac,  # => fort: n_elt_mats, #f2py figures out from arrays
            elastic_props.c_IJ,
            elastic_props.rho,
            fm.ac_mesh_from_em,
            fm.mesh_mail_fname,
            #fm.n_msh_pts,           #f2py figures out from arrays
            #fm.n_msh_el,            #f2py figures out from arrays
            fm.node_physindex,  # => fort: type_nod
            fm.elnd_to_mshpt,       # => fort: elnd_to_mshpt
            fm.v_el_2_mat_idx,  # => fort: type_el
            fm.v_nd_xy,         # => fort: v_nd_xy
        )

        (
            elnd_to_mshpt_out,
            type_el_out,
            v_nd_xy_out,
            self.eigs_nu,
            self.fem_evecs,
            self.mode_pol,
        ) = process_fortran_return(resm, "solving for acoustic modes")

        self.fem_mesh.store_fortran_ac_mesh_properties(type_el_out, elnd_to_mshpt_out, v_nd_xy_out)

                # Retrieve the material properties of each mesh point.
        self.ls_material = nb_fortran.array_material_ac(
            fm.n_msh_el,
            elastic_props.n_mats_ac,
            fm.v_el_2_mat_idx,
            elastic_props.rho,
            elastic_props.c_IJ,
            elastic_props.p_ijkl,
            elastic_props.eta_ijkl,
        )

        # TODO: ls_material is isotropic parts of all elastic properties (repeated 6 times) by mesh elt: [10 x 6 xn_msh_el]
        # rho, c11, c12, c44, p11, p12, p44, eta11, eta12, eta44
        # doesn't seem very useful. May as well turn off.

    def calc_field_powers(self):

        # FEM Eigenvalue is frequency, rather than angular frequency Omega
        Omega_AC = self.eigs_nu * twopi  # DELETE ME

        # Calc unnormalised power in each AC mode - PRA Eq. 18.
        # This quantity is of order 10^12
        # P ~= -2i Omega |c_ijkl| |u|^2/w_x A
        #   with Omega=1e10, |c_ijkl|=1e9 u=1, w_x=1e-6 A=1e-12

        fm = self.fem_mesh
        tstruc = self.structure
        elastic_props = self.structure.elastic_props

        if True or self.calc_AC_mode_power:
            print("doing AC mode power")
            if tstruc.using_linear_elements():
                # Semi-analytic integration following KD 9/9/16 notes. Fastest!
                resm = nb_fortran.ac_mode_power_analytic(
                    self.n_modes,
                    fm.n_msh_el,
                    fm.n_msh_pts,
                    fm.v_nd_xy,
                    fm.elnd_to_mshpt,
                    elastic_props.n_mats_ac,
                    fm.v_el_2_mat_idx,
                    elastic_props.c_IJ,
                    self.q_AC,
                    Omega_AC,
                    self.fem_evecs,
                )
                (self.AC_mode_power,) = process_fortran_return(resm, "finding ac mode power analytic")
                #print("AC mode powers analytic", self.AC_mode_power)
            else:
                if not tstruc.using_curvilinear_elements():
                    print(
                        "Warning: ac_mode_power_int - not sure if mesh contains curvi-linear elements",
                        "\n using slow quadrature integration by default.\n\n",
                    )
                # Integration by quadrature. Slowest.
                resm = nb_fortran.ac_mode_power_quadrature(
                    self.n_modes,
                    fm.n_msh_el,
                    fm.n_msh_pts,
                    fm.v_nd_xy,
                    fm.elnd_to_mshpt,
                    elastic_props.n_mats_ac,
                    fm.v_el_2_mat_idx,
                    #elastic_props.c_IJ,
                    #elastic_props.c_ijkz,
                    elastic_props.c_zjkl,
                    self.q_AC,
                    Omega_AC,
                    self.fem_evecs,
                )

                (self.AC_mode_power,) = process_fortran_return(resm, "finding ac mode power quadrature")
                #print("AC mode powers quadrature", self.AC_mode_power)

    def calc_field_energies(self):
# Calc unnormalised elastic energy in each AC mode - PRA Eq. 16.
        print("Doing AC mode energy")
        fm = self.fem_mesh
        tstruc = self.structure
        elastic_props = self.structure.elastic_props

        Omega_AC = self.eigs_nu * twopi  # DELETE ME

        # This quantity is of order 10^12 since the integration units of (microns)^2 are not accounted for
        if tstruc.using_linear_elements():
            # Semi-analytic integration. Fastest!
            resm = nb_fortran.ac_mode_energy_analytic(
                self.n_modes,
                fm.n_msh_el,
                fm.n_msh_pts,
                #fm.n_nodes,
                fm.elnd_to_mshpt,
                fm.v_nd_xy,
                elastic_props.n_mats_ac,
                fm.v_el_2_mat_idx,
                elastic_props.rho,
                Omega_AC,
                self.fem_evecs,
            )

            (self.AC_mode_energy,) = process_fortran_return(resm, "finding ac mode energy analytic")

        else:
            if not tstruc.using_curvilinear_elements():
                print(
                    "Warning: ac_mode_elastic_energy_int - not sure if mesh contains curvi-linear elements",
                    "\n using slow quadrature integration by default.\n\n",
                )
            # Integration by quadrature. Slowest.
            resm = nb_fortran.ac_mode_energy_quadrature(
                self.n_modes,
                fm.n_msh_el,
                fm.n_msh_pts,
                fm.v_nd_xy,
                #fm.n_nodes,
                fm.elnd_to_mshpt,
                elastic_props.n_mats_ac,
                fm.v_el_2_mat_idx,
                elastic_props.rho,
                Omega_AC,
                self.fem_evecs,
            )
            (self.AC_mode_energy,) = process_fortran_return(resm, "finding ac mode energy quadrature")


    def calc_modes(self, bcs=None):
        """Run a Fortran FEM calculation to find the acoustic modes.

        Returns a ``Simulation`` object that has these key values:

        eigs_nu: a 1D array of Eigenvalues (frequencies) in [1/s]

        fem_evecs: the associated Eigenvectors, ie. the fields, stored as
               [field comp, node num on element, Eig value, el num]

        AC_mode_energy: the elastic power in the acoutic modes.
        """

        print("\n\nCalculating elastic modes")


        self.do_main_eigensolve(bcs)

        self.calc_field_powers()

        self.calc_field_energies()

        self.calc_acoustic_losses()

        self.make_result()

    def choose_eigensolver_frequency_shift(self):
        # Calculate where to center the Eigenmode solver around. (ARPACK Shift and invert FEM method)
        # If no user preference, aim for a frequency a little below the first longitudinal mode

        if self.shift_Hz is None:
            elastic_props = self.structure.elastic_props
            # For AC problem shift is a frequency; [shift] = s^-1.
            v_list = []
            for el in range(elastic_props.n_mats_ac):
                # Using acoustic velocity of longitudinal mode pg 215 Auld vol 1.
                v_list.append(
                    np.sqrt(elastic_props.c_IJ[0, 0][el] / elastic_props.rho[el])
                )
                # # Using acoustic velocity of shear mode pg 215 Auld vol 1.
                # v_list.append(np.sqrt(self.structure.actens_c_IJ[3,3][el]/self.structure.rho[el]))
            AC_velocity = np.real(v_list).min()
            shift = np.real(AC_velocity * self.q_AC / (2.0 * np.pi))

            shift_nu = 0.9 * shift

        else:
            shift_nu = self.shift_Hz
        return shift_nu

    # TODO: make sure this is not done more than once on the same Simulation
    def calc_acoustic_losses(self, fixed_Q=None):
        alpha = None

        tstruc = self.structure
        elastic_props = self.structure.elastic_props

        fm = self.fem_mesh
        Omega_AC = self.eigs_nu * twopi

        if fixed_Q is None:
            self.Q_method = QAcMethod.Intrinsic

            # Calc alpha (loss) Eq. 45
            print("Acoustic loss calculation with intrinsic losses")

            if tstruc.using_linear_elements():
                resm = nb_fortran.ac_alpha_analytic(
                    self.n_modes,
                    fm.n_msh_el,
                    fm.n_msh_pts,
                    fm.elnd_to_mshpt,
                    fm.v_nd_xy,
                    elastic_props.n_mats_ac,
                    fm.v_el_2_mat_idx,
                    elastic_props.eta_ijkl,
                    self.q_AC,
                    Omega_AC,
                    self.fem_evecs,
                    # sim_AC.AC_mode_power) # appropriate for alpha in [1/m]
                    self.AC_mode_energy,
                )  # appropriate for alpha in [1/s]
                (alpha,) = process_fortran_return(resm, "finding ac alpha analytic")
            else:
                if not tstruc.using_curvilinear_elements():
                    print(
                        "Warning: ac_alpha_int - not sure if mesh contains curvi-linear elements",
                        "\n using slow quadrature integration by default.\n\n",
                    )

                # not sure why this is needed by ac_alpha_int
                # overlap = np.zeros(self.n_modes, dtype=complex)
                resm = nb_fortran.ac_alpha_quadrature(
                    self.n_modes,
                    fm.n_msh_el,
                    fm.n_msh_pts,
                    fm.v_nd_xy,
                    fm.elnd_to_mshpt,
                    elastic_props.n_mats_ac,
                    fm.v_el_2_mat_idx,
                    elastic_props.eta_ijkl,
                    self.q_AC,
                    Omega_AC,
                    self.fem_evecs,
                    # sim_AC.AC_mode_power, Fortran_debug) # appropriate for alpha in [1/m]
                    self.AC_mode_energy,
                )  # appropriate for alpha in [1/s]
                (alpha,) = process_fortran_return(resm, "finding ac alpha quadrature")

            self.ac_alpha_t = np.real(alpha)
            # Q_factors = 0.5*(q_AC/alpha)*np.ones(n_modes) # appropriate for alpha in [1/m]
            # appropriate for alpha in [1/s]

            self.ac_Qmech =  0.5 * (np.real(Omega_AC) / self.ac_alpha_t) * np.ones(self.n_modes)

        else:
            print("Acoustic losses calculated from specified fixed Q.")
            self.Q_method = QAcMethod.Fixed
            # factor of a 1/2 because alpha is for power!
            # alpha [1/m] = Omega_AC/(2*vg*fixed_Q) = q_AC/fixed_Q
            # alpha [1/s] = vg * alpha [1/m]
            # alpha [1/s] = Omega_AC/(2*fixed_Q)
            # alpha = 0.5*(q_AC/fixed_Q)*np.ones(n_modes) # appropriate for alpha in [1/m]
            self.ac_Qmech = fixed_Q * np.ones(self.n_modes)
            # appropriate for alpha in [1/s]
            self.ac_alpha_t =  0.5 * (np.real(Omega_AC) / fixed_Q) * np.ones(self.n_modes)

        # SBS linewidth of each resonance in [Hz]   #TODO: not sure about the 1/pi.
        self.ac_linewidth = self.ac_alpha_t / np.pi
        # If linewdith should be amplitude rate in Hz, wouldn't it be
        # alpha/(2 * 2pi)  since alpha is a power decay rate


        #TODO: this Sim/SimResult ownership is a mess.
        # find a cleaner way to deal with the call of this function from gain_calculation
        if self.sim_result is not None:
            self.sim_result.ac_alpha_t = self.ac_alpha_t
            self.sim_result.ac_Qmech = self.ac_Qmech
            self.sim_result.ac_linewidth = self.ac_linewidth

    def make_result(self):
        self.sim_result = ACSimResult(self)






def em_mode_calculation(wg, num_modes, wl_nm, n_eff, Stokes, debug, **args):

    sim = EMSimulation(wg, num_modes=num_modes, wl_nm=wl_nm,
                        n_eff_target=n_eff, Stokes=Stokes, debug=debug, **args)


    if args.get('shortrun', False ):
        print('Got short for EM')
    else:
        args['shortrun'] = False

    sim.calc_modes(args['shortrun'])


    return sim.get_sim_result()

def ac_mode_calculation(wg, num_modes, q_AC, shift_Hz, EM_sim, bcs, debug, **args):

    sim = ACSimulation(wg, num_modes=num_modes, q_AC=q_AC,
                    shift_Hz=shift_Hz, simres_EM=EM_sim, debug=debug, **args)

    sim.calc_modes(bcs)

    return sim.get_sim_result()

# deprecated
def bkwd_Stokes_modes(sim):
    return sim.bkwd_Stokes_modes()
