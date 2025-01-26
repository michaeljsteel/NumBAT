# mode_calcs.py is a subroutine of NumBAT that contains methods to
# calculate the EM and Acoustic modes of a structure.

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
from pathlib import Path


import numbat
from modes import ModeAC, ModeEM, ModePlotHelper
from numbattools import process_fortran_return
from plottools import progressBar
import plotmodes

from nbtypes import (
    FieldType,
    FieldCode,
    PointGroup,
    QAcMethod,
    SI_nm,
    SimType,
    SymRep,
    SI_speed_c,
    twopi,
)

from femmesh import FemMesh

import integration
from fortran import nb_fortran


class SimResult:

    def __init__(self, sim):
        self._structure = ( sim.structure )
        # TODO: limit and ultimately remove access to these
        self._sim = sim
        self._mode_plot_helper = None
        self.sim_type = None  # Unknown at this point

        self.n_modes = sim.n_modes
        self.mode_set = []
        self.r0_offset = [0, 0]  # passed to modes when created

        self.sym_reps = None
        self.point_group = PointGroup.Unknown

    def is_EM(self):
        return False

    def is_AC(self):
        return False

    def make_H_fields(self):  # keep linter quiet, should never be called
        raise NotImplementedError

    def clean_for_pickle(self):
        if self._mode_plot_helper is not None:
            self._mode_plot_helper.cleanup()

        if self.is_AC():
            # Acoustic sims can contain EM sims which must also be clean for saving
            if self._sim.simres_EM is not None:  # TODO: to clean_for_asve9
                self._sim.simres_EM.clean_for_pickle()

    def get_modes_on_fem_mesh(self, md, field_type):

        fem_evecs = self.fem_evecs_for_ft(field_type)
        n_msh_el = self.fem_mesh.n_msh_el

        v_Fx6p = np.zeros(6*n_msh_el, dtype=np.complex128)
        v_Fy6p = np.zeros(6*n_msh_el, dtype=np.complex128)
        v_Fz6p = np.zeros(6*n_msh_el, dtype=np.complex128)

        i = 0
        for i_el in range(n_msh_el):
            for i_node in range(6):  # TODO: make this one xyz array so we can broadcast
                v_Fx6p[i] = fem_evecs[0, i_node, md, i_el]
                v_Fy6p[i] = fem_evecs[1, i_node, md, i_el]
                v_Fz6p[i] = fem_evecs[2, i_node, md, i_el]
                i += 1

        v_Fa6p = np.sqrt(np.abs(v_Fx6p)**2 +
                        np.abs(v_Fy6p)**2 + np.abs(v_Fz6p)**2)

        return (v_Fx6p, v_Fy6p, v_Fz6p, v_Fa6p)

    def save_simulation(self, prefix): # must be overidden

        self.clean_for_pickle()
        np.savez(prefix, simulation=self)

    # TODO: make this a property?
    def get_mode_helper(self):
        if self._mode_plot_helper is None:
            self._mode_plot_helper = ModePlotHelper(self)
        return self._mode_plot_helper

    def get_xyshift(self):
        if self.is_EM():
            return self._structure.shift_em_x, self._structure.shift_em_y
        else:
            return self._structure.shift_ac_x, self._structure.shift_ac_y

    # this is clumsy and only works if called after modes have been calced.
    def set_r0_offset(self, rx, ry):
        # print('set_r0_offset is depreacated')
        self.r0_offset = [rx, ry]
        for m in self.get_all_modes():
            m.set_r0_offset(rx, ry)

    def analyse_all_modes(self, n_points=501):
        """Perform modal property analysis on complete set of eigenmodes."""
        modes = self.get_all_modes()
        for m in modes:
            m.analyse_mode(n_points=n_points)

    def _build_modes(self):
        for m in range(self.n_modes):
            if self.is_EM():
                mode = ModeEM(self, m)
            else:
                mode = ModeAC(self, m)
            # awkward and specific to do this here, but might have already been set in the Simulation object befores modes are created
            mode.set_r0_offset(self.r0_offset[0], self.r0_offset[1])
            self.mode_set.append(mode)

    def get_all_modes(self):
        """Returns an array of class `Mode` containing the solved electromagnetic or acoustic modes.

        :rtype: numarray(Mode)
        """
        if not self.mode_set:
            self._build_modes()
        return self.mode_set

    def get_mode(self, m):
        if not self.mode_set:
            self._build_modes()

        return self.mode_set[m]

    def symmetry_classification(self, m):
        """If the point group of the structure has been specified, returns the symmetry class of the given mode.

        :param int m: Index of the mode of interest.
        :rtype: PointGroup
        """
        if self.point_group == PointGroup.Unknown:
            return ""
        return f"{self.point_group.name}:{self.sym_reps[m].name}"

    def analyse_symmetries(self, ptgrp):

        print("Analyse symmetries is out of action")
        return
        self.point_group = ptgrp
        symlist = integration.symmetries(self)
        self.sym_reps = []
        if ptgrp == PointGroup.C2V:
            for m, sym in enumerate(symlist):
                if sym == (1, 1, 1):
                    self.sym_reps.append(SymRep.A)
                elif sym == (-1, 1, -1):
                    self.sym_reps.append(SymRep.B1)
                elif sym == (1, -1, -1):
                    self.sym_reps.append(SymRep.B2)
                elif sym == (-1, -1, 1):
                    self.sym_reps.append(SymRep.B3)
                else:
                    print("Warning: Unknown symmetry pattern", sym)
                    self.sym_reps.append(SymRep.Unknown)

        else:
            print("unknown symmetry properties in mode_calcs")

    def plot_modes(
        self,
        ivals=None,
        n_points=501,
        quiver_points=30,
        xlim_min=0,
        xlim_max=0,
        ylim_min=0,
        ylim_max=0,
        aspect = 1.0,
        field_type="EM_E",
        hide_vector_field=False,
        num_ticks=None,
        colorbar=True,
        contours=False,
        contour_lst=None,
        prefix="",
        suffix="",
        ticks=True,
        comps=(),
        decorator=None,
        suppress_imimre=True,
    ):
        """Plot E or H fields of EM mode, or the AC modes displacement fields.

        Args:
            sim_result : A ``Struct`` instance that has had calc_modes calculated

        Keyword Args:
            ivals  (list): mode numbers of modes you wish to plot

            n_points  (int): The number of points across unitcell to
                interpolate the field onto

            xlim_min  (float): Limit plotted xrange to xlim_min:(1-xlim_max) of unitcell

            xlim_max  (float): Limit plotted xrange to xlim_min:(1-xlim_max) of unitcell

            ylim_min  (float): Limit plotted yrange to ylim_min:(1-ylim_max) of unitcell

            ylim_max  (float): Limit plotted yrange to ylim_min:(1-ylim_max) of unitcell

            field_type  (str): Either 'EM' or 'AC' modes

            num_ticks  (int): Number of tick marks

            contours  (bool): Controls contours being overlaid on fields

            contour_lst  (list): Specify contour values

            stress_fields  (bool): Calculate acoustic stress fields

            pdf_png  (str): File type to save, either 'png' or 'pdf'

            prefix  (str): Add a string to start of file name

            suffix  (str): Add a string to end of file name.

            modal_gains (float array): Pre-calculated gain for each acoustic mode given chosen optical fields.
        """

        field_code = FieldCode(field_type, self.is_AC())

        if field_code.is_EM_H():
            self.make_H_fields()

        modetype = field_code.mode_type_as_str()
        ival_range = ivals if ivals is not None else range(self.n_modes)

        ntoplot = len(ival_range)

        if ntoplot > 1:
            print(f"Plotting {ntoplot} {modetype} modes in range m=[{ival_range[0]},{ival_range[-1]}]:")
        else:
            print(f"Plotting {modetype} mode m={ival_range[0]}.")

        mode_helper = self.get_mode_helper()
        mode_helper.define_plot_grid_2D(n_pts=n_points)


        nbapp = numbat.NumBATApp()
        pf = Path(nbapp.outpath_fields(prefix=prefix))

        if not pf.exists():
            pf.mkdir()

        plotparams = plotmodes.PlotParams2D()

        #mode_helper.update_plot_params(
        plotparams.update(
            {
                'xlim_min': xlim_min,
                'xlim_max': xlim_max,
                'ylim_min': ylim_min,
                'ylim_max': ylim_max,
                'aspect': aspect,
                #'field_type': field_type,
                'field_type': field_code.as_field_type(),
                'hide_vector_field': hide_vector_field,
                'quiver_points': quiver_points,
                'num_ticks': num_ticks,
                'colorbar': colorbar,
                'contours': contours,
                'contour_lst': contour_lst,
                'prefix': prefix,
                'suffix': suffix,
                'ticks': ticks,
                'decorator': decorator,
                'suppress_imimre': suppress_imimre,
            }
        )

        for m in progressBar(ival_range, prefix="  Progress:", length=20):
                self.get_mode(m).plot_mode(comps, field_code.as_field_type(), plot_params = plotparams)

    def plot_modes_1D(self,
        scut,
        val1,
        val2=None,
        ivals=None,
        n_points=501,
        field_type="EM_E",
        num_ticks=None,
        prefix="",
        suffix="",
        ticks=True,
        comps=[],
        decorator=None,
        suppress_imimre=True,
    ):



        field_code = FieldCode(field_type, self.is_AC())

        if field_code.is_EM_H():
            self.make_H_fields()

        modetype = field_code.mode_type_as_str() #"acoustic" if field_type == FieldType.AC else "em"

        ival_range = ivals if ivals is not None else range(self.n_modes)

        ntoplot = len(ival_range)

        if ntoplot > 1:
            print(f"Plotting 1D cuts of {ntoplot} {modetype} modes in range m=[{ival_range[0]},{ival_range[-1]}]:")
        else:
            print(f"Plotting 1D cut of {modetype} mode m={ival_range[0]}.")

        mode_helper = self.get_mode_helper()
        mode_helper.define_plot_grid_2D(n_pts=n_points)

        nbapp = numbat.NumBATApp()
        pf = Path(nbapp.outpath_fields(prefix=prefix))

        if not pf.exists():
            pf.mkdir()


        plotparams = plotmodes.PlotParams1D()
        plotparams.update(
            {
                'field_type': field_code.as_field_type(),
                'num_ticks': num_ticks,
                'prefix': prefix,
                'suffix': suffix,
                'ticks': ticks,
                'decorator': decorator,
                'suppress_imimre': suppress_imimre,
            }
        )

        for m in progressBar(ival_range, prefix="  Progress:", length=20):
            self.get_mode(m).plot_mode_1D_cut(scut, val1, val2,
                                                comps, field_code.as_field_type(), plot_params=plotparams)


class EMSimResult(SimResult):
    def __init__(self, sim):
        # grant temp access to sim to take what we need
        # ultimately get these things out of sim for good

        super().__init__(sim)

        self.lambda_m = sim.lambda_m
        self.omega_EM = (
            twopi * SI_speed_c / self.lambda_m
        )  # Angular freq in units of rad/s
        self.k_0 = twopi / self.lambda_m

        self.fem_mesh = sim.fem_mesh

        # Eigen solutions
        self.fem_evecs = sim.fem_evecs
        self.fem_evecs_H = None

        self.eigs_kz = sim.eigs_kz

        # Additional quantities
        self.EM_mode_energy = sim.EM_mode_energy
        self.EM_mode_power = sim.EM_mode_power


    # def clean_for_pickle(self):
    #     if self._mode_plot_helper is not None:
    #         self._mode_plot_helper.cleanup()

    # def save_simulation(self, prefix):
    #     self.clean_for_pickle()
    #     np.savez(prefix, simulation=self)



    def fem_evecs_for_ft(self, ft):
        return self.fem_evecs_H if ft == FieldType.EM_H else self.fem_evecs

    def make_H_fields(self):
        n_modes = len(self.eigs_kz)
        fm = self.fem_mesh
        self.fem_evecs_H = nb_fortran.h_mode_field_ez(
            self.k_0,
            n_modes,
            fm.n_msh_el,
            fm.n_msh_pts,
            fm.n_nodes,
            fm.elnd_to_mshpt,
            fm.v_nd_xy,
            self.eigs_kz,
            self.fem_evecs,
        )

    def is_EM(self):
        return True

    def neff(self, m):
        """Returns the effective index of EM mode `m`.

        :param int m: Index of the mode of interest.
        :rtype: float
        """

        return np.real(self.eigs_kz[m] * self.lambda_m / (twopi))

    def neff_all(self):
        """Return an array of the effective index of all electromagnetic modes.

        :return: numpy array of effective indices
        :rtype: array(float)
        """

        return np.real(self.eigs_kz * self.lambda_m / (twopi))

    def ngroup_EM_available(self):
        """Returns true if a measure of the electromagnetic group index is available."""
        return not (self.EM_mode_energy is None or self.EM_mode_power is None)

    def ngroup_EM(self, m):
        """Returns the group index of electromagnetic mode `m`, if available, otherwise returns zero with a warning message.

        :param int m: Index of the mode of interest.
        :return: Group index of the mode.
        :rtype: float
        """

        if not self.ngroup_EM_available():
            print(
                """EM group index requires calculation of mode energy and mode power when calculating EM modes.
               Set calc_EM_mode_energy=True and calc_AC_mode_power=True in call to Simulation"""
            )
            ng = 0

        if abs(self.EM_mode_energy[m]) > 0.0:
            vg = np.real(self.EM_mode_power[m] / self.EM_mode_energy[m])
            ng = SI_speed_c / vg
        else:
            ng = 0

        return ng

    def ngroup_EM_all(self):
        """Returns a numpy array of the group index of all electromagnetic modes, if available,
        otherwise returns a zero numarray with a warning message.

        :return: numpy array of  index of the mode.
        :rtype: array(float)
        """
        if not self.ngroup_EM_available():
            print(
                """EM group index requires calculation of mode energy and mode power when calculating EM modes.
               Set calc_EM_mode_energy=True in call to calc_EM_modes"""
            )
            return np.zeros(len(self.eigs_kz), dtype=float)
        vg = np.real(self.EM_mode_power / self.EM_mode_energy)
        ng = SI_speed_c / vg
        return ng

    def kz_EM(self, m):
        """Returns the wavevector in 1/m of electromagnetic mode `m`.

        :param int m: Index of the mode of interest.
        :return: Wavevector k in 1/m.
        :rtype: float
        """

        return self.eigs_kz[m]

    def kz_EM_all(self):
        """Return an array of the wavevector in 1/m of all electromagnetic modes.

        :return: numpy array of wavevectors in 1/m
        :rtype: array(float)
        """
        return self.eigs_kz


class ACSimResult(SimResult):
    def __init__(self, sim):

        super().__init__(sim)

        self.q_AC = sim.q_AC  # input wavevector

        self.fem_mesh = sim.fem_mesh

        # eigen solution
        self.eigs_nu = sim.eigs_nu  # eigenvalues are frequency nu in Hz
        self.fem_evecs = sim.fem_evecs

        # additional properties
        self.Omega_AC = self.eigs_nu * twopi  # DELETE ME

        self.AC_mode_energy = sim.AC_mode_energy
        self.AC_mode_power = sim.AC_mode_power

        self.Q_method = sim.Q_method
        self.ac_alpha_t = sim.ac_alpha_t
        self.ac_Qmech = sim.ac_Qmech
        self.ac_linewidth = sim.ac_linewidth

    # def clean_for_pickle(self):
    #     if self._mode_plot_helper is not None:
    #         self._mode_plot_helper.cleanup()

    # def save_simulation(self, prefix):
    #     self.clean_for_pickle()
    #     np.savez(prefix, simulation=self)



    def fem_evecs_for_ft(self, ft):
        return self.fem_evecs

    def is_AC(self):
        return True

    def nu_AC(self, m):
        """Returns the frequency in Hz of acoustic mode `m`.

        :param int m: Index of the mode of interest.
        :return: Frequency :math:`\\nu` in Hz
        :rtype: float
        """

        return self.eigs_nu[m]

    def nu_AC_all(self):
        """Return an array of the frequency in Hz of all acoustic modes.

        :return: numpy array of frequencies in Hz
        :rtype: array(float)
        """
        return self.eigs_nu

    def Omega_AC_all(self):
        """Return an array of the angular frequency in 1/s of all acoustic modes.

        :return: numpy array of angular frequencies in 1/s
        :rtype: array(float)
        """
        return self.eigs_nu * twopi

    def vp_AC(self, m):
        """Return the phase velocity in m/s of acoustic mode `m`.

        :return: Phase velocity of acoustic mode `m` in m/s
        :rtype: float
        """
        return np.pi * 2 * np.real(self.eigs_nu[m]) / self.q_AC

    def vp_AC_all(self):
        """Return an array of the phase velocity in m/s of all acoustic modes.

        :return: numpy array of elastic phase velocities in m/s
        :rtype: array(float)
        """
        return np.pi * 2 * np.real(self.eigs_nu) / self.q_AC

    def vgroup_AC_available(self):
        """Returns true if a measure of the acoustic group velocity is available."""
        return not (self.AC_mode_energy is None or self.AC_mode_power is None)

    def vg_AC(self, m):
        """Return group velocity of AC mode m in m/s"""
        if self.AC_mode_energy is None or self.AC_mode_power is None:
            print(
                """AC group velocity requires calculation of mode energy and mode power when calculating AC modes.
               Set calc_AC_mode_power=True in call to calc_AC_modes"""
            )
            return 0
        vg = np.real(self.AC_mode_power[m] / self.AC_mode_energy[m])
        return vg

    def vg_AC_all(self):
        """Return group velocity of all AC modes in m/s"""

        if self.AC_mode_energy is None or self.AC_mode_power is None:
            print(
                """AC group velocity requires calculation of mode energy and mode power when calculating AC modes.
               Set calc_AC_mode_power=True in call to calc_AC_modes"""
            )
            return np.zeros(len(self.eigs_nu), dtype=float)
        vg = np.real(self.AC_mode_power / self.AC_mode_energy)
        return vg

    def alpha_t_AC(self, m):
        return self.ac_alpha_t[m]

    def alpha_t_AC_all(self):
        return self.ac_alpha_t

    # spatial loss [1/m]  #TODO: surely this should be the group velocity for scaling between spatial and temporal decays #Which is a problem because it requires knowing vg
    def alpha_s_AC(self, m):
        return self.alpha_t_AC(m) / self.vg_AC(m)

    def alpha_s_AC_all(self):  # spatial loss [1/m]
        return self.alpha_t_AC_all / self.vg_AC_all

    def Qmech_AC(self, m):
        return self.ac_Qmech[m]

    def Qmech_AC_all(self):
        return self.ac_Qmech

    def linewidth_AC(self, m):
        return self.ac_linewidth[m]

    def linewidth_AC_all(self):
        return self.ac_linewidth


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

        self.sim_type = None

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

    def is_EM(self):
        """Returns true if the solver is setup for an electromagnetic problem."""
        return self.sim_type == SimType.EM

    def is_AC(self):
        """Returns true if the solver is setup for an acoustic problem."""
        return self.sim_type == SimType.AC


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
    ):

        super().__init__(structure, num_modes, debug)

        self.sim_type = SimType.EM

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

    def make_result(self):
        self.sim_result = EMSimResult(self)

    def calc_modes(self):
        """Run a Fortran FEM calculation to find the optical modes.

        Returns a ``Simulation`` object that has these key values:

        eigs_kz: a 1D array of Eigenvalues (propagation constants) in [1/m]

        fem_evecs: the associated Eigenvectors, ie. the fields, stored as [field comp, node nu on element, Eig value, el nu]

        EM_mode_power: the power in the optical modes. Note this power is negative for modes travelling in the negative
                       z-direction, eg the Stokes wave in backward SBS.
        """

        print("\n\nCalculating EM modes:")

        tstruc = self.structure

        if self.n_modes < 20:
            self.n_modes = 20
            print("Warning: ARPACK needs >= 20 modes so set num_modes=20.")

        self.E_H_field = 1  # Selected formulation (1=E-Field, 2=H-Field)
        # bnd_cdn_i = 2       # Boundary conditions (0=Dirichlet,1=Neumann,2=domain_x)
        bnd_cdn_i = 0  # Boundary conditions (0=Dirichlet,1=Neumann,2=domain_x)  TODO: this has been set to periodic for some time?!

        itermax = 30  # Maximum number of iterations for convergence
        EM_FEM_debug = self.debug  # Fortran routines will display & save add. info

        print(
            " Boundary conditions:",
            str({0: "Dirichlet", 1: "Neumann", 2: "Periodic"}[bnd_cdn_i]),
        )

        # Calculate where to center the Eigenmode solver around.
        # (Shift and invert FEM method)
        shift_ksqr = self.n_eff_target**2 * self.k_0**2

        EM_FEM_debug = 0

        fm = self.fem_mesh
        opt_props = tstruc.optical_props

        #print('optprops', fm.n_msh_pts, fm.n_msh_el, opt_props.n_mats_em, opt_props.v_refindexn)
        #print(f'bloch:', self.k_perp)
        #print(f'ksqr: {np.real(shift_ksqr):.4e} {np.imag(shift_ksqr):.4e}')


        resm = nb_fortran.calc_em_modes(
            self.n_modes,
            self.lambda_m,
            self.d_in_m,
            self.k_perp,
            shift_ksqr,     # passing single complex numbers with intel is fraught
            self.E_H_field,
            bnd_cdn_i,
            itermax,
            EM_FEM_debug,
            fm.mesh_mail_fname,
            fm.n_msh_pts,
            fm.n_msh_el,
            #opt_props.n_mats_em,   # f2py figures this out from v_refindex, probably because not using f2py(depends)
            opt_props.v_refindexn,
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

        #print("modepol", self.mode_pol)
        #print("ls material: n", self.ls_material, self.ls_material.shape)

        self.fem_mesh.store_fortran_em_mesh_properties(type_el, node_physindex, elnd_to_mshpt, v_nd_xy)

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



        # Calc linear energy density (not power) in each EM mode - PRA Eq. 6.

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

        # This group velocity calc is not accurate in the presence of dispersion!
        # self.group_velocity_EM = self.EM_mode_power/self.EM_mode_power_energy

        # If considering a the backwards propagating Stokes field.
        if self.Stokes:
            self.eigs_kz = -1 * self.eigs_kz
            self.fem_evecs = np.conj(self.fem_evecs)

        self.make_result()

        # Not necessary because EM FEM mesh always normalised in area to unity.
        # print area
        # x_tmp = []
        # y_tmp = []
        # for i in np.arange(self.n_msh_pts):
        #     x_tmp.append(self.v_nd_xy[0,i])
        #     y_tmp.append(self.v_nd_xy[1,i])
        # x_min = np.min(x_tmp); x_max=np.max(x_tmp)
        # y_min = np.min(y_tmp); y_max=np.max(y_tmp)
        # area = abs((x_max-x_min)*(y_max-y_min))
        # print area
        # self.EM_mode_power = self.EM_mode_power*area


# def calc_EM_mode_energy(self):  # these require extraction of numerical props from Fortran. Clean that up first.
#      assert(self.is_EM())
#      if not self.EM_mode_energy is None: return  # nothing to do

#    def calc_EM_mode_power(self):
#      assert(self.is_EM())
#      if not self.EM_mode_power is None: return  # nothing to do

#    def calc_EM_mode_power_and_energy(self):
#      assert(self.is_EM())
#      self.calc_EM_mode_power()
#      self.calc_EM_mode_energy()


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

        self.sim_type = SimType.AC
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

    def calc_modes(self, bcs=None):
        """Run a Fortran FEM calculation to find the acoustic modes.

        Returns a ``Simulation`` object that has these key values:

        eigs_nu: a 1D array of Eigenvalues (frequencies) in [1/s]

        fem_evecs: the associated Eigenvectors, ie. the fields, stored as
               [field comp, node num on element, Eig value, el num]

        AC_mode_energy: the elastic power in the acoutic modes.
        """

        print("\n\nCalculating elastic modes")

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

        show_mem_est = False

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
            show_mem_est,
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

        # FEM Eigenvalue is frequency, rather than angular frequency Omega
        Omega_AC = self.eigs_nu * twopi  # DELETE ME

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

        # Calc unnormalised power in each AC mode - PRA Eq. 18.
        # This quantity is of order 10^12
        # P ~= -2i Omega |c_ijkl| |u|^2/w_x A
        #   with Omega=1e10, |c_ijkl|=1e9 u=1, w_x=1e-6 A=1e-12

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

        # Calc unnormalised elastic energy in each AC mode - PRA Eq. 16.
        print("Doing AC mode energy")

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
            print("Acoustic loss calc")

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


def bkwd_Stokes_modes(EM_sim):  # TODO: make a member function
    """Defines the backward travelling Stokes waves as the conjugate
        of the forward travelling pump waves.

    Returns a ``Simulation`` object that has these key values:

    Eig_values: a 1D array of Eigenvalues (propagation constants) in [1/m]

    fem_evecs: the associated Eigenvectors, ie. the fields, stored as
           [field comp, node nu on element, Eig value, el nu]

    EM_mode_power: the power in the Stokes modes. Note this power is negative because the modes
                   are travelling in the negative z-direction.
    """

    # EM_sim.clean_for_pickle()

    Stokes_modes = copy.deepcopy(EM_sim)
    Stokes_modes.fem_evecs = np.conj(Stokes_modes.fem_evecs)
    Stokes_modes.eigs_kz = -1.0 * Stokes_modes.eigs_kz
    Stokes_modes.EM_mode_power = -1.0 * Stokes_modes.EM_mode_power

    return Stokes_modes


def fwd_Stokes_modes(EM_sim):  # TODO: make a member function
    """Defines the forward travelling Stokes waves as a copy
        of the forward travelling pump waves.

    Returns a ``Simulation`` object that has these key values:

    """

    EM_sim.clean_for_pickle()

    Stokes_modes = copy.deepcopy(EM_sim)
    return Stokes_modes



def em_mode_calculation(wg, num_modes, wl_nm, n_eff, Stokes, debug, **args):

    sim = EMSimulation(wg, num_modes=num_modes, wl_nm=wl_nm,
                        n_eff_target=n_eff, Stokes=Stokes, debug=debug, **args)

    sim.calc_modes()

    return sim.get_sim_result()

def ac_mode_calculation(wg, num_modes, q_AC, shift_Hz, EM_sim, bcs, debug, **args):

    sim = ACSimulation(wg, num_modes=num_modes, q_AC=q_AC,
                    shift_Hz=shift_Hz, simres_EM=EM_sim, debug=debug, **args)

    sim.calc_modes(bcs)

    return sim.get_sim_result()
