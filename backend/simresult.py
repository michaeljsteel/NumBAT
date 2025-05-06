# modecalcs.py is a subroutine of NumBAT that contains methods to
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

import numpy as np
#from pathlib import Path

import copy

from modes import ModeAC, ModeEM, ModeInterpolator

from plottools import progressBar
import plotmodes
import reporting

from nbtypes import (
    FieldType,
    FieldTag,
    PointGroup,
    SymRep,
    SI_speed_c,
    twopi,
)

import integration
from fortran import nb_fortran

class SimResult:

    def __init__(self, sim):
        self._structure = sim.structure
        # TODO: limit and ultimately remove access to these
        self._sim = sim
        self._mode_interpolator = None

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
        if self._mode_interpolator is not None:
            self._mode_interpolator.cleanup()

        if self.is_AC():
            # Acoustic sims can contain EM sims which must also be clean for saving
            if self._sim.simres_EM is not None:  # TODO: to clean_for_asve9
                self._sim.simres_EM.clean_for_pickle()

    def get_mode_fields_on_fem_mesh(self, md, field_type):

        fem_evecs = self.fem_evecs_for_ft(field_type)
        #n_msh_el = self.fem_mesh.n_msh_el

        # v_Fx6p = np.zeros(6*n_msh_el, dtype=np.complex128)
        # v_Fy6p = np.zeros(6*n_msh_el, dtype=np.complex128)
        # v_Fz6p = np.zeros(6*n_msh_el, dtype=np.complex128)

        # # i = 0
        # for i_el in range(n_msh_el):
        #     for i_node in range(6):  # TODO: make this one xyz array so we can broadcast
        #         v_Fx6p[i] = fem_evecs[0, i_node, md, i_el]
        #         v_Fy6p[i] = fem_evecs[1, i_node, md, i_el]
        #         v_Fz6p[i] = fem_evecs[2, i_node, md, i_el]
        #         i += 1

        # This is slice equiv of above loops
        # The ordering matches the xy coords defind at XXXXXX

        v_Fx6p = np.array(fem_evecs[0, :6, md, :].T.flat)
        v_Fy6p = np.array(fem_evecs[1, :6, md, :].T.flat)
        v_Fz6p = np.array(fem_evecs[2, :6, md, :].T.flat)

        v_Fa6p = np.sqrt(np.abs(v_Fx6p)**2 +
                        np.abs(v_Fy6p)**2 + np.abs(v_Fz6p)**2)

        return (v_Fx6p, v_Fy6p, v_Fz6p, v_Fa6p)

    def save_simulation(self, prefix): # must be overidden

        self.clean_for_pickle()
        np.savez(prefix, simulation=self)

    # TODO: make this a property?
    def get_mode_interpolator(self):
        if self._mode_interpolator is None:
            self._mode_interpolator = ModeInterpolator(self)
        return self._mode_interpolator

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
        mode_indices=None,
        n_points=501,
        field_type=FieldType.EM_E,
        comps=(),
        prefix="",
        suffix="",
        xlim_min=None,
        xlim_max=None,
        ylim_min=None,
        ylim_max=None,
        aspect = 1.0,
        hide_vector_field=False,
        quiver_points=30,
        ticks=True,
        num_ticks=None,
        title=True,
        frame=True,
        colorbar=True,
        contours=False,
        contour_lst=None,
        suppress_imimre=True,
        logamp=False,
        decorator=None,
        **kwargs
    ):
        """Plot E or H fields of EM mode, or the AC modes displacement fields.

        Args:
            sim_result : A ``Struct`` instance that has had calc_modes calculated

        Keyword Args:
            mode_indices  (list): mode numbers of modes you wish to plot

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

        if 'ivals' in kwargs:
            reporting.deprecated_parameter_exit('ivals', 'mode_indices',
                                                'SimResult.plot_modes')
        d_args = locals()
        del d_args['self']

        ftag = FieldTag.make_from_field(field_type, self.is_AC())
        if ftag.is_EM_H():
            self.make_H_fields()

        domain_s = ftag.domain_type_as_str()
        field_s = ftag.field_type_as_str()


        md_range = mode_indices if mode_indices is not None else range(self.n_modes)

        ntoplot = len(md_range)

        if ntoplot > 1:
            print(f"Plotting {field_s} fields of {ntoplot} {domain_s} modes in range m=[{md_range[0]},{md_range[-1]}]:")
        else:
            print(f"Plotting {field_s} field of {domain_s} mode m={md_range[0]}.")

        md_interp = self.get_mode_interpolator()
        md_interp.define_plot_grid_2D(n_pts=n_points)

        plot_params = plotmodes.PlotParams2D()

        plot_params.update(d_args)
        plot_params['field_type'] = ftag.as_field_type()


        # plot_params.update(
        #     {
        #         'xlim_min': xlim_min,
        #         'xlim_max': xlim_max,
        #         'ylim_min': ylim_min,
        #         'ylim_max': ylim_max,
        #         'aspect': aspect,
        #         #'field_type': field_type,
        #         'field_type': ftag.as_field_type(),
        #         'hide_vector_field': hide_vector_field,
        #         'quiver_points': quiver_points,
        #         'num_ticks': num_ticks,
        #         'colorbar': colorbar,
        #         'contours': contours,
        #         'contour_lst': contour_lst,
        #         'prefix': prefix,
        #         'suffix': suffix,
        #         'ticks': ticks,
        #         'decorator': decorator,
        #         'suppress_imimre': suppress_imimre,
        #     }
        # )

        for m in progressBar(md_range, prefix="  Progress:", length=20):
                self.get_mode(m).plot_mode(comps,
                                           ftag.as_field_type(), plot_params = plot_params)

    def plot_modes_1D(self,
        scut,
        val1,
        val2=None,
        mode_indices=None,
        n_points=501,
        field_type="EM_E",
        num_ticks=None,
        prefix="",
        suffix="",
        ticks=True,
        comps=[],
        decorator=None,
        suppress_imimre=True,
        logx=False,
        logy=False,
    ):

        #field_code = FieldCode(field_type, self.is_AC())

        ftag = FieldTag.make_from_field(field_type, self.is_AC())
        if ftag.is_EM_H():
            self.make_H_fields()

        #if field_code.is_EM_H():
        #   self.make_H_fields()

      #  modetype = field_code.mode_type_as_str() #"acoustic" if field_type == FieldType.AC else "em"

        domain_s = ftag.domain_type_as_str()
        field_s = ftag.field_type_as_str()


        md_range = mode_indices if mode_indices is not None else range(self.n_modes)

        ntoplot = len(md_range)

        if ntoplot > 1:
            print(f"Plotting 1D cuts of {field_s} fields of {ntoplot} {domain_s} modes in range m=[{md_range[0]},{md_range[-1]}]:")
        else:
            print(f"Plotting 1D cut of {field_s} field of {domain_s} mode m={md_range[0]}.")

        md_interp = self.get_mode_interpolator()
        md_interp.define_plot_grid_2D(n_pts=n_points)

        plot_params = plotmodes.PlotParams1D()
        plot_params.update(
            {
                'field_type': ftag.as_field_type(),
                'num_ticks': num_ticks,
                'prefix': prefix,
                'suffix': suffix,
                'ticks': ticks,
                'decorator': decorator,
                'suppress_imimre': suppress_imimre,
                'logx': logx,
                'logy': logy,
            }
        )

        for m in progressBar(md_range, prefix="  Progress:", length=20):
            self.get_mode(m).plot_mode_1D(scut, val1, val2, comps,
                                          ftag.as_field_type(),
                                          plot_params=plot_params)

    def write_modes(self, prefix='', field_type = FieldType.EM_E, mode_indices=None, n_points=501):
        """Writes a set of mode profiles as ascii files on a rectangular grid containing the FEM structure to the current output fields director.

        All components Fx, Fy, Fz are written with separate fields for each component and the real and imaginary parts.
            Args:
                prefix: adjust the output directory
                mode_indices: list of mode indices or None for all
                field_type: for EM fields, write E or H fields
                n_points: dimension of rectangular grid
        """

        ftag = FieldTag.make_from_field(field_type, self.is_AC())

        if ftag.is_EM_H():
            self.make_H_fields()


        md_range = mode_indices if mode_indices is not None else range(self.n_modes)

        for m in md_range:
            self.get_mode(m).write_mode(prefix, n_points, ftag.as_field_type())

    def write_modes_1D(self, s_cut, val1, val2=None, mode_indices=None, prefix='',
                       n_points=501, field_type = FieldType.EM_E):


        #ft = FieldType.AC if self.is_AC() else field_type

        ftag = FieldTag.make_from_field(field_type, self.is_AC())

        if ftag.is_EM_H():
            self.make_H_fields()


        md_range = mode_indices if mode_indices is not None else range(self.n_modes)

        for m in md_range:
            self.get_mode(m).write_mode_1D(s_cut, val1, val2,
                                           n_points, prefix, ftag.as_field_type())


class EMSimResult(SimResult):
    def __init__(self, sim):
        # grant temp access to sim to take what we need
        # ultimately get these things out of sim for good

        super().__init__(sim)

        # Pull EM properties from Simulation

        self.lambda_m = sim.lambda_m

        # Angular freq in units of rad/s
        self.omega_EM = twopi * SI_speed_c / self.lambda_m

        self.k_0 = twopi / self.lambda_m
        self.eigs_kz = sim.eigs_kz
        self.EM_mode_energy = sim.EM_mode_energy
        self.EM_mode_power = sim.EM_mode_power

        self.fem_mesh = sim.fem_mesh
        self.fem_evecs = sim.fem_evecs
        self.fem_evecs_H = None


    def fem_evecs_for_ft(self, ft):
        return self.fem_evecs_H if ft == FieldType.EM_H else self.fem_evecs

    def make_H_fields(self):
        if self.fem_evecs_H is not None:
            return

        n_modes = len(self.eigs_kz)
        fm = self.fem_mesh

        self.fem_evecs_H = nb_fortran.h_mode_field_ez(
            self.k_0,
            n_modes,
            fm.n_msh_el,
            fm.n_msh_pts,
            fm.m_elnd_to_mshpt,
            fm.v_mshpt_xy,
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

    def bkwd_Stokes_modes(self):
        reporting.deprecated_function('numbat.bkwd_backwd_modes()',
                                      'SimResult.clone_as_backward_modes()')
        self.clean_for_pickle()
        return self.clone_as_backward_modes()

    def fwd_Stokes_modes(self):
        reporting.deprecated_function('numbat.fwd_Stokes_modes()',
                                      'SimResult.clone_as_forward_modes()')

        return self.clone_as_forward_modes()

    def clone_as_backward_modes(self):
        """Defines the backward travelling Stokes waves as the conjugate
            of the forward travelling pump waves.

        Returns a ``Simulation`` object that has these key values:

        Eig_values: a 1D array of Eigenvalues (propagation constants) in [1/m]

        fem_evecs: the associated Eigenvectors, ie. the fields, stored as
            [field comp, node nu on element, Eig value, el nu]

        EM_mode_power: the power in the Stokes modes. Note this power is negative because the modes
                    are travelling in the negative z-direction.
        """


        backwd_modes = copy.deepcopy(self)
        backwd_modes.fem_evecs = np.conj(backwd_modes.fem_evecs)
        backwd_modes.eigs_kz = -1.0 * backwd_modes.eigs_kz
        backwd_modes.EM_mode_power = -1.0 * backwd_modes.EM_mode_power

        return backwd_modes


    def clone_as_forward_modes(self):
        """Defines the forward travelling Stokes waves as a copy
            of the forward travelling pump waves.

        Returns a ``Simulation`` object that has these key values:

        """


        self.clean_for_pickle()

        fwd_modes = copy.deepcopy(self)
        return fwd_modes



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
    #     if self._mode_interpolator is not None:
    #         self._mode_interpolator.cleanup()

    # def save_simulation(self, prefix):
    #     self.clean_for_pickle()
    #     np.savez(prefix, simulation=self)


    def is_AC(self):
        return True


    def fem_evecs_for_ft(self, ft):
        return self.fem_evecs

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
