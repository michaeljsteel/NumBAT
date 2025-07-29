"""
Integration and gain calculation routines for NumBAT.

This module provides classes and functions for calculating SBS gain, Q-factors, and related quantities
from electromagnetic and acoustic simulation results. It includes the GainProps class for storing and
accessing gain and loss data, as well as functions for performing the necessary integrals and interpolations.
"""


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



import csv
import numpy as np
from scipy import interpolate
import matplotlib

from nbtypes import SI_permittivity_eps0
from numbattools import np_min_max, process_fortran_return
import reporting
import plotgain
from fortran import nb_fortran


class GainProps(object):
    """
    Stores and manages SBS gain, Q-factor, and related properties for a set of EM and acoustic modes.

    Attributes:
        _max_pumps_m (int): Maximum number of pump EM modes.
        _max_Stokes_m (int): Maximum number of Stokes EM modes.
        _max_ac_m (int): Maximum number of acoustic modes.
        _allowed_pumps_m (list): Indices of allowed pump EM modes.
        _allowed_Stokes_m (list): Indices of allowed Stokes EM modes.
        _allowed_ac_m (list): Indices of allowed acoustic modes.
        _gain_tot (np.ndarray): Total SBS gain array.
        _gain_PE (np.ndarray): Photoelastic SBS gain array.
        _gain_MB (np.ndarray): Moving boundary SBS gain array.
        def_m_pump (int): Default pump mode index.
        def_m_Stokes (int): Default Stokes mode index.
        linewidth_Hz (float or np.ndarray): Acoustic linewidth(s) in Hz.
        alpha (float or np.ndarray): Acoustic loss(es).
        Q_factor (float or np.ndarray): Acoustic Q-factor(s).
        sim_AC: Acoustic simulation object.
    """
    @staticmethod
    def _set_allowed_ms(mlist, m_allow, maxm):
        """
        Set allowed mode indices in mlist based on m_allow and maxm.

        Args:
            mlist (list): List to update with allowed indices.
            m_allow (int, list, or 'All'): Allowed mode(s).
            maxm (int): Maximum number of modes.
        """
        if m_allow == "All":
            mlist[:] = list(range(maxm))
        elif isinstance(m_allow, int):
            mlist[:] = [m_allow]
        else:
            mlist[:] = m_allow
        if max(mlist) >= maxm:
            reporting.report_and_exit(
                "Requested mode range too large in GainProps object: " + str(m_allow)
            )

    def __init__(self):
        """
        Initialize GainProps with default values and all modes allowed.
        """
        self._max_pumps_m = 1
        self._max_Stokes_m = 1
        self._max_ac_m = 1
        self._allowed_pumps_m = []  # Pumps for which gains have been calced
        self._allowed_Stokes_m = []  # Stokes for which gains have been calced
        self._allowed_ac_m = []  # Acoustics for which gains have been calced
        self._gain_tot = None
        self._gain_PE = None
        self._gain_MB = None
        self.def_m_pump = 0  # must be one of allowed_pumps_m TODO: needs to be checked
        self.def_m_Stokes = 0  # must be one of allowed_Stokes_m TODO: needs to be checked
        self.linewidth_Hz = None
        self.alpha = None
        self.Q_factor = None
        self.sim_AC = None
        # these choices guaranteed to work
        self.set_allowed_EM_pumps('All')
        self.set_allowed_EM_Stokes('All')
        self.set_allowed_AC('All')
        self.set_EM_modes(0,0)

    def _set_sim_AC(self, sac):
        """Set the acoustic simulation object."""
        self.sim_AC = sac

    def set_allowed_EM_pumps(self, m_allow):
        """Set allowed pump EM modes."""
        self._set_allowed_ms(self._allowed_pumps_m, m_allow, self._max_pumps_m)

    def set_allowed_EM_Stokes(self, m_allow):
        """Set allowed Stokes EM modes."""
        self._set_allowed_ms(self._allowed_Stokes_m, m_allow, self._max_Stokes_m)

    def set_allowed_AC(self, m_allow):
        """Set allowed acoustic modes."""
        self._set_allowed_ms(self._allowed_ac_m, m_allow, self._max_ac_m)

    def set_EM_modes(self, mP, mS):
        """Set the default pump and Stokes EM mode indices."""
        self.def_m_pump = mP
        self.def_m_Stokes = mS

    def gain_total(self, m_AC):
        """Return the total SBS gain for the current pump, Stokes, and given acoustic mode index."""
        return self._gain_tot[self.def_m_pump, self.def_m_Stokes, m_AC]

    def gain_PE(self, m_AC):
        """Return the photoelastic SBS gain for the current pump, Stokes, and given acoustic mode index."""
        return self._gain_PE[self.def_m_pump, self.def_m_Stokes, m_AC]

    def gain_MB(self, m_AC):
        """Return the moving boundary SBS gain for the current pump, Stokes, and given acoustic mode index."""
        return self._gain_MB[self.def_m_pump, self.def_m_Stokes, m_AC]

    def gain_total_all(self):
        """Return the total SBS gain for the current pump and Stokes over all acoustic modes."""
        return self._gain_tot[self.def_m_pump, self.def_m_Stokes, :]

    def gain_PE_all(self):
        """Return the photoelastic SBS gain for the current pump and Stokes over all acoustic modes."""
        return self._gain_PE[self.def_m_pump, self.def_m_Stokes, :]

    def gain_MB_all(self):
        """Return the moving boundary SBS gain for the current pump and Stokes over all acoustic modes."""
        return self._gain_MB[self.def_m_pump, self.def_m_Stokes, :]

    def gain_total_all_by_em_modes(self, m_pump, m_Stokes):
        """Return the total SBS gain for given pump and Stokes indices over all acoustic modes."""
        return self._gain_tot[m_pump, m_Stokes, :]

    def gain_PE_all_by_em_modes(self, m_pump, m_Stokes):
        """Return the photoelastic SBS gain for given pump and Stokes indices over all acoustic modes."""
        return self._gain_PE[m_pump, m_Stokes, :]

    def gain_MB_all_by_em_modes(self, m_pump, m_Stokes):
        """Return the moving boundary SBS gain for given pump and Stokes indices over all acoustic modes."""
        return self._gain_MB[m_pump, m_Stokes, :]

    def alpha_all(self):
        """Return the acoustic loss (alpha) for all modes."""
        return self.alpha

    def Q_factor_all(self):
        """Return the Q-factor for all modes."""
        return self.Q_factor

    def linewidth_Hz_all(self):
        """Return the linewidth (Hz) for all modes."""
        return self.linewidth_Hz

    def gain_total_raw(self):
        """Return the raw total SBS gain array."""
        return self._gain_tot

    def gain_PE_raw(self):
        """Return the raw photoelastic SBS gain array."""
        return self._gain_PE

    def gain_MB_raw(self):
        """Return the raw moving boundary SBS gain array."""
        return self._gain_MB

    def _set_gain_tot(self, g):
        """Set the total SBS gain array and update mode counts and allowed modes."""
        self._gain_tot = g
        (self._max_pumps_m, self._max_Stokes_m, self._max_ac_m) = self._gain_tot.shape
        # some reasonable default choices
        self.set_allowed_EM_pumps(0)
        self.set_allowed_EM_Stokes(0)
        self.set_allowed_AC(range(self._max_ac_m))

    def _set_gain_PE(self, g):
        """Set the photoelastic SBS gain array."""
        self._gain_PE = g

    def _set_gain_MB(self, g):
        """Set the moving boundary SBS gain array."""
        self._gain_MB = g

    def _set_alpha(self, a):
        """Set the acoustic loss (alpha)."""
        self.alpha = a

    def _set_linewidth_Hz(self, lwhz):
        """Set the linewidth (Hz)."""
        self.linewidth_Hz = lwhz

    def _set_Q_factor(self, qf):
        """Set the Q-factor."""
        self.Q_factor = qf

    def check_acoustic_expansion_size(self):
        """
        Check if the maximum gain occurs near the upper end of the acoustic mode range.
        Warn if so, as this may indicate the need for more acoustic modes.
        """
        for mP in self._allowed_pumps_m:
            for mS in self._allowed_Stokes_m:
                t_gains = self._gain_tot[mP, mS, :]
                num_AC = len(t_gains)
                imaxg = np.argmax(np.abs(t_gains))
                if imaxg > num_AC * 0.75:
                    maxg = np.abs(t_gains[imaxg])
                    reporting.register_warning(
                        f"""For pump and Stokes indices {mP} and {mS}, the maximum total SBS gain of {maxg:.3e} was found for acoustic mode {imaxg} which is in the upper """
                        + r"25% of the number of acoustic modes in the calculation."
                        + "\nYou should probably check the consistency of the calculation with a larger number of acoustic modes."
                    )

    def plot_spectra(
        self,
        freq_min=0.0,
        freq_max=50e9,
        num_interp_pts=0,
        dB=False,
        dB_peak_amp=10,
        mode_comps=False,
        logy=False,
        save_txt=False,
        prefix="",
        suffix="",
        decorator=None,
        show_gains="All",
        mark_modes_thresh=0.02,
    ):
        """
        Plot SBS gain spectra using the stored gain arrays and simulation data.

        Args:
            freq_min (float): Minimum frequency for plot (Hz).
            freq_max (float): Maximum frequency for plot (Hz).
            num_interp_pts (int): Number of interpolation points.
            dB (bool): Plot in dB scale if True.
            dB_peak_amp (float): Peak amplitude for dB plot.
            mode_comps (bool): Plot mode components if True.
            logy (bool): Logarithmic y-axis if True.
            save_txt (bool): Save data to text file if True.
            prefix (str): Filename prefix for output.
            suffix (str): Filename suffix for output.
            decorator (callable): Function to decorate the plot.
            show_gains (str): Which gains to show ('All' or specific).
            mark_modes_thresh (float): Threshold for marking modes.
        Returns:
            Output of plotgain.plot_gain_spectra().
        """
        # TODO: this needs work
        em_pump_m = self._allowed_pumps_m[0]
        em_Stokes_m = self._allowed_Stokes_m[0]
        return plotgain.plot_gain_spectra(
            self.sim_AC,
            self._gain_tot,
            self._gain_PE,
            self._gain_MB,
            self.linewidth_Hz,
            em_pump_m,
            em_Stokes_m,
            self._allowed_ac_m,
            freq_min,
            freq_max,
            num_interp_pts,
            dB,
            dB_peak_amp,
            mode_comps,
            logy,
            save_txt,
            prefix,
            suffix,
            decorator,
            show_gains,
            mark_modes_thresh
        )


def get_gains_and_qs(
    sim_EM_pump,
    sim_EM_Stokes,
    sim_AC,
    q_AC,
    EM_mode_index_pump=0,
    EM_mode_index_Stokes=0,
    AC_mode_index=0,
    fixed_Q=None,
    typ_select_out=None,
    **kwargs
):

    if 'EM_ival_pump' in kwargs:

        reporting.report_and_exit('The parameter EM_ival_pump is now called EM_mode_index_pump')

    if 'EM_ival_Stokes' in kwargs:
        reporting.report_and_exit('The parameter EM_ival_Stokes is now called EM_mode_index_Stokes')

    if 'AC_ival' in kwargs:
        reporting.report_and_exit('The parameter AC_ival is now called AC_mode_index')

    # TODO: get rid of this old backend
    SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = gain_and_qs(
        sim_EM_pump,
        sim_EM_Stokes,
        sim_AC,
        q_AC,
        EM_mode_index_pump,
        EM_mode_index_Stokes,
        AC_mode_index,
        fixed_Q,
        typ_select_out,
        new_call_format=True,

    )

    gain = GainProps()
    gain._set_sim_AC(sim_AC)
    gain._set_gain_tot(SBS_gain)
    gain._set_gain_PE(SBS_gain_PE)
    gain._set_gain_MB(SBS_gain_MB)
    gain._set_alpha(alpha)
    gain._set_linewidth_Hz(linewidth_Hz)
    gain._set_Q_factor(Q_factors)

    gain.set_allowed_EM_pumps(EM_mode_index_pump)
    gain.set_allowed_EM_Stokes(EM_mode_index_Stokes)
    gain.set_allowed_AC(AC_mode_index)
    gain.set_EM_modes(EM_mode_index_pump, EM_mode_index_Stokes)

    gain.check_acoustic_expansion_size()

    return gain


def gain_and_qs(
    simres_EM_pump,
    simres_EM_Stokes,
    simres_AC,
    q_AC,
    EM_mode_index_pump=0,
    EM_mode_index_Stokes=0,
    AC_mode_index=0,
    fixed_Q=None,
    typ_select_out=None,
    new_call_format=False,
):
    r""" Calculate interaction integrals and SBS gain.

        Implements Eqs. 33, 41, 45, 91 of
        Wolff et al. PRA 92, 013836 (2015) doi/10.1103/PhysRevA.92.013836
        These are for Q_photoelastic, Q_moving_boundary, the Acoustic loss "alpha",
        and the SBS gain respectively.

        Note there is a sign error in published Eq. 41. Also, in implementing Eq. 45 we use integration by parts, with a
        boundary integral term set to zero on physical grounds, and filled in some missing subscripts. We prefer to express
        Eq. 91 with the Lorentzian explicitly visible, which makes it clear how to transform to frequency space.
        The final integrals are

        .. math::
            Q^{\rm PE} = -\varepsilon_0 \int_A {\rm d}^2r \sum_{ijkl} \varepsilon^2_r e^{(s)\star}_i e^{(p)}_j p_{ijkl} \partial_k u_l^{*},\\

            Q^{\rm MB} =  \int_C {\rm d \mathbf{r} (\mathbf{u}^{*} \cdot \hat n}) \big[ (\varepsilon_a - \varepsilon_b)
            \varepsilon_0 ({\rm \hat n \times \mathbf{e}}) \cdot ({\rm \hat n \times \mathbf{e}}) -
            (\varepsilon_a^{-1} - \varepsilon_b^{-1})  \varepsilon_0^{-1} ({\rm \hat n \cdot \mathbf{d}})
            \cdot ({\rm \hat n \cdot \mathbf{d}}) \big],\\

            \alpha = \frac{\Omega^2}{\mathcal{E}_{ac}} \int {\rm d}^2r \sum_{ijkl} \partial_i u_j^{*} \eta_{ijkl} \partial_k u_l,\\

            \Gamma =  \frac{2 \omega \Omega {\rm Re} (Q_1 Q_1^*)}{P_p P_s \mathcal{E}_{ac}} \frac{1}{\alpha} \frac{\alpha^2}{\alpha^2 + \kappa^2}.


        Args:
            sim_EM_pump  (``Simulation`` object): Contains all info on pump EM modes

            sim_EM_Stokes  (``Simulation`` object): Contains all info on Stokes EM modes

            sim_AC  (``Simulation`` object): Contains all info on AC modes

            q_AC  (float): Propagation constant of acoustic modes.

        Keyword Args:
            EM_mode_index_pump  (int/string): Specify mode number of EM mode 1 (pump mode)
                to calculate interactions for.
                Numbering is python index so runs from 0 to num_EM_modes-1,
                with 0 being fundamental mode (largest prop constant).
                Can also set to 'All' to include all modes.

            EM_mode_index_Stokes  (int/string): Specify mode number of EM mode 2 (stokes mode)
                to calculate interactions for.
                Numbering is python index so runs from 0 to num_EM_modes-1,
                with 0 being fundamental mode (largest prop constant).
                Can also set to 'All' to include all modes.

            AC_mode_index  (int/string): Specify mode number of AC mode
                to calculate interactions for.
                Numbering is python index so runs from 0 to num_AC_modes-1,
                with 0 being fundamental mode (largest prop constant).
                Can also set to 'All' to include all modes.

            fixed_Q  (int): Specify a fixed Q-factor for the AC modes, rather than
                calculating the acoustic loss (alpha).

        Returns:
            SBS_gain  : The SBS gain including both photoelastic and moving boundary contributions.
                        Note this will be negative for backwards SBS because gain is expressed as
                        gain in power as move along z-axis in positive direction, but the Stokes
                        waves experience gain as they propagate in the negative z-direction.
                        Dimensions = [n_modes_EM_Stokes,n_modes_EM_pump,n_modes_AC].

            SBS_gain_PE  : The SBS gain for only the photoelastic effect.
                           The comment about negative gain (see SBS_gain above) holds here also.
                           Dimensions = [n_modes_EM_Stokes,n_modes_EM_pump,n_modes_AC].

            SBS_gain_MB  : The SBS gain for only the moving boundary effect.
                           The comment about negative gain (see SBS_gain above) holds here also.
                           Dimensions = [n_modes_EM_Stokes,n_modes_EM_pump,n_modes_AC].

            alpha  : The acoustic power loss for each mode in [1/s]. Dimensions = [n_modes_AC].
    """

    # Notes about internals of fortran integration
    # Calc overlap of basis functions (and PE tensor and epsilon)
    # Then use this multiple times for calc of each mode field values

    # phi is values of Lagrange polynomials (1-6) at that node.
    # grad is value of gradient of Lagrange polynomials (1-6) at that node.
    # i variables refer to E field
    # j variables refer to H field
    # ww weight function
    # coeff numerical integration

    if not new_call_format:
        print(
            "Note:, integration.gain_and_qs() is deprecated.  You are encouraraged to switch to the Gain() interface provided by integration.get_gains_and_qs()"
        )

    if EM_mode_index_pump == "All":
        EM_mode_index_pump_fortran = -1
    else:
        EM_mode_index_pump_fortran = EM_mode_index_pump + 1  # convert back to Fortran indexing

    if EM_mode_index_Stokes == "All":
        EM_mode_index_Stokes_fortran = -1
    else:
        EM_mode_index_Stokes_fortran = EM_mode_index_Stokes + 1  # convert back to Fortran indexing

    if AC_mode_index == "All":
        AC_mode_index_fortran = -1
    else:
        AC_mode_index_fortran = AC_mode_index + 1  # convert back to Fortran indexing

    # TODO : bad !
    sim_EM_pump = simres_EM_pump._sim
    sim_EM_Stokes = simres_EM_Stokes._sim
    sim_AC = simres_AC._sim

    sim_AC.calc_acoustic_losses(fixed_Q)

    ncomps = 3

    n_modes_EM_pump = sim_EM_pump.n_modes
    n_modes_EM_Stokes = sim_EM_Stokes.n_modes
    n_modes_AC = sim_AC.n_modes

    fem_ac = sim_AC.fem_mesh
    struc = sim_EM_pump.structure
    #opt_props = struc.optical_props
    el_props = struc.elastic_props

    nnodes = fem_ac.n_nodes

    trimmed_EM_pump_field = np.zeros(
        (ncomps, nnodes, n_modes_EM_pump, fem_ac.n_msh_el), dtype=complex
    )
    trimmed_EM_Stokes_field = np.zeros(
        (ncomps, nnodes, n_modes_EM_Stokes, fem_ac.n_msh_el), dtype=complex
    )

    # We want EM fields on the AC grid and only the first six FEM nodes per elt
    # of the full 13 FEM
    for el in range(fem_ac.n_msh_el):
        new_el = fem_ac.el_convert_tbl[el]
        trimmed_EM_pump_field[:, :, :, el] = sim_EM_pump.fem_evecs[:, :6, :, new_el]
        trimmed_EM_Stokes_field[:, :, :, el] = sim_EM_Stokes.fem_evecs[:, :6, :, new_el]

    acoustic_eps_effs = el_props.v_acoustic_eps_eff

    alpha = simres_AC.alpha_t_AC_all()
    elastic_props = sim_AC.structure.elastic_props
    # Calc Q_photoelastic Eq. 33

    is_curvi = False
    if struc.using_linear_elements():
        print("\n Photoelastic calc: linear elements")
    else:
        print("\n Photoelastic calc: curvilinear elements")
        is_curvi = True

    resm = nb_fortran.photoelastic_int_common(
        is_curvi,
        sim_EM_pump.n_modes,
        sim_EM_Stokes.n_modes,
        sim_AC.n_modes,
        EM_mode_index_pump_fortran,
        EM_mode_index_Stokes_fortran,
        AC_mode_index_fortran,
        fem_ac.n_msh_el,
        fem_ac.n_msh_pts,
        fem_ac.m_elnd_to_mshpt,
        fem_ac.v_mshpt_xy,
        elastic_props.n_mats_ac,
        fem_ac.v_el_2_mat_idx,
        elastic_props.p_ijkl,
        q_AC,
        trimmed_EM_pump_field,
        trimmed_EM_Stokes_field,
        sim_AC.fem_evecs,
        acoustic_eps_effs,
    )

    (Q_PE,) = process_fortran_return(
        resm, "finding linear element photoelastic couplings"
    )

 #    for ii in range(3):
 #        for jj in range(3):
 #            for kk in range(3):
 #                for ll in range(3):
 #                    pijkl_comp = copy.deepcopy(elastic_props.p_ijkl)
 #                    pijkl_comp[ii,jj,kk,ll]=0.0
 #                    pijkl_switch = elastic_props.p_ijkl - pijkl_comp
 #
 #                    resm = nb_fortran.photoelastic_int_common(
 #                        True or is_curvi,
 #                        sim_EM_pump.n_modes,
 #                        sim_EM_Stokes.n_modes,
 #                        sim_AC.n_modes,
 #                        EM_mode_index_pump_fortran,
 #                        EM_mode_index_Stokes_fortran,
 #                        AC_mode_index_fortran,
 #                        fem_ac.n_msh_el,
 #                        fem_ac.n_msh_pts,
 #                        fem_ac.m_elnd_to_mshpt,
 #                        fem_ac.v_mshpt_xy,
 #                        elastic_props.n_mats_ac,
 #                        fem_ac.v_el_2_mat_idx,
 #                        pijkl_switch,
 #                        q_AC,
 #                        trimmed_EM_pump_field,
 #                        trimmed_EM_Stokes_field,
 #                        sim_AC.fem_evecs,
 #                        acoustic_eps_effs,)
 #
 #                    (t_Qpe,) = process_fortran_return( resm, "finding linear element photoelastic couplings")
 #                    II = voigt.to_Voigt[ii,jj]
 #                    JJ = voigt.to_Voigt[kk,ll]
 #                    print('got Qpe', II, JJ, t_Qpe[0,0,1])
 #



    # Calc Q_moving_boundary Eq. 41
    #TODO: Needs major fixing to get multiple boundaries right
    typ_select_in = 1  # first element in acoustic_eps_effs list, in fortan indexing
    if len(acoustic_eps_effs) == 2 and acoustic_eps_effs[0]!= acoustic_eps_effs[1]: # This check needed in case two regions actually have the same index and are not the interesting boundary
        typ_select_out = 2
    elif typ_select_out is None:
        typ_select_out = -1
    print("\n Moving boundary calc")


    resm = nb_fortran.moving_boundary(
        sim_EM_pump.n_modes,
        sim_EM_Stokes.n_modes,
        sim_AC.n_modes,
        EM_mode_index_pump_fortran,
        EM_mode_index_Stokes_fortran,
        AC_mode_index_fortran,
        fem_ac.n_msh_el,
        fem_ac.n_msh_pts,
        fem_ac.m_elnd_to_mshpt,
        fem_ac.v_mshpt_xy,
        elastic_props.n_mats_ac,
        fem_ac.v_el_2_mat_idx,
        typ_select_in,
        typ_select_out,
        trimmed_EM_pump_field,
        trimmed_EM_Stokes_field,
        sim_AC.fem_evecs,
        acoustic_eps_effs,
    )

    (Q_MB,) = process_fortran_return(resm, "finding moving boundary coupling")

    Q = Q_PE + Q_MB  # TODO: the Q couplings come out as non trivially complex. Why?

    omEM = simres_EM_pump.omega_EM
    OmAC = simres_AC.Omega_AC

    # first find the numerators of gain. Sturmberg Eq (12)
    gain = 2 * omEM * OmAC * np.real(Q * np.conj(Q))
    gain_PE = 2 * omEM * OmAC * np.real(Q_PE * np.conj(Q_PE))
    gain_MB = 2 * omEM * OmAC * np.real(Q_MB * np.conj(Q_MB))

    normal_fact = np.zeros(
        (n_modes_EM_Stokes, n_modes_EM_pump, n_modes_AC), dtype=complex
    )
    for i in range(
        n_modes_EM_Stokes
    ):  # TODO: express this as some one line outer product?
        P1 = sim_EM_Stokes.EM_mode_power[i]
        for j in range(n_modes_EM_pump):
            P2 = sim_EM_pump.EM_mode_power[j]
            for k in range(n_modes_AC):
                # P3 = sim_AC.AC_mode_power[k]
                P3 = simres_AC.AC_mode_energy[k]
                normal_fact[i, j, k] = P1 * P2 * P3 * alpha[k]

    # now normalise to find the final gain. Sturmberg Eq (12)
    SBS_gain = np.real(gain / normal_fact)
    SBS_gain_PE = np.real(gain_PE / normal_fact)
    SBS_gain_MB = np.real(gain_MB / normal_fact)


    return (
        SBS_gain,
        SBS_gain_PE,
        SBS_gain_MB,
        simres_AC.linewidth_AC_all(),
        simres_AC.Qmech_AC_all(),
        simres_AC.alpha_t_AC_all(),
    )


#### Categorise modes by their symmetries #############################################
def symmetries(simres, n_points=10, negligible_threshold=1e-5):
    """ Plot EM mode fields.

        Args:
            simres : A ``Struct`` instance that has had calc_modes calculated

        Keyword Args:
            n_points  (int): The number of points across unitcell to \
                interpolate the field onto.
    """

    mode_fields = simres.fem_evecs

    # field mapping
    x_tmp = []
    y_tmp = []
    for i in np.arange(simres.n_msh_pts):
        x_tmp.append(simres.v_mshpt_xy[0, i])
        y_tmp.append(simres.v_mshpt_xy[1, i])
    x_min, x_max = np_min_max(x_tmp)
    y_min, y_max = np_min_max(y_tmp)

    area = abs((x_max - x_min) * (y_max - y_min))
    n_pts_x = int(n_points * abs(x_max - x_min) / np.sqrt(area))
    n_pts_y = int(n_points * abs(y_max - y_min) / np.sqrt(area))
    v_x = np.zeros(n_pts_x * n_pts_y)
    v_y = np.zeros(n_pts_x * n_pts_y)
    i = 0
    for x in np.linspace(x_min, x_max, n_pts_x):
        for y in np.linspace(y_min, y_max, n_pts_y):
            v_x[i] = x
            v_y[i] = y
            i += 1
    v_x = np.array(v_x)
    v_y = np.array(v_y)

    # unrolling data for the interpolators
    m_elnd_to_mshpt = simres.m_elnd_to_mshpt.T
    v_mshpt_xy = simres.v_mshpt_xy.T

    sym_list = []

    for mode_index in range(len(simres.eigs_kz)):
        # dense triangulation with multiple points
        v_x6p = np.zeros(6 * simres.n_msh_el)
        v_y6p = np.zeros(6 * simres.n_msh_el)
        v_Ex6p = np.zeros(6 * simres.n_msh_el, dtype=np.complex128)
        v_Ey6p = np.zeros(6 * simres.n_msh_el, dtype=np.complex128)
        v_triang6p = []

        i = 0
        for i_el in np.arange(simres.n_msh_el):
            # triangles
            idx = np.arange(6 * i_el, 6 * (i_el + 1))
            triangles = [
                [idx[0], idx[3], idx[5]],
                [idx[1], idx[4], idx[3]],
                [idx[2], idx[5], idx[4]],
                [idx[3], idx[4], idx[5]],
            ]
            v_triang6p.extend(triangles)

            for i_node in np.arange(6):
                # index for the coordinates
                i_ex = m_elnd_to_mshpt[i_el, i_node] - 1
                # values
                v_x6p[i] = v_mshpt_xy[i_ex, 0]
                v_y6p[i] = v_mshpt_xy[i_ex, 1]
                v_Ex6p[i] = mode_fields[0, i_node, mode_index, i_el]
                v_Ey6p[i] = mode_fields[1, i_node, mode_index, i_el]
                i += 1

        # dense triangulation with unique points
        v_triang1p = []
        for i_el in np.arange(simres.n_msh_el):
            # triangles
            triangles = [
                [
                    m_elnd_to_mshpt[i_el, 0] - 1,
                    m_elnd_to_mshpt[i_el, 3] - 1,
                    m_elnd_to_mshpt[i_el, 5] - 1,
                ],
                [
                    m_elnd_to_mshpt[i_el, 1] - 1,
                    m_elnd_to_mshpt[i_el, 4] - 1,
                    m_elnd_to_mshpt[i_el, 3] - 1,
                ],
                [
                    m_elnd_to_mshpt[i_el, 2] - 1,
                    m_elnd_to_mshpt[i_el, 5] - 1,
                    m_elnd_to_mshpt[i_el, 4] - 1,
                ],
                [
                    m_elnd_to_mshpt[i_el, 3] - 1,
                    m_elnd_to_mshpt[i_el, 4] - 1,
                    m_elnd_to_mshpt[i_el, 5] - 1,
                ],
            ]
            v_triang1p.extend(triangles)

        # triangulations
        triang6p = matplotlib.tri.Triangulation(v_x6p, v_y6p, v_triang6p)
        triang1p = matplotlib.tri.Triangulation(
            v_mshpt_xy[:, 0], v_mshpt_xy[:, 1], v_triang1p
        )

        # building interpolators: triang1p for the finder, triang6p for the values
        finder = matplotlib.tri.TrapezoidMapTriFinder(triang1p)
        ReEx = matplotlib.tri.LinearTriInterpolator(
            triang6p, v_Ex6p.real, trifinder=finder
        )
        ImEx = matplotlib.tri.LinearTriInterpolator(
            triang6p, v_Ex6p.imag, trifinder=finder
        )
        ReEy = matplotlib.tri.LinearTriInterpolator(
            triang6p, v_Ey6p.real, trifinder=finder
        )
        ImEy = matplotlib.tri.LinearTriInterpolator(
            triang6p, v_Ey6p.imag, trifinder=finder
        )

        ### plotting
        # interpolated fields
        m_ReEx = ReEx(v_x, v_y).reshape(n_pts_x, n_pts_y)
        m_ReEy = ReEy(v_x, v_y).reshape(n_pts_x, n_pts_y)
        m_ImEx = ImEx(v_x, v_y).reshape(n_pts_x, n_pts_y)
        m_ImEy = ImEy(v_x, v_y).reshape(n_pts_x, n_pts_y)
        m_Ex = m_ReEx + 1j * m_ImEx
        m_Ey = m_ReEy + 1j * m_ImEy

        if np.max(np.abs(m_Ex[~np.isnan(m_Ex)])) < negligible_threshold:
            m_Ex = np.zeros(np.shape(m_Ex))
        if np.max(np.abs(m_Ey[~np.isnan(m_Ey)])) < negligible_threshold:
            m_Ey = np.zeros(np.shape(m_Ey))

        m_Ex_ymirror = np.zeros((n_pts_x, n_pts_y), dtype=np.complex128)
        m_Ex_xmirror = np.zeros((n_pts_x, n_pts_y), dtype=np.complex128)
        m_Ex_rotated = np.zeros((n_pts_x, n_pts_y), dtype=np.complex128)
        m_Ey_ymirror = np.zeros((n_pts_x, n_pts_y), dtype=np.complex128)
        m_Ey_xmirror = np.zeros((n_pts_x, n_pts_y), dtype=np.complex128)
        m_Ey_rotated = np.zeros((n_pts_x, n_pts_y), dtype=np.complex128)
        Ex_sigma_y = 0
        Ey_sigma_y = 0
        Ex_sigma_x = 0
        Ey_sigma_x = 0
        Ex_C_2 = 0
        Ey_C_2 = 0

        for ix in range(n_pts_x):
            for iy in range(n_pts_y):
                m_Ex_ymirror[ix, iy] = m_Ex[ix, n_pts_y - iy - 1]
                m_Ey_ymirror[ix, iy] = -1 * (m_Ey[ix, n_pts_y - iy - 1])
                m_Ex_xmirror[ix, iy] = -1 * (m_Ex[n_pts_x - ix - 1, iy])
                m_Ey_xmirror[ix, iy] = m_Ey[n_pts_x - ix - 1, iy]
                m_Ex_rotated[ix, iy] = -1 * (m_Ex[n_pts_x - ix - 1, n_pts_y - iy - 1])
                m_Ey_rotated[ix, iy] = -1 * (m_Ey[n_pts_x - ix - 1, n_pts_y - iy - 1])

        Ex_sigma_y = np.sum(np.abs(m_Ex - m_Ex_ymirror))
        Ey_sigma_y = np.sum(np.abs(m_Ey - m_Ey_ymirror))
        Ex_sigma_x = np.sum(np.abs(m_Ex - m_Ex_xmirror))
        Ey_sigma_x = np.sum(np.abs(m_Ey - m_Ey_xmirror))
        Ex_C_2 = np.sum(np.abs(m_Ex - m_Ex_rotated))
        Ey_C_2 = np.sum(np.abs(m_Ey - m_Ey_rotated))
        sigma_y = (Ex_sigma_y + Ey_sigma_y) / (n_pts_x * n_pts_y)
        sigma_x = (Ex_sigma_x + Ey_sigma_x) / (n_pts_x * n_pts_y)
        C_2 = (Ex_C_2 + Ey_C_2) / (n_pts_x * n_pts_y)

        if abs(C_2) > 0.2:
            C_2_print = -1
        else:
            C_2_print = 1
        if abs(sigma_y) > 0.1:
            sigma_y_print = -1
        else:
            sigma_y_print = 1
        if abs(sigma_x) > 0.1:
            sigma_x_print = -1
        else:
            sigma_x_print = 1
        sym_list.append((C_2_print, sigma_y_print, sigma_x_print))

        # v_plots = [np.real(m_Ex_ymirror),np.real(m_Ey_ymirror),np.real(m_Ez_ymirror),
        #     np.imag(m_Ex_ymirror),np.imag(m_Ey_ymirror),np.imag(m_Ez_ymirror)]
        # # field plots
        # plt.clf()
        # plt.figure(figsize=(13,13))
        # for i_p,plot in enumerate(v_plots):
        #     ax = plt.subplot(3,3,i_p+1)
        #     im = plt.imshow(plot.T,cmap='inferno');
        #     # colorbar
        #     divider = make_axes_locatable(ax)
        #     cax = divider.append_axes("right", size="5%", pad=0.1)
        #     cbar = plt.colorbar(im, cax=cax)
        # plt.savefig('fields/field_%(i)i-ymirror.pdf' %
        #         {'i' : mode_index}, bbox_inches='tight')
        # plt.close()
        # v_plots = [np.real(m_Ex_xmirror),np.real(m_Ey_xmirror),np.real(m_Ez_xmirror),
        #     np.imag(m_Ex_xmirror),np.imag(m_Ey_xmirror),np.imag(m_Ez_xmirror)]
        # # field plots
        # plt.clf()
        # plt.figure(figsize=(13,13))
        # for i_p,plot in enumerate(v_plots):
        #     ax = plt.subplot(3,3,i_p+1)
        #     im = plt.imshow(plot.T,cmap='inferno');
        #     # colorbar
        #     divider = make_axes_locatable(ax)
        #     cax = divider.append_axes("right", size="5%", pad=0.1)
        #     cbar = plt.colorbar(im, cax=cax)
        # plt.savefig('fields/field_%(i)i-xmirror.pdf' %
        #         {'i' : mode_index}, bbox_inches='tight')
        # plt.close()
        # v_plots = [np.real(m_Ex_rotated),np.real(m_Ey_rotated),np.real(m_Ez_rotated),
        #     np.imag(m_Ex_rotated),np.imag(m_Ey_rotated),np.imag(m_Ez_rotated)]
        # # field plots
        # plt.clf()
        # plt.figure(figsize=(13,13))
        # for i_p,plot in enumerate(v_plots):
        #     ax = plt.subplot(3,3,i_p+1)
        #     im = plt.imshow(plot.T,cmap='inferno');
        #     # colorbar
        #     divider = make_axes_locatable(ax)
        #     cax = divider.append_axes("right", size="5%", pad=0.1)
        #     cbar = plt.colorbar(im, cax=cax)
        # plt.savefig('fields/field_%(i)i-rotated.pdf' %
        #         {'i' : mode_index}, bbox_inches='tight')
        # plt.close()
        # v_plots = [m_ReEx,m_ReEy,m_ReEz,
        #     m_ImEx,m_ImEy,m_ImEz]
        # # field plots
        # plt.clf()
        # plt.figure(figsize=(13,13))
        # for i_p,plot in enumerate(v_plots):
        #     ax = plt.subplot(3,3,i_p+1)
        #     im = plt.imshow(plot.T,cmap='inferno');
        #     # colorbar
        #     divider = make_axes_locatable(ax)
        #     cax = divider.append_axes("right", size="5%", pad=0.1)
        #     cbar = plt.colorbar(im, cax=cax)
        # plt.savefig('fields/field_%(i)i.pdf' %
        #         {'i' : mode_index}, bbox_inches='tight')
        # plt.close()

    return sym_list


