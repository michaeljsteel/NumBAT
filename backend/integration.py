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
    @staticmethod
    def _set_allowed_ms(mlist, m_allow, maxm):
        if m_allow == "All":
            mlist = range(maxm)
        elif isinstance(m_allow, int):
            mlist[:] = [m_allow]
        else:
            mlist[:] = m_allow
        if max(mlist) >= maxm:
            reporting.report_and_exit(
                "Requested mode range too large in GainProps object: " + str(m_allow)
            )

    def __init__(self):
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

        # must be one of allowed_Stokes_m TODO: needs to be checked
        self.def_m_Stokes = 0

        self.linewidth_Hz = None
        self.alpha = None
        self.Q_factor = None
        self.sim_AC = None

        # these choices guaranteed to work
        self.set_allowed_EM_pumps(0)
        self.set_allowed_EM_Stokes(0)
        self.set_allowed_AC(0)

    def _set_sim_AC(self, sac):
        self.sim_AC = sac

    def set_allowed_EM_pumps(self, m_allow):
        self._set_allowed_ms(self._allowed_pumps_m, m_allow, self._max_pumps_m)

    def set_allowed_EM_Stokes(self, m_allow):
        self._set_allowed_ms(self._allowed_Stokes_m, m_allow, self._max_Stokes_m)

    def set_allowed_AC(self, m_allow):
        self._set_allowed_ms(self._allowed_ac_m, m_allow, self._max_ac_m)

    def set_EM_modes(self, mP, mS):
        self.def_m_pump = mP
        self.def_m_Stokes = mS

    def gain_total(self, m_AC):
        return self._gain_tot[self.def_m_pump, self.def_m_Stokes, m_AC]

    def gain_PE(self, m_AC):
        return self._gain_PE[self.def_m_pump, self.def_m_Stokes, m_AC]

    def gain_MB(self, m_AC):
        return self._gain_MB[self.def_m_pump, self.def_m_Stokes, m_AC]

    def gain_total_all(self):
        return self._gain_tot[self.def_m_pump, self.def_m_Stokes, :]

    def gain_PE_all(self):
        return self._gain_PE[self.def_m_pump, self.def_m_Stokes, :]

    def gain_MB_all(self):
        return self._gain_MB[self.def_m_pump, self.def_m_Stokes, :]

    def gain_total_all_by_em_modes(self, m_pump, m_Stokes):
        return self._gain_tot[m_pump, m_Stokes, :]

    def gain_PE_all_by_em_modes(self, m_pump, m_Stokes):
        return self._gain_PE[m_pump, m_Stokes, :]

    def gain_MB_all_by_em_modes(self, m_pump, m_Stokes):
        return self._gain_MB[m_pump, m_Stokes, :]

    def alpha_all(self):
        return self.alpha

    def Q_factor_all(self):
        return self.Q_factor

    def linewidth_Hz_all(self):
        return self.linewidth_Hz

    def gain_total_raw(self):
        return self._gain_tot

    def gain_PE_raw(self):
        return self._gain_PE

    def gain_MB_raw(self):
        return self._gain_MB

    # def alpha_raw(self):        return self.alpha
    # def Q_factor_raw(self):     return self.Q_factor
    # def linewidth_Hz_raw(self): return self.linewidth_Hz

    def _set_gain_tot(self, g):
        self._gain_tot = g
        (self._max_pumps_m, self._max_Stokes_m, self._max_ac_m) = self._gain_tot.shape

        # some reasonable default choices
        self.set_allowed_EM_pumps(0)
        self.set_allowed_EM_Stokes(0)
        self.set_allowed_AC(range(self._max_ac_m))

    def _set_gain_PE(self, g):
        self._gain_PE = g

    def _set_gain_MB(self, g):
        self._gain_MB = g

    def _set_alpha(self, a):
        self.alpha = a

    def _set_linewidth_Hz(self, lwhz):
        self.linewidth_Hz = lwhz

    def _set_Q_factor(self, qf):
        self.Q_factor = qf

    def check_acoustic_expansion_size(self):
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
        num_interp_pts=3000,
        dB=False,
        dB_peak_amp=10,
        mode_comps=False,
        logy=False,
        pdf_png="png",
        save_txt=False,
        prefix="",
        suffix="",
        decorator=None,
        show_gains="All",
        mark_modes_thresh=0.02,
    ):
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
            pdf_png,
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
    EM_ival_pump=0,
    EM_ival_Stokes=0,
    AC_ival=0,
    fixed_Q=None,
    typ_select_out=None,
):

    # TODO: get rid of this old backend
    SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = gain_and_qs(
        sim_EM_pump,
        sim_EM_Stokes,
        sim_AC,
        q_AC,
        EM_ival_pump,
        EM_ival_Stokes,
        AC_ival,
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

    gain.set_allowed_EM_pumps(EM_ival_pump)
    gain.set_allowed_EM_Stokes(EM_ival_Stokes)
    gain.set_allowed_AC(AC_ival)
    gain.set_EM_modes(EM_ival_pump, EM_ival_Stokes)

    gain.check_acoustic_expansion_size()

    return gain


def gain_and_qs(
    simres_EM_pump,
    simres_EM_Stokes,
    simres_AC,
    q_AC,
    EM_ival_pump=0,
    EM_ival_Stokes=0,
    AC_ival=0,
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
            EM_ival_pump  (int/string): Specify mode number of EM mode 1 (pump mode)
                to calculate interactions for.
                Numbering is python index so runs from 0 to num_EM_modes-1,
                with 0 being fundamental mode (largest prop constant).
                Can also set to 'All' to include all modes.

            EM_ival_Stokes  (int/string): Specify mode number of EM mode 2 (stokes mode)
                to calculate interactions for.
                Numbering is python index so runs from 0 to num_EM_modes-1,
                with 0 being fundamental mode (largest prop constant).
                Can also set to 'All' to include all modes.

            AC_ival  (int/string): Specify mode number of AC mode
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

    if EM_ival_pump == "All":
        EM_ival_pump_fortran = -1
    else:
        EM_ival_pump_fortran = EM_ival_pump + 1  # convert back to Fortran indexing
    if EM_ival_Stokes == "All":
        EM_ival_Stokes_fortran = -1
    else:
        EM_ival_Stokes_fortran = EM_ival_Stokes + 1  # convert back to Fortran indexing
    if AC_ival == "All":
        AC_ival_fortran = -1
    else:
        AC_ival_fortran = AC_ival + 1  # convert back to Fortran indexing

    # TODO : bad !
    sim_EM_pump = simres_EM_pump._sim
    sim_EM_Stokes = simres_EM_Stokes._sim
    sim_AC = simres_AC._sim

    ncomps = 3

    n_modes_EM_pump = sim_EM_pump.n_modes
    n_modes_EM_Stokes = sim_EM_Stokes.n_modes
    n_modes_AC = sim_AC.n_modes

    fem_ac = sim_AC.fem_mesh
    struc = sim_EM_pump.structure
    opt_props = struc.optical_props
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

        # for n in range(nnodes):
        #     for x in range(ncomps):
        #         for ival in range(n_modes_EM_pump):
        #             trimmed_EM_pump_field[x, n, ival, el] = sim_EM_pump.fem_evecs[
        #                 x, n, ival, new_el
        #             ]

        #         for ival in range(n_modes_EM_Stokes):
        #             trimmed_EM_Stokes_field[x, n, ival, el] = sim_EM_Stokes.fem_evecs[
        #                 x, n, ival, new_el
        #             ]

    relevant_eps_effs = []  # TODO: move this into ElProps
    for el_typ in range(opt_props.n_mats_em):
        if el_typ + 1 in el_props.typ_el_AC:
            relevant_eps_effs.append(opt_props.v_refindexn[el_typ] ** 2)

    sim_AC.calc_acoustic_losses(fixed_Q)

    # print("\n-----------------------------------------------")
    #    if fixed_Q is None:
    #        # Calc alpha (loss) Eq. 45
    #        print("Acoustic loss calc")
    #        start = time.time()
    #        try:
    #            if sim_EM_pump.structure.inc_shape in sim_EM_pump.structure.linear_element_shapes:
    #                alpha = nb_fortran.ac_alpha_int_v2(sim_AC.n_modes,
    #                    sim_AC.n_msh_el, sim_AC.n_msh_pts, nnodes,
    #                    sim_AC.elnd_to_mshpt, sim_AC.v_el_2_mat_idx, sim_AC.v_nd_xy,
    #                    sim_AC.structure.n_mats_ac, sim_AC.structure.elastic_props.eta_ijkl,
    #                    q_AC, sim_AC.Omega_AC, sim_AC.fem_evecs,
    #                    # sim_AC.AC_mode_power) # appropriate for alpha in [1/m]
    #                    sim_AC.AC_mode_energy) # appropriate for alpha in [1/s]
    #            else:
    #                if sim_EM_pump.structure.inc_shape not in sim_EM_pump.structure.curvilinear_element_shapes:
    #                    print("Warning: ac_alpha_int - not sure if mesh contains curvi-linear elements",
    #                        "\n using slow quadrature integration by default.\n\n")
    #                alpha = nb_fortran.ac_alpha_int(sim_AC.n_modes,
    #                    sim_AC.n_msh_el, sim_AC.n_msh_pts, nnodes,
    #                    sim_AC.elnd_to_mshpt, sim_AC.v_el_2_mat_idx, sim_AC.v_nd_xy,
    #                    sim_AC.structure.n_mats_ac, sim_AC.structure.elastic_props.eta_ijkl,
    #                    q_AC, sim_AC.Omega_AC, sim_AC.fem_evecs,
    #                    # sim_AC.AC_mode_power, Fortran_debug) # appropriate for alpha in [1/m]
    #                    sim_AC.AC_mode_energy, Fortran_debug) # appropriate for alpha in [1/s]
    #        except KeyboardInterrupt:
    #            print("\n\n Routine ac_alpha_int interrupted by keyboard.\n\n")
    #        alpha = np.real(alpha)
    #        # Q_factors = 0.5*(q_AC/alpha)*np.ones(n_modes_AC) # appropriate for alpha in [1/m]
    #        Q_factors = 0.5*(sim_AC.Omega_AC/alpha)*np.ones(n_modes_AC) # appropriate for alpha in [1/s]
    #        end = time.time()
    #        print("     time (sec.)", (end - start))
    #        print("remove me alpha", alpha)
    #        print("remove me Qs", Q_factors)
    #    else:
    #        # factor of a 1/2 because alpha is for power!
    #        # alpha [1/m] = Omega_AC/(2*vg*fixed_Q) = q_AC/fixed_Q
    #        # alpha [1/s] = vg * alpha [1/m]
    #        # alpha [1/s] = Omega_AC/(2*fixed_Q)
    #        # alpha = 0.5*(q_AC/fixed_Q)*np.ones(n_modes_AC) # appropriate for alpha in [1/m]
    #        alpha = 0.5*(sim_AC.Omega_AC/fixed_Q)*np.ones(n_modes_AC) # appropriate for alpha in [1/s]
    #        Q_factors = fixed_Q*np.ones(n_modes_AC)
    #        print("fixed q:", alpha, Q_factors)
    #
    #    linewidth_Hz = alpha/np.pi # SBS linewidth of each resonance in [Hz]

    # Calc Q_photoelastic Eq. 33

    alpha = simres_AC.alpha_t_AC_all()
    elastic_props = sim_AC.structure.elastic_props

    sim_AC.fem_evecs[2, :, :, :] = 0

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
        EM_ival_pump_fortran,
        EM_ival_Stokes_fortran,
        AC_ival_fortran,
        fem_ac.n_msh_el,
        fem_ac.n_msh_pts,
        fem_ac.elnd_to_mshpt,
        fem_ac.v_nd_xy,
        elastic_props.n_mats_ac,
        fem_ac.v_el_2_mat_idx,
        elastic_props.p_ijkl,
        q_AC,
        trimmed_EM_pump_field,
        trimmed_EM_Stokes_field,
        sim_AC.fem_evecs,
        relevant_eps_effs,
    )

    (Q_PE,) = process_fortran_return(
        resm, "finding linear element photoelastic couplings"
    )

    # Calc Q_moving_boundary Eq. 41
    typ_select_in = 1  # first element in relevant_eps_effs list, in fortan indexing
    if len(relevant_eps_effs) == 2:
        typ_select_out = 2
    elif typ_select_out is None:
        typ_select_out = -1
    print("\n Moving boundary calc")

    resm = nb_fortran.moving_boundary(
        sim_EM_pump.n_modes,
        sim_EM_Stokes.n_modes,
        sim_AC.n_modes,
        EM_ival_pump_fortran,
        EM_ival_Stokes_fortran,
        AC_ival_fortran,
        fem_ac.n_msh_el,
        fem_ac.n_msh_pts,
        fem_ac.elnd_to_mshpt,
        fem_ac.v_nd_xy,
        elastic_props.n_mats_ac,
        fem_ac.v_el_2_mat_idx,
        typ_select_in,
        typ_select_out,
        trimmed_EM_pump_field,
        trimmed_EM_Stokes_field,
        sim_AC.fem_evecs,
        relevant_eps_effs,
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
        x_tmp.append(simres.v_nd_xy[0, i])
        y_tmp.append(simres.v_nd_xy[1, i])
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
    elnd_to_mshpt = simres.elnd_to_mshpt.T
    v_nd_xy = simres.v_nd_xy.T

    sym_list = []

    for ival in range(len(simres.eigs_kz)):
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
                i_ex = elnd_to_mshpt[i_el, i_node] - 1
                # values
                v_x6p[i] = v_nd_xy[i_ex, 0]
                v_y6p[i] = v_nd_xy[i_ex, 1]
                v_Ex6p[i] = mode_fields[0, i_node, ival, i_el]
                v_Ey6p[i] = mode_fields[1, i_node, ival, i_el]
                i += 1

        # dense triangulation with unique points
        v_triang1p = []
        for i_el in np.arange(simres.n_msh_el):
            # triangles
            triangles = [
                [
                    elnd_to_mshpt[i_el, 0] - 1,
                    elnd_to_mshpt[i_el, 3] - 1,
                    elnd_to_mshpt[i_el, 5] - 1,
                ],
                [
                    elnd_to_mshpt[i_el, 1] - 1,
                    elnd_to_mshpt[i_el, 4] - 1,
                    elnd_to_mshpt[i_el, 3] - 1,
                ],
                [
                    elnd_to_mshpt[i_el, 2] - 1,
                    elnd_to_mshpt[i_el, 5] - 1,
                    elnd_to_mshpt[i_el, 4] - 1,
                ],
                [
                    elnd_to_mshpt[i_el, 3] - 1,
                    elnd_to_mshpt[i_el, 4] - 1,
                    elnd_to_mshpt[i_el, 5] - 1,
                ],
            ]
            v_triang1p.extend(triangles)

        # triangulations
        triang6p = matplotlib.tri.Triangulation(v_x6p, v_y6p, v_triang6p)
        triang1p = matplotlib.tri.Triangulation(
            v_nd_xy[:, 0], v_nd_xy[:, 1], v_triang1p
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
        #         {'i' : ival}, bbox_inches='tight')
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
        #         {'i' : ival}, bbox_inches='tight')
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
        #         {'i' : ival}, bbox_inches='tight')
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
        #         {'i' : ival}, bbox_inches='tight')
        # plt.close()

    return sym_list


def grad_u(dx, dy, u_mat, q_AC):
    """Take the gradient of field as well as of conjugate of field."""

    m_ux = u_mat[0]
    m_uy = u_mat[1]
    m_uz = u_mat[2]
    del_x_ux = np.gradient(m_ux, dx, axis=0)
    del_y_ux = np.gradient(m_ux, dy, axis=1)
    del_x_uy = np.gradient(m_uy, dx, axis=0)
    del_y_uy = np.gradient(m_uy, dy, axis=1)
    del_x_uz = np.gradient(m_uz, dx, axis=0)
    del_y_uz = np.gradient(m_uz, dy, axis=1)
    del_z_ux = 1j * q_AC * m_ux
    del_z_uy = 1j * q_AC * m_uy
    del_z_uz = 1j * q_AC * m_uz
    del_x_ux_star = np.gradient(np.conj(m_ux), dx, axis=0)
    del_y_ux_star = np.gradient(np.conj(m_ux), dy, axis=1)
    del_x_uy_star = np.gradient(np.conj(m_uy), dx, axis=0)
    del_y_uy_star = np.gradient(np.conj(m_uy), dy, axis=1)
    del_x_uz_star = np.gradient(np.conj(m_uz), dx, axis=0)
    del_y_uz_star = np.gradient(np.conj(m_uz), dy, axis=1)
    del_z_ux_star = -1j * q_AC * np.conj(m_ux)
    del_z_uy_star = -1j * q_AC * np.conj(m_uy)
    del_z_uz_star = -1j * q_AC * np.conj(m_uz)

    del_u_mat = np.array(
        [
            [del_x_ux, del_x_uy, del_x_uz],
            [del_y_ux, del_y_uy, del_y_uz],
            [del_z_ux, del_z_uy, del_z_uz],
        ]
    )
    del_u_mat_star = np.array(
        [
            [del_x_ux_star, del_x_uy_star, del_x_uz_star],
            [del_y_ux_star, del_y_uy_star, del_y_uz_star],
            [del_z_ux_star, del_z_uy_star, del_z_uz_star],
        ]
    )

    return del_u_mat, del_u_mat_star


def comsol_fields(data_file, n_points, ival=0):
    """Load Comsol field data on (assumed) grid mesh."""

    with open(data_file, "rt", encoding="ascii") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=" ")  # , quotechar='|')
        for _ in range(9):  # skip header
            next(spamreader)
        x_coord = []
        y_coord = []
        u_x = []
        u_y = []
        u_z = []
        for row in spamreader:
            row = [_f for _f in row if _f]
            row = [float(x) for x in row]
            x_coord.append(row[0])
            y_coord.append(row[1])
            u_x.append(row[(ival * 6) + 2] + 1j * row[(ival * 6) + 3])
            u_y.append(row[(ival * 6) + 4] + 1j * row[(ival * 6) + 5])
            u_z.append(row[(ival * 6) + 6] + 1j * row[(ival * 6) + 7])

    x_coord = np.array(x_coord).reshape(n_points, n_points)
    y_coord = np.array(y_coord).reshape(n_points, n_points)
    x_coord = np.swapaxes(x_coord, 0, 1)
    y_coord = np.swapaxes(y_coord, 0, 1)
    u_x = np.array(u_x).reshape(n_points, n_points)
    u_y = np.array(u_y).reshape(n_points, n_points)
    u_z = np.array(u_z).reshape(n_points, n_points)
    u_x = np.swapaxes(u_x, 0, 1)
    u_y = np.swapaxes(u_y, 0, 1)
    u_z = np.swapaxes(u_z, 0, 1)
    field_mat = np.array([u_x, u_y, u_z])

    return x_coord, y_coord, field_mat


def interp_py_fields(
    sim_EM_pump,
    sim_EM_Stokes,
    sim_AC,
    q_AC,
    n_points,
    EM_ival_pump=0,
    EM_ival_Stokes=0,
    AC_ival=0,
):
    """Interpolate fields from FEM mesh to square grid."""

    # Trim EM fields to non-vacuum area where AC modes are defined
    # n_modes_EM = sim_EM_pump.n_modes
    # n_modes_AC = sim_AC.n_modes
    n_msh_el_AC = sim_AC.n_msh_el
    ncomps = 3
    nnodes = 6
    trimmed_EM_field_p = np.zeros((ncomps, nnodes, n_msh_el_AC), dtype=complex)
    trimmed_EM_field_S = np.zeros((ncomps, nnodes, n_msh_el_AC), dtype=complex)
    trimmed_EM_n = np.zeros((1, nnodes, n_msh_el_AC), dtype=complex)
    for el in range(n_msh_el_AC):
        new_el = sim_AC.el_convert_tbl[el]
        for n in range(nnodes):
            for x in range(ncomps):
                trimmed_EM_field_p[x, n, el] = sim_EM_pump.fem_evecs[
                    x, n, EM_ival_pump, new_el
                ]
                trimmed_EM_field_S[x, n, el] = sim_EM_Stokes.fem_evecs[
                    x, n, EM_ival_Stokes, new_el
                ]
            trimmed_EM_n[0, n, el] = sim_EM_pump.ls_material[0, n, new_el]

    # field mapping
    x_tmp = []
    y_tmp = []
    for i in np.arange(sim_AC.n_msh_pts):
        x_tmp.append(sim_AC.v_nd_xy[0, i])
        y_tmp.append(sim_AC.v_nd_xy[1, i])
    x_min = np.min(x_tmp)
    x_max = np.max(x_tmp)
    y_min = np.min(y_tmp)
    y_max = np.max(y_tmp)
    # area = abs((x_max-x_min)*(y_max-y_min))
    n_pts_x = n_points
    n_pts_y = n_points
    v_x = np.zeros(n_pts_x * n_pts_y)
    v_y = np.zeros(n_pts_x * n_pts_y)
    i = 0
    for x in np.linspace(x_min, x_max, n_pts_x):
        for y in np.linspace(y_max, y_min, n_pts_y):
            v_x[i] = x
            v_y[i] = y
            i += 1
    v_x = np.array(v_x)
    v_y = np.array(v_y)

    # unrolling data for the interpolators
    elnd_to_mshpt = sim_AC.elnd_to_mshpt.T
    v_nd_xy = sim_AC.v_nd_xy.T

    # dense triangulation with multiple points
    v_x6p = np.zeros(6 * sim_AC.n_msh_el)
    v_y6p = np.zeros(6 * sim_AC.n_msh_el)
    v_ux6p = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_uy6p = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_uz6p = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_Ex6p_E_p = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_Ey6p_E_p = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_Ez6p_E_p = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_Ex6p_E_S = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_Ey6p_E_S = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_Ez6p_E_S = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    v_n = np.zeros(6 * sim_AC.n_msh_el, dtype=np.complex128)
    # v_triang6p = []

    i = 0
    for i_el in np.arange(sim_AC.n_msh_el):
        for i_node in np.arange(6):
            # index for the coordinates
            i_ex = elnd_to_mshpt[i_el, i_node] - 1
            # values
            v_x6p[i] = v_nd_xy[i_ex, 0]
            v_y6p[i] = v_nd_xy[i_ex, 1]
            v_ux6p[i] = sim_AC.fem_evecs[0, i_node, AC_ival, i_el]
            v_uy6p[i] = sim_AC.fem_evecs[1, i_node, AC_ival, i_el]
            v_uz6p[i] = sim_AC.fem_evecs[2, i_node, AC_ival, i_el]
            v_Ex6p_E_p[i] = trimmed_EM_field_p[0, i_node, i_el]
            v_Ey6p_E_p[i] = trimmed_EM_field_p[1, i_node, i_el]
            v_Ez6p_E_p[i] = trimmed_EM_field_p[2, i_node, i_el]
            v_Ex6p_E_S[i] = trimmed_EM_field_S[0, i_node, i_el]
            v_Ey6p_E_S[i] = trimmed_EM_field_S[1, i_node, i_el]
            v_Ez6p_E_S[i] = trimmed_EM_field_S[2, i_node, i_el]
            v_n[i] = trimmed_EM_n[0, i_node, i_el]
            i += 1

    xy = list(zip(v_x6p, v_y6p))
    grid_x, grid_y = np.mgrid[
        x_min : x_max : n_pts_x * 1j, y_min : y_max : n_pts_y * 1j
    ]
    # pump mode
    m_ReEx_E = interpolate.griddata(
        xy, v_Ex6p_E_p.real, (grid_x, grid_y), method="cubic"
    )
    m_ReEy_E = interpolate.griddata(
        xy, v_Ey6p_E_p.real, (grid_x, grid_y), method="cubic"
    )
    m_ReEz_E = interpolate.griddata(
        xy, v_Ez6p_E_p.real, (grid_x, grid_y), method="cubic"
    )
    m_ImEx_E = interpolate.griddata(
        xy, v_Ex6p_E_p.imag, (grid_x, grid_y), method="cubic"
    )
    m_ImEy_E = interpolate.griddata(
        xy, v_Ey6p_E_p.imag, (grid_x, grid_y), method="cubic"
    )
    m_ImEz_E = interpolate.griddata(
        xy, v_Ez6p_E_p.imag, (grid_x, grid_y), method="cubic"
    )
    m_Ex_E = m_ReEx_E + 1j * m_ImEx_E
    m_Ey_E = m_ReEy_E + 1j * m_ImEy_E
    m_Ez_E = m_ReEz_E + 1j * m_ImEz_E
    m_Ex_E = m_Ex_E.reshape(n_pts_x, n_pts_y)
    m_Ey_E = m_Ey_E.reshape(n_pts_x, n_pts_y)
    m_Ez_E = m_Ez_E.reshape(n_pts_x, n_pts_y)
    E_mat_p = np.array([m_Ex_E, m_Ey_E, m_Ez_E])
    # Stokes mode
    m_ReEx_E = interpolate.griddata(
        xy, v_Ex6p_E_S.real, (grid_x, grid_y), method="cubic"
    )
    m_ReEy_E = interpolate.griddata(
        xy, v_Ey6p_E_S.real, (grid_x, grid_y), method="cubic"
    )
    m_ReEz_E = interpolate.griddata(
        xy, v_Ez6p_E_S.real, (grid_x, grid_y), method="cubic"
    )
    m_ImEx_E = interpolate.griddata(
        xy, v_Ex6p_E_S.imag, (grid_x, grid_y), method="cubic"
    )
    m_ImEy_E = interpolate.griddata(
        xy, v_Ey6p_E_S.imag, (grid_x, grid_y), method="cubic"
    )
    m_ImEz_E = interpolate.griddata(
        xy, v_Ez6p_E_S.imag, (grid_x, grid_y), method="cubic"
    )
    m_Ex_E = m_ReEx_E + 1j * m_ImEx_E
    m_Ey_E = m_ReEy_E + 1j * m_ImEy_E
    m_Ez_E = m_ReEz_E + 1j * m_ImEz_E
    m_Ex_E = m_Ex_E.reshape(n_pts_x, n_pts_y)
    m_Ey_E = m_Ey_E.reshape(n_pts_x, n_pts_y)
    m_Ez_E = m_Ez_E.reshape(n_pts_x, n_pts_y)
    E_mat_S = np.array([m_Ex_E, m_Ey_E, m_Ez_E])
    # AC mode
    m_Reux = interpolate.griddata(xy, v_ux6p.real, (grid_x, grid_y), method="cubic")
    m_Reuy = interpolate.griddata(xy, v_uy6p.real, (grid_x, grid_y), method="cubic")
    m_Reuz = interpolate.griddata(xy, v_uz6p.real, (grid_x, grid_y), method="cubic")
    m_Imux = interpolate.griddata(xy, v_ux6p.imag, (grid_x, grid_y), method="cubic")
    m_Imuy = interpolate.griddata(xy, v_uy6p.imag, (grid_x, grid_y), method="cubic")
    m_Imuz = interpolate.griddata(xy, v_uz6p.imag, (grid_x, grid_y), method="cubic")
    m_ux = m_Reux + 1j * m_Imux
    m_uy = m_Reuy + 1j * m_Imuy
    m_uz = m_Reuz + 1j * m_Imuz
    m_ux = m_ux.reshape(n_pts_x, n_pts_y)
    m_uy = m_uy.reshape(n_pts_x, n_pts_y)
    m_uz = m_uz.reshape(n_pts_x, n_pts_y)
    u_mat = np.array([m_ux, m_uy, m_uz])

    dx = grid_x[-1, 0] - grid_x[-2, 0]
    dy = grid_y[0, -1] - grid_y[0, -2]
    del_u_mat, del_u_mat_star = grad_u(dx, dy, u_mat, q_AC)

    m_Ren = interpolate.griddata(xy, v_n.real, (grid_x, grid_y), method="cubic")
    m_Imn = interpolate.griddata(xy, v_n.imag, (grid_x, grid_y), method="cubic")
    m_n = m_Ren + 1j * m_Imn
    m_n = m_n.reshape(n_pts_x, n_pts_y)

    return (
        n_pts_x,
        n_pts_y,
        dx,
        dy,
        E_mat_p,
        E_mat_S,
        u_mat,
        del_u_mat,
        del_u_mat_star,
        m_n,
    )


def grid_integral(
    m_n,
    sim_AC_structure,
    sim_AC_Omega_AC,
    n_pts_x,
    n_pts_y,
    dx,
    dy,
    E_mat_p,
    E_mat_S,
    u_mat,
    del_u_mat,
    del_u_mat_star,
    AC_ival,
):
    """Quadrature integration of AC energy density, AC loss (alpha), and PE gain."""

    # AC energy density integral
    F_AC_energy = 0
    for i in range(3):
        integrand_AC = np.conj(u_mat[i]) * u_mat[i] * sim_AC_structure.rho
        # do a 1-D integral over every row
        I_en = np.zeros(n_pts_x)
        for r in range(n_pts_x):
            I_en[r] = np.trapz(np.real(integrand_AC[r, :]), dx=dy)
        # then an integral over the result
        F_AC_energy += np.trapz(I_en, dx=dx)
        # Adding imag comp
        I_en = np.zeros(n_pts_x)
        for r in range(n_pts_x):
            I_en[r] = np.trapz(np.imag(integrand_AC[r, :]), dx=dy)
        F_AC_energy += 1j * np.trapz(I_en, dx=dx)
    energy_py = 2 * F_AC_energy * sim_AC_Omega_AC[AC_ival] ** 2

    # AC loss (alpha) integral
    F_alpha = 0
    for i in range(3):
        for k in range(3):
            for l in range(3):
                for j in range(3):
                    integrand = (
                        del_u_mat[i, j]
                        * del_u_mat_star[k, l]
                        * sim_AC_structure.elastic_props.eta_ijkl[i, j, k, l]
                    )
                    I_en = np.zeros(n_pts_x)
                    for r in range(n_pts_x):
                        I_en[r] = np.trapz(np.real(integrand[r, :]), dx=dy)
                    F_alpha += np.trapz(I_en, dx=dx)
                    I_en = np.zeros(n_pts_x)
                    for r in range(n_pts_x):
                        I_en[r] = np.trapz(np.imag(integrand[r, :]), dx=dy)
                    F_alpha += 1j * np.trapz(I_en, dx=dx)
    alpha_py = np.real(F_alpha * sim_AC_Omega_AC[AC_ival] ** 2 / energy_py)

    # PE gain integral

    F_PE = 0
    for i in range(3):
        for k in range(3):
            for l in range(3):
                for j in range(3):
                    # integrand_PE = relevant_eps_effs[0]**2 * E_mat_p[j]*np.conj(E_mat_S[i])*sim_AC_structure.elastic_props.p_ijkl[i,j,k,l]*del_u_mat_star[k,l]
                    integrand_PE = (
                        m_n**4
                        * E_mat_p[j]
                        * np.conj(E_mat_S[i])
                        * sim_AC_structure.elastic_props.p_ijkl[i, j, k, l]
                        * del_u_mat_star[k, l]
                    )
                    I_en = np.zeros(n_pts_x)
                    for r in range(n_pts_x):
                        I_en[r] = np.trapz(np.real(integrand_PE[r, :]), dx=dy)
                    F_PE += np.trapz(I_en, dx=dx)
                    I_en = np.zeros(n_pts_x)
                    for r in range(n_pts_x):
                        I_en[r] = np.trapz(np.imag(integrand_PE[r, :]), dx=dy)
                    F_PE += 1j * np.trapz(I_en, dx=dx)
    Q_PE_py = F_PE * SI_permittivity_eps0

    return energy_py, alpha_py, Q_PE_py


def gain_python(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC, comsol_data_file, comsol_ivals=1
):
    """Calculate interaction integrals and SBS gain in python.
    Load in acoustic mode displacement and calculate gain from this also.
    """

    print("gain python is out of action")
    return
    n_modes_EM = sim_EM_pump.n_modes
    # n_modes_AC = sim_AC.n_modes
    EM_ival_pump = 0
    EM_ival_Stokes = 0

    n_points = 100
    n_points_comsol_data = 100

    # relevant_eps_effs =[]
    # for el_typ in range(sim_EM_pump.structure.n_mats_em):
    #     if el_typ+1 in sim_AC.typ_el_AC:
    #         relevant_eps_effs.append(sim_EM_pump.v_refindexn[el_typ]**2)

    energy_py = np.zeros(comsol_ivals, dtype=np.complex128)
    alpha_py = np.zeros(comsol_ivals)
    Q_PE_py = np.zeros(
        (len(sim_EM_pump.eigs_kz), len(sim_EM_Stokes.eigs_kz), comsol_ivals),
        dtype=np.complex128,
    )
    energy_comsol = np.zeros(comsol_ivals, dtype=np.complex128)
    alpha_comsol = np.zeros(comsol_ivals)
    Q_PE_comsol = np.zeros(
        (len(sim_EM_pump.eigs_kz), len(sim_EM_Stokes.eigs_kz), comsol_ivals),
        dtype=np.complex128,
    )

    for AC_ival in range(comsol_ivals):  # Comsol data only contains some AC modes
        # Interpolate NumBAT FEM fields onto grid
        (
            n_pts_x,
            n_pts_y,
            dx,
            dy,
            E_mat_p,
            E_mat_S,
            u_mat,
            del_u_mat,
            del_u_mat_star,
            m_n,
        ) = interp_py_fields(
            sim_EM_pump,
            sim_EM_Stokes,
            sim_AC,
            q_AC,
            n_points,
            EM_ival_pump=EM_ival_pump,
            EM_ival_Stokes=EM_ival_Stokes,
            AC_ival=AC_ival,
        )

        # Carry out integration
        (
            energy_py[AC_ival],
            alpha_py[AC_ival],
            Q_PE_py[EM_ival_pump, EM_ival_Stokes, AC_ival],
        ) = grid_integral(
            m_n,
            sim_AC.structure,
            sim_AC.Omega_AC,
            n_pts_x,
            n_pts_y,
            dx,
            dy,
            E_mat_p,
            E_mat_S,
            u_mat,
            del_u_mat,
            del_u_mat_star,
            AC_ival,
        )

        # Load Comsol FEM fields onto grid - acoustic displacement fields
        x_coord, y_coord, u_mat_comsol = comsol_fields(
            comsol_data_file, n_points_comsol_data, ival=AC_ival
        )
        dx_comsol = x_coord[-1, 0] - x_coord[-2, 0]
        dy_comsol = y_coord[0, -1] - y_coord[0, -2]
        del_u_mat_comsol, del_u_mat_star_comsol = grad_u(dx, dy, u_mat_comsol, q_AC)

        # Carry out integration
        n_pts_x_comsol = n_points_comsol_data
        n_pts_y_comsol = n_points_comsol_data
        (
            energy_comsol[AC_ival],
            alpha_comsol[AC_ival],
            Q_PE_comsol[EM_ival_pump, EM_ival_Stokes, AC_ival],
        ) = grid_integral(
            m_n,
            sim_AC.structure,
            sim_AC.Omega_AC,
            n_pts_x_comsol,
            n_pts_y_comsol,
            dx_comsol,
            dy_comsol,
            E_mat_p,
            E_mat_S,
            u_mat_comsol,
            del_u_mat_comsol,
            del_u_mat_star_comsol,
            AC_ival,
        )

    # Note this is only the PE contribution to gain.
    gain_PE_py = (
        2
        * sim_EM_pump.omega_EM
        * sim_AC.Omega_AC[:comsol_ivals]
        * np.real(Q_PE_py * np.conj(Q_PE_py))
    )
    normal_fact_py = np.zeros((n_modes_EM, n_modes_EM, comsol_ivals), dtype=complex)
    gain_PE_comsol = (
        2
        * sim_EM_pump.omega_EM
        * sim_AC.Omega_AC[:comsol_ivals]
        * np.real(Q_PE_comsol * np.conj(Q_PE_comsol))
    )
    normal_fact_comsol = np.zeros((n_modes_EM, n_modes_EM, comsol_ivals), dtype=complex)
    for i in range(n_modes_EM):
        P1 = sim_EM_pump.EM_mode_power[i]
        for j in range(n_modes_EM):
            P2 = sim_EM_Stokes.EM_mode_power[j]
            for k in range(comsol_ivals):
                P3_py = energy_py[k]
                normal_fact_py[i, j, k] = P1 * P2 * P3_py * alpha_py[k]
                P3_comsol = energy_comsol[k]
                normal_fact_comsol[i, j, k] = P1 * P2 * P3_comsol * alpha_comsol[k]
    SBS_gain_PE_py = np.real(gain_PE_py / normal_fact_py)
    SBS_gain_PE_comsol = np.real(gain_PE_comsol / normal_fact_comsol)

    return SBS_gain_PE_py, alpha_py, SBS_gain_PE_comsol, alpha_comsol
