# Copyright (C) 2017-2025  Michael Steel

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
import sys
from typing import Any, Callable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
from matplotlib.axes import Axes as MplAxes

import numpy as np
from numpy.typing import NDArray

import scipy.optimize as sciopt

from plotting.plottools import get_rgb_for_poln
import numtools.rootbracket as ntrb
from nbtypes import SI_um, NormalisationConstants, SI_GHz

from nbanalytic.elasticmodeplots import _field_plot_1d, _field_plot_2d, ModeFunction1D



#######################################################################################

g_norms = NormalisationConstants()

twopi = 2 * math.pi

def elasticfreeslab_Lamb_chareq_even(Omega: float, Vs: float, Vl: float, wid: float, q: float) -> float:
    """Dispersion relation for TE slab optical waveguide in normalised units.


    Agrees with Mathematica notes in christoffel_2.nb
    """

    ONE = 1.0 + 0.0j
    kappa_s = np.sqrt(ONE * (Omega / Vs) ** 2 - q**2)
    kappa_l = np.sqrt(ONE * (Omega / Vl) ** 2 - q**2)

    # lhs = np.tan(kappa_l * wid/2)/np.tan(kappa_s * wid/2)
    # rhs = - (q**2-kappa_s**2)**2/(4 *q**2 *  kappa_l * kappa_s)

    lhs = np.tan(kappa_l * wid / 2) * (4 * q**2 * kappa_l * kappa_s)
    rhs = -((q**2 - kappa_s**2) ** 2) * np.tan(kappa_s * wid / 2)

    eq = wid**4 * (lhs - rhs)
    if abs(np.imag(eq)) > abs(np.real(eq)):
        return np.imag(eq)
    else:
        return np.real(eq)


def elasticfreeslab_Lamb_chareq_odd(Omega: float, Vs: float, Vl: float, wid: float, q: float) -> float:
    """Dispersion relation for TE slab optical waveguide in normalised units."""
    ONE = 1.0 + 0.0j
    kappa_s = np.sqrt(ONE * (Omega / Vs) ** 2 - q**2)
    kappa_l = np.sqrt(ONE * (Omega / Vl) ** 2 - q**2)

    # lhs = np.tan(kappa_s * wid/2)/np.tan(kappa_l * wid/2)
    # rhs = - (q**2-kappa_s**2)**2/(4 *q**2 * kappa_l * kappa_s)
    lhs = np.tan(kappa_s * wid / 2) * (4 * q**2 * kappa_l * kappa_s)
    rhs = -((q**2 - kappa_s**2) ** 2) * np.tan(kappa_l * wid / 2)

    eq = wid**4 * (lhs - rhs)
    if abs(np.imag(eq)) > abs(np.real(eq)):
        return np.imag(eq)
    else:
        return np.real(eq)



def ombrak_is_at_tan_resonance(ombrak: Sequence[float], q: float, Vs: float, Vl: float, width: float) -> bool:
    """Check for bracketed frequencies where tan(kappa_i w/2) = tan((Omega/V_i)^2 - q^2) blow ups.

    Looking for Om brackets that surround (Omega/V_i)^2 - q^2)w/2 = (2n+1) pi/2, so (Omega/V_i)^2 - q^2)w/pi = (2n+1)

    """
    widonpi = width/(np.pi)
    qsq = q**2
    if (ombrak[0] / Vs) ** 2 - qsq > 0:
        sarg0 = np.sqrt((ombrak[0] / Vs) ** 2 - qsq) * widonpi

        sarg1 = np.sqrt((ombrak[1] / Vs) ** 2 - qsq) * widonpi

        if math.ceil(sarg0) == math.floor(sarg1) and math.ceil(sarg0) % 2 == 1:  # an odd multiple of pi/2
                return True

    if (ombrak[0] / Vl) ** 2 - qsq > 0:
        larg0 = np.sqrt((ombrak[0] / Vl) ** 2 - qsq) * widonpi

        larg1 = np.sqrt((ombrak[1] / Vl) ** 2 - qsq) * widonpi

        if math.ceil(larg0) == math.floor(larg1) and math.ceil(larg0) % 2 == 1:  # an odd multiple of pi/2
            return True
    return False

def qbrak_is_at_tan_resonance(Omega: float, qbrak: Sequence[float], Vs: float, Vl: float, width: float) -> bool:
    """Check for bracketed wavenumbers where tan(kappa_i w/2) = tan((Omega/V_i)^2 - q^2) blow ups.

    Looking for q brackets that surround (Omega/V_i)^2 - q^2)w/2 = (2n+1) pi/2, so  (Omega/V_i)^2 - q^2)w/pi = (2n+1)
    """

    widonpi = width/(np.pi)
    OmonVs_sq = (Omega / Vs) ** 2

    if OmonVs_sq - qbrak[1] ** 2 > 0:
        sarg0 = np.sqrt(OmonVs_sq - qbrak[0] ** 2) * widonpi
        sarg1 =  np.sqrt(OmonVs_sq - qbrak[1] ** 2) * widonpi

        if math.ceil(sarg0) == math.floor(sarg1) and math.ceil(sarg0) % 2 == 1 :  # an odd multiple of pi/2
            return True

    OmonVl_sq = (Omega / Vl) ** 2
    if OmonVl_sq - qbrak[1] ** 2 > 0:
        larg0 = np.sqrt(OmonVl_sq - qbrak[0] ** 2) * widonpi
        larg1 = np.sqrt(OmonVl_sq - qbrak[1] ** 2) * widonpi

        if math.ceil(larg0) == math.floor(larg1) and math.ceil(larg0) % 2 == 1 :  # an odd multiple of pi/2
            return True
    return False


def _find_Lamb_q_brackets_smart(self, Omega: float, Vs, Vl, width,
                                    drfunc: Callable[[float], float]) -> List[Tuple[float, float]]:
        # attempts to predict the crossings but doesn't get them all, becuse not all are about tan blowups

        # qlo = 1e-9
        # qhi=Omega/self._Vbulk_shear * 2

        # need to avoid blowups of tan functions
        # kappa_s w/2 = (2n+1) pi/2
        # kappa_l w/2 = (2n+1) pi/2
        hi_ns = np.ceil(width * Omega / (twopi * Vs) - 0.5)
        hi_nl = np.ceil(width * Omega / (twopi * Vl) - 0.5)

        bad_qs = (
            (Omega / Vs) ** 2
            - ((np.arange(0, hi_ns) + 0.5) * twopi / width) ** 2
        ) ** 0.5
        bad_ql = (
            (Omega / Vl) ** 2
            - ((np.arange(0, hi_nl) + 0.5) * twopi / width) ** 2
        ) ** 0.5

        # print(' hins', Omega, self._Vbulk_shear,hi_ns, hi_nl, self.width*Omega/(2*np.pi*self._Vbulk_long)-.5, bad_qs, bad_ql)

        bad_q = np.append(bad_qs, bad_ql)
        bad_q = np.sort(bad_q)
        qbraks = []
        # tollo = 1-1e-7
        # tolhi = 1+1e-7

        # print('Doing Omega', Omega)
        # qsteps=2000
        # qlo = 1e-9
        # qhi=Omega/self._Vbulk_shear * 2
        # v_q = np.linspace(qlo, qhi,  qsteps)
        # v_disprel = list(map(drfunc, v_q))
        # fig,axs=plt.subplots()
        # axs.plot(v_q*1e-6, v_disprel)
        # axs.plot(v_q*1e-6,0*v_q,':')
        # axs.set_ylim(-1000,1000)
        # #axs.set_xlim(0,20)

        # # carve out the bad points
        # for bq in bad_q:
        #     qbraks.append([qlo, bq*tollo])
        #     qlo = bq*tolhi
        #     axs.plot(qlo*self.width,0,'x')
        # qbraks.append([qlo, qhi])

        # fig.savefig('qscan.png')
        # sys.exit(0)

        return qbraks


def _find_Lamb_q_brackets(Omega: float, Vs, Vl, width,
                          drfunc: Callable[[float], float],
                          q_bracket_steps: int = 1000) -> List[Tuple[float, float]]:
        # brute force bracketing
        qsteps = q_bracket_steps

        qlo = 1e-9
        qhi = Omega / Vs * 4
        v_q = np.linspace(qlo, qhi, qsteps)

        v_disprel = np.array(list(map(drfunc, v_q)))

        qbraks = []

        is_small = 0.1 * max(abs(v_disprel))
        fi_m1 = v_disprel[0]
        for i in range(1, len(v_disprel) - 1):
            fi = v_disprel[i]
            if fi * fi_m1 < 0:  # a possible crossing?
                if abs(fi) < is_small and abs(fi_m1) < is_small:  # legit
                    qbrak = [v_q[i - 1], v_q[i]]
                    # remove erroneous  brackets due to tan functions
                    if not qbrak_is_at_tan_resonance(Omega, qbrak, Vs, Vl, width ):
                        qbraks.append(qbrak)
                # else: # a tan function infinity, update but don't keep
                #    pass
            fi_m1 = fi


        return qbraks


def _find_Lamb_Om_brackets_at_q(q: float, Vs: float, Vl: float, width: float,
                               drfunc: Callable[[float], float],
                               om_max: float = 200e9, om_bracket_steps: int = 1000) -> List[Tuple[float, float]]:

        # brute force bracketing
        omsteps = om_bracket_steps

        omlo = 1e-9
        omhi = om_max

        use_ntrb = False
        if use_ntrb:
            d_rbr = {
            "threshold_factor": 100,
            "min_abs_value_for_sign_change": 1e-6,
            "ignore_large_jumps": True,
            }
            ombraks, removed = ntrb.robust_bracket_roots(drfunc, omlo, omhi, om_bracket_steps, **d_rbr)
            return ombraks


        v_om = np.linspace(omlo, omhi, omsteps)
        v_disprel = np.array(list(map(drfunc, v_om)))
        np.set_printoptions(precision=15)

        ombraks = []
        removeds = []

        is_small = 0.1 * max(abs(v_disprel))

        # We skip over the very first possible bracket.
        # If drfunc(om=\eps) \approx 0, then if we get drfunc(om=\eps)=\pm 1e-15,
        # we can get a spurious root. So just don't look there.
        i_st = 1

        fi_m1 = v_disprel[i_st]
        for i in range(i_st + 1, len(v_disprel) - 1):
            fi = v_disprel[i]
            if fi * fi_m1 < 0:  # a possible crossing?
                if abs(fi) < is_small and abs(fi_m1) < is_small:  # legit
                    ombrak = [v_om[i - 1], v_om[i]]
                    if not ombrak_is_at_tan_resonance(ombrak, q, Vs, Vl, width):
                        ombraks.append(ombrak)
                    else:
                        removeds.append(ombrak)
                # else: # a tan function infinity, update but don't keep
                #    pass
            fi_m1 = fi

        # plotscan=True
        plotscan = False

        # omnormfac = self.width/self._Vbulk_shear/np.pi   # genuine normalised angular frequency
        # omnormlab = r'$\Omega w/(\pi V_s)$'

        omnormfac = 1 / (twopi * SI_GHz)  # actual frequency in GHz
        omnormlab = r"$\Omega/(2\pi)$ [GHz]"

        if plotscan:
            # print(ombraks)
            fig, axs = plt.subplots()
            v_omnorm = v_om * omnormfac
            # axs.plot(v_omnorm, np.log10(abs(v_disprel)))
            axs.plot(v_omnorm, v_disprel)

            axs.plot(v_omnorm, 0 * v_om, ":")
            # axs.set_ylim(-1*is_small, 1*is_small)
            axs.set_ylim(-1e3, 1e3)
            # axs.set_xlim(1.3,1.5)
            for iob, ob in enumerate(ombraks):
                obnorm = ob[0] * omnormfac
                lab = "" if iob else "genuine roots"
                axs.plot(obnorm, 0, "x", ms=20, label=lab)
            for iob, ob in enumerate(removeds):
                obnorm = ob[0] * omnormfac
                lab = "" if iob else "false roots"
                axs.plot(obnorm, 0, "o", mfc="none", ms=20, label=lab)

            axs.legend()
            axs.set_xlabel(omnormlab)
            axs.set_ylabel(r"dispersion relation $(\Omega, q)$")

            fig.savefig("omscan.png")
            sys.exit(0)

        return ombraks



class ElasticIsotropicFreeSlab:
    """Elastic slab waveguide solver for isotropic materials.
    Finds the dispersion of the Lamb modes of isolated flat plate
    All units are in SI.

    Dispersion relation has even and odd parts

    Even:
      tan(kappa_l w/2)/tan(kappa_s w/2) = - (q^2-kappa_s^2)/(4 q^2 kappa_l kappa_s)

    Odd:
      tan(kappa_s w/2)/tan(kappa_l w/2) = - (q^2-kappa_s^2)/(4 q^2 kappa_l kappa_s)

    where kappa_s^2 = (Omega/V_s)^2 - q^2
    where kappa_l^2 = (Omega/V_l)^2 - q^2

    The two sides switch from real and imag at the same time so all makes sense at all values of q
    """

    def __init__(
        self,
        material: Any,
        thickness_SI: float,
        max_omega: float = 100e9,
        omega_bracketing_resolution: int = 1000,
        q_bracketing_resolution: int = 1000,
    ) -> None:
        self._material = material  # slab material
        self._Vbulk_shear = material.Vac_shear()
        self._Vbulk_long = material.Vac_longitudinal()
        self._Vl_norm = self._Vbulk_long/self._Vbulk_shear

        self.width = thickness_SI

        self._max_omega = max_omega
        self.omega_bracketing_resolution = omega_bracketing_resolution
        self.q_bracketing_resolution = q_bracketing_resolution




    def find_Lamb_Omegas_at_q(self, q: float, max_modes: int, even_modes: bool = True) -> NDArray:
        """Find disperison of Lamb"""
        # q must be smaller than Omega/V_l

        def dreven(om: float) -> float:
            return elasticfreeslab_Lamb_chareq_even(om, self._Vbulk_shear, self._Vbulk_long, self.width, q)

        def drodd(om: float) -> float:
            return elasticfreeslab_Lamb_chareq_odd(om, self._Vbulk_shear, self._Vbulk_long, self.width, q)

        if even_modes:
            drfunc = dreven
        else:
            drfunc = drodd

        omsols = []

        ombraks = _find_Lamb_Om_brackets_at_q(q, self._Vbulk_shear, self._Vbulk_long, self.width,
                                              drfunc, self._max_omega,
                                              self.omega_bracketing_resolution)

        for ombrak in ombraks:
            if drfunc(ombrak[0]) * drfunc(ombrak[1]) > 0:
                print("False bracket", ombrak, drfunc(ombrak[0]), drfunc(ombrak[1]))
                continue

            om_sol = sciopt.root_scalar(drfunc, bracket=ombrak)
            if not om_sol.converged:
                raise ValueError(om_sol.flag)
            else:
                omsols.append(om_sol.root)

            if len(omsols) == max_modes:
                break

        rv = np.array(omsols)
        return rv

    def find_Lamb_qs_at_Omega(self, Omega: float, max_modes: int, even_modes: bool = True) -> NDArray:
        """Finds up to max_modes values of q satisfying the dispersion relation D_Omega(q).

        Solutions are found in decreasing order from largest to smallest q, since for fixed Omega, larger q means lower mode.

        """
        # q must be smaller than Omega/V_l

        if even_modes:
            def drfunc(q: float) -> float:
                return elasticfreeslab_Lamb_chareq_even(Omega, self._Vbulk_shear, self._Vbulk_long, self.width, q)
        else:
            def drfunc(q: float) -> float:
                return elasticfreeslab_Lamb_chareq_odd(Omega, self._Vbulk_shear, self._Vbulk_long, self.width, q)

        qsols = []

        qbraks = self._find_Lamb_q_brackets(Omega, self._Vbulk_shear, self._Vbulk_long, self.width,
                                            drfunc, self.q_bracketing_resolution)

        #if not len(qbraks):
        #    print("No brackets for Omega ", Omega)
        # Search in reverse order to find highest qs first
        qbraks.reverse()
        for qbrak in qbraks:
            if drfunc(qbrak[0]) * drfunc(qbrak[1]) > 0:
                print("False bracket", qbrak, drfunc(qbrak[0]), drfunc(qbrak[1]))
                continue

            q_sol = sciopt.root_scalar(drfunc, bracket=qbrak)
            if not q_sol.converged:
                raise ValueError(q_sol.flag)
            else:
                qsols.append(q_sol.root)

            if len(qsols) == max_modes:
                break


        return np.array(qsols)


    def find_Lamb_dispersion_for_q_bands_at_om(self, v_Omega: Sequence[float], max_modes: int, even_modes: bool = True) -> NDArray:
        """Returns q(Omega) for each band"""

        m_qs = np.zeros([len(v_Omega), max_modes])

        for iOm, Om in enumerate(v_Omega):
            qsols = self.find_Lamb_qs_at_Omega(Om, max_modes, even_modes)

            # We might have less than max_modes, or even none
            if len(qsols):
                m_qs[iOm, : len(qsols)] = qsols

        return m_qs

    def find_Lamb_dispersion_for_om_bands_at_q(self, v_q: Sequence[float], max_modes: int, even_modes: bool = True) -> NDArray:
        """Returns Omega(q) for each band"""

        m_Om = np.zeros([len(v_q), max_modes])

        for iq, q in enumerate(v_q):
            omsols = self.find_Lamb_Omegas_at_q(q, max_modes, even_modes)

            for iom, om in enumerate(omsols):
                m_Om[iq, iom] = om

        return m_Om

    def find_Lamb_q_at_Omega_oddmode_zero_asymptotic_small_qw(self, Omega: float) -> float:
        om = Omega * self.width/self._Vbulk_shear
        rv = self._Vbulk_long_norm
        vsq = om * np.sqrt((1-1/rv**2)/3) + om*om * (1/(6*rv**2)-9/40)
        V_asy = self._Vbulk_shear*np.sqrt(vsq)
        q = Omega/ V_asy
        return q

    def find_Lamb_Omega_at_q_oddmode_zero_asymptotic_small_qw(self, q: float, even_modes: bool = True) -> float:
        # returns q(Omega) for small qw
        if even_modes:
                return 0  # don't know answer yet
        else:
            X = q*self.width
            rv = self._Vbulk_long_norm
            V_asy = self._Vbulk_shear*np.sqrt( X**2*(1-rv**2)/3 )
            Omega = q * V_asy
            return Omega

    def find_SH_dispersion_for_q_bands_at_om(self, v_Omega: Sequence[float], max_modes: int, col_array: bool = False) -> Tuple[NDArray, Optional[NDArray]]:
        # Returns q(Omega) for each band
        # Bands are numbered from 0?
        m_qs = np.zeros([len(v_Omega), max_modes])
        m_col = None

        for m in range(max_modes):
            # Only take square root of positive numbers
            qsq = (v_Omega / self._Vbulk_shear) ** 2 - (m * np.pi / self.width) ** 2
            qsq = qsq * (qsq > 0)

            # m_qs[:,m] = np.where(qsq>0, np.sqrt(qsq), 0)
            m_qs[:, m] = np.sqrt(qsq)

        if col_array:  # All SH modes are purely x polarised by definition
            m_col = np.zeros([len(v_Omega), max_modes, 3])
            rgb = get_rgb_for_poln(1, 0, 0)
            m_col[:, :] = rgb

        return m_qs, m_col


    def find_Rayleigh_dispersion_for_v(self, v_Omega: Sequence[float], col_array: bool = False) -> Tuple[NDArray, Optional[NDArray]]:
        # find Rayleigh wavenumber v_q of Rayleigh mode for np array v_Omega
        vR = self._material.Vac_Rayleigh()
        v_qR = v_Omega/vR

        v_col = None

        if col_array:
            # Is Rayleigh mode polarisation frequency dependent?
            v_col = np.zeros([len(v_Omega), 3])
            rgb = self._get_rgb_for_poln(1, 0, 0)
            v_col[:] = rgb

        return v_qR, v_col

    def find_Rayleigh_profile_1d(self, Omega: float, depth: float, npts: int = 200) -> Tuple[NDArray, NDArray]:

        Vs = self._Vbulk_shear
        Vl = self._Vbulk_long
        Vr = self._material.Vac_Rayleigh()

        qR = Omega/Vr

        v_x = np.linspace(-depth, 0, npts)
        m_uxyz = np.zeros([npts,3], dtype=np.complex64)

        C_ONE = 1+0j
        gam_s = np.sqrt(C_ONE*(qR**2-(Omega/Vs)**2))
        gam_l = np.sqrt(C_ONE*(qR**2-(Omega/Vl)**2))
        exp_glx = np.exp(gam_l*v_x)
        exp_gsx = np.exp(gam_s*v_x)

        v_uy = ((exp_glx - 2*exp_gsx)*qR**2 + gam_s**2 * exp_glx)/(2*gam_s)
        v_uz = 1j*qR*( -2*gam_l*gam_s * exp_gsx + exp_glx*(qR**2+gam_s**2) )/(2*gam_s*gam_l)

        m_uxyz[:,1] = v_uy
        m_uxyz[:,2] = v_uz

        max_u = np.max(np.abs(m_uxyz))
        if max_u > 0:
            m_uxyz /= np.max(np.abs(m_uxyz))

        mode_prof = ModeFunction1D(Omega, Vr, v_x, m_uxyz)

        return mode_prof


    def plot_Rayleigh_profile_1d(self, Omega: float, depth: float, npts: int = 200, ax: Optional[MplAxes] = None, legend: bool = False) -> None:
        mode_prof = self.find_Rayleigh_profile_1d(Omega, depth, npts)
        fname = mode_prof.plot_mode_profile_1d(ax=ax, legend=legend)
        return fname


    def plot_Rayleigh_profile_2d(self, Omega: float, depth: float, zperiods: int = 1, npts: int = 20, ax: Optional[MplAxes] = None,
                                 displacement_scale: float = 0.05, use_arrows: bool = False) -> None:

        mode_prof = self.find_Rayleigh_profile_1d(Omega, depth, npts)
        fname = mode_prof.plot_mode_profile_2d(ax=ax, npts=npts, zperiods=zperiods,
                                      displacement_scale=displacement_scale  )
        return fname


    def plot_slab_mode_profile_1d(self, mode: int, Omega: float, q: float, is_even: bool, npts: int = 200, ax: Optional[MplAxes] = None, legend: bool = False) -> float:

        mode_prof = self.find_mode_profile_1d(Omega, q, npts=npts, even_mode=is_even)
        fname = mode_prof.plot_mode_profile_1d(ax=ax, npts=npts, legend=legend)
        return fname


    def plot_slab_mode_profile_2d(self, mode: int, Omega: float, q: float, is_even: bool, npts: int = 30, zperiods: int = 1,
                                  ax: Optional[MplAxes] = None, displacement_scale: float = 0.1) -> None:

        mode_prof = self.find_mode_profile_1d(Omega, q, npts=npts, even_mode=is_even)
        fname = mode_prof.plot_mode_profile_2d(ax=ax, npts=npts, zperiods=zperiods,
                                      displacement_scale=displacement_scale )
        return fname


    def find_mode_profile_1d(self, Omega: float, q: float, npts: int = 200, even_mode: bool = True) -> Tuple[NDArray, NDArray, float]:
        # Solutions from Auld Vol 2, 10.22

        hwid = self.width/2

        v_x = np.linspace(-self.width/2, self.width/2, npts)
        m_uxyz = np.zeros([npts,3], dtype=np.complex64)

        C_ONE = 1+0j
        kap_l = np.sqrt(C_ONE*((Omega/self._Vbulk_long)**2 - q**2))
        kap_s = np.sqrt(C_ONE*((Omega/self._Vbulk_shear)**2 - q**2))

        delq2kaps2 = q**2 - kap_s**2

        #qrat_y= (q**2-kap_s**2)/(2*kap_s)
        #qrat_z= (q**2-kap_s**2)/(2*q)


        cos_klw = np.cos(kap_l*hwid)
        cos_ksw = np.cos(kap_s*hwid)

        cos_klx = np.cos(kap_l*v_x)
        sin_klx = np.sin(kap_l*v_x)
        cos_ksx = np.cos(kap_s*v_x)
        sin_ksx = np.sin(kap_s*v_x)


        if even_mode:
           # vely =   - (kap_l     * np.sin(kap_l*v_x) * np.cos(kap_s*hwid)
           #             + qrat_y * np.sin(kap_s*v_x) * np.cos(kap_l*hwid) )
           # velz = -1j*( q       * np.cos(kap_l*v_x) * np.cos(kap_s*hwid)
           #             - qrat_z * np.cos(kap_s*v_x) * np.cos(kap_l*hwid) )
           vely = (q/kap_s) * ( (2*kap_s*kap_l/delq2kaps2) * cos_ksw*sin_klx + cos_klw*sin_ksx)

           velz = -1j*( (2*q**2/delq2kaps2) * cos_ksw * cos_klx - cos_klw * cos_ksx )



        else:
            #vely =   (kap_l     * np.cos(kap_l*v_x) * np.sin(kap_s*hwid)
            #          + qrat_y * np.cos(kap_s*v_x) * np.sin(kap_l*hwid) )
            #velz = -1j*( q       * np.sin(kap_l*v_x) * np.sin(kap_s*hwid)
            #            - qrat_z * np.sin(kap_s*v_x) * np.sin(kap_l*hwid) )

           vely = (1j / (2* q * kap_s)) * ( 2*q**2 * cos_klw*cos_ksx - delq2kaps2 * cos_ksw*cos_klx)

           velz =  cos_klw*sin_ksx + cos_ksw*sin_klx * delq2kaps2/(2*kap_s*kap_l)

        # branch cut has changed phases, flip them back
        if np.max(np.abs(np.real(vely))) < np.max(np.abs(np.imag(vely))):
            vely *= 1.0j
            velz *= 1.0j

        m_uxyz[:,1] = vely
        m_uxyz[:,2] = velz

        max_u = np.max(np.abs(m_uxyz))
        if max_u > 0:
            m_uxyz /= np.max(np.abs(m_uxyz))

        # en_y = np.sum(np.abs(vely)**2)
        # en_z = np.sum(np.abs(velz)**2)
        # frac_trans = en_y/(en_y+en_z)

        Vp = Omega/q
        mode_prof = ModeFunction1D(Omega, Vp, v_x, m_uxyz)

        return mode_prof

