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
from typing import Optional, Tuple

import numpy as np
import scipy.optimize as sciopt

import matplotlib.pyplot as plt
import scipy.special as sp
from numpy.typing import NDArray

from nbanalytic.emconstants import EMPoln
import reporting

class TwoLayerFiberEM(object):
    """Analytic modes of a two-layer (step-index) circular optical fibre.

    This class provides characteristic equations and simple solvers for
    TE/TM (m=0) and hybrid (HE/EH, m>=1) modes of a cylindrical step-index
    fibre using Bessel/modified-Bessel functions. It also includes helpers to
    compute the fibre V-number.

    All inputs/outputs are SI units; refractive indices are dimensionless.
    """

    def __init__(self, n_core: float, n_clad: float, arad: float) -> None:
        """Construct a two-layer step-index fibre model.

        Parameters
        ----------
        n_core : float
            Core refractive index (dimensionless).
        n_clad : float
            Cladding refractive index (dimensionless). Typically n_core > n_clad.
        arad : float
            Core radius in metres (> 0).

        Notes
        -----
        No validation is enforced here; callers should ensure physical values.
        """
        self._n_core: float = n_core
        self._n_clad: float = n_clad
        self._rad_core: float = arad

        self._nco_sq: float = self._n_core**2
        self._ncl_sq: float = self._n_clad**2
        self._nsq_ratio: float = self._ncl_sq / self._nco_sq

    def plot_em_chareq_at_k(self, kvac: float, rad_core: float, n_core: float, n_clad: float) -> None:
        """Plot characteristic equations versus effective index for a given k.

        This debugging/visualisation aid samples the characteristic equations
        for TE, TM, and hybrid families across a range of neff between the
        cladding and core indices at a fixed wavenumber k.

        Parameters
        ----------
        k : float
            Free-space wavenumber (rad/m), k = 2π/λ.
        rad_core : float
            Core radius (m). Used only for plotting; does not change state.
        n_core : float
            Core refractive index (dimensionless).
        n_clad : float
            Cladding refractive index (dimensionless).
        """
        nbrack = 500
        dn = n_core - n_clad
        v_neff = np.linspace(n_core - dn / nbrack, n_clad + dn / nbrack, nbrack)
        v_dr_TE = np.zeros(nbrack, float)
        v_dr_TM = np.zeros(nbrack, float)
        v_dr_hy = np.zeros(nbrack, float)

        fig, axs = plt.subplots(1, 3, figsize=(20, 6))
        for m in range(5):
            for i, neff in enumerate(v_neff):
                v_dr_TE[i] = self.chareq_em_fib2_TE_m(neff, m, kvac)
                v_dr_TM[i] = self.chareq_em_fib2_TM_m(neff, m, kvac)
                v_dr_hy[i] = self.chareq_em_fib2_hy_m(neff, m, kvac)
            axs[0].plot(v_neff, v_dr_TE)
            axs[1].plot(v_neff, v_dr_TM)
            axs[2].plot(v_neff, v_dr_hy)
        for i in range(3):
            axs[i].set_ylim(-20, 20)
        fig.savefig("tut_12-em_chareq.png")


    def chareq_em_fib2_TE_m(
        self,
        neff: float,
        m: int,
        kvac: float,
    ) -> float:
        """Characteristic equation residual for TE (m=0) modes.

        Parameters
        ----------
        neff : float
            Trial effective index (n_clad < neff < n_core).
        m : int
            Azimuthal order. Unused; TE is strictly m = 0.
        k : float
            Free-space wavenumber (rad/m), k = 2π/λ.

        Returns
        -------
        float
            Residual of the TE characteristic equation. Zeros indicate modes.
        """
        # m is un-used (TE is m=0)
        u = self._rad_core * kvac * np.sqrt(self._nco_sq - neff**2)
        w = self._rad_core * kvac * np.sqrt(neff**2 - self._ncl_sq)

        return sp.jv(1, u) / (sp.jv(0, u) * u) + sp.kv(1, w) / (sp.kv(0, w) * w)


    def chareq_em_fib2_TM_m(
        self,
        neff: float,
        m: int,
        kvac: float,
    ) -> float:
        """Characteristic equation residual for TM (m=0) modes.

        Parameters
        ----------
        neff : float
            Trial effective index (n_clad < neff < n_core).
        m : int
            Azimuthal order. Unused; TM is strictly m = 0.
        kvac : float
            Free-space wavenumber (rad/m), k = 2π/λ.

        Returns
        -------
        float
            Residual of the TM characteristic equation. Zeros indicate modes.
        """
        # m is un-used (TM is m=0)
        u = self._rad_core * kvac * np.sqrt(self._nco_sq - neff**2)
        w = self._rad_core * kvac * np.sqrt(neff**2 - self._ncl_sq)

        fac1 = np.real(sp.jv(1, u) / (sp.jv(0, u) * u))
        fac2 = np.real(sp.kv(1, w) / (sp.kv(0, w) * w))

        return fac1 + self._nsq_ratio * fac2


    def chareq_em_fib2_hy_m(
        self,
        neff: float,
        m: int,
        kvac: float,
    ) -> float:
        """Characteristic equation residual for hybrid (HE/EH) modes.

        Parameters
        ----------
        neff : float
            Trial effective index (n_clad < neff < n_core).
        m : int
            Azimuthal order (m ≥ 1) for hybrid modes.
        kvac : float
            Free-space wavenumber (rad/m), k = 2π/λ.

        Returns
        -------
        float
            Residual of the hybrid characteristic equation. Zeros are modes.
        """
        u = self._rad_core * kvac * np.sqrt(self._nco_sq - neff**2)
        w = self._rad_core * kvac * np.sqrt(neff**2 - self._ncl_sq)

        n2rat = self._nsq_ratio

        invu2 = 1.0 / u**2
        invw2 = 1.0 / w**2
        jrat = sp.jvp(m, u) / (sp.jv(m, u) * u)
        krat = sp.kvp(m, w) / (sp.kv(m, w) * w)
        fac1 = jrat + krat
        fac2 = jrat + n2rat * krat
        fac3 = m * m * (invu2 + invw2) * (invu2 + n2rat * invw2)
        return fac1 * fac2 - fac3


# TODO Seems unused
# def del_solve_chareq_em_fib2_disprel_multim(poln, k, nmodes, mlo, mhi, r_core, n_core, n_clad):
#     # solves for modes with azimuthal order in [mlo, mhi] inclusive

#     nbrack = 500
#     dn = n_core-n_clad

#     # look for up to nmodes in [n_clad, n_core]
#     neff = np.linspace(n_core-dn/nbrack, n_clad+dn/nbrack, nbrack)
#     sol_neff = np.zeros(nmodes, dtype=np.float64)
#     imode = 0
#     for m in range(mlo, mhi+1):

#         last_neff = neff[0]
#         last_res = f_disprel(last_neff, m, k, r_core, n_core, n_clad)

#         nobrack = True  # if we find no roots with a given m, there will be no higher m roots and we can give up
#         ineff = 1
#         while imode < nmodes and ineff < nbrack:
#             t_neff = neff[ineff]
#             t_res = f_disprel(t_neff, m, k, r_core, n_core, n_clad)

#             if ((family == EMPoln.HY and last_res * t_res < 0)  # Hybrid eig curves are smooth
#                     or (last_res < 0 and t_res > 0)):  # TE and TM curves have +inf to -inf breaks which are not brackets
#                 # a bracket! go find the root
#                 nobrack = False
#                 # root, rootres
#                 root = sciopt.brentq(f_disprel, last_neff, t_neff,
#                                      args=(m, k, r_core, n_core, n_clad))
#                 sol_neff[imode] = root
#                 imode += 1

#             last_neff = t_neff
#             last_res = t_res
#             ineff += 1

#         if nobrack:
#             break  # exhausted all roots at this k

#     return sol_neff



    def _solve_chareq_em_fib2_disprel(
        self,
        poln: EMPoln,
        kvac: float,
        nmodes: int,
        azi_lo: int,
        azi_hi: int,
    ) -> Tuple[int, NDArray[np.float64]]:
        """Internal: scan for neff roots of the characteristic equation.

        Brackets roots over a uniform grid of trial neff between n_clad and
        n_core for each azimuthal order in [azi_lo, azi_hi], then refines each
        bracket with Brent's method.

        Parameters
        ----------
        poln : EMPoln
            Mode family: TE, TM, or HY (hybrid).
        kvac : float
            Free-space wavenumber (rad/m), k = 2π/λ.
        nmodes : int
            Maximum number of modes to return.
        azi_lo : int
            Lowest azimuthal order to search (inclusive).
        azi_hi : int
            Highest azimuthal order to search (inclusive).

        Returns
        -------
        Tuple[int, NDArray[np.float64]]
            A pair (nfound, neffs) where nfound is the number of roots found
            and neffs contains the effective indices in ascending order in the
            first nfound entries.
        """
        # solves for modes with azimuthal order in [mlo, mhi] inclusive

        nbrack = 500
        dn = self._n_core - self._n_clad

        # look for up to nmodes in [n_clad, n_core]
        neff = np.linspace(self._n_core - dn / nbrack, self._n_clad + dn / nbrack, nbrack)
        v_neff: NDArray[np.float64] = np.zeros(nmodes, dtype=np.float64)
        imode = 0
        match poln:
            case EMPoln.TE:
                f_disprel = self.chareq_em_fib2_TE_m
            case EMPoln.TM:
                f_disprel = self.chareq_em_fib2_TM_m
            case EMPoln.HY:
                f_disprel = self.chareq_em_fib2_hy_m
            case _:
                reporting.report_and_exit(f"Unknown fiber polarisation type: {poln}")

        for m in range(azi_lo, azi_hi + 1):

            last_neff = neff[0]
            last_res = f_disprel(last_neff, m, kvac)

            nobrack = True  # if we find no roots with a given m, there will be no higher m roots and we can give up
            ineff = 1
            while imode < nmodes and ineff < nbrack:
                t_neff = neff[ineff]
                t_res = f_disprel(t_neff, m, kvac)

                if (
                    poln == EMPoln.HY and last_res * t_res < 0
                ) or (  # Hybrid eig curves are smooth
                    last_res < 0 and t_res > 0
                ):  # TE and TM curves have +inf to -inf breaks which are not brackets
                    # a bracket! go find the root
                    nobrack = False
                    # root, rootres
                    root = sciopt.brentq(
                        f_disprel, last_neff, t_neff, args=(m, kvac)
                    )
                    v_neff[imode] = root
                    imode += 1

                last_neff = t_neff
                last_res = t_res
                ineff += 1

            if nobrack:
                break  # exhausted all roots at this k

        return (imode, v_neff)



    def Vnum_k(self, k: float) -> float:
        """Return the fibre V-number at wavenumber k.

        V = k a sqrt(n_core^2 - n_clad^2)

        Parameters
        ----------
        k : float
            Free-space wavenumber (rad/m), k = 2π/λ.

        Returns
        -------
        float
            The dimensionless V-number.
        """
        return k * self._rad_core * math.sqrt(self._nco_sq - self._ncl_sq)

    def Vnumb_lam(self, lam: float) -> float:
        """Return the fibre V-number at wavelength lam.

        Parameters
        ----------
        lam : float
            Free-space wavelength (m).

        Returns
        -------
        float
            The dimensionless V-number.
        """
        return self.Vnum_k(2 * math.pi / lam)

    def find_neffs_for_k(
        self,
        kvac: float,
        poln: EMPoln,
        azi_m_lo: int,
        azi_m_hi: int,
        nmax_modes: int,
    ) -> Tuple[int, NDArray[np.float64]]:
        """Find effective indices for modes at wavenumber k.

        Parameters
        ----------
        kvac : float
            Free-space wavenumber (rad/m), k = 2π/λ.
        poln : EMPoln
            Mode family: TE, TM, or HY (hybrid).
        azi_m_lo : int
            Lowest azimuthal order to search (inclusive).
        azi_m_hi : int
            Highest azimuthal order to search (inclusive).
        nmax_modes : int
            Maximum number of modes to return.

        Returns
        -------
        Tuple[int, NDArray[np.float64]]
            (nfound, neffs) where neffs holds up to nfound effective indices.
        """
        (nmodes, v_neff) = self._solve_chareq_em_fib2_disprel(
            poln,
            kvac,
            nmax_modes,
            azi_m_lo,
            azi_m_hi,
            #self._rad_core,
            #self._n_core,
            #self._n_clad,
        )
        return (nmodes, v_neff)

    def find_neff_HE11_for_k(self, kvac: float) -> Optional[float]:
        """Return the HE11 effective index at wavenumber k, if found.

        Parameters
        ----------
        kvac : float
            Free-space wavenumber (rad/m), k = 2π/λ.

        Returns
        -------
        Optional[float]
            neff of the fundamental HE11 mode if found; otherwise None.
        """
        (nmodes, v_neff) = self.find_neffs_for_k(kvac, EMPoln.HY, 1, 1, 1)

        if nmodes == 1:
            return v_neff[0]
        else:
            return None