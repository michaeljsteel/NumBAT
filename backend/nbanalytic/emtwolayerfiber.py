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

import numpy as np
import scipy.optimize as sciopt

import matplotlib.pyplot as plt
import scipy.special as sp

from nbanalytic.emconstants import EMPoln
import reporting

def plot_em_chareq_at_k(k, rcore, ncore, nclad):
    nbrack = 500
    dn = ncore - nclad
    v_neff = np.linspace(ncore - dn / nbrack, nclad + dn / nbrack, nbrack)
    v_dr_TE = np.zeros(nbrack, float)
    v_dr_TM = np.zeros(nbrack, float)
    v_dr_hy = np.zeros(nbrack, float)

    fig, axs = plt.subplots(1, 3, figsize=(20, 6))
    for m in range(5):
        for i, neff in enumerate(v_neff):
            v_dr_TE[i] = chareq_em_fib2_TE_m(neff, m, k, rcore, ncore, nclad)
            v_dr_TM[i] = chareq_em_fib2_TM_m(neff, m, k, rcore, ncore, nclad)
            v_dr_hy[i] = chareq_em_fib2_hy_m(neff, m, k, rcore, ncore, nclad)
        axs[0].plot(v_neff, v_dr_TE)
        axs[1].plot(v_neff, v_dr_TM)
        axs[2].plot(v_neff, v_dr_hy)
    for i in range(3):
        axs[i].set_ylim(-20, 20)
    fig.savefig("tut_12-em_chareq.png")


def chareq_em_fib2_TE_m(neff, m, k, rco, ncore, nclad):
    # m is un-used (TE is m=0)
    u = rco * k * np.sqrt(ncore**2 - neff**2)
    w = rco * k * np.sqrt(neff**2 - nclad**2)

    return sp.jv(1, u) / (sp.jv(0, u) * u) + sp.kv(1, w) / (sp.kv(0, w) * w)


def chareq_em_fib2_TM_m(neff, m, k, rco, ncore, nclad):
    # m is un-used (TM is m=0)
    u = rco * k * np.sqrt(ncore**2 - neff**2)
    w = rco * k * np.sqrt(neff**2 - nclad**2)

    n2rat = (nclad / ncore) ** 2

    fac1 = np.real(sp.jv(1, u) / (sp.jv(0, u) * u))
    fac2 = np.real(sp.kv(1, w) / (sp.kv(0, w) * w))

    return fac1 + n2rat * fac2


def chareq_em_fib2_hy_m(neff, m, k, rco, ncore, nclad):
    u = rco * k * np.sqrt(ncore**2 - neff**2)
    w = rco * k * np.sqrt(neff**2 - nclad**2)

    n2rat = (nclad / ncore) ** 2

    invu2 = 1.0 / u**2
    invw2 = 1.0 / w**2
    jrat = sp.jvp(m, u) / (sp.jv(m, u) * u)
    krat = sp.kvp(m, w) / (sp.kv(m, w) * w)
    fac1 = jrat + krat
    fac2 = jrat + n2rat * krat
    fac3 = m * m * (invu2 + invw2) * (invu2 + n2rat * invw2)
    return fac1 * fac2 - fac3


# TODO Seems unused
# def del_solve_chareq_em_fib2_disprel_multim(poln, k, nmodes, mlo, mhi, rco, ncore, nclad):
#     # solves for modes with azimuthal order in [mlo, mhi] inclusive

#     nbrack = 500
#     dn = ncore-nclad

#     # look for up to nmodes in [nclad, ncore]
#     neff = np.linspace(ncore-dn/nbrack, nclad+dn/nbrack, nbrack)
#     sol_neff = np.zeros(nmodes, dtype=np.float64)
#     imode = 0
#     for m in range(mlo, mhi+1):

#         last_neff = neff[0]
#         last_res = f_disprel(last_neff, m, k, rco, ncore, nclad)

#         nobrack = True  # if we find no roots with a given m, there will be no higher m roots and we can give up
#         ineff = 1
#         while imode < nmodes and ineff < nbrack:
#             t_neff = neff[ineff]
#             t_res = f_disprel(t_neff, m, k, rco, ncore, nclad)

#             if ((family == EMPoln.HY and last_res * t_res < 0)  # Hybrid eig curves are smooth
#                     or (last_res < 0 and t_res > 0)):  # TE and TM curves have +inf to -inf breaks which are not brackets
#                 # a bracket! go find the root
#                 nobrack = False
#                 # root, rootres
#                 root = sciopt.brentq(f_disprel, last_neff, t_neff,
#                                      args=(m, k, rco, ncore, nclad))
#                 sol_neff[imode] = root
#                 imode += 1

#             last_neff = t_neff
#             last_res = t_res
#             ineff += 1

#         if nobrack:
#             break  # exhausted all roots at this k

#     return sol_neff


def _solve_chareq_em_fib2_disprel(
    poln, kvec, nmodes, azi_lo, azi_hi, rco, ncore, nclad
):
    # solves for modes with azimuthal order in [mlo, mhi] inclusive

    nbrack = 500
    dn = ncore - nclad

    # look for up to nmodes in [nclad, ncore]
    neff = np.linspace(ncore - dn / nbrack, nclad + dn / nbrack, nbrack)
    v_neff = np.zeros(nmodes, dtype=np.float64)
    imode = 0
    match poln:
        case EMPoln.TE:
            f_disprel = chareq_em_fib2_TE_m
        case EMPoln.TM:
            f_disprel = chareq_em_fib2_TM_m
        case EMPoln.HY:
            f_disprel = chareq_em_fib2_hy_m
        case _:
            reporting.report_and_exit(f"Unknown fiber polarisation type: {poln}")

    for m in range(azi_lo, azi_hi + 1):

        last_neff = neff[0]
        last_res = f_disprel(last_neff, m, kvec, rco, ncore, nclad)

        nobrack = True  # if we find no roots with a given m, there will be no higher m roots and we can give up
        ineff = 1
        while imode < nmodes and ineff < nbrack:
            t_neff = neff[ineff]
            t_res = f_disprel(t_neff, m, kvec, rco, ncore, nclad)

            if (
                poln == EMPoln.HY and last_res * t_res < 0
            ) or (  # Hybrid eig curves are smooth
                last_res < 0 and t_res > 0
            ):  # TE and TM curves have +inf to -inf breaks which are not brackets
                # a bracket! go find the root
                nobrack = False
                # root, rootres
                root = sciopt.brentq(
                    f_disprel, last_neff, t_neff, args=(m, kvec, rco, ncore, nclad)
                )
                v_neff[imode] = root
                imode += 1

            last_neff = t_neff
            last_res = t_res
            ineff += 1

        if nobrack:
            break  # exhausted all roots at this k

    return (imode, v_neff)


class TwoLayerFiberEM(object):
    """Exact analytic solutions of the two-layer step-index fibre.

    All incoming and outgoing parameters are in base SI units.
    """

    def __init__(self, ncore, nclad, arad):
        self._ncore = ncore
        self._nclad = nclad
        self._arad = arad

    def Vnum_k(self, k):
        return k * self._arad * math.sqrt(self._ncore**2 - self._nclad**2)

    def Vnumb_lam(self, lam):
        return self.Vnum_k(2 * math.pi / lam)

    def find_neffs_for_k(self, kvec, poln, azi_m_lo, azi_m_hi, nmax_modes):
        """
        Solves fiber exact eigenproblem at angular frequency omega in 1/s.

        Returns tuple (nfound, v_k) containing number of modes found and vector
        of eigen wavenumbers k_i[0:nfound-1].
        0 <= nfound <= nmax
        """
        (nmodes, v_neff) = _solve_chareq_em_fib2_disprel(
            poln,
            kvec,
            nmax_modes,
            azi_m_lo,
            azi_m_hi,
            self._arad,
            self._ncore,
            self._nclad,
        )
        return (nmodes, v_neff)

    def find_neff_HE11_for_k(self, kvec):
        (nmodes, v_neff) = self.find_neffs_for_k(kvec, EMPoln.HY, 1, 1, 1)

        if nmodes == 1:
            return v_neff[0]
        else:
            return None