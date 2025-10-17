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


"""elasticcylinder.py
----------------------
Analytic routines for elastic wave dispersion in cylindrical geometry.

This module contains closed-form characteristic equations for an isotropic
elastic cylinder in vacuum, helper functions to compute wavenumbers, and a
small solver to locate roots of the characteristic equations (mode frequencies).

All inputs and outputs are in SI units unless explicitly noted in a function
docstring. The functions are primarily used by higher-level analytical and
testing code in the NumBAT project.

References
----------
Many expressions are based on standard textbooks on elastodynamics (e.g. Auld)
and on analytic results for Pochhammer-Chree/Rayleigh-Lamb type modes in
cylindrical structures.
"""

import math
from typing import Tuple

import numpy as np

import scipy.optimize as sciopt
import scipy.signal
import scipy.special as sp

twopi = 2 * math.pi


#######################################################################################


def getqt(Om: float, q: float, V: float) -> complex:
    """Compute transverse/longitudinal radial wavenumber component.

    Parameters
    ----------
    Om : float
        Angular frequency (rad/s).
    q : float
        Axial wavenumber (1/m).
    V : float
        Bulk wave speed (m/s).

    Returns
    -------
    complex
        Radial wavenumber sqrt((Om/V)**2 - q**2). If the argument is negative
        the function returns a purely imaginary value (decaying radial behaviour).
    """
    dsq = (Om / V) ** 2 - q**2
    if dsq > 0:
        return np.sqrt(dsq)

    if dsq < 0:
        return complex(0, np.sqrt(-dsq))

    return 0


# Note that these disp rels continue to work fine when qts, qtl go imaginary
# The imaginary terms always come in product pairs giving a real answer

# p=0 torsional


def chareq_elastic_cylinder_p0_1D(Om: float, p: int, q: float, rho: float, c11: float, c12: float, c44: float, a_nm: float) -> float:
    """Characteristic function for p=0 torsional (1D) modes.

    This returns a real scalar constructed from the complex characteristic
    expression; callers typically search for zeros of this function to find
    resonant frequencies.

    Parameters
    ----------
    Om : float
        Angular frequency (rad/s).
    p : int
        Azimuthal order (unused for p=0 case but kept for signature parity).
    q : float
        Axial wavenumber (1/m).
    rho, c11, c12, c44 : float
        Material density and stiffness constants in SI units.
    a_nm : float
        Cylinder radius in nanometres (function converts to metres internally).

    Returns
    -------
    float
        A real-valued scalar made from the real and imaginary parts of the
        analytic characteristic expression; roots correspond to modal frequencies.
    """

    Vs = np.sqrt(c44 / rho)
    qts = getqt(Om, q, Vs)
    a = a_nm * 1e-9
    qtsa = qts * a

    J0s = sp.jv(0, qtsa)
    J1s = sp.jv(1, qtsa)

    m22 = -(qts**2) * J0s + 2 * qts / a * J1s

    chareq = m22 / Om**2

    return np.real(chareq) + np.imag(chareq)

    # p=0 Pochammer


def chareq_elastic_cylinder_p0_2D(Om: float, p: int, q: float, rho: float, c11: float, c12: float, c44: float, a_nm: float) -> float:
    """Characteristic function for p=0 (2D Pochhammer) modes.

    Parameters
    ----------
    Om, p, q, rho, c11, c12, c44, a_nm
        See :func:`chareq_elastic_cylinder_p0_1D` for parameter meanings.

    Returns
    -------
    float
        Real scalar formed from characteristic determinant used for root finding.
    """
    Vl = np.sqrt(c11 / rho)
    Vs = np.sqrt(c44 / rho)

    qts = getqt(Om, q, Vs)
    qtl = getqt(Om, q, Vl)

    a = a_nm * 1e-9

    qtsa = qts * a
    qtla = qtl * a

    J0l = sp.jv(p, qtla)
    J0pl = sp.jvp(p, qtla)
    J0ppl = 0.5 * (sp.jvp(p - 1, qtla) - sp.jvp(p + 1, qtla))

    J1s = sp.jv(1, qtsa)
    J1ps = sp.jvp(1, qtsa)

    m00 = 2 * c44 * qtl**2 * J0ppl - c12 * (qtl**2 + q**2) * J0l
    m01 = 2j * c44 * q * qts * J1ps

    m10 = -2j * q * qtl * J0pl
    m11 = (-(qts**2) + q**2) * J1s

    d1 = m00 * m11
    d2 = m01 * m10

    disc = d1 - d2

    # Because of catastrophic cancellation, can't reasonably hope for closer to zero than this
    tol = 1e-12
    if np.abs(disc) < tol * max(abs(d1), abs(d2)):
        chareq = 0
    else:
        chareq = disc / (Om**4)  # Om**4 makes numbers nice sie

    # Branch cuts behave so that it Works nicely to just add real and imag parts.
    return chareq.real + chareq.imag


def chareq_elastic_cylinder_ppos(Om: float, p: int, q: float, rho: float, c11: float, c12: float, c44: float, a_nm: float) -> float:  # p is azimuthal order

    Vl = np.sqrt(c11 / rho)  # = sqrt((lambda + 2 mu)/rho)
    Vs = np.sqrt(c44 / rho)  # = sqrt(mu/rho)

    qts = getqt(Om, q, Vs)
    qtl = getqt(Om, q, Vl)

    a = a_nm * 1e-9

    qtsa = qts * a
    qtla = qtl * a

    J0l = sp.jv(p, qtla)
    J0pl = sp.jvp(p, qtla)
    J0ppl = 0.5 * (sp.jvp(p - 1, qtla) - sp.jvp(p + 1, qtla))

    J0s = sp.jv(p, qtsa)
    J0ps = sp.jvp(p, qtsa)
    J0pps = 0.5 * (sp.jvp(p - 1, qtsa) - sp.jvp(p + 1, qtsa))

    J1s = sp.jv(p + 1, qtsa)
    J1ps = sp.jvp(p + 1, qtsa)

    m00 = 2 * c44 * qtl**2 * J0ppl - c12 * (qtl**2 + q**2) * J0l
    m01 = 1j * 2 * c44 * q * qts * J1ps
    m02 = 2 * c44 * p / a * (J0ps * qts - J0s / a)

    m10 = 1j * 2 * q * qtl * J0pl
    m11 = J1s * ((qts**2 - q**2) - (p + p * p) / (a * a)) - J1ps * p * qts / a
    m12 = 1j * q * p * J0s / a

    m20 = 2 * p / a * (qtl * J0pl - J0l / a)
    m21 = 1j * J1s * (q / a + p * q / a) - 1j * J1ps * q * qts
    m22 = qts**2 * (2 * J0pps + J0s)

    d0a = m00 * m11 * m22
    d0b = -m00 * m12 * m21
    d1a = m01 * m12 * m20
    d1b = -m01 * m10 * m22
    d2a = m02 * m10 * m21
    d2b = -m02 * m11 * m20

    bigd = np.max(np.abs(np.array([d0a, d0b, d1a, d1b, d2a, d2b])))
    disc = d0a + d0b + d1a + d1b + d2a + d2b

    tol = 1e-12  # Because of catastrophic cancellation, can't hope for closer to zero than this
    if np.abs(disc) < tol * bigd:
        chareq = 0
    else:
        chareq = disc / (Om**5)  # Om**5 makes numbers nice sie

    # Branch cuts behave so that it Works nicely to just add real and imag parts.
    return chareq.real + chareq.imag


def chareq_elastic_cylinder(Om: float, m: int, q: float, rho: float, c11: float, c12: float, c44: float, arad: float) -> float:  # m is azimuthal order
    a_nm = arad * 1e9
    if m == -1:  # actually denotes a p=0 1D state
        return chareq_elastic_cylinder_p0_1D(Om, 0, q, rho, c11, c12, c44, a_nm)

    if m == 0:
        return chareq_elastic_cylinder_p0_2D(Om, m, q, rho, c11, c12, c44, a_nm)

    return chareq_elastic_cylinder_ppos(Om, m, q, rho, c11, c12, c44, a_nm)


def _findroots_elastic_cylinder_chareq(q: float, nu_hi: float, m: int, nmodes: int, Vl: float, rho: float, c11: float, c12: float, c44: float, arad: float) -> Tuple[int, np.ndarray]:

    num_Om = 500
    Omlo = 0.01e6
    Omhi = twopi * nu_hi

    # look for up to nmodes in [Omlo, Omhi]
    Om = np.linspace(Omlo, Omhi, num_Om)
    drval = np.zeros(len(Om))

    v_Om = np.zeros(nmodes, dtype=np.float64)
    imode = 0

    t_Om = Om[0]
    t_res = chareq_elastic_cylinder(t_Om, m, q, rho, c11, c12, c44, arad)
    drval[0] = t_res

    for iOm in range(1, num_Om):
        last_Om = t_Om
        last_res = t_res

        t_Om = Om[iOm]
        t_res = chareq_elastic_cylinder(t_Om, m, q, rho, c11, c12, c44, arad)
        drval[iOm] = t_res

        if imode < nmodes and last_res * t_res < 0:  # a bracket! go find the root
            root, root_res = sciopt.brentq(
                chareq_elastic_cylinder,
                last_Om,
                t_Om,
                args=(m, q, rho, c11, c12, c44, arad),
                full_output=True,
                disp=False,
            )

            if not root_res.converged:

                # Must be at a nasty point. Try bisection which can't really fail.
                root, root_res = sciopt.bisect(
                    chareq_elastic_cylinder,
                    last_Om,
                    t_Om,
                    args=(m, q, rho, c11, c12, c44, arad),
                    full_output=True,
                    disp=False,
                )
            if root_res.converged:
                v_Om[imode] = root
                imode += 1
            else:
                print(
                    "Brent and bisection unable to converge for p=%d, q=%.4f/um, nu in [%.6f, %.6f] GHz, err in [%.4e, %.4e]. "
                    % (
                        m,
                        q * 1e-6,
                        last_Om / (2e9 * math.pi),
                        t_Om / (2e9 * math.pi),
                        last_res,
                        t_res,
                    )
                )

    # This dispersion relation tends to exhibit double roots that may not trigger a bracketing.
    # Instead we can look for any points which are local min/max and very close to zero
    # Picking the right tolerances here is rather sensitive
    bigdr = np.max(np.abs(drval))
    smalldr = np.min(np.abs(drval))
    if (
        smalldr < 0.001 * bigdr
    ):  # There is a significant degree of smallness worth checking
        # Find all local maxima
        quasi_roots = scipy.signal.find_peaks(-np.abs(drval))[0]
        # print(f'potential extra roots for p={p}, q={q*1e-6:.4f}\n',
        #      quasi_roots, '\n  ', Om[quasi_roots]/(1e9*twopi), '\n   ',
        #      drval[quasi_roots]/(1e9*twopi))
        for qrt in quasi_roots:
            # Found all required modes
            if imode == nmodes:
                break

            # Don't bother with end points
            if qrt in [0, len(Om) - 1]:
                continue

            # Bracketed root, so will have been picked up above by brentq
            if (drval[qrt - 1] * drval[qrt] < 0) or (drval[qrt] * drval[qrt + 1] < 0):
                continue

            droot = drval[qrt]
            if np.abs(droot) < 1e-8 and np.abs(droot) < 1e-9 * bigdr:
                # found a likely double root, fit a parabola to get accurate location
                a = (drval[qrt + 1] + drval[qrt - 1] - 2 * drval[qrt]) / 2
                b = (drval[qrt + 1] - drval[qrt - 1]) / 2
                v_Om[imode] = Om[qrt] - b / (2 * a) * (Om[1] - Om[0])
                print(
                    f"Found an additional double root for p={m} at nu={v_Om[imode]/(1e9*twopi):.4e}."
                )
                imode += 1

    return (imode, v_Om)


class ElasticIsotropicCylinder(object):
    """Exact analytic solutions for an isotropic elastic cylinder in vacuum.

    This convenience class wraps the characteristic-equation routines and
    provides search helpers to find modal frequencies at given axial wavenumbers.

    Parameters
    ----------
    rho : float
        Material mass density (kg/m^3).
    c11, c12, c44 : float
        Elastic stiffness constants (SI units, Pa or derived units consistent with rho).
    arad : float
        Cylinder radius in metres.
    """

    def __init__(self, rho, c11, c12, c44, arad):
        self._rho = rho
        self._c11 = c11
        self._c12 = c12
        self._c44 = c44
        self._arad = arad
        self._Vl = np.sqrt(c11 / rho)
        self._Vs = np.sqrt(c44 / rho)

    def find_Omega_at_q(self, q, Om_hi, m, nmax_modes):
        """Find discrete angular frequencies for a given axial wavenumber.

        Parameters
        ----------
        q : float
            Axial wavenumber (1/m).
        Om_hi : float
            Upper search frequency (angular frequency, rad/s).
        m : int
            Azimuthal order.
        nmax_modes : int
            Maximum number of modes to find.

        Returns
        -------
        tuple
            (nmodes_found, omegas) where omegas is an array of angular frequencies.
        """

        (nmodes, v_Om) = _findroots_elastic_cylinder_chareq(
            q,
            Om_hi/twopi,
            m,
            nmax_modes,
            self._Vl,
            self._rho,
            self._c11,
            self._c12,
            self._c44,
            self._arad,
        )
        return (nmodes, v_Om)

    def dispersion_relation_at_q_nu(self, q, nu, m):
        """Evaluate the characteristic equation at given axial wavenumber and frequency.

        Parameters
        ----------
        q : float
            Axial wavenumber (1/m).
        nu : float
            Frequency in Hz.
        m : int
            Azimuthal order.

        Returns
        -------
        float
            The (real) value of the characteristic function at the supplied inputs.
        """

        return chareq_elastic_cylinder(
            np.pi * 2 * nu, m, q, self._rho, self._c11, self._c12, self._c44, self._arad
        )  # m is azimuthal order
