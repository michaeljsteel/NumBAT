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
from numpy.typing import NDArray

import numpy as np
from functools import partial

import scipy.optimize as sciopt
import scipy.signal
import scipy.special as sp


twopi = 2 * math.pi

# Root finding tolerance constants
CANCELLATION_TOL = 1e-12  # Relative tolerance for catastrophic cancellation in determinants
DOUBLE_ROOT_ABS_TOL = 1e-8  # Absolute tolerance for double root detection
DOUBLE_ROOT_REL_TOL = 1e-9  # Relative tolerance (vs max) for double root detection
DOUBLE_ROOT_SMALLNESS_RATIO = 0.001  # Threshold for considering double root search


#######################################################################################


def get_q_transverse(Om: float, q: float, V: float) -> complex:
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
        return complex(np.sqrt(dsq), 0.0)

    if dsq < 0:
        return complex(0, np.sqrt(-dsq))

    return complex(0.0, 0.0)



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
    rad_core : float
        Cylinder radius in metres.
    """

    # Note that these disp rels continue to work fine when qts, qtl go imaginary
    # The imaginary terms always come in product pairs giving a real answer

    # p=0 torsional


    def chareq_elastic_cylinder_p0_1D(self, Om: float, p: int, q: float) -> float:
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

        Returns
        -------
        float
            A real-valued scalar made from the real and imaginary parts of the
            analytic characteristic expression; roots correspond to modal frequencies.
        """

        qts = get_q_transverse(Om, q, self._Vs)
        qtsa = qts * self._rad_core

        J0s = sp.jv(0, qtsa)
        J1s = sp.jv(1, qtsa)

        m22 = -(qts**2) * J0s + 2 * qts / self._rad_core * J1s

        chareq = m22 / Om**2

        return np.real(chareq) + np.imag(chareq)

        # p=0 Pochammer


    def chareq_elastic_cylinder_p0_2D(self, Om: float, p: int, q: float) -> float:
        """Characteristic function for p=0 (2D Pochhammer) modes.

        Parameters
        ----------
        Om : float
            Angular frequency (rad/s).
        p : int
            Azimuthal order (unused for p=0 case but kept for signature parity).
        q : float
            Axial wavenumber (1/m).

        Returns
        -------
        float
            Real scalar formed from characteristic determinant used for root finding.
        """

        qts = get_q_transverse(Om, q, self._Vs)
        qtl = get_q_transverse(Om, q, self._Vl)

        qtsa = qts * self._rad_core
        qtla = qtl * self._rad_core

        J0l = sp.jv(p, qtla)
        J0pl = sp.jvp(p, qtla)
        J0ppl = 0.5 * (sp.jvp(p - 1, qtla) - sp.jvp(p + 1, qtla))

        J1s = sp.jv(1, qtsa)
        J1ps = sp.jvp(1, qtsa)

        m00 = 2 * self._c44 * qtl**2 * J0ppl - self._c12 * (qtl**2 + q**2) * J0l
        m01 = 2j * self._c44 * q * qts * J1ps

        m10 = -2j * q * qtl * J0pl
        m11 = (-(qts**2) + q**2) * J1s

        d1 = m00 * m11
        d2 = m01 * m10

        disc = d1 - d2

        # Because of catastrophic cancellation, can't reasonably hope for closer to zero than this
        if np.abs(disc) < CANCELLATION_TOL * max(abs(d1), abs(d2)):
            chareq = 0.0
        else:
            chareq = disc / (Om**4)  # Om**4 makes numbers nice size

        # Branch cuts behave so that it Works nicely to just add real and imag parts.
        return chareq.real + chareq.imag


    def chareq_elastic_cylinder_ppos(self, Om: float, p: int, q: float) -> float:
        """Characteristic function for p>0 (general azimuthal order) modes.

        This computes the full 3x3 determinant for modes with positive azimuthal
        order p, including flexural and other complex mode families.

        Parameters
        ----------
        Om : float
            Angular frequency (rad/s).
        p : int
            Azimuthal order (must be positive).
        q : float
            Axial wavenumber (1/m).

        Returns
        -------
        float
            Real scalar formed from the 3x3 characteristic determinant used for root finding.
        """

        a = self._rad_core

        qts = get_q_transverse(Om, q, self._Vs)
        qtl = get_q_transverse(Om, q, self._Vl)

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

        m00 = 2 * self._c44 * qtl**2 * J0ppl - self._c12 * (qtl**2 + q**2) * J0l
        m01 = 1j * 2 * self._c44 * q * qts * J1ps
        m02 = 2 * self._c44 * p / a * (J0ps * qts - J0s / a)

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

        # Because of catastrophic cancellation, can't hope for closer to zero than this
        if np.abs(disc) < CANCELLATION_TOL * bigd:
            chareq = 0.0
        else:
            chareq = disc / (Om**5)  # Om**5 makes numbers nice size

        # Branch cuts behave so that it Works nicely to just add real and imag parts.
        return chareq.real + chareq.imag


    def chareq_elastic_cylinder(self, Om: float, m: int, q: float) -> float:
        """Dispatcher for characteristic equation based on azimuthal order.

        Routes to the appropriate characteristic function based on the azimuthal
        order m: torsional (m=-1), axisymmetric Pochhammer (m=0), or general (m>0).

        Parameters
        ----------
        Om : float
            Angular frequency (rad/s).
        m : int
            Azimuthal order. Use -1 for torsional modes, 0 for axisymmetric, >0 for flexural.
        q : float
            Axial wavenumber (1/m).
        rho, c11, c12, c44 : float
            Material density and stiffness constants in SI units.
        rad_core : float
            Cylinder radius in metres.

        Returns
        -------
        float
            Value of the characteristic equation at the given parameters.
        """
        if m == -1:  # actually denotes a p=0 1D state
            return self.chareq_elastic_cylinder_p0_1D(Om, 0, q)

        elif m == 0:
            return self.chareq_elastic_cylinder_p0_2D(Om, 0, q)

        else:
            return self.chareq_elastic_cylinder_ppos(Om, m, q)


    def _findroots_elastic_cylinder_chareq(self, q: float, nu_hi: float, m: int, nmodes: int) -> Tuple[int, NDArray[np.float64]]:
        """Find roots of the characteristic equation (mode frequencies) at given q.

        Scans the frequency range [Omlo, Omhi] for up to `nmodes` roots using
        bracket-and-refine with Brent's method. Also detects potential double roots
        via local extrema analysis.

        Parameters
        ----------
        q : float
            Axial wavenumber (1/m).
        nu_hi : float
            Upper search frequency in Hz.
        m : int
            Azimuthal order.
        nmodes : int
            Maximum number of modes to find.
        rho, c11, c12, c44 : float
            Material density and stiffness constants in SI units.
        rad_core : float
            Cylinder radius in metres.

        Returns
        -------
        Tuple[int, NDArray[np.float64]]
            (number_of_modes_found, array_of_angular_frequencies)
            The array has length `nmodes`; only the first `number_of_modes_found`
            entries are valid roots.
        """

        num_Om = 500
        Omlo = 0.01e6
        Omhi = twopi * nu_hi

        # look for up to nmodes in [Omlo, Omhi]
        Om = np.linspace(Omlo, Omhi, num_Om)
        drval = np.zeros(len(Om))

        v_Om = np.zeros(nmodes, dtype=np.float64)
        imode = 0

        chareq = partial(self.chareq_elastic_cylinder, m=m, q=q)
        t_Om = Om[0]
        t_res = chareq(t_Om)
        drval[0] = t_res

        for iOm in range(1, num_Om):
            last_Om = t_Om
            last_res = t_res

            t_Om = Om[iOm]
            t_res = chareq(t_Om)
            drval[iOm] = t_res

            if imode < nmodes and last_res * t_res < 0:  # a bracket! go find the root
                root, root_res = sciopt.brentq(
                    chareq,
                    last_Om,
                    t_Om,
                    full_output=True,
                    disp=False,
                )

                if not root_res.converged:

                    # Must be at a nasty point. Try bisection which can't really fail.
                    root, root_res = sciopt.bisect(
                        chareq,
                        last_Om,
                        t_Om,
                        full_output=True,
                        disp=False,
                    )
                if root_res.converged:
                    v_Om[imode] = root
                    imode += 1
                else:
                    print(
                        f"Brent and bisection unable to converge for p={m}, "
                        f"q={q*1e-6:.4f}/um, "
                        f"nu in [{last_Om/(2e9*math.pi):.6f}, {t_Om/(2e9*math.pi):.6f}] GHz, "
                        f"err in [{last_res:.4e}, {t_res:.4e}]. "
                    )

        # This dispersion relation tends to exhibit double roots that may not trigger a bracketing.
        # Instead we can look for any points which are local min/max and very close to zero
        # Picking the right tolerances here is rather sensitive
        bigdr = np.max(np.abs(drval))
        smalldr = np.min(np.abs(drval))
        if smalldr < DOUBLE_ROOT_SMALLNESS_RATIO * bigdr:
            # There is a significant degree of smallness worth checking
            # Find all local maxima
            quasi_roots = scipy.signal.find_peaks(-np.abs(drval))[0]

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
                if np.abs(droot) < DOUBLE_ROOT_ABS_TOL and np.abs(droot) < DOUBLE_ROOT_REL_TOL * bigdr:
                    # found a likely double root, fit a parabola to get accurate location
                    a = (drval[qrt + 1] + drval[qrt - 1] - 2 * drval[qrt]) / 2
                    b = (drval[qrt + 1] - drval[qrt - 1]) / 2

                    # Guard against near-zero curvature (would cause divide-by-zero)
                    if np.abs(a) > 1e-15:
                        v_Om[imode] = Om[qrt] - b / (2 * a) * (Om[1] - Om[0])
                        print(
                            f"Found an additional double root for p={m} at nu={v_Om[imode]/(1e9*twopi):.4e}."
                        )
                    imode += 1

        return (imode, v_Om)



    def __init__(self, rho: float, c11: float, c12: float, c44: float, rad_core: float) -> None:
        self._rho = rho
        self._c11 = c11
        self._c12 = c12
        self._c44 = c44
        self._rad_core = rad_core
        self._Vl = np.sqrt(c11 / rho)
        self._Vs = np.sqrt(c44 / rho)

    def find_Omega_at_q(self, q: float, Om_hi: float,
                        m: int, nmax_modes: int) -> Tuple[int, NDArray[np.float64]]:
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
        Tuple[int, NDArray[np.float64]]
            (nmodes_found, omegas) where omegas is an array of angular frequencies.
        """

        (nmodes, v_Om) = self._findroots_elastic_cylinder_chareq(
            q,
            Om_hi/twopi,
            m,
            nmax_modes,
        )
        return (nmodes, v_Om)

    def dispersion_relation_at_q_nu(self, q: float, nu: float, m: int) -> float:
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

        return self.chareq_elastic_cylinder(
            np.pi * 2 * nu, m, q
        )
