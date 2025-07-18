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

from math import sqrt, atan, floor
import numpy as np
import scipy.optimize as sciopt


from nbanalytic.emconstants import EMPoln


def _emslab_chareq_TE(V, gamma, m, b):
    """Dispersion relation for TE slab optical waveguide in normalised units."""
    return (
        V * sqrt(1 - b)
        - m * np.pi
        - atan(sqrt(b / (1 - b)))
        - atan(sqrt((b + gamma) / (1 - b)))
    )


def _emslab_chareq_TM(V, gamma, ns, nf, nc, m, b):
    """Dispersion relation for TE slab optical waveguide in normalised units."""
    tmfac1 = (nf / ns) ** 2
    tmfac2 = (nf / nc) ** 2
    return (
        V * sqrt(1 - b)
        - m * np.pi
        - atan(tmfac1 * sqrt(b / (1 - b)))
        - atan(tmfac2 * sqrt((b + gamma) / (1 - b))))


class EMSlab:
    """EM slab waveguide  solver.

    All units are in SI_ums.
    """

    def __init__(self, ns, nf, nc, wid):
        # Allows incoming complex value refractive indices but makes them real

        self._ns = ns.real  # substrate index
        self._nf = nf.real  # film index
        self._nc = nc.real  # cover index
        self.width = wid.real  # width

        if not self._nc <= self._ns < self._nf:
            raise ValueError("Indices do not specify a guiding structure.")
        self._gamma = (self._ns**2 - self._nc**2) / (self._nf**2 - self._ns**2)

        self._tmfac_s = (self._nf / self._ns) ** 2
        self._tmfac_c = (self._nf / self._nc) ** 2

    def lambda_to_V(self, lam):
        return 2 * np.pi / lam * self.width * sqrt(self._nf**2 - self._ns**2)

    def V_to_lambda(self, V):
        return 2 * np.pi / V * self.width * sqrt(self._nf**2 - self._ns**2)

    def b_to_neff(self, b):
        return sqrt(self._ns**2 + (self._nf**2 - self._ns**2) * b)

    def cutoff_V(self, m, poln):
        match poln:
            case EMPoln.TE:
                return m * np.pi + atan(sqrt(self._gamma))

            case EMPoln.TM:
                return m * np.pi + atan(self._tmfac_c * sqrt(self._gamma))

    def num_modes(self, V, poln):
        match poln:

            case EMPoln.TE:
                Vcorr = V - atan(sqrt(self._gamma))
                return floor(Vcorr / np.pi) + 1

            case EMPoln.TM:
                tmfac = (self._nf / self._nc) ** 2
                Vcorr = V - atan(tmfac * sqrt(self._gamma))
                return floor(Vcorr / np.pi) + 1

    def find_b_for_V(self, V, m, poln):
        """Assumes that cutoffs have been checked and there is a solution."""

        match poln:

            case EMPoln.TE:
                bfunc = lambda b: _emslab_chareq_TE(V, self._gamma, m, b)

            case EMPoln.TM:
                bfunc = lambda b: _emslab_chareq_TM(
                    V, self._gamma, self._ns, self._nf, self._nc, m, b
                )

        tol = 1e-10
        b_sol = sciopt.root_scalar(bfunc, bracket=[tol, 1 - tol])

        if not b_sol.converged:
            raise ValueError(b_sol.flag)
        else:
            return b_sol.root

    def find_dispersion_for_band_m(self, v_V, m, poln):
        v_b = np.zeros(v_V.shape)

        for iv, V in enumerate(v_V):
            if V > self.cutoff_V(m, poln):
                v_b[iv] = self.find_b_for_V(V, m, poln)
        return v_b

    def find_dispersion_for_bands_at_V(self, v_V, num_bands, poln):
        m_b = np.zeros([len(v_V), num_bands])

        for m in range(num_bands):
            m_b[:, m] = self.find_dispersion_for_band_m(v_V, m, poln)

        return m_b
