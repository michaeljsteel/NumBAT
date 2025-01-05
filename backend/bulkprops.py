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
import scipy.linalg

import voigt


def power_flux_christoffel(kapv, v_p, evec, c_stiff):
    r"""Evaluates the power flux P=-v^* \cdot T for a given unit wavevector kapv, eigenvector evec and implicit wavenumber k and frequency omega.

    Factors of 2 seem to be correct here, but good to write out in full in docs.
    """

    # S_I = \nabla_Ij (uj e^i(k kap . r)) = i k (\nabla_Ij  e^i(k kap . r)) . uj
    # Evaluate the S_I 6x1 vector, Auld 1.50, 1.49
    S_I = 1j * np.matmul(voigt.kvec_to_symmetric_gradient(kapv), evec)

    # Auld 3.20   # Indices off by 1 from zero count
    T_I = np.matmul(c_stiff.value(), S_I)

    T_ij = voigt.stress_6col_to_3mat(T_I)

    # Pcomp = - 1/2 v^* . T ,   Auld 2.30       , 5.77
    #       =  -1/2 (i \omega u^*) .  T
    # om = 1.0 # unit k and omega
    # vsx, vsy, vsz = 1j*np.conj(evec)
    # Pcomp = - 0.5 * np.array([
    #    vsx*T_I[0] + vsy*T_I[5] + vsz*T_I[4],
    #    vsx*T_I[5] + vsy*T_I[1] + vsz*T_I[3],
    #    vsx*T_I[4] + vsy*T_I[3] + vsz*T_I[2] ])

    Pcomp = -0.5 * 1j * np.matmul(np.conj(evec), T_ij)

    # u_s = 1/2  k^2 S_I c_IJ SJ  # Auld, 5.35
    # real vs complex fields?
    u_s = 0.5 * np.matmul(np.matmul(S_I, c_stiff.value()), S_I)

    # vg = Pcomp/u_s ->  Pcomp/us   (omega k)/(k^2) = v_p Pcomp/us
    v_g = -np.real(v_p * Pcomp / u_s)

    return v_g


def Gamma_christoffel(vkap, c_stiff, rho):
    """Returns Gamma_ij = 1/V0^2   mD.Cij.md^T/rho  in units of (km/s)^2
    vkap is unit wavevector.
    V0=1km/s

    See Auld V1. Sec 7.D
    """

    mD = voigt.kvec_to_symmetric_gradient(vkap).T
    v0sq = 1.0e6
    m_Gamma = np.matmul(np.matmul(mD, c_stiff.value()), mD.T) / (v0sq * rho)

    return m_Gamma


def chareq_christoffel(vkap, c_stiff, rho, v_p):
    r"""Returns Om=|Gamma_ij -\omega^2 I|, the characteristic function of the Christoffel equation.

    v_p is in km/s
    See Auld V1. Sec 7.D
    """

    op_chris = Gamma_christoffel(vkap, c_stiff, rho) - v_p**2 * np.eye(3)

    return scipy.linalg.det(op_chris)


def solve_christoffel(vkap, c_stiff, rho):
    """Solve eigenproblem of Christoffel equation in the direction vkap (a 2D unit vector). Returns for each of 3 modes:
        phase velocity                v_phase[m]
        polarisation eigenvetors      evecs[:,m]
        group velocity vectors.       v_group[m:x/y/z]

    Modes are sorted by decreasing phase velocity.

    See Auld V1. Sec 7.D
    """

    m_Gamma = Gamma_christoffel(vkap, c_stiff, rho)

    # Solve and normalise
    evals, evecs = scipy.linalg.eig(m_Gamma)
    for i in range(3):
        evecs[:, i] /= np.linalg.norm(evecs[:, i])
        # TODO: make a oneliner:"
        # evecs *= 1/np.sqrt(np.diag(np.real(evecs.T @ evecs)))

    vphases = np.sqrt(np.real(evals))  # result is in km/s

    # Sort according to velocity

    ivs = np.argsort(-vphases)  # most likely get pwave first

    v_vphase = np.sqrt(np.real(evals[ivs])).copy()
    v_evecs = evecs[:, ivs].copy()

    # now look for vg here
    # vg = - nabla_k Om/ dOm/dom = nabla_kappa Om/ dOm/dvp =

    v_vgroup = np.zeros([3, 3])  # first index is mode, second is component of \vec{v}_g

    for m in range(3):  # for each mode at this vkap
        v_p = v_vphase[m]
        v_vgroup[m, :] = power_flux_christoffel(vkap, v_p, v_evecs[:, m], c_stiff)

    return v_vphase, v_evecs, v_vgroup
