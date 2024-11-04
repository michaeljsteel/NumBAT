import math
from enum import Enum
import sys

import matplotlib.pyplot as plt

import numpy as np
import scipy.optimize as sciopt
import scipy.signal
import scipy.special as sp

import reporting

twopi = 2*math.pi
#cvac = nbtypes.SI_speed_c
from math import sqrt, atan


class EMPoln(Enum):
    TE = 'TE'
    TM = 'TM'
    HY = 'HY'


def plot_em_chareq_at_k(k, rcore, ncore, nclad):
    nbrack = 500
    dn = ncore-nclad
    v_neff = np.linspace(ncore-dn/nbrack, nclad+dn/nbrack, nbrack)
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
    fig.savefig('tut_12-em_chareq.png')


def chareq_em_fib2_TE_m(neff, m, k, rco, ncore, nclad):
    # m is un-used (TE is m=0)
    u = rco*k*np.sqrt(ncore**2-neff**2)
    w = rco*k*np.sqrt(neff**2-nclad**2)

    return sp.jv(1, u)/(sp.jv(0, u)*u) + sp.kv(1, w)/(sp.kv(0, w)*w)


def chareq_em_fib2_TM_m(neff, m, k, rco, ncore, nclad):
    # m is un-used (TM is m=0)
    u = rco*k*np.sqrt(ncore**2-neff**2)
    w = rco*k*np.sqrt(neff**2-nclad**2)

    n2rat = (nclad/ncore)**2

    fac1 = np.real(sp.jv(1, u)/(sp.jv(0, u)*u))
    fac2 = np.real(sp.kv(1, w)/(sp.kv(0, w)*w))

    return fac1 + n2rat * fac2


def chareq_em_fib2_hy_m(neff, m, k, rco, ncore, nclad):
    u = rco*k*np.sqrt(ncore**2-neff**2)
    w = rco*k*np.sqrt(neff**2-nclad**2)

    n2rat = (nclad/ncore)**2

    invu2 = 1.0/u**2
    invw2 = 1.0/w**2
    jrat = sp.jvp(m, u)/(sp.jv(m, u)*u)
    krat = sp.kvp(m, w)/(sp.kv(m, w)*w)
    fac1 = jrat + krat
    fac2 = jrat + n2rat * krat
    fac3 = m*m*(invu2 + invw2)*(invu2 + n2rat*invw2)
    return fac1*fac2-fac3


#TODO Seems unused
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


def _solve_chareq_em_fib2_disprel(poln, kvec, nmodes, azi_lo, azi_hi, rco, ncore, nclad):
    # solves for modes with azimuthal order in [mlo, mhi] inclusive

    nbrack = 500
    dn = ncore-nclad

    # look for up to nmodes in [nclad, ncore]
    neff = np.linspace(ncore-dn/nbrack, nclad+dn/nbrack, nbrack)
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
            reporting.report_and_exit(f'Unknown fiber polarisation type: {poln}')

    for m in range(azi_lo, azi_hi+1):

        last_neff = neff[0]
        last_res = f_disprel(last_neff, m, kvec, rco, ncore, nclad)

        nobrack = True  # if we find no roots with a given m, there will be no higher m roots and we can give up
        ineff = 1
        while imode < nmodes and ineff < nbrack:
            t_neff = neff[ineff]
            t_res = f_disprel(t_neff, m, kvec, rco, ncore, nclad)

            if ((poln == EMPoln.HY and last_res * t_res < 0)  # Hybrid eig curves are smooth
                    or (last_res < 0 and t_res > 0)):  # TE and TM curves have +inf to -inf breaks which are not brackets
                # a bracket! go find the root
                nobrack = False
                # root, rootres
                root = sciopt.brentq(f_disprel, last_neff, t_neff,
                                     args=(m, kvec, rco, ncore, nclad))
                v_neff[imode] = root
                imode += 1

            last_neff = t_neff
            last_res = t_res
            ineff += 1

        if nobrack:
            break  # exhausted all roots at this k

    return (imode, v_neff)


class TwoLayerFiberEM(object):
    '''Exact analytic solutions of the two-layer step-index fibre.

       All incoming and outgoing parameters are in base SI units.
    '''

    def __init__(self, ncore, nclad, arad):
        self._ncore = ncore
        self._nclad = nclad
        self._arad = arad

    def Vnum_k(self, k):
        return k*self._arad*math.sqrt(self._ncore**2-self._nclad**2)

    def Vnumb_lam(self, lam):
        return self.Vnum_k(twopi/lam)

    def find_neffs_for_k(self, kvec, poln, azi_m_lo, azi_m_hi, nmax_modes):
        '''
        Solves fiber exact eigenproblem at angular frequency omega in 1/s.

        Returns tuple (nfound, v_k) containing number of modes found and vector
        of eigen wavenumbers k_i[0:nfound-1].
        0 <= nfound <= nmax
        '''
        (nmodes, v_neff) = _solve_chareq_em_fib2_disprel(poln, kvec,
                                                         nmax_modes, azi_m_lo, azi_m_hi,
                                                         self._arad, self._ncore, self._nclad)
        return (nmodes, v_neff)

    def find_neff_HE11_for_k(self, kvec):
        (nmodes, v_neff) = self.find_neffs_for_k(kvec, EMPoln.HY, 1, 1, 1)

        if nmodes == 1:
            return v_neff[0]
        else:
            return None

#######################################################################################


def getqt(Om, q, V):
    dsq = (Om/V)**2 - q**2
    if dsq > 0:
        return np.sqrt(dsq)

    if dsq < 0:
        return complex(0, np.sqrt(-dsq))

    return 0

# Note that these disp rels continue to work fine when qts, qtl go imaginary
# The imaginary terms always come in product pairs giving a real answer

    # p=0 torsional


def chareq_elastic_rod_p0_1d(Om, p, q, rho, c11, c12, c44, a_nm):

    Vs = np.sqrt(c44/rho)
    qts = getqt(Om, q, Vs)
    a = a_nm*1e-9
    qtsa = qts*a

    J0s = sp.jv(0, qtsa)
    J1s = sp.jv(1, qtsa)

    m22 = -qts**2 * J0s + 2*qts/a*J1s

    chareq = m22/Om**2

    return np.real(chareq) + np.imag(chareq)

    # p=0 Pochammer


def chareq_elastic_rod_p0_2d(Om, p, q, rho, c11, c12, c44, a_nm):
    Vl = np.sqrt(c11/rho)
    Vs = np.sqrt(c44/rho)

    qts = getqt(Om, q, Vs)
    qtl = getqt(Om, q, Vl)

    a = a_nm*1e-9

    qtsa = qts*a
    qtla = qtl*a

    J0l = sp.jv(p, qtla)
    J0pl = sp.jvp(p, qtla)
    J0ppl = 0.5*(sp.jvp(p-1, qtla) - sp.jvp(p+1, qtla))

    J1s = sp.jv(1, qtsa)
    J1ps = sp.jvp(1, qtsa)

    m00 = 2*c44 * qtl**2 * J0ppl - c12*(qtl**2 + q**2)*J0l
    m01 = 2j*c44 * q * qts * J1ps

    m10 = -2j * q * qtl * J0pl
    m11 = (-qts**2+q**2) * J1s

    d1 = m00*m11
    d2 = m01*m10

    disc = d1-d2

    # Because of catastrophic cancellation, can't reasonably hope for closer to zero than this
    tol = 1e-12
    if np.abs(disc) < tol*max(abs(d1), abs(d2)):
        chareq = 0
    else:
        chareq = disc/(Om**4)  # Om**4 makes numbers nice sie

    # Branch cuts behave so that it Works nicely to just add real and imag parts.
    return chareq.real + chareq.imag


def chareq_elastic_rod_ppos(Om, p, q, rho, c11, c12, c44, a_nm):  # p is azimuthal order

    Vl = np.sqrt(c11/rho)   # = sqrt((lambda + 2 mu)/rho)
    Vs = np.sqrt(c44/rho)   # = sqrt(mu/rho)

    qts = getqt(Om, q, Vs)
    qtl = getqt(Om, q, Vl)

    a = a_nm*1e-9

    qtsa = qts*a
    qtla = qtl*a

    J0l = sp.jv(p, qtla)
    J0pl = sp.jvp(p, qtla)
    J0ppl = 0.5*(sp.jvp(p-1, qtla) - sp.jvp(p+1, qtla))

    J0s = sp.jv(p, qtsa)
    J0ps = sp.jvp(p, qtsa)
    J0pps = 0.5*(sp.jvp(p-1, qtsa) - sp.jvp(p+1, qtsa))

    J1s = sp.jv(p+1, qtsa)
    J1ps = sp.jvp(p+1, qtsa)

    m00 = 2*c44 * qtl**2 * J0ppl - c12*(qtl**2 + q**2)*J0l
    m01 = 1j*2*c44 * q * qts * J1ps
    m02 = 2*c44*p/a * (J0ps*qts-J0s/a)

    m10 = 1j*2 * q * qtl * J0pl
    m11 = J1s*((qts**2-q**2) - (p+p*p)/(a*a)) - J1ps*p*qts/a
    m12 = 1j * q * p * J0s/a

    m20 = 2*p/a*(qtl * J0pl-J0l/a)
    m21 = 1j*J1s*(q/a+p*q/a)-1j*J1ps*q*qts
    m22 = qts**2 * (2*J0pps+J0s)

    d0a = m00*m11*m22
    d0b = -m00*m12*m21
    d1a = m01*m12*m20
    d1b = -m01*m10*m22
    d2a = m02*m10*m21
    d2b = -m02*m11*m20

    bigd = np.max(np.abs(np.array([d0a, d0b, d1a, d1b, d2a, d2b])))
    disc = d0a+d0b+d1a+d1b+d2a+d2b

    tol = 1e-12  # Because of catastrophic cancellation, can't hope for closer to zero than this
    if np.abs(disc) < tol*bigd:
        chareq = 0
    else:
        chareq = disc/(Om**5)  # Om**5 makes numbers nice sie

    # Branch cuts behave so that it Works nicely to just add real and imag parts.
    return chareq.real + chareq.imag


def chareq_elastic_rod(Om, m, q, rho, c11, c12, c44, arad):  # m is azimuthal order
    a_nm = arad*1e9
    if m == -1:  # actually denotes a p=0 1D state
        return chareq_elastic_rod_p0_1d(Om, 0, q, rho, c11, c12, c44, a_nm)

    if m == 0:
        return chareq_elastic_rod_p0_2d(Om, m, q, rho, c11, c12, c44, a_nm)

    return chareq_elastic_rod_ppos(Om, m, q, rho, c11, c12, c44, a_nm)


def _findroots_elastic_rod_chareq(q, nu_hi, m, nmodes, Vl, rho, c11, c12, c44, arad):

    num_Om = 500
    Omlo = 0.01e6
    Omhi = twopi*nu_hi

    # look for up to nmodes in [Omlo, Omhi]
    Om = np.linspace(Omlo, Omhi, num_Om)
    drval = np.zeros(len(Om))

    v_Om = np.zeros(nmodes, dtype=np.float64)
    imode = 0

    t_Om = Om[0]
    t_res = chareq_elastic_rod(t_Om, m, q, rho, c11, c12, c44, arad)
    drval[0] = t_res

    for iOm in range(1, num_Om):
        last_Om = t_Om
        last_res = t_res

        t_Om = Om[iOm]
        t_res = chareq_elastic_rod(t_Om, m, q, rho, c11, c12, c44, arad)
        drval[iOm] = t_res

        if (imode < nmodes and last_res * t_res < 0):  # a bracket! go find the root
            root, root_res = sciopt.brentq(chareq_elastic_rod, last_Om, t_Om,
                                           args=(m, q, rho, c11, c12, c44, arad), full_output=True, disp=False)

            if not root_res.converged:

                # Must be at a nasty point. Try bisection which can't really fail.
                root, root_res = sciopt.bisect(chareq_elastic_rod, last_Om, t_Om,
                                               args=(m, q, rho, c11, c12, c44, arad), full_output=True, disp=False)
            if root_res.converged:
                v_Om[imode] = root
                imode += 1
            else:
                print('Brent and bisection unable to converge for p=%d, q=%.4f/um, nu in [%.6f, %.6f] GHz, err in [%.4e, %.4e]. ' % (
                    m, q*1e-6, last_Om/(2e9*math.pi), t_Om/(2e9*math.pi), last_res, t_res))

    # This dispersion relation tends to exhibit double roots that may not trigger a bracketing.
    # Instead we can look for any points which are local min/max and very close to zero
    # Picking the right tolerances here is rather sensitive
    bigdr = np.max(np.abs(drval))
    smalldr = np.min(np.abs(drval))
    if smalldr < 0.001 * bigdr:  # There is a significant degree of smallness worth checking
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
            if qrt in [0, len(Om)-1]:
                continue

            # Bracketed root, so will have been picked up above by brentq
            if (drval[qrt-1] * drval[qrt] < 0) or (drval[qrt] * drval[qrt+1] < 0):
                continue

            droot = drval[qrt]
            if np.abs(droot) < 1e-8 and np.abs(droot) < 1e-9 * bigdr:
                # found a likely double root, fit a parabola to get accurate location
                a = (drval[qrt+1]+drval[qrt-1]-2*drval[qrt])/2
                b = (drval[qrt+1]-drval[qrt-1])/2
                v_Om[imode] = Om[qrt] - b/(2*a)*(Om[1]-Om[0])
                print(f'Found an additional double root for p={m} at nu={v_Om[imode]/(1e9*twopi):.4e}.')
                imode += 1

    return (imode, v_Om)


class ElasticRod(object):
    '''Exact analytic solutions of the elastic problem for an isotropic cylinder in vacuum.

       All incoming and outgoing parameters are in base SI units.
    '''

    def __init__(self, rho, c11, c12, c44, arad):
        self._rho = rho
        self._c11 = c11
        self._c12 = c12
        self._c44 = c44
        self._arad = arad
        self._Vl = np.sqrt(c11/rho)
        self._Vs = np.sqrt(c44/rho)

    def find_nu_for_q(self, q, nu_hi, m, nmax_modes):
        (nmodes, v_Om) = _findroots_elastic_rod_chareq(q, nu_hi, m, nmax_modes, self._Vl,
                                                       self._rho, self._c11, self._c12, self._c44,
                                                       self._arad)
        return (nmodes, v_Om)

    def dispersion_relation_at_q_nu(self, q, nu, m):

        return chareq_elastic_rod(np.pi*2*nu, m, q,
                                  self._rho, self._c11, self._c12, self._c44, self._arad)  # m is azimuthal order




def emslab_chareq_TE(V, gamma, m, b):
    '''Dispersion relation for TE slab optical waveguide in normalised units.'''
    return V * sqrt(1-b) - m*np.pi - atan(sqrt(b/(1-b))) - atan(sqrt((b+gamma)/(1-b)))

def emslab_chareq_TM(V, gamma, ns, nf, nc, m, b):
    '''Dispersion relation for TE slab optical waveguide in normalised units.'''
    tmfac1 = (nf/ns)**2
    tmfac2 = (nf/nc)**2
    return V * sqrt(1-b) - m*np.pi - atan(tmfac1*sqrt(b/(1-b))) - atan(tmfac2*sqrt((b+gamma)/(1-b)))


class EMSlab:
    '''EM slab waveguide  solver.

       All units are in microns.
    '''

    def __init__(self, ns, nf, nc, wid):
        # Allows incoming complex value refractive indices but makes them real

        self._ns = ns.real  # substrate index
        self._nf = nf.real  # film index
        self._nc = nc.real  # cover index
        self._wid = wid.real   # width

        if not self._nc <= self._ns < self._nf:
            raise ValueError('Indices do not specify a guiding structure.')
        self._gamma = (self._ns**2-self._nc**2)/(self._nf**2-self._ns**2)

        self._tmfac_s = (self._nf/self._ns)**2
        self._tmfac_c = (self._nf/self._nc)**2

    def lambda_to_V(self, lam):
        return 2*np.pi/lam * self._wid * sqrt(self._nf**2-self._ns**2)

    def V_to_lambda(self, V):
        return 2*np.pi/V * self._wid * sqrt(self._nf**2-self._ns**2)


    def b_to_neff(self, b):
        return sqrt(self._ns**2 + (self._nf**2 - self._ns**2) * b)

    def cutoff_V(self, m, poln):
        match poln:
            case EMPoln.TE:
                return m*np.pi + atan(sqrt(self._gamma))

            case EMPoln.TM:
                return m*np.pi + atan(self._tmfac_c * sqrt(self._gamma))

    def num_modes(self, V, poln):
        match poln:

            case EMPoln.TE:
                Vcorr = V - atan(sqrt(self._gamma))
                return math.floor(Vcorr/np.pi) + 1

            case EMPoln.TM:
                tmfac =  (self._nf/self._nc)**2
                Vcorr = V - atan(tmfac * sqrt(self._gamma))
                return math.floor(Vcorr/np.pi) + 1

    def find_b_for_V(self, V, m, poln):
        '''Assumes that cutoffs have been checked and there is a solution.'''

        match poln:

            case EMPoln.TE:
                bfunc = lambda b: emslab_chareq_TE(V, self._gamma, m, b)

            case EMPoln.TM:
                bfunc = lambda b: emslab_chareq_TM(V, self._gamma, self._ns, self._nf, self._nc, m, b)

        tol = 1e-10
        b_sol = sciopt.root_scalar(bfunc, bracket=[tol,1-tol])

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

    def find_dispersion_for_bands(self, v_V, num_bands, poln):
        m_b = np.zeros([len(v_V), num_bands])

        for m in range(num_bands):
            m_b[:, m] = self.find_dispersion_for_band_m(v_V, m, poln)

        return m_b


def elasticfreeslab_Lamb_chareq_even(Omega, Vs, Vl, wid, q):
    '''Dispersion relation for TE slab optical waveguide in normalised units.'''

    'Dispersion relation makes sense if both kappas are real or both imaginary but there seem to be no solutions'
    ONE = 1.0+0.0j
    kappa_s = np.sqrt(ONE*(Omega/Vs)**2 - q**2)
    kappa_l = np.sqrt(ONE*(Omega/Vl)**2 - q**2)

    #lhs = np.tan(kappa_l * wid/2)/np.tan(kappa_s * wid/2)
    #rhs = - (q**2-kappa_s**2)**2/(4 *q**2 *  kappa_l * kappa_s)

    lhs = np.tan(kappa_l * wid/2)*(4 *q**2 *  kappa_l * kappa_s)
    rhs = - (q**2-kappa_s**2)**2 * np.tan(kappa_s * wid/2)


    eq = wid**4 * (lhs - rhs)
    if abs(np.imag(eq)) > abs(np.real(eq)):
        return np.imag(eq)
    else:
        return np.real(eq)

def elasticfreeslab_Lamb_chareq_odd(Omega, Vs, Vl, wid, q):
    '''Dispersion relation for TE slab optical waveguide in normalised units.'''
    ONE = 1.0+0.0j
    kappa_s = np.sqrt(ONE*(Omega/Vs)**2 - q**2)
    kappa_l = np.sqrt(ONE*(Omega/Vl)**2 - q**2)

    #lhs = np.tan(kappa_s * wid/2)/np.tan(kappa_l * wid/2)
    #rhs = - (q**2-kappa_s**2)**2/(4 *q**2 * kappa_l * kappa_s)
    lhs = np.tan(kappa_s * wid/2)*(4 *q**2 * kappa_l * kappa_s)
    rhs = - (q**2-kappa_s**2)**2 * np.tan(kappa_l * wid/2)

    eq = wid**4  * (lhs - rhs)
    if abs(np.imag(eq)) > abs(np.real(eq)):
        return np.imag(eq)
    else:
        return np.real(eq)



class ElasticFreeSlab:
    '''Elastic slab waveguide solver for isotropic materials.
       Finds the dispersion of the Lamb modes of isolated flat plate
       All units are in microns.

       Dispersion relation has even and odd parts

       Even:
         tan(kappa_l w/2)/tan(kappa_s w/2) = - (q^2-kappa_s^2)/(4 q^2 kappa_l kappa_s)

       Odd:
         tan(kappa_s w/2)/tan(kappa_l w/2) = - (q^2-kappa_s^2)/(4 q^2 kappa_l kappa_s)

       where kappa_s^2 = (Omega/V_s)^2 - q^2
       where kappa_l^2 = (Omega/V_l)^2 - q^2

       The two sides switch from real and imag at the same time so all makes sense at all values of q
    '''

    def __init__(self, mat_s, wid):
        self._mats = mat_s  # slab material
        self._Vs = mat_s.Vac_shear()

        #self._Vl = mat_s.Vac_longitudinal()
        self._Vl = self._Vs*1.9056
        self._wid = wid


    def _find_Lamb_q_brackets_smart(self, Omega, drfunc):
        # attempts to predict the crossings but doesn't get them all, becuse not all are about tan blowups

        qlo = 1e-9
        qhi=Omega/self._Vs * 2
        # need to avoid blowsups of tan functions
        # kappa_s w/2 = (2n+1) pi/2
        # kappa_l w/2 = (2n+1) pi/2
        hi_ns = np.ceil(self._wid*Omega/(2*np.pi*self._Vs)-.5)
        hi_nl = np.ceil(self._wid*Omega/(2*np.pi*self._Vl)-.5)

        bad_qs = ((Omega/self._Vs)**2 - ( (np.arange(0,hi_ns) + .5)*2*np.pi/self._wid)**2)**0.5
        bad_ql = ((Omega/self._Vl)**2 - ((
            np.arange(0,hi_nl)+.5)*2*np.pi/self._wid)**2)**0.5

        #print(' hins', Omega, self._Vs,hi_ns, hi_nl, self._wid*Omega/(2*np.pi*self._Vl)-.5, bad_qs, bad_ql)

        bad_q = np.append(bad_qs, bad_ql)
        bad_q = np.sort(bad_q)
        qbraks=[]
        tollo = 1-1e-7
        tolhi = 1+1e-7


        # print('Doing Omega', Omega)
        # qsteps=2000
        # qlo = 1e-9
        # qhi=Omega/self._Vs * 2
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
        #     axs.plot(qlo*self._wid,0,'x')
        # qbraks.append([qlo, qhi])



        # fig.savefig('qscan.png')
        # sys.exit(0)

        return qbraks


    def at_tan_resonance_om(self, ombrak, q):
        if (ombrak[0]/self._Vs)**2-q**2 >0:
            sarg0 =  np.sqrt((ombrak[0]/self._Vs)**2-q**2)*self._wid/2/(np.pi/2)
            sarg1 =  np.sqrt((ombrak[1]/self._Vs)**2-q**2)*self._wid/2/(np.pi/2)
            if math.ceil(sarg0) == math.floor(sarg1) and math.ceil(sarg0)%2 ==1:  # an odd multiple of pi/2
                return True

        if (ombrak[0]/self._Vl)**2-q**2 >0:
            larg0 = np.sqrt((ombrak[0]/self._Vl)**2-q**2)*self._wid/2/(np.pi/2)
            larg1 =  np.sqrt((ombrak[1]/self._Vl)**2-q**2)*self._wid/2/(np.pi/2)
            if math.ceil(larg0) == math.floor(larg1) and math.ceil(larg0)%2 ==1:  # an odd multiple of pi/2
                return True
        return False

    def at_tan_resonance(self, Omega, qbrak):
        if (Omega/self._Vs)**2-qbrak[1]**2 >0:
            sarg0 =  np.sqrt((Omega/self._Vs)**2-qbrak[0]**2)*self._wid/2/(np.pi/2)
            sarg1 =  np.sqrt((Omega/self._Vs)**2-qbrak[1]**2)*self._wid/2/(np.pi/2)
            if math.ceil(sarg0) == math.floor(sarg1) and math.ceil(sarg0)%2 ==1:  # an odd multiple of pi/2
                return True

        if (Omega/self._Vl)**2-qbrak[1]**2 >0:
            larg0 = np.sqrt((Omega/self._Vl)**2-qbrak[0]**2)*self._wid/2/(np.pi/2)
            larg1 =  np.sqrt((Omega/self._Vl)**2-qbrak[1]**2)*self._wid/2/(np.pi/2)
            if math.ceil(larg0) == math.floor(larg1) and math.ceil(larg0)%2 ==1:  # an odd multiple of pi/2
                return True
        return False

    def _find_Lamb_q_brackets(self, Omega, drfunc):
        # brute force bracketing
        qsteps=1000

        qlo = 1e-9
        qhi=Omega/self._Vs * 2
        v_q = np.linspace(qlo, qhi,  qsteps)


        v_disprel = np.array(list(map(drfunc, v_q)))


        qbraks = []

        is_small = 0.1*max(abs(v_disprel))
        fim1 = v_disprel[0]
        for i in range(1,len(v_disprel)-1):
            fi = v_disprel[i]
            if fi * fim1 < 0 :  # a possible crossing?
                if abs(fi) < is_small and abs(fim1) < is_small: # legit
                    qbrak = [v_q[i-1], v_q[i]]
                    # remove erroneous  brackets due to tan functions
                    if not self.at_tan_resonance(Omega, qbrak):
                        qbraks.append(qbrak)
                #else: # a tan function infinity, update but don't keep
                #    pass
            fim1 = fi
        # print(qbraks)
        # fig,axs=plt.subplots()
        # axs.plot(v_q*1e-6, v_disprel)
        # axs.plot(v_q*1e-6,0*v_q,':')
        # axs.set_ylim(-1e6,1e6)
        # #axs.set_xlim(0,20)
        # for qb in qbraks:
        #     axs.plot(qb[0]*self._wid,0,'x', ms=20)
        # fig.savefig('qscan.png')
        # sys.exit(0)

        return qbraks

    def _find_Lamb_om_brackets(self, q, drfunc):
        # brute force bracketing
        omsteps=5000

        omlo = 1e-9
        omhi=100e9
        v_om = np.linspace(omlo, omhi,  omsteps)


        v_disprel = np.array(list(map(drfunc, v_om)))


        ombraks = []
        removeds =[]

        is_small = 0.1*max(abs(v_disprel))
        fim1 = v_disprel[0]
        for i in range(1,len(v_disprel)-1):
            fi = v_disprel[i]
            if fi * fim1 < 0 :  # a possible crossing?
                if abs(fi) < is_small and abs(fim1) < is_small: # legit
                    ombrak = [v_om[i-1], v_om[i]]
                    if not self.at_tan_resonance_om(ombrak, q):
                        ombraks.append(ombrak)
                    else:
                        removeds.append(ombrak)
                #else: # a tan function infinity, update but don't keep
                #    pass
            fim1 = fi
        #plotscan=True
        plotscan=False
        if plotscan:
            print(ombraks)
            fig,axs=plt.subplots()
            v_omnorm = v_om*self._wid/self._Vs/np.pi
            #axs.plot(v_omnorm, np.log10(abs(v_disprel)))
            axs.plot(v_omnorm,  v_disprel)

            axs.plot(v_omnorm,0*v_om,':')
            #axs.set_ylim(-1*is_small, 1*is_small)
            axs.set_ylim(-1e3,1e3)
            axs.set_xlim(1.3,1.5)
            for ob in ombraks:
                obnorm = ob[0] * self._wid/self._Vs/np.pi
                axs.plot(obnorm,0,'x', ms=20)
            for ob in removeds:
                obnorm = ob[0] * self._wid/self._Vs/np.pi
                axs.plot(obnorm,0,'o',mfc='none', ms=20)


            fig.savefig('omscan.png')
            sys.exit(0)

        return ombraks


    def find_Lamb_Oms_for_q(self, q, max_modes, even_modes=True):
        # q must be smaller than Omega/V_l

        if even_modes:
            drfunc = lambda om: elasticfreeslab_Lamb_chareq_even(om, self._Vs, self._Vl, self._wid, q)
        else:
            drfunc = lambda om: elasticfreeslab_Lamb_chareq_odd(om, self._Vs, self._Vl, self._wid, q)

        omsols = []

        #qbraks = self._find_Lamb_q_brackets_smart(Omega, drfunc)
        ombraks = self._find_Lamb_om_brackets(q, drfunc)

        #print('got braks', qbraks)
        for ombrak in ombraks:
            if drfunc(ombrak[0]) * drfunc(ombrak[1])>0:
                print('False bracket', ombrak, drfunc(ombrak[0]), drfunc(ombrak[1]))
                continue

            om_sol = sciopt.root_scalar(drfunc, bracket=ombrak)
            if not om_sol.converged:
                raise ValueError(om_sol.flag)
            else:
                omsols.append(om_sol.root)

            if len(omsols) == max_modes: break

        return omsols

    def find_Lamb_qs_for_Omega(self, Omega, max_modes, even_modes=True):
        # q must be smaller than Omega/V_l

        if even_modes:
            drfunc = lambda q: elasticfreeslab_Lamb_chareq_even(Omega, self._Vs, self._Vl, self._wid, q)
        else:
            drfunc = lambda q: elasticfreeslab_Lamb_chareq_odd(Omega, self._Vs, self._Vl, self._wid, q)

        qsols = []

        #qbraks = self._find_Lamb_q_brackets_smart(Omega, drfunc)
        qbraks = self._find_Lamb_q_brackets(Omega, drfunc)

        #print('got braks', qbraks)
        for qbrak in qbraks:
         #   print('doing brak ', Omega, qbrak)
            if drfunc(qbrak[0]) * drfunc(qbrak[1])>0:
                print('False bracket', qbrak, drfunc(qbrak[0]), drfunc(qbrak[1]))
                continue

            q_sol = sciopt.root_scalar(drfunc, bracket=qbrak)
            if not q_sol.converged:
                raise ValueError(q_sol.flag)
            else:
                qsols.append(q_sol.root)

            if len(qsols) == max_modes: break

        return qsols

    def find_Lamb_dispersion_for_bands(self, v_Omega, max_modes, even_modes=True):
        m_qs = np.zeros([len(v_Omega), max_modes])

        for iOm, Om in enumerate(v_Omega):
            qsols = self.find_Lamb_qs_for_Omega(Om, max_modes, even_modes)

            for iq, q in enumerate(qsols):
                m_qs[iOm, iq] = q

        return m_qs

    def find_Lamb_dispersion_for_bands_of_q(self, v_q, max_modes, even_modes=True):
        m_Om = np.zeros([len(v_q), max_modes])

        for iq, q in enumerate(v_q):
            omsols = self.find_Lamb_Oms_for_q(q, max_modes, even_modes)

            for iom, om in enumerate(omsols):
                m_Om[iq, iom] = om

        return m_Om

    def disprel_rayleigh(self,  vR):
        vl = self._Vl/self._Vs

        return vR**6 - 8 * vR**4 + vR**2*(24-16/vl**2 )+16*(1/vl**2-1)

    def find_Rayleigh_dispersion(self, v_Omega):
        # find v_q of Rayleigh mode for v_Omega
        m_qs = np.zeros(len(v_Omega))


        dr_rayleigh = lambda v: self.disprel_rayleigh(v)

        # Calculation works in units of self._Vs
        vRlo = 0.001
        vRhi = 1

        vv = np.linspace(vRlo, vRhi,100)


        vres = sciopt.root_scalar(dr_rayleigh, bracket=(vRlo, vRhi))
        if not vres.converged:
                raise ValueError(vres.flag)
        else:
            vR = vres.root


        # put back units of self._Vs
        return v_Omega/(vR* self._Vs)


class ElasticSlab:
    '''Elastic slab waveguide solver for isotropic materials.
       Finds the dispersion of the Lamb modes of isolated flat plate
       All units are in microns.
    '''

    def __init__(self, mat_s, mat_f, mat_c, wid):
        self._mats = mats  # substrate material
        self._matf = matf  # film material
        self._matc = matc  # cover material
        self._wid = wid    # width


