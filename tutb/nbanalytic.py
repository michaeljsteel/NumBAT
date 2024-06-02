from enum import Enum
import math
import numpy as np
import scipy.special as sp
import scipy.optimize as sciopt
import scipy.signal
import matplotlib.pyplot as plt


import reporting
import nbtypes


twopi = 2*math.pi
cvac = nbtypes.speed_c


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
            v_dr_TE[i] = chareq_em_fib2_TE_m(
                neff, m, k, rcore, ncore, nclad)+.1
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
