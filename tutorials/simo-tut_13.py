""" Example showing how the 'onion' geometry template can be used
    to calculate the dispersion profile of standard SMF-28 optical fibre.
"""

import os
import time
import sys
import queue
import copy
import threading
import math
import scipy.signal
import scipy.optimize as sciopt
import scipy.special as sp

import matplotlib.pyplot as plt
import numpy as np


sys.path.append("../backend/")

import materials
import objects
import plotting
from numbattools import launch_worker_threads_and_wait

import starter

twopi = math.pi*2.0

# a helper function to make consistent pretty plots
def decorate_and_save_plot (fig, ax, fname, title, xlab, ylab, xlims, ylims,
                            legend_loc=None, legend_ncol=1):
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_xlim(xlims[0], xlims[1])
    ax.set_ylim(ylims[0], ylims[1])
    ax.set_title(title)
    if legend_loc is not None:
        ax.legend(loc=legend_loc, fontsize='xx-small', ncol=legend_ncol)
    else:
        ax.legend(fontsize='xx-small', ncol=legend_ncol)
    fig.savefig(fname)

def plot_and_label (ax, vx, mat, sty, lab, col=None):
    if col is None:
        ax.plot(vx, mat, sty, ms=5)
        ax.plot(vx, mat[:,0], sty, label=lab, ms=5)
    else:
        ax.plot(vx, mat, sty, color=col, ms=5)
        ax.plot(vx, mat[:,0], sty, label=lab, color=col, ms=5)

#
#
#Elastic analytic calculations
#
#

def getqt(Om, q, V):
    dsq = (Om/V)**2 - q**2
    if dsq>0:
        return np.sqrt(dsq)

    if dsq<0:
        return complex(0, np.sqrt(-dsq))

    return 0

#Note that these disp rels continue to work fine when qts, qtl go imaginary
#The imaginary terms always come in product pairs giving a real answer

    #p=0 torsional
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

    #p=0 Pochammer
def chareq_elastic_rod_p0_2d(Om, p, q, rho, c11, c12, c44, a_nm):
    Vl = np.sqrt(c11/rho)
    Vs = np.sqrt(c44/rho)

    qts = getqt(Om, q, Vs)
    qtl = getqt(Om, q, Vl)

    a = a_nm*1e-9

    qtsa=qts*a
    qtla=qtl*a


    J0l = sp.jv(p, qtla)
    J0pl =sp.jvp(p, qtla)
    J0ppl = 0.5*( sp.jvp(p-1, qtla) - sp.jvp(p+1, qtla))

    J1s = sp.jv(1, qtsa)
    J1ps = sp.jvp(1, qtsa)

    m00 = 2*c44 * qtl**2 * J0ppl - c12*(qtl**2 + q**2)*J0l
    m01 = 2j*c44 * q * qts * J1ps

    m10 = -2j * q * qtl * J0pl
    m11 = (-qts**2+q**2) * J1s

    d1=m00*m11
    d2=m01*m10

    disc = d1-d2

    # Because of catastrophic cancellation, can't reasonably hope for closer to zero than this
    tol=1e-12
    if np.abs(disc) < tol*max(abs(d1),abs(d2)):
        chareq=0
    else:
        chareq = disc/(Om**4) # Om**4 makes numbers nice sie

    # Branch cuts behave so that it Works nicely to just add real and imag parts.
    return chareq.real +chareq.imag

def chareq_elastic_rod_ppos(Om, p, q, rho, c11, c12, c44, a_nm): #p is azimuthal order

    Vl = np.sqrt(c11/rho)   # = sqrt((lambda + 2 mu)/rho)
    Vs = np.sqrt(c44/rho)   # = sqrt(mu/rho)

    qts = getqt(Om, q, Vs)
    qtl = getqt(Om, q, Vl)

    a = a_nm*1e-9

    qtsa=qts*a
    qtla=qtl*a

    J0l = sp.jv(p, qtla)
    J0pl = sp.jvp(p, qtla)
    J0ppl = 0.5*( sp.jvp(p-1, qtla) - sp.jvp(p+1, qtla))

    J0s = sp.jv(p, qtsa)
    J0ps = sp.jvp(p, qtsa)
    J0pps = 0.5*( sp.jvp(p-1, qtsa) - sp.jvp(p+1, qtsa))

    J1s = sp.jv(p+1, qtsa)
    J1ps = sp.jvp(p+1, qtsa)


    m00 = 2*c44 * qtl**2 * J0ppl - c12*(qtl**2 + q**2)*J0l
    m01 = 1j*2*c44 * q * qts * J1ps
    m02 = 2*c44*p/a *( J0ps*qts-J0s/a)

    m10 = 1j*2 * q * qtl * J0pl
    m11 = J1s*((qts**2-q**2) - (p+p*p)/(a*a) ) -J1ps*p*qts/a
    m12 = 1j* q * p * J0s/a

    m20 = 2*p/a*(qtl * J0pl-J0l/a)
    m21 = 1j*J1s*(q/a+p*q/a)-1j*J1ps*q*qts
    m22 = qts**2 * (2*J0pps+J0s)

    d0a =  m00*m11*m22
    d0b = -m00*m12*m21
    d1a =  m01*m12*m20
    d1b = -m01*m10*m22
    d2a =  m02*m10*m21
    d2b = -m02*m11*m20

    bigd= np.max(np.abs(np.array([d0a,d0b, d1a,d1b, d2a,d2b])))
    disc = d0a+d0b+d1a+d1b+d2a+d2b

    tol=1e-12  # Because of catastrophic cancellation, can't hope for closer to zero than this
    if np.abs(disc) < tol*bigd:
        chareq=0
    else:
        chareq = disc/(Om**5) # Om**5 makes numbers nice sie

    # Branch cuts behave so that it Works nicely to just add real and imag parts.
    return chareq.real +chareq.imag



def chareq_elastic_rod(Om, p, q, rho, c11, c12, c44, a_nm): #p is azimuthal order
    if p==-1:  # actually denotes a p=0 1D state
        return chareq_elastic_rod_p0_1d(Om, 0, q, rho, c11, c12, c44, a_nm)

    if p==0:
        return chareq_elastic_rod_p0_2d(Om, p, q, rho, c11, c12, c44, a_nm)

    return chareq_elastic_rod_ppos(Om, p, q, rho, c11, c12, c44, a_nm)

def findroots_elastic_rod_chareq(Vl, nu_hi, p, q, nmodes, rho, c11, c12, c44, arad):

    num_Om = 500
    Omlo = 0.01e6
    Omhi = twopi*nu_hi

    # look for up to nmodes in [Omlo, Omhi]
    Om = np.linspace(Omlo, Omhi, num_Om)
    drval = np.zeros(len(Om))

    sol_Om = np.zeros(nmodes, dtype=np.float64)
    imode = 0

    t_Om = Om[0]
    t_res =chareq_elastic_rod(t_Om, p, q, rho, c11, c12, c44, arad)
    drval[0] = t_res

    for iOm in range(1,num_Om):
        last_Om = t_Om
        last_res  = t_res

        t_Om = Om[iOm]
        t_res =  chareq_elastic_rod(t_Om, p, q, rho, c11, c12, c44, arad)
        drval[iOm] = t_res

        if (imode < nmodes and last_res * t_res < 0): # a bracket! go find the root
            root, root_res  = sciopt.brentq(chareq_elastic_rod, last_Om, t_Om,
                    args=(p, q, rho, c11, c12, c44, arad), full_output=True, disp=False)

            if not root_res.converged:

                #Must be at a nasty point. Try bisection which can't really fail.
                root, root_res  = sciopt.bisect(chareq_elastic_rod, last_Om, t_Om,
                    args=(p, q, rho, c11, c12, c44, arad), full_output=True, disp=False)
            if root_res.converged:
                sol_Om[imode] = root
                imode += 1
            else:
                print('Brent and bisection unable to converge for p=%d, q=%.4f/um, nu in [%.6f, %.6f] GHz, err in [%.4e, %.4e]. '%(
                    p,q*1e-6, last_Om/(2e9*math.pi), t_Om/(2e9*math.pi), last_res, t_res))


    # This dispersion relation tends to exhibit double roots that may not trigger a bracketing.
    # Instead we can look for any points which are local min/max and very close to zero
    # Picking the right tolerances here is rather sensitive
    bigdr = np.max(np.abs(drval))
    smalldr = np.min(np.abs(drval))
    if smalldr < 0.001 * bigdr:  # There is a significant degree of smallness worth checking
        quasi_roots = scipy.signal.find_peaks(-np.abs(drval))[0] # Find all local maxima
        #print(f'potential extra roots for p={p}, q={q*1e-6:.4f}\n', 
        #      quasi_roots, '\n  ', Om[quasi_roots]/(1e9*twopi), '\n   ',
        #      drval[quasi_roots]/(1e9*twopi))
        for qrt in quasi_roots:  
            # Found all required modes
            if imode == nmodes:  break
 
            # Don't bother with end points
            if qrt in [0,len(Om)-1]: continue
 
            # Bracketed root, so will have been picked up above by brentq
            if (drval[qrt-1] * drval[qrt] < 0) or (drval[qrt] * drval[qrt+1] < 0): continue
 
            droot=drval[qrt]
            if np.abs(droot)<1e-8 and np.abs(droot)<1e-9 * bigdr:
                # found a likely double root, fit a parabola to get accurate location
                a=(drval[qrt+1]+drval[qrt-1]-2*drval[qrt])/2
                b=(drval[qrt+1]-drval[qrt-1])/2
                sol_Om[imode]=Om[qrt]- b/(2*a)*(Om[1]-Om[0]) 
                print('Found an additional double root for p=%d at nu=%.4e.'%(p, 
                                                                 sol_Om[imode]/(1e9*twopi)))
                imode+=1


    return sol_Om

def solve_elastic_rod_analytical(prefix, qvec, nmodes, coremat, arad):
    cmat = coremat.c_tensor
    c11 = cmat[1,1].real
    c12 = cmat[1,2].real
    c44 = cmat[4,4].real
    rho = coremat.rho
    Vl = coremat.Vac_longitudinal()
    Vs = coremat.Vac_shear()
    nu_hi = 1.75*Vs*qvec[-1]/twopi  # Gives a suitable number of modes


    print('Material properties:')
    print(f'c11: {c11/1e9} GPa, c12: {c12/1e9} GPa, c44: {c44/1e9} GPa, Vs {Vs:.2f} m/s, Vl {Vl:.2f} m/s')

    #chareq_elastic_rod_p0_2d(2*pi*3.5e9, 0, 3.1e6, rho, c11, c12, c44, arad)
    #sys.exit(0)

    # The worker function passed to CalcThread to do one task
    def solelrod_caller(args):
        (p, iq, q) = args # Matches queues passed to CalcThread
        v_Om_modes = findroots_elastic_rod_chareq(Vl, nu_hi, p, q, nmodes, rho, c11, c12, c44, arad)
        return (p, iq, v_Om_modes)


    # Two simplified char eqs apply for p=0. We use p=-1 for one of them.
    plo = -1
    phi = 7

    num_cores = os.cpu_count()
    q_work = queue.Queue()        # for assigning the work
    q_result = queue.Queue()      # for returning the answers

    # build all the tasks
    for p in range(plo, phi+1):
        for iq,q in enumerate(qvec):
            q_work.put((p, iq, q))

    print(f'Assigned {q_work.qsize()} tasks to {num_cores} threads.')

    # Make plots of the char equation at q as func of Omega (to check if it behaves nicely for root-finding)
    plot_chareqs=False
    if plot_chareqs:
        for p in range(plo, phi+1):
            figp,axs = plt.subplots(2,1)
            for iq,q in enumerate(qvec):
                if iq%5: continue

                Omlo = 1.e5
                Omvec = np.linspace(Omlo, twopi*2e9, 200)
                drvec = np.zeros(len(Omvec))

                for iOm, Om in enumerate(Omvec):
                    drvec[iOm] = chareq_elastic_rod(Om, p, q, rho, c11, c12, c44, arad)
                axs[0].plot(Omvec/(twopi*1e9), drvec, '.', label=r'$q={0:.2f}/$ μm'.format(q*1e-6))
                axs[1].plot(Omvec/(twopi*1e9), np.log10(np.abs(drvec)), '.', label=r'$q={0:.2f}/$ μm'.format(q*1e-6))
            ysc = {-1: 5e-9, 0:1e-5, 1:2e-3, 2:1e3, 3:.1e-9, 4:1e-0, 5:.1e-9, 6:1e-8, 7:1}[p]
            axs[0].set_ylim(-ysc, ysc) # Now the drvec=1 values representing no sols will be skipped
            for iax in range(2):
                axs[iax].set_xlabel(r'$\Omega/2\pi$  [GHz]')
                axs[iax].axhline(0, color='gray', linestyle='dotted')
                axs[iax].legend()
            figp.savefig('omdr_p%d.png'%p)


    # Run the main dispersion relation calculation
    # No plotting is done in this calc so no plot queue is required.
    launch_worker_threads_and_wait(num_cores, solelrod_caller, q_result, q_work)


    #Collect all the data back into one matrix
    nu_an=np.zeros([phi-plo+1, len(qvec), nmodes], dtype=float)
    while not q_result.empty():  #only one thread now, so no blocking to worry about
        (p, iq, v_Om_modes) = q_result.get()
        nu_an[p+1, iq,:len(v_Om_modes)] = v_Om_modes/twopi

    #for p in range(plo, phi+1):
    #    for iq in range(len(qvec)):
    #        print(p, iq, len(nu_an[p,iq,:]))


    # Plot of the analytic solution alone
    fig,ax = plt.subplots()
    fign,axn = plt.subplots()
    cmap = plt.get_cmap("tab10")
    Vref=Vs

    # broadcast qvec against the middle (q) index of nu_an
    # (numpy broadcasting matches indices or length 1 from the right leftwards)
    neff_an = Vref*qvec[:,np.newaxis]/(nu_an*twopi+1e-5) # small bit to avoid div by zero

    for p in range(plo, phi+1):
        if p==-1: lab = r'$p=0$, Torsional'
        elif p==0:lab = r'$p=0$, Poch.'
        else: lab = r'$p=%d$'%p


        col = cmap(p-plo)

        plot_and_label (ax, qvec*1e-6, nu_an[p-plo,:,:]*1e-9, '.', lab, col)
        plot_and_label (axn, qvec*1e-6, neff_an[p-plo,:,:], '.', lab, col)

        if False:  # individual dispersion plots if desired
        #Analytic dispersion for given p
            figb,axb=plt.subplots()
            plot_and_label (axb, qvec*1e-6, nu_an[p-plo,:,:]*1e-9, '.', lab, col)
            axb.plot(qvec*1e-6, qvec*Vl/(twopi*1e9), ':', color='gray')
            axb.plot(qvec*1e-6, qvec*Vs/(twopi*1e9), ':', color='gray')

            decorate_and_save_plot(figb, axb, prefix + '-acdisp_qnu-p%d.png'%p,
                                   'Dispersion p=%d'%p, r'$q\, [\mathrm{μm}^{-1}]$ ',
                                   r'$\nu$ [GHz]', (0, qvec[-1]*1e-6), (0,nu_hi*1e-9))


    ax.plot(qvec*1e-6, qvec*Vl/(twopi*1e9), ':', color='gray', lw=.5)
    ax.plot(qvec*1e-6, qvec*Vs/(twopi*1e9), ':', color='gray', lw=.5)
    decorate_and_save_plot(fig, ax, prefix + '-acdisp_qnu_exact.png',
                           'Analytic acoustic dispersion', r'$q\, [\mathrm{μm}^{-1}]$ ',
                           r'$\nu$ [GHz]', (0, qvec[-1]*1e-6), (0,nu_hi*1e-9), legend_ncol=2)

    axn.axhline(1, color='gray', linestyle='dotted')
    axn.axhline(Vs/Vl, color='gray', linestyle='dotted')
    decorate_and_save_plot(fign, axn, prefix + '-acdisp_qneff_exact.png', 'Analytic acoustic dispersion',
              r'$q\, [\mathrm{μm}^{-1}]$ ', r'$\bar{n}=v_s/v$',
              (0, qvec[-1]*1e-6), (.4,1.2), legend_ncol=2)

    return nu_an

def solve_elastic_rod_numerical(prefix, qvec, nmodes, wguide, sim_EM, cmat):

    Vl = cmat.Vac_longitudinal()
    Vs = cmat.Vac_shear()
    nu_hi = 1.75*Vs*qvec[-1]/twopi  # Gives a suitable number of modes


    nbasis_AC = min(80, nmodes*2) # 80 ish

    n_field_outs = 2  # how many ksteps to write full set of fields
    field_out_skip = max(int(len(qvec)/n_field_outs), 1)

    # Worker function that can be passed to CalcThread
    def accalc_caller(args):
        (iq, tq, doplot, wg) = args # Matches queues passed to CalcThread

        shift_Hz=0.8*tq*Vs/twopi  # look for modes not to far from the shear sound line
        sim_AC= wg.calc_AC_modes(nbasis_AC, tq, shift_Hz=shift_Hz, EM_sim=sim_EM) #, bcs='Open')
        v_nu_num= np.sort(np.real(sim_AC.nu_AC_all()))[:nmodes]

        if doplot: # Only worker 1 will ever do this
            print('{0} is plotting elastic modes at iq = {1:d} of [0..{2:d}].'.format(
                threading.current_thread().name, iq, len(qvec)-1))
            plotting.plot_mode_fields(sim_AC, ivals=range(nmodes), prefix=prefix+'_%d'%iq)

        return (iq, tq, v_nu_num)


    # Multithread-safe queues for communicating with workers
    q_work_plot = queue.Queue()   # for work requiring plotting
    q_work = queue.Queue()        # everything else
    q_result = queue.Queue()      # for returning the answers

    # Create work queues
    for iq, tq in enumerate(qvec):
        doplot =  (iq % field_out_skip == 0) # Time for some field plots!
        wg = copy.deepcopy(wguide)  # wguide and sim are not thread-safe when we plot mode profiles
        if doplot:
            q_work_plot.put((iq, tq, doplot, wg))
        else:
            q_work.put((iq, tq, doplot, wg))


    num_cores = os.cpu_count()
    launch_worker_threads_and_wait(num_cores, accalc_caller, q_result, q_work, q_work_plot)

    #Collect all the data back into one matrix
    nu_num = np.zeros([len(qvec), nmodes], dtype=float)
    while not q_result.empty():  #only one thread now, so no blocking to worry about
        (iq, tq, nu_modes) = q_result.get()
        nu_num[iq,:] = nu_modes

    fig,ax = plt.subplots()

    plot_and_label(ax, qvec*1e-6, nu_num*1e-9, 'xk', 'FEM')
    ax.plot(qvec*1e-6, qvec*Vl/(twopi*1e9), ':', color='gray')
    ax.plot(qvec*1e-6, qvec*Vs/(twopi*1e9), ':', color='gray')

    # Finish (q,nu) plot
    decorate_and_save_plot(fig, ax, prefix + '-acdisp_qnu_num.png',
          'Elastic dispersion of single rod',
          r'$q\, [\mathrm{μm}^{-1}]$ ', r'$\nu$ [GHz]',
           (0, qvec[-1]*1e-6), (0,nu_hi*1e-9))

    Vref=Vs
    neff_num = Vref*qvec[:,np.newaxis]/(nu_num*twopi+1e-5) # small bit to avoid div by zero

    fign,axn = plt.subplots()
    plot_and_label (axn, qvec*1e-6, neff_num, 'xk', 'FEM')

    axn.axhline(1, color='gray', linestyle='dotted')
    axn.axhline(Vs/Vl, color='gray', linestyle='dotted')
    decorate_and_save_plot(fign, axn, prefix + '-acdisp_qneff_num.png',
              'Elastic dispersion of single rod',
              r'$q\, [\mathrm{μm}^{-1}]$ ', r'$\bar{n}=v_s/v$',
              (0, qvec[-1]*1e-6), (.4,1.2), 'upper right')


    return nu_num


def solve_elastic_dispersion(prefix, ssys, wguide, sim_EM, rcore, mat_core):
    print('\n\nAcoustic dispersion problem')
    print(    '---------------------------')

    Vl = mat_core.Vac_longitudinal()
    Vs = mat_core.Vac_shear()

    qsteps_num=51
    qsteps_analytic=201

    #qsteps_num=3
    #qsteps_analytic=21

    qlo=.0001e6
    qhi=5e6
    qvec_an=np.linspace(qlo, qhi, qsteps_analytic)
    qvec_num=np.linspace(qlo, qhi, qsteps_num)
    nu_hi = 1.75*Vs*qhi/twopi  # Gives a suitable number of modes


    nmodes = 40


    print('Doing analytic problem')
    nu_an = solve_elastic_rod_analytical(prefix, qvec_an, nmodes, mat_core, rcore)

    print('Doing numerical problem')
    nu_num = solve_elastic_rod_numerical(prefix, qvec_num, nmodes, wguide, sim_EM, mat_core)
    #nu_num = np.zeros([len(qvec_num), nmodes], dtype=float)


    #Joint comparison plot
    fig, ax = plt.subplots()
    fign, axn = plt.subplots()
    cmap = plt.get_cmap("tab10")

    Vref=Vs
    neff_num = Vref*qvec_num[:,np.newaxis]/(nu_num*twopi+1e-5)
    neff_an  = Vref*qvec_an[:,np.newaxis]/(nu_an*twopi+1e-5) # avoid div by zero

    for ip in range(nu_an.shape[0]):
        if ip==0: lab = r'$p=0$, Torsional'  #Indices different because of 2 p=0 cases pushing columns to right
        elif ip==1: lab = r'$p=0$, Poch.'
        else: lab = r'$p=%d$'%(ip-1)


        col = cmap(ip)
        plot_and_label(ax, qvec_an*1e-6, nu_an[ip,:,:]*1e-9, '.', lab, col)
        plot_and_label(axn, qvec_an*1e-6, neff_an[ip,:,:], '.', lab, col)

        #write data to text file with qvec as 0th column
        mdata = np.hstack((qvec_an[:,np.newaxis], nu_an[ip,:,:])) 
        np.savetxt(prefix+'-acdisp_qnu_p%d.txt'%(ip-1), mdata, fmt='%.8e',
              header='q [/m]    '+'nu_j [Hz]   '*nmodes)

    # Numerical part of joint (q,nu) plot
    plot_and_label(ax, qvec_num*1e-6, nu_num*1e-9, 'xk', 'FEM')
    ax.plot(qvec_num*1e-6, qvec_num*Vl/(twopi*1e9), ':', color='gray', lw=.5)
    ax.plot(qvec_num*1e-6, qvec_num*Vs/(twopi*1e9), ':', color='gray', lw=.5)

    #write numerical dispersion data with qvec as 0th column
    mdata = np.hstack((qvec_num[:,np.newaxis], nu_num)) 
    np.savetxt(prefix+'-acdisp_qnu_num.txt', mdata, fmt='%.8e', 
              header='q [/m]    '+'nu_j [Hz]   '*nmodes)

    # Finish joint (q,nu) plot
    decorate_and_save_plot(fig, ax, prefix + '-acdisp_qnu.png',
          'Elastic dispersion of single rod',
          r'$q\, [\mathrm{μm}^{-1}]$ ', r'$\nu$ [GHz]',
           (0, qvec_num[-1]*1e-6), (0,nu_hi*1e-9), legend_ncol=2)

    # Joint (nu, neff) plot
    plot_and_label (axn, qvec_num*1e-6, neff_num, 'xk', 'FEM')
    axn.axhline(1, color='gray', linestyle='dotted')
    axn.axhline(Vs/Vl, color='gray', linestyle='dotted')
    decorate_and_save_plot(fign, axn, prefix + '-acdisp_qneff.png',
                           'Acoustic effective index of single rod',
              r'$q\, [\mathrm{μm}^{-1}]$ ', r'$\bar{n}=v_s/v$',
              (0, qvec_num[-1]*1e-6), (.4,1.2), legend_loc='upper left',legend_ncol=2)



def do_main():

    start = time.time()

    pref0, refine_fac = starter.read_args(13, sys.argv, refine=4)

    # Geometric Parameters - all in nm.

    smf28 = True  # which system to study: SMF-28 or silicon rod in silica

    mat_SiO2_GeO2= materials.make_material("SiO2GeO2_smf28")
    mat_SiO2= materials.make_material("SiO2_smf28")
    mat_vac= materials.make_material("Vacuum")
    #mat_Si= materials.make_material("Si_2021_Poulton")
    mat_As2S3= materials.make_material("As2S3_2023_Steel")

    if smf28:
        ssys = 'SMF-28 optical fibre'
        dcore = 8200    # 8 micron core diameter
        dclad = 125000  # 125 micron cladding diameter
        mat_core = mat_SiO2_GeO2
        mat_clad = mat_SiO2
        mat_bkg  = mat_vac
        prefix =pref0+'-em-smf28'
    else:
        ssys = 'chalc rod in silica'
        dcore = 300      # 150 nm core diameter
        dclad = 800      # 500 nm cladding diameter
        mat_core = mat_As2S3
        mat_clad = mat_SiO2
        mat_bkg  = mat_vac
        prefix =pref0+'-em-sil'


    # Elastic problem of a single rod in vacuum
    rcore = 1000
    acore = 2*rcore
    inc_shape = 'circ_onion1'
    unitcell_x = rcore*2.5  # system size in nm
    unitcell_y = unitcell_x
    #mat_PMMA = materials.make_material('PMMA')  #

    mat_core = mat_As2S3
    mat_bkg  = mat_vac

    wguide = objects.Structure(unitcell_x, acore, inc_shape=inc_shape,
                            unitcell_y=unitcell_y, inc_b_x  =rcore*.1,
                            material_bkg=mat_bkg,
                            material_a=mat_core,
                            lc_bkg=.1, lc_refine_1=3.0*refine_fac, lc_refine_2=3*refine_fac)

    prefix =pref0
    #wguide.plot_mesh(prefix)
    sim_EM= wguide.calc_EM_modes(40, 1550,1.5) # solve one EM step to prep the waveguide meshing

    solve_elastic_dispersion(prefix, ssys, wguide, sim_EM, rcore, mat_core)




if __name__ == '__main__':
    do_main()
