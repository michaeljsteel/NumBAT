""" Example showing how the 'onion' geometry template can be used
    to calculate the dispersion profile of standard SMF-28 optical fibre.
"""

import os
import time
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.optimize as sciopt
import scipy.special as sp
import queue 
import copy
import traceback
import threading

sys.path.append("../backend/")
import materials
import objects
import mode_calcs
import integration
import plotting
from fortran import NumBAT
from math import *
from numbattools import launch_worker_threads_and_wait

import starter

# a helper function to make consistent pretty plots
def decorate_and_save_plot (fig, ax, fname, title, xlab, ylab, xlims, ylims):
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_xlim(xlims[0], xlims[1])
    ax.set_ylim(ylims[0], ylims[1])
    ax.set_title(title)
    ax.legend()
    fig.savefig(fname)

def plot_and_label (ax, vx, mat, sty, lab):
    ax.plot(vx, mat, sty)
    ax.plot(vx, mat[:,0], sty, label=lab)

#
#
#Elastic analytic calculations
#
#

#Note that these disp rels continue to work fine when qts, qtl go imaginary
#The imaginary terms always come in product pairs giving a real answer

    #p=0 torsional 
def chareq_elastic_rod_p0_1d(Om, p, q, rho, c11, c12, c44, a_nm): 
    Vs = np.sqrt(c44/rho)

    qts = np.sqrt((Om/Vs)**2-q**2+0j)  


    a = a_nm*1e-9

    Js = sp.jv(p, qts*a)
    Jps = sp.jvp(p, qts*a)
    m22 = qts**2 * Js + 2*qts/a*Jps 

    chareq = m22/Om**2

    #print('p01d {0:.3e} {1:.3e} {2:.3e} {3:.3e} '.format( a, Om, Vs, q), qts, Js, Jps, m22,chareq)
    return np.real(chareq)

    #p=0 Pochammer 
def chareq_elastic_rod_p0_2d(Om, p, q, rho, c11, c12, c44, a_nm): 
    Vl = np.sqrt(c11/rho)
    Vs = np.sqrt(c44/rho)

    qts = np.sqrt((Om/Vs)**2-q**2+0j)
    qtl = np.sqrt((Om/Vl)**2-q**2+0j) 


    a = a_nm*1e-9


    Jl = sp.jv(p, qtl*a)
    Jpl =sp.jvp(p, qtl*a)
    Jppl = 0.5*( sp.jvp(p-1, qtl*a) - sp.jvp(p+1, qtl*a))

    Js = sp.jv(p, qts*a)
    Jps = sp.jvp(p, qts*a)
    Jpps = 0.5*( sp.jvp(p-1, qts*a) - sp.jvp(p+1, qts*a))

    m00 = 2*c44 * qtl**2 * Jppl - c12*(qtl**2 + q**2)*Jl
    m01 = 1j*2*c44 * q * qts * Jpps

    m10 = 1j*2 * q * qtl * Jpl
    m11 = (qts**2-q**2) * Jps 

    chareq = (m00*m11-m01*m10)/(Om**4) # Om**4 makes numbers nice sie

    return chareq.real


def chareq_elastic_rod_ppos(Om, p, q, rho, c11, c12, c44, a_nm): #p is azimuthal order

    Vl = np.sqrt(c11/rho)
    Vs = np.sqrt(c44/rho)

    qts = np.sqrt((Om/Vs)**2-q**2+0j)
    qtl = np.sqrt((Om/Vl)**2-q**2+0j) 

    a = a_nm*1e-9


    Jl = sp.jv(p, qtl*a)
    Jpl = sp.jvp(p, qtl*a)
    Jppl = 0.5*( sp.jvp(p-1, qtl*a) - sp.jvp(p+1, qtl*a))

    Js = sp.jv(p, qts*a)
    Jps = sp.jvp(p, qts*a)
    Jpps = 0.5*( sp.jvp(p-1, qts*a) - sp.jvp(p+1, qts*a))

    m00 = 2*c44 * qtl**2 * Jppl - c12*(qtl**2 + q**2)*Jl
    m01 = 1j*2*c44 * q * qts * Jpps
    m02 = 2*c44*p/a *( qts * Jps - Js/a)  + 2*c12*p/a *qts *Jps

    m10 = 1j*2 * q * qtl * Jpl
    m11 = (qts**2-q**2) * Jps 
    m12 = 1j* q * p * Js/a

    m20 = 2*p/a*(qtl * Jpl-Jl/a)
    m21 = 2*1j*q*p/(qts*a)*(qts * Jps-Js/a)
    m22 = qts**2 * Js + 2*qts/a*Jps 

    chareq = (m00*(m11*m22-m12*m21) + m01*(m12*m20-m10*m22) + m02*(m10*m21-m11*m20))/(Om**5) # Om**5 makes numbers nice sie

    return chareq.real

def chareq_elastic_rod(Om, p, q, rho, c11, c12, c44, a_nm): #p is azimuthal order
    if p==-1:
        return chareq_elastic_rod_p0_1d(Om, p, q, rho, c11, c12, c44, a_nm)
    elif p==0:
        return chareq_elastic_rod_p0_2d(Om, p, q, rho, c11, c12, c44, a_nm)
    else:
        return chareq_elastic_rod_ppos(Om, p, q, rho, c11, c12, c44, a_nm)

def findroots_elastic_rod_chareq(p, q, nmodes, rho, c11, c12, c44, arad):

    nbrack = 500

    Vl = np.sqrt(c11/rho).real
    Vs = np.sqrt(c44/rho).real

    if p==-1:       #ensure transverse wave numbers are real
        Omlo = q*Vs
    else:
        Omlo = q*Vl

    Omhi = Omlo + 25e9
    dOm = Omhi - Omlo

    # look for up to nmodes in [Omlo, Omhi]
    Om = np.linspace(Omlo+dOm/nbrack, Omhi-dOm/nbrack, nbrack)

    sol_Om = np.zeros(nmodes, dtype=np.float64)
    imode = 0

    last_Om = Om[0]
    last_res =chareq_elastic_rod(last_Om, p, q, rho, c11, c12, c44, arad) 

    #nobrack = True # if we find no roots with a given p, there will be no higher p roots and we can give up
    iOm = 1
    while imode < nmodes and iOm < nbrack:
        t_Om = Om[iOm]
        t_res =  chareq_elastic_rod(t_Om, p, q, rho, c11, c12, c44, arad)
        if (last_res * t_res < 0):
            # a bracket! go find the root
    #        nobrack = False
            #print('Bracketed:', p, q, last_Om, t_Om, last_res, t_res)
            root  = sciopt.brentq(chareq_elastic_rod, last_Om, t_Om,
                        args=(p, q, rho, c11, c12, c44, arad))
            sol_Om[imode] = root
            imode += 1
              
        last_Om = t_Om
        last_res  = t_res
        iOm+=1
 
    return sol_Om

def solve_elastic_rod_analytical(prefix, qvec, nmodes, coremat, arad):
    cmat = coremat.c_tensor
    c11 = cmat[1,1].real
    c12 = cmat[1,2].real
    c44 = cmat[4,4].real
    rho = coremat.rho
    Vl = np.sqrt(c11/rho).real
    Vs = np.sqrt(c44/rho).real

    # The worker function passed to CalcThread to do one task
    def solelrod_caller(args):
        (iq, p, q) = args # Matches queues passed to CalcThread
        v_Om_modes = findroots_elastic_rod_chareq(p, q, nmodes, rho, c11, c12, c44, arad)
        return (iq, p, v_Om_modes)


    # Two simplified char eqs apply for p=0. We use p=-1 for one of them. 
    plo = -1
    phi = 6

    q_work = queue.Queue()        # for assigning the work
    q_result = queue.Queue()      # for returning the answers

    # build all the tasks
    for p in range(plo, phi+1):  
        for iq,q in enumerate(qvec):
            q_work.put((iq, p, q))

    plot_chareqs=True
    # Make plots of the char equation at q as func of Omega (to check if it behaves)
    if plot_chareqs:
        for p in range(plo, phi+1):
            figp,axp = plt.subplots()   
            for iq,q in enumerate(qvec):
   
                if p == -1:
                    Omlo = Vs*q
                else:
                    Omlo = Vl*q
                Omlo = 1.e5
                Omvec = np.linspace(Omlo, 2*pi*10e9, 1000)
                drvec = np.zeros(len(Omvec))
    
                for iOm, Om in enumerate(Omvec):
                    #drvec[iOm] = chareq_elastic_rod(Om, p, q, rho, c11, c12, c44, arad)
                    drvec[iOm] = log10(abs(chareq_elastic_rod(Om, p, q, rho, c11, c12, c44, arad)))
                if iq%20==0:
                    axp.plot(Omvec/(2*pi*1e9), drvec, '.', label=r'$q={0:.2f}/$ μm'.format(q*1e-6))
                else:
                    axp.plot(Omvec/(2*pi*1e9), drvec, '.')
            ymin=abs(np.min(drvec)) # find largest negative value and take that as representative range
            #axp.set_ylim(-1.5*ymin, 5*ymin) # Now the drvec=1 values representing no sols will be skipped
            axp.set_ylim(-2e-7,2e-7)
            axp.set_xlabel(r'$\Omega/2\pi$  [GHz]')
            axp.legend()
            figp.savefig('omdr_p%d.png'%p)

        
    num_cores = os.cpu_count()
    num_cores = 1
    launch_worker_threads_and_wait(num_cores, solelrod_caller, q_result, q_work) # no plot queue in this case

    nu_modes_an=np.zeros([phi-plo+1, len(qvec), nmodes], dtype=float)

    #Collect all the data back into one matrix
    while not q_result.empty():  #only one thread now, so no blocking to worry about
        (iq, p, v_Om_modes) = q_result.get()
        nu_modes_an[p+1, iq,:] = v_Om_modes/(2*pi)

    
    # Plot of the analytic solution alone
    fig,ax = plt.subplots()
    for p in range(plo, phi+1): 
        if p==-1: lab = r'$p=0$, Torsional'
        elif p==0:lab = r'$p=0$, Poch.'
        else: lab = r'$p=%d$'%p
        plot_and_label (ax, qvec*1e-6, nu_modes_an[p,:,:]*1e-9, '.', lab) 

    decorate_and_save_plot(fig, ax, prefix + '-disp_qnu-exact.png', 'Analytic acoustic dispersion',
              r'$q\, [\mathrm{μm}^{-1}]$ ', r'$\nu$ [GHz]',
              (0, qvec[-1]*1e-6), (0,10))

    return nu_modes_an

def solve_elastic_rod_numerical(prefix, qvec_num, nmodes, wguide, sim_EM):

    nu_num = np.zeros([len(qvec_num), nmodes], dtype=float)

    nbasis_AC = nmodes* 2 # 80 ish

    n_field_outs = 2  # how many ksteps to write full set of fields
    field_out_skip = max(int(len(qvec_num)/n_field_outs), 1)

    # Worker function that can be passed to CalcThread
    def accalc_caller(args):
        (iq, tq, doplot, wg) = args # Matches queues passed to CalcThread

        sim_AC= wg.calc_AC_modes(nbasis_AC, tq, shift_Hz=5e9, EM_sim=sim_EM) #, bcs='Open')
        nu_num= np.sort(np.real(sim_AC.nu_AC_all()))[:nmodes]

        if doplot: # Only worker 1 will ever do this
            print('{0} is plotting elastic modes at iq = {1:d} of [0..{2:d}].'.format(
                threading.current_thread().name, iq, len(qvec_num)-1))
            plotting.plot_mode_fields(sim_AC, ivals=range(nmodes), prefix=prefix+'_%d'%iq)

        return (iq, tq, nu_num)


    # Multithread-safe queues for communicating with workers
    q_work_plot = queue.Queue()   # for work requiring plotting
    q_work = queue.Queue()        # everything else
    q_result = queue.Queue()      # for returning the answers

    # Create work queues
    for iq, tq in enumerate(qvec_num):
        doplot =  (iq % field_out_skip == 0) # Time for some field plots!
        wg = copy.deepcopy(wguide)  # wguide and sim are not thread-safe when we plot mode profiles
        if doplot:
            q_work_plot.put((iq, tq, doplot, wg))
        else:
            q_work.put((iq, tq, doplot, wg))
    

    num_cores = os.cpu_count()
    launch_worker_threads_and_wait(num_cores, accalc_caller, q_result, q_work, q_work_plot)

    #Collect all the data back into one matrix
    while not q_result.empty():  #only one thread now, so no blocking to worry about
        (iq, tq, nu_modes) = q_result.get()
        nu_num[iq,:] = nu_modes

    fig,ax = plt.subplots()
    plot_and_label(ax, qvec_num*1e-6, nu_num*1e-9, 'xk', 'FEM')
    decorate_and_save_plot(fig, ax, prefix + '-disp_qnu-num.png', 
              'Elastic dispersion of single rod',
              r'$q\, [\mathrm{μm}^{-1}]$ ', r'$\nu$ [GHz]',
              (0, qvec_num[-1]*1e-6), (0,10))


    return nu_num


def solve_elastic_dispersion(prefix, ssys, wguide, sim_EM, rcore, mat_core):
    print('\n\nAcoustic dispersion problem')
    print(    '---------------------------')


    qsteps_num=10
    qsteps_analytic=100
    
    qlo=.1e6    
    qhi=1e7
    qvec_an=np.linspace(qlo, qhi, qsteps_analytic)
    qvec_num=np.linspace(qlo, qhi, qsteps_num)


    nmodes = 40


    print('Doing analytic problem')
    nu_an = solve_elastic_rod_analytical(prefix, qvec_an, nmodes, mat_core, rcore)

    sys.exit(0)

    print('Doing numerical problem')
    nu_num = solve_elastic_rod_numerical(prefix, qvec_num, nmodes, wguide, sim_EM)


    #Joint comparison plot
    fig, ax = plt.subplots()
    for ip in range(nu_an.shape[0]):
        if ip==0: lab = r'$p=0$, Torsional'  #Indices different because of 2 p=0 cases pushing columns to right
        elif ip==1: lab = r'$p=0$, Poch.'
        else: lab = r'$p=%d$'%(ip-1)
        plot_and_label(ax, qvec_an*1e-6, nu_an[ip,:,:]*1e-9, '.', lab)

    plot_and_label(ax, qvec_num*1e-6, nu_num*1e-9, 'xk', 'FEM')
    
    decorate_and_save_plot(fig, ax, prefix + '-disp_qnu.png', 
          'Numerical elastic dispersion of single rod',
          r'$q\, [\mathrm{μm}^{-1}]$ ', r'$\nu$ [GHz]',
           (0, qvec_num[-1]*1e-6), (0,10))


#
#
#EM analytic calculations
#
#


def chareq_em_fib2_TE_m(neff, m, k, rco, ncore, nclad):
    # m is un-used (TE is m=0)
    u = rco*k*np.sqrt(ncore**2-neff**2)
    w = rco*k*np.sqrt(neff**2-nclad**2)

    return sp.jv(1,u)/(sp.jv(0,u)*u) + sp.kv(1,w)/(sp.kv(0,w)*w)   

def chareq_em_fib2_TM_m(neff, m, k, rco, ncore, nclad):
    # m is un-used (TM is m=0)
    u = rco*k*np.sqrt(ncore**2-neff**2)
    w = rco*k*np.sqrt(neff**2-nclad**2)

    n2rat = (nclad/ncore)**2

    fac1 = np.real(sp.jv(1,u)/(sp.jv(0,u)*u) )
    fac2 = np.real(sp.kv(1,w)/(sp.kv(0,w)*w) ) 

    return fac1 + n2rat * fac2

def chareq_em_fib2_hy_m(neff, m, k, rco, ncore, nclad):
    u = rco*k*np.sqrt(ncore**2-neff**2)
    w = rco*k*np.sqrt(neff**2-nclad**2)

    n2rat = (nclad/ncore)**2

    invu2 = 1.0/u**2
    invw2 = 1.0/w**2
    jrat =  sp.jvp(m,u)/(sp.jv(m,u)*u) 
    krat =  sp.kvp(m,w)/(sp.kv(m,w)*w)
    fac1 =  jrat + krat
    fac2 =  jrat + n2rat * krat
    fac3 = m*m*(invu2 + invw2)*(invu2 + n2rat*invw2)
    return fac1*fac2-fac3

    
def plot_em_chareq_at_k(k, rcore, ncore, nclad):
    nbrack = 500
    dn = ncore-nclad
    v_neff = np.linspace(ncore-dn/nbrack, nclad+dn/nbrack, nbrack)
    v_dr_TE = np.zeros(nbrack, float)
    v_dr_TM = np.zeros(nbrack, float)
    v_dr_hy = np.zeros(nbrack, float)

    fig, axs = plt.subplots(1,3, figsize=(20,6))
    for m in range(5):
        for i,neff in enumerate(v_neff):
            v_dr_TE[i] = chareq_em_fib2_TE_m(neff, m, k, rcore, ncore, nclad)+.1
            v_dr_TM[i] = chareq_em_fib2_TM_m(neff, m, k, rcore, ncore, nclad)
            v_dr_hy[i] = chareq_em_fib2_hy_m(neff, m, k, rcore, ncore, nclad)
        axs[0].plot(v_neff, v_dr_TE)
        axs[1].plot(v_neff, v_dr_TM)
        axs[2].plot(v_neff, v_dr_hy)
    for i in range(3):
        axs[i].set_ylim(-20,20)
    fig.savefig('tut_12-em_chareq.png')

def solve_chareq_em_fib2_disprel(f_disprel, family, k, nmodes, mlo, mhi, rco, ncore, nclad):
    # solves for modes with azimuthal order in [mlo, mhi] inclusive

    nbrack = 500
    dn = ncore-nclad

    # look for up to nmodes in [nclad, ncore]
    neff = np.linspace(ncore-dn/nbrack, nclad+dn/nbrack, nbrack)
    sol_neff = np.zeros(nmodes, dtype=np.float64)
    imode = 0
    for m in range(mlo, mhi+1):

        last_neff = neff[0]
        last_res = f_disprel(last_neff, m, k, rco, ncore, nclad)
 
        nobrack = True # if we find no roots with a given m, there will be no higher m roots and we can give up
        ineff = 1
        while imode < nmodes and ineff < nbrack: 
            t_neff = neff[ineff]
            t_res = f_disprel(t_neff, m, k, rco, ncore, nclad)
 
            if ((family == 'Hy' and last_res * t_res < 0) # Hybrid eig curves are smooth
                 or (last_res<0 and t_res>0)): # TE and TM curves have +inf to -inf breaks which are not brackets
                # a bracket! go find the root
                nobrack = False
                #root, rootres 
                root  = sciopt.brentq(f_disprel, last_neff, t_neff,
                                     args=(m, k, rco, ncore, nclad))
                sol_neff[imode] = root
                imode += 1
             
            last_neff = t_neff
            last_res  = t_res
            ineff+=1
 
        if nobrack: break # exhausted all roots at this k


    return sol_neff
    

def solve_analytic_em_two_layer_fiber(kvec, nmodes, rco, ncore, nclad):

    fake = False

    neff_TE_an=np.zeros([len(kvec), nmodes], dtype=float)
    neff_TM_an=np.zeros([len(kvec), nmodes], dtype=float)
    an_neff_hy=np.zeros([len(kvec), nmodes], dtype=float)

    # for skipping over
    if fake:
        neff_TE_an += (ncore+nclad)/2
        neff_TM_an += (ncore+nclad)/2
        an_neff_hy += (ncore+nclad)/2
        return neff_TE_an, neff_TM_an, an_neff_hy

    v_neff_TE = np.zeros(nmodes, dtype=float)
    v_neff_TM = np.zeros(nmodes, dtype=float)
    v_neff_hy = np.zeros(nmodes, dtype=float)
    v_neff = np.zeros(nmodes, dtype=float)

    for ik, k in enumerate(kvec):
        mhy_lo = 1
        mhy_hi = 5
        v_neff_TE = solve_chareq_em_fib2_disprel(chareq_em_fib2_TE_m, 'TE', k, nmodes, 0, 0, rco, ncore, nclad)
        v_neff_TM = solve_chareq_em_fib2_disprel(chareq_em_fib2_TM_m, 'TM', k, nmodes, 0, 0, rco, ncore, nclad)
        v_neff_hy = solve_chareq_em_fib2_disprel(chareq_em_fib2_hy_m, 'Hy', k, nmodes, mhy_lo, mhy_hi, rco, ncore, nclad) # m in [1,5)
        neff_TE_an[ik,:] = np.sort(v_neff_TE)[-nmodes:]
        neff_TM_an[ik,:] = np.sort(v_neff_TM)[-nmodes:]
        an_neff_hy[ik,:] = np.sort(v_neff_hy)[-nmodes:]

    return neff_TE_an, neff_TM_an, an_neff_hy

    
def make_em_plots(prefix, ssys, Vvec_an, neff_TE_an, neff_TM_an, neff_Hy_an, 
                  Vvec_num, neff_num, ncore, nclad):

    Vlo, Vhi =Vvec_an[0], Vvec_an[-1]  

    fig, ax = plt.subplots()
    plot_and_label(ax, Vvec_an, neff_TE_an, '-b', 'TE')
    plot_and_label(ax, Vvec_an, neff_TM_an, '-g', 'TM')
    plot_and_label(ax, Vvec_an, neff_Hy_an, '-r', 'Hy')
    decorate_and_save_plot(fig, ax, prefix + '-em_disp_Vneff_exact.png',
            'Analytic dispersion of '+ssys,
            r'$V$ number ', r'$\bar{n}$', (Vlo, Vhi), (nclad, ncore))
      
    fig, ax = plt.subplots()
    ax.plot(Vvec_num, neff_num, 'xk')
    decorate_and_save_plot(fig, ax, prefix + '-em_disp_Vneff_num.png',
            'Numerical dispersion of '+ssys,
            r'$V$ number ', r'$\bar{n}$', (Vlo, Vhi), (nclad, ncore))
    
    
    fig, ax = plt.subplots()
    plot_and_label(ax, Vvec_an, neff_TE_an, '-b', 'TE')
    plot_and_label(ax, Vvec_an, neff_TM_an, '-g', 'TM')
    plot_and_label(ax, Vvec_an, neff_Hy_an, '-r', 'Hy')
    ax.plot(Vvec_num, neff_num, 'xk')
    decorate_and_save_plot(fig, ax, prefix + '-em_disp_Vneff.png',
            'Dispersion of '+ssys,
            r'$V$ number ', r'$\bar{n}$', (Vlo, Vhi), (nclad, ncore))

    # Plots of normalised waveguide parameter b
    b_TE_an = (neff_TE_an**2-nclad**2)/(ncore**2-nclad**2)
    b_TM_an = (neff_TM_an**2-nclad**2)/(ncore**2-nclad**2)
    b_Hy_an = (neff_Hy_an**2-nclad**2)/(ncore**2-nclad**2)
    b_num =  (neff_num**2-nclad**2)/(ncore**2-nclad**2)

    fig, ax = plt.subplots()
    plot_and_label(ax, Vvec_an, b_TE_an, '-b', 'TE')
    plot_and_label(ax, Vvec_an, b_TM_an, '-g', 'TM')
    plot_and_label(ax, Vvec_an, b_Hy_an, '-r', 'Hy')
    ax.plot(Vvec_num, b_num, 'xk')
    decorate_and_save_plot(fig, ax, prefix + '-em_disp_Vb.png',
            'Dispersion of '+ssys, r'$V$ number ', r'$b$', (Vlo, Vhi), (0, 1))

    # Find group index by numerical differentiation
    dV=Vvec_an[1]-Vvec_an[0]
    dndV_TE_an = (neff_TE_an[2:,:]-neff_TE_an[:-2,:])/(2*dV)
    dndV_TM_an = (neff_TM_an[2:,:]-neff_TM_an[:-2,:])/(2*dV)
    dndV_Hy_an = (neff_Hy_an[2:,:]-neff_Hy_an[:-2,:])/(2*dV)

    dV=Vvec_num[1]-Vvec_num[0]
    dndV_num = (neff_num[2:,:]-neff_num[:-2,:])/(2*dV)

    ng_TE_an = dndV_TE_an + neff_TE_an[1:-1,:]
    ng_TM_an = dndV_TM_an + neff_TM_an[1:-1,:]
    ng_Hy_an = dndV_Hy_an + neff_Hy_an[1:-1,:]
    ng_num = dndV_num + neff_num[1:-1,:]

    ng_num *= neff_num[1:-1,:]>nclad   # only plot for V>V_cutoff

    fig, ax = plt.subplots()
    plot_and_label(ax, Vvec_an[1:-1], ng_TE_an, '.b', 'TE')
    plot_and_label(ax, Vvec_an[1:-1], ng_TM_an, '.g', 'TM')
    plot_and_label(ax, Vvec_an[1:-1], ng_Hy_an, '.r', 'Hy')
    ax.plot(Vvec_num[1:-1], ng_num, 'xk') 

    decorate_and_save_plot(fig, ax, prefix + '-em_disp_ng.png',
            'Dispersion of '+ssys,
            r'$V$ number ', r'Group index $n_g$', (Vlo, Vhi), (nclad, ncore))



def solve_em_dispersion(prefix, ssys, wguide, rcore, ncore, nclad):
    print('\n\nElectromagnetic dispersion problem')
    print(    '----------------------------------')
    nmodes = 20
    nbasis_EM = 40
    EM_ival_pump = range(nmodes)

    ksteps_num=40
    ksteps_analytic=200
    
    lam_hi=2000    # wavelength range in nm
    lam_lo=500    # wavelength range in nm
    klo=2*pi/lam_hi  # this is free space wavenumber in inverse nm
    khi=2*pi/lam_lo
    kvec_an=np.linspace(klo, khi, ksteps_analytic)
    kvec_num=np.linspace(klo, khi, ksteps_num)

    cvac = 299_792_458*1e9  # speed of light in nm/s
    Vvec_an = kvec_an*rcore*np.sqrt(ncore**2-nclad**2)
    Vvec_num = kvec_num*rcore*np.sqrt(ncore**2-nclad**2)
    Vlo = Vvec_an[0]
    Vhi = Vvec_an[-1]

    
    
    print('Doing analytic problem')
    
    #plot_em_chareq_at_k(kvec[-1], rcore, ncore, nclad) # for checking zero finding behaviour

    neff_TE_an, neff_TM_an, neff_Hy_an = solve_analytic_em_two_layer_fiber(
            kvec_an, nmodes, rcore, ncore, nclad)

    
    print('Doing numerical problem')

    neff_num = np.zeros([len(kvec_num), nmodes], dtype=float)

    # Expected effective index of fundamental guided mode.
    n_eff = (ncore+nclad)/2.0
    n_field_outs = 5  # how many ksteps to write full set of fields
    field_out_skip = max(int(len(kvec_num)/n_field_outs),1)

    for ik, tk in enumerate(kvec_num):
    
      t_lambda_nm=2*pi/tk
      print('\n Doing wavelength {0:d} of {1:d}: λ = {2:.2f} nm, k = 2π/λ = {3:.3f}/um'.format(
          ik+1, len(kvec_num), t_lambda_nm, tk*1000))
    
      sim_EM= wguide.calc_EM_modes(nbasis_EM, t_lambda_nm, n_eff)
      neff_num[ik]= np.sort(np.real(sim_EM.neff_all()))[-nmodes:]

      if ik % field_out_skip == 0: # Time for some field plots!
          plotting.plot_mode_fields(sim_EM, EM_AC='EM_E', ivals=range(nmodes),
          prefix=prefix+'_%d'%ik, ticks=True)
    

    make_em_plots(prefix, ssys, Vvec_an, neff_TE_an, neff_TM_an, neff_Hy_an, 
                  Vvec_num, neff_num, ncore, nclad)
    return sim_EM

    
def do_main():
    
    start = time.time()

    pref0, refine_fac = starter.read_args(12, sys.argv)
    
    # Geometric Parameters - all in nm.

    smf28 = True  # which system to study: SMF-28 or silicon rod in silica

    mat_SiO2_GeO2= materials.make_material("SiO2GeO2_smf28")
    mat_SiO2= materials.make_material("SiO2_smf28")
    mat_vac= materials.make_material("Vacuum")
    mat_Si= materials.make_material("Si_2021_Poulton")

    if smf28:
        ssys = 'SMF-28 optical fibre'
        dcore = 8200    # 8 micron core diameter
        dclad = 125000  # 125 micron cladding diameter
        mat_core = mat_SiO2_GeO2
        mat_clad = mat_SiO2
        mat_bkg  = mat_vac
        prefix =pref0+'-em-smf28'
    else:
        ssys = 'silicon rod in silica'
        dcore = 300      # 150 nm core diameter
        dclad = 800      # 500 nm cladding diameter
        mat_core = mat_Si
        mat_clad = mat_SiO2
        mat_bkg  = mat_vac
        prefix =pref0+'-em-sil'


    rcore = dcore/2.0    
    rclad = dclad/2.0    
    ncore = np.real(mat_core.refindex_n)
    nclad = np.real(mat_clad.refindex_n)
    acore = rcore*2           #annular thickness of each layer (first is diam)
    aclad = rclad - rcore     #annular thickness of each layer

    onion=1
    if onion == 1: 
        inc_shape = 'circ_onion1'
        unitcell_x = rcore*10  # system size in nm
        unitcell_y = unitcell_x
        mat_bkg = mat_clad  # because the background is now the second layer, ie the cladding
    elif onion == 2: 
        inc_shape = 'onion2'
        unitcell_x = rcore*10  # system size in nm
        unitcell_y = unitcell_x
    else:
        #inc_shape = 'circ_onion'
        inc_shape = 'onion'
        unitcell_x = rcore*35  # system size in nm
        unitcell_y = unitcell_x


    wguide = objects.Structure(unitcell_x, acore, inc_shape=inc_shape,  # remove these factors of 2
                               inc_b_x=aclad,
                            unitcell_y=unitcell_y,
                            material_bkg=mat_bkg, 
                            material_a=mat_core, 
                            material_b=mat_clad,
                            lc_bkg=.1, lc_refine_1=4.0*refine_fac, lc_refine_2=4*refine_fac)
    
    #wguide.plot_mesh(prefix)

    #sim_EM = solve_em_dispersion(prefix, ssys, wguide, rcore, ncore, nclad)


    # Elastic problem of a single rod in vacuum
    rcore = 1000
    acore = 2*rcore
    inc_shape = 'circ_onion1'
    unitcell_x = rcore*2.5  # system size in nm
    unitcell_y = unitcell_x
    mat_PMMA = materials.make_material('PMMA')  # 

    mat_core = mat_Si
    mat_bkg  = mat_vac
    
    wguide = objects.Structure(unitcell_x, acore, inc_shape=inc_shape,  
                            unitcell_y=unitcell_y, inc_b_x  =rcore*.1,
                            material_bkg=mat_bkg,
                            material_a=mat_core, 
                            lc_bkg=.1, lc_refine_1=4.0*refine_fac, lc_refine_2=4*refine_fac)

    prefix =pref0+'-ac'
    #wguide.plot_mesh(prefix)
    #sim_EM= wguide.calc_EM_modes(40, 1550,1.5) # solve one EM step to prop the waveguide meshing
    #plotting.plot_mode_fields(sim_EM, EM_AC='EM_E', ivals=range(5), prefix=prefix)
    sim_EM=None

    solve_elastic_dispersion(prefix, ssys, wguide, sim_EM, rcore, mat_core)



    
if __name__ == '__main__': 
    do_main()
