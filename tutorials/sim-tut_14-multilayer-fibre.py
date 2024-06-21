""" Example showing how the 'onion' geometry template can be used
    to calculate the dispersion profile of a multilayer Bragg fibre
    as described in C. M. de Sterke et al, "Differential losses in Bragg fibers",
    J. Appl. Phys. 76, 680 (1994).
"""

import os
import time
import sys
import copy
import math
import multiprocessing
import scipy.optimize as sciopt
import scipy.special as sp

import matplotlib.pyplot as plt
import numpy as np


from pathlib import Path
sys.path.append(str(Path('../backend')))


import numbat
import materials
from numbattools import launch_worker_processes_and_wait

import starter

twopi = math.pi*2.0

# a helper function to make consistent pretty plots
def decorate_and_save_plot (fig, ax, fname, title, xlab, ylab, xlims = None, ylims = None,
                            legloc=None):
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if xlims is not None: ax.set_xlim(xlims[0], xlims[1])
    if ylims is not None: ax.set_ylim(ylims[0], ylims[1])
    ax.set_title(title)
    if legloc is not None:
        ax.legend(loc=legloc, fontsize='x-small')
    else:
        ax.legend(fontsize='x-small')
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


# Solve the analytic dispersion relation using worker threads
def solve_em_multilayer_fiber_analytic(kvec, nmodes, rco, ncore, nclad):

    # The worker function passed to CalcThread to do one task
    def solemrod_caller(args):
        (ik, k) = args # Matches queues passed to CalcThread
        mhy_lo = 1
        mhy_hi = 5
        v_neff_TE = solve_chareq_em_fib2_disprel(chareq_em_fib2_TE_m, 'TE', k, nmodes, 0, 0, rco, ncore, nclad)
        v_neff_TM = solve_chareq_em_fib2_disprel(chareq_em_fib2_TM_m, 'TM', k, nmodes, 0, 0, rco, ncore, nclad)
        v_neff_hy = solve_chareq_em_fib2_disprel(chareq_em_fib2_hy_m, 'Hy', k, nmodes, mhy_lo, mhy_hi, rco, ncore, nclad) # m in [1,5)
        return (ik, v_neff_TE, v_neff_TM, v_neff_hy)


    num_workers = os.cpu_count()
    manager = multiprocessing.Manager()
    q_work = multiprocessing.JoinableQueue()      # for assigning the work
    q_result = manager.Queue()                    # for collecting the work

    # build all the tasks
    for ik,k in enumerate(kvec): q_work.put((ik, k))

    launch_worker_processes_and_wait(num_workers, solemrod_caller, q_result, q_work)

    #Collect all the data back into one matrix per polarisation
    neff_TE_an=np.zeros([len(kvec), nmodes], dtype=float)
    neff_TM_an=np.zeros([len(kvec), nmodes], dtype=float)
    an_neff_hy=np.zeros([len(kvec), nmodes], dtype=float)
    while not q_result.empty():  #only one thread now, so no blocking to worry about
        (ik, v_neff_TE, v_neff_TM, v_neff_hy) = q_result.get()
        neff_TE_an[ik,:] = np.sort(v_neff_TE)[-nmodes:]
        neff_TM_an[ik,:] = np.sort(v_neff_TM)[-nmodes:]
        an_neff_hy[ik,:] = np.sort(v_neff_hy)[-nmodes:]

    return neff_TE_an, neff_TM_an, an_neff_hy

def solve_em_multilayer_fiber_numerical(prefix, wguide, kvec, nmodes, nbasis, rcore, ncore, nclad):

    Vvec= kvec*rcore*np.sqrt(np.real(ncore**2-nclad**2))
    print('\n\nKvec and Vnumber:', kvec, Vvec, '\n\n')
    neff_num = np.zeros([len(kvec), nmodes], dtype=float)

    # Expected effective index of fundamental guided mode.
    n_eff = (ncore+nclad)/2.0
    n_field_outs = 3  # how many ksteps to write full set of fields
    field_out_skip = max(int(len(kvec)/n_field_outs),1)

    # Multiprocess-safe queues for communicating with workers
    # The Manager.Queue is needed for technical reasons
    manager = multiprocessing.Manager()
    q_work = multiprocessing.JoinableQueue()      # for assigning the work
    q_result = manager.Queue()                    # for collecting the work

    # Create work queues
    for ik, tk in enumerate(kvec):
        doplot =  (ik % field_out_skip == 0) # Time for some field plots!
        wg = copy.deepcopy(wguide)  # wguide and sim are not thread-safe when we plot mode profiles
        q_work.put((ik, tk, doplot, wg))


    def emcalc_caller(args):
        (ik, tk, doplot, wg) = args # Matches queues passed to launch_worker...

        t_lambda_nm = twopi/tk

        print('Starting mode solve', ik)
        simres_EM= wguide.calc_EM_modes(nbasis, t_lambda_nm, n_eff)
        print('Done mode solve', ik)
        neff_k= np.sort(np.real(simres_EM.neff_all()))[-nmodes:]

        if doplot: # Only worker 1 will ever do this
            print('Doing plot of modes', ik)
            simres_EM.plot_modes(ivals=range(nmodes), prefix=f'{prefix}_{ik}', ticks=True)
            print('Done plot of modes', ik)

        return (ik, tk, neff_k)

    #num_workers = os.cpu_count()
    #num_workers = 2
    num_workers = 0   # just do direct in the main process and single thread

    launch_worker_processes_and_wait(num_workers, emcalc_caller, q_result, q_work)


    #Collect all the data back into one matrix
    while not q_result.empty():  #only one thread now, so no blocking to worry about
        (ik, tk, neff_k) = q_result.get()
        neff_num[ik,:] = neff_k

    #fig, ax = plt.subplots()
    #ax.plot(Vvec, neff_num, 'xk')
    #decorate_and_save_plot(fig, ax, prefix + '-emdisp_Vneff_num.png',
    #        'Numerical dispersion of multilayer fibre',
    #        r'$V$ number ', r'$\bar{n}$', (Vvec[0], Vvec[-1]), (nclad, ncore))

    print('V and ns', Vvec, neff_num)

    return (Vvec, neff_num)


def make_em_plots(prefix, Vvec_num, neff_num, ncore, nclad):

    Vlo, Vhi =Vvec_num[0], Vvec_num[-1]

    #fig, ax = plt.subplots()
    #plot_and_label(ax, Vvec_an, neff_TE_an, '-b', 'TE')
    #plot_and_label(ax, Vvec_an, neff_TM_an, '-g', 'TM')
    #plot_and_label(ax, Vvec_an, neff_Hy_an, '-r', 'Hy')
    #decorate_and_save_plot(fig, ax, prefix + '-emdisp_Vneff_exact.png',
    #        'Analytic dispersion of '+ssys,
    #        r'$V$ number ', r'$\bar{n}$', (Vlo, Vhi), (nclad, ncore))

    fig, ax = plt.subplots()
    ax.plot(Vvec_num, neff_num, 'xk')
    decorate_and_save_plot(fig, ax, prefix + '-emdisp_Vneff.png',
            'Numerical dispersion of multilayer',
            r'$V$ number ', r'$\bar{n}$', (Vlo, Vhi)) # , (nclad, ncore))


    #fig, ax = plt.subplots()
    #plot_and_label(ax, Vvec_an, neff_TE_an, '-b', 'TE')
    #plot_and_label(ax, Vvec_an, neff_TM_an, '-g', 'TM')
    #plot_and_label(ax, Vvec_an, neff_Hy_an, '-r', 'Hy')
    #ax.plot(Vvec_num, neff_num, 'xk')
    #decorate_and_save_plot(fig, ax, prefix + '-emdisp_Vneff.png',
    #        'Dispersion of '+ssys,
    #        r'$V$ number ', r'$\bar{n}$', (Vlo, Vhi), (nclad, ncore))

    # Plots of normalised waveguide parameter b
    #b_TE_an = (neff_TE_an**2-nclad**2)/(ncore**2-nclad**2)
    #b_TM_an = (neff_TM_an**2-nclad**2)/(ncore**2-nclad**2)
    #b_Hy_an = (neff_Hy_an**2-nclad**2)/(ncore**2-nclad**2)
    #b_num =  (neff_num**2-nclad**2)/(ncore**2-nclad**2)

    #fig, ax = plt.subplots()
    #plot_and_label(ax, Vvec_an, b_TE_an, '-b', 'TE')
    #plot_and_label(ax, Vvec_an, b_TM_an, '-g', 'TM')
    #plot_and_label(ax, Vvec_an, b_Hy_an, '-r', 'Hy')
    #ax.plot(Vvec_num, b_num, 'xk')
    #decorate_and_save_plot(fig, ax, prefix + '-emdisp_Vb.png',
    #        'Dispersion of multilayer fiber', r'$V$ number ', r'$b$', (Vlo, Vhi), (0, 1))

    # Find group index by numerical differentiation
    #dV=Vvec_an[1]-Vvec_an[0]
    #dndV_TE_an = (neff_TE_an[2:,:]-neff_TE_an[:-2,:])/(2*dV)
    #dndV_TM_an = (neff_TM_an[2:,:]-neff_TM_an[:-2,:])/(2*dV)
    #dndV_Hy_an = (neff_Hy_an[2:,:]-neff_Hy_an[:-2,:])/(2*dV)

    dV=Vvec_num[1]-Vvec_num[0]
    dndV_num = (neff_num[2:,:]-neff_num[:-2,:])/(2*dV)

    #ng_TE_an = dndV_TE_an + neff_TE_an[1:-1,:]
    #ng_TM_an = dndV_TM_an + neff_TM_an[1:-1,:]
    #ng_Hy_an = dndV_Hy_an + neff_Hy_an[1:-1,:]
    ng_num = dndV_num + neff_num[1:-1,:]

    ng_num *= neff_num[1:-1,:]>nclad   # only plot for V>V_cutoff

    fig, ax = plt.subplots()
    #plot_and_label(ax, Vvec_an[1:-1], ng_TE_an, '.b', 'TE')
    #plot_and_label(ax, Vvec_an[1:-1], ng_TM_an, '.g', 'TM')
    #plot_and_label(ax, Vvec_an[1:-1], ng_Hy_an, '.r', 'Hy')
    ax.plot(Vvec_num[1:-1], ng_num, 'xk')

    decorate_and_save_plot(fig, ax, prefix + '-emdisp_ng.png',
            'Dispersion of multilayer fibre',
            r'$V$ number ', r'Group index $n_g$', (Vlo, Vhi), (nclad, ncore))



def solve_em_dispersion(prefix, wguide, ncore, nclad, rcore):
    print('\n\nElectromagnetic dispersion problem')
    print(    '----------------------------------')

    nmodes = 20
    nbasis = 40

    ksteps_num=31

    lam_hi=5000    # wavelength range in nm
    lam_lo=500    # wavelength range in nm
    klo=twopi/lam_hi  # this is free space wavenumber in inverse nm
    khi=twopi/lam_lo
    kvec=np.linspace(klo, khi, ksteps_num)

    (Vvec_num, neff_num) = solve_em_multilayer_fiber_numerical(prefix, wguide, kvec, nmodes, nbasis, rcore, ncore, nclad)

    make_em_plots(prefix, Vvec_num, neff_num, ncore, nclad)

def do_main():

    start = time.time()

    prefix, refine_fac = starter.read_args(14, sys.argv, refine=4)
    prefix = 'tt'
    refine_fac = 2

    nbapp=numbat.NumBATApp(prefix)

    # Geometric Parameters - all in nm.



    # TODOL This doesn't work because material tensors aren't reconstructed after initial setup
    # set_refractive_index() needs to trigger a rebuild
#    mat_a = materials.make_material("SiO2_smf28")
#    mat_b = copy.deepcopy(mat_a)
#    mat_vac= materials.make_material("Vacuum")
#
#    deltan = 0.05
#    ncore = 1.50 + deltan
#    nclad = 1.50 - deltan
#    mat_a.set_refractive_index(nclad)
#    mat_b.set_refractive_index(ncore)

    mat_a = materials.make_material("Si_2021_Poulton")
    mat_b = materials.make_material("SiO2_smf28")
    mat_vac= materials.make_material("Vacuum")
    ncore = np.real(mat_a.refindex_n)
    nclad = np.real(mat_b.refindex_n)


    rcore = 3250  # radius 3.25 micron
    rn = 1000  # annular layer thickness
    acore = 2*rcore  # central layer diameter

    inc_shape = 'circ_onion'
    unitcell_x = rcore*7  # system size in nm

    wguide = nbapp.make_structure(inc_shape, unitcell_x,  inc_a_x = acore, # remove these factors of 2
                               inc_b_x=rn, inc_c_x=rn, inc_d_x=rn, inc_e_x=rn,
                               inc_f_x=rn, inc_g_x=rn, inc_h_x=rn, inc_i_x=rn,
                               inc_j_x=rn, inc_k_x=rn, inc_l_x=rn, inc_m_x=rn,
                               inc_n_x=rn, inc_o_x=rn,
                               material_bkg=mat_vac, material_a=mat_a,
                               material_b=mat_b, material_c=mat_a, material_d=mat_b, material_e=mat_a,
                               material_f=mat_b, material_g=mat_a, material_h=mat_b, material_i=mat_a,
                               material_j=mat_b, material_k=mat_a, material_l=mat_b, material_m=mat_a,
                               material_n=mat_b, material_o=mat_a,
                            lc_bkg=.1, lc_refine_1=3.0*refine_fac, lc_refine_2=3*refine_fac)

    #wguide.plot_mesh(prefix)
    #wguide.check_mesh()

    solve_em_dispersion(prefix, wguide, ncore, nclad, rcore)

    print(nbapp.final_report())


if __name__ == '__main__':
    do_main()
