""" Example showing how the 'onion' geometry template can be used
    to calculate the dispersion profile of standard SMF-28 optical fibre.
"""

import os
import time
import sys
import copy
import multiprocessing
import math
#import scipy.signal
#import scipy.optimize as sciopt
#import scipy.special as sp

import matplotlib.pyplot as plt
import numpy as np


sys.path.append("../backend/")

import numbat
import materials
import plotting
from numbattools import launch_worker_processes_and_wait
from nbanalytic import TwoLayerFiberEM, EMPoln

import starter

twopi = math.pi*2.0

# a helper function to make consistent pretty plots
def decorate_and_save_plot (fig, ax, fname, title, xlab, ylab, xlims, ylims,
                            legloc=None):
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_xlim(xlims[0], xlims[1])
    ax.set_ylim(ylims[0], ylims[1])
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


# Solve the analytic dispersion relation using worker processes
def solve_em_two_layer_fiber_analytic(kvec, nmodes, rco, ncore, nclad):

    fib_em = TwoLayerFiberEM(ncore, nclad, rco)

    # The worker function passed to CalcThread to do one task
    def solemrod_caller(args):
        (ik, k) = args # Matches queues passed to CalcThread

        mhy_lo = 1
        mhy_hi = 5
        #v_neff_TE = solve_chareq_em_fib2_disprel(chareq_em_fib2_TE_m, 'TE', k, nmodes, 0, 0, rco, ncore, nclad)
        #v_neff_TM = solve_chareq_em_fib2_disprel(chareq_em_fib2_TM_m, 'TM', k, nmodes, 0, 0, rco, ncore, nclad)
        #v_neff_hy = solve_chareq_em_fib2_disprel(chareq_em_fib2_hy_m, 'Hy', k, nmodes, mhy_lo, mhy_hi, rco, ncore, nclad) # m in [1,5)
        v_neff_TE = fib_em.find_neffs_for_k(k, EMPoln.TE, 0, 0, nmodes)[1]
        v_neff_TM = fib_em.find_neffs_for_k(k, EMPoln.TM, 0, 0, nmodes)[1]
        v_neff_hy = fib_em.find_neffs_for_k(k, EMPoln.HY, mhy_lo, mhy_hi, nmodes)[1]

        return (ik, v_neff_TE, v_neff_TM, v_neff_hy)


    # Multiprocessing here is for the demonstration of how to do it
    # It's probably not actually improving performance

    # Choose a reasonable number of processes to test the system
    # without grinding the computer to a halt
    num_cores = max(2,int(os.cpu_count()/4))

    manager = multiprocessing.Manager()
    q_work = multiprocessing.JoinableQueue()        # for assigning the work
    q_result = manager.Queue()      # for returning the answers

    # build all the tasks
    for ik,k in enumerate(kvec): q_work.put((ik, k))

    launch_worker_processes_and_wait(num_cores, solemrod_caller, q_result, q_work)

    #Collect all the data back into one matrix per polarisation
    neff_TE_an=np.zeros([len(kvec), nmodes], dtype=float)
    neff_TM_an=np.zeros([len(kvec), nmodes], dtype=float)
    an_neff_hy=np.zeros([len(kvec), nmodes], dtype=float)
    while not q_result.empty():  #only one process now, so no blocking to worry about
        (ik, v_neff_TE, v_neff_TM, v_neff_hy) = q_result.get()
        neff_TE_an[ik,:] = np.sort(v_neff_TE)[-nmodes:]
        neff_TM_an[ik,:] = np.sort(v_neff_TM)[-nmodes:]
        an_neff_hy[ik,:] = np.sort(v_neff_hy)[-nmodes:]

    return neff_TE_an, neff_TM_an, an_neff_hy

def solve_em_two_layer_fiber_numerical(prefix, wguide, kvec, nmodes, nbasis, rco, ncore, nclad):
    neff_num = np.zeros([len(kvec), nmodes], dtype=float)

    # Expected effective index of fundamental guided mode.
    n_eff = (ncore+nclad)/2.0
    n_field_outs = 5  # how many ksteps to write full set of fields
    field_out_skip = max(int(len(kvec)/n_field_outs),1)

    # Multiprocess-safe queues for communicating with workers
    q_work = multiprocessing.JoinableQueue()        # everything else
    q_result = multiprocessing.JoinableQueue()      # for returning the answers

    # Create work queues
    for ik, tk in enumerate(kvec):
        doplot =  (ik % field_out_skip == 0) # Time for some field plots!
        wg = copy.deepcopy(wguide)  # wguide and sim are not thread-safe when we plot mode profiles
        q_work.put((ik, tk, doplot, wg))


    # function to perform one worker FEM task
    def emcalc_caller(args):
        (ik, tk, doplot, wg) = args # Matches queues passed to launch_worker...

        t_lambda_nm = twopi/tk
        sim_EM= wguide.calc_EM_modes(nbasis, t_lambda_nm, n_eff)   # main work happens here!
        neff_k= np.sort(np.real(sim_EM.neff_all()))[-nmodes:]

        if doplot: # Only worker 1 will ever do this
            print('{0} is plotting elastic modes at iq = {1:d} of [0..{2:d}].'.format(
                multiprocessing.current_process().name, ik, len(kvec)-1))
            plotting.plot_mode_fields(sim_EM, EM_AC='EM_E', ivals=range(nmodes),
                                      prefix=prefix+'_%d'%ik, ticks=True)

        return (ik, tk, neff_k)

    num_cores = max(2,int(os.cpu_count()/4))
    num_cores=1
    launch_worker_processes_and_wait(num_cores, emcalc_caller, q_result, q_work)


    #Collect all the data back into one matrix
    while not q_result.empty():  #only one process now, so no blocking to worry about
        (ik, tk, neff_k) = q_result.get()
        neff_num[ik,:] = neff_k

    return neff_num


def make_em_plots(prefix, ssys, Vvec_an, neff_TE_an, neff_TM_an, neff_Hy_an,
                  Vvec_num, neff_num, ncore, nclad):

    Vlo, Vhi =Vvec_an[0], Vvec_an[-1]

    # Analytical only plot
    #fig, ax = plt.subplots()
    #plot_and_label(ax, Vvec_an, neff_TE_an, '-b', 'TE')
    #plot_and_label(ax, Vvec_an, neff_TM_an, '-g', 'TM')
    #plot_and_label(ax, Vvec_an, neff_Hy_an, '-r', 'Hy')
    #decorate_and_save_plot(fig, ax, prefix + '-emdisp_Vneff_exact.png',
    #        'Analytic dispersion of '+ssys,
    #        r'$V$ number ', r'$\bar{n}$', (Vlo, Vhi), (nclad, ncore))


    # Numerical only plot
    #fig, ax = plt.subplots()
    #ax.plot(Vvec_num, neff_num, 'xk')
    #decorate_and_save_plot(fig, ax, prefix + '-emdisp_Vneff_num.png',
    #        'Numerical dispersion of '+ssys,
    #        r'$V$ number ', r'$\bar{n}$', (Vlo, Vhi), (nclad, ncore))


    # General plot for neff
    fig, ax = plt.subplots()
    plot_and_label(ax, Vvec_an, neff_TE_an, '-b', 'TE')
    plot_and_label(ax, Vvec_an, neff_TM_an, '-g', 'TM')
    plot_and_label(ax, Vvec_an, neff_Hy_an, '-r', 'Hy')
    ax.plot(Vvec_num, neff_num, 'xk')
    decorate_and_save_plot(fig, ax, prefix + '-emdisp_Vneff.png',
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
    decorate_and_save_plot(fig, ax, prefix + '-emdisp_Vb.png',
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

    decorate_and_save_plot(fig, ax, prefix + '-emdisp_ng.png',
            'Dispersion of '+ssys,
            r'$V$ number ', r'Group index $n_g$', (Vlo, Vhi), (nclad, ncore))



def solve_em_dispersion(prefix, ssys, wguide, rcore, ncore, nclad):
    print('\n\nElectromagnetic dispersion problem')
    print(    '----------------------------------')
    nmodes = 20
    nbasis = 40

    ksteps_num=50
    ksteps_analytic=200

    lam_hi=3000    # wavelength range in nm
    lam_lo=500    # wavelength range in nm
    klo=twopi/lam_hi  # this is free space wavenumber in inverse nm
    khi=twopi/lam_lo
    kvec_an=np.linspace(klo, khi, ksteps_analytic)
    kvec_num=np.linspace(klo, khi, ksteps_num)

    Vvec_an = kvec_an*rcore*np.sqrt(ncore**2-nclad**2)
    Vvec_num = kvec_num*rcore*np.sqrt(ncore**2-nclad**2)


    print('Doing analytic problem')

    #plot_em_chareq_at_k(kvec[-1], rcore, ncore, nclad) # for checking zero finding behaviour

    neff_TE_an, neff_TM_an, neff_Hy_an = solve_em_two_layer_fiber_analytic(kvec_an, nmodes, rcore, ncore, nclad)

    print('Doing numerical problem')
    neff_num = solve_em_two_layer_fiber_numerical(prefix, wguide, kvec_num, nmodes, nbasis, rcore, ncore, nclad)

    make_em_plots(prefix, ssys, Vvec_an, neff_TE_an, neff_TM_an, neff_Hy_an,
                  Vvec_num, neff_num, ncore, nclad)

def do_main():

    start = time.time()

    pref0, refine_fac = starter.read_args(12, sys.argv)


    # Geometric Parameters - all in nm.

    smf28 = False  # which system to study: SMF-28 or chalc rod in silica

    mat_SiO2_GeO2= materials.make_material("SiO2GeO2_smf28")
    mat_SiO2= materials.make_material("SiO2_smf28")
    mat_vac= materials.make_material("Vacuum")
    mat_Si= materials.make_material("Si_2021_Poulton")
    mat_As2S3= materials.make_material("As2S3_2023_Steel")

    if smf28:
        ssys = 'SMF-28 optical fibre'
        dcore = 8200    # 8 micron core diameter
        dclad = 125000  # 125 micron cladding diameter
        mat_core = mat_SiO2_GeO2
        mat_clad = mat_SiO2
        mat_bkg  = mat_vac
        prefix =pref0+'-smf28'
    else:
        ssys = 'chalc rod in silica'
        dcore = 300      # 150 nm core diameter
        dclad = 800      # 500 nm cladding diameter
        mat_core = mat_As2S3
        mat_clad = mat_SiO2
        mat_bkg  = mat_vac
        prefix =pref0+'-sil'

    nbapp=numbat.NumBATApp(prefix)

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
        inc_shape = 'circ_onion2'
        unitcell_x = rcore*10  # system size in nm
        unitcell_y = unitcell_x


    wguide = nbapp.make_structure(unitcell_x, acore, inc_shape=inc_shape,  # remove these factors of 2
                               inc_b_x=aclad,
                            unitcell_y=unitcell_y,
                            material_bkg=mat_bkg,
                            material_a=mat_core,
                            material_b=mat_clad,
                            lc_bkg=.1, lc_refine_1=3.0*refine_fac, lc_refine_2=3*refine_fac)

    #wguide.plot_mesh(prefix)

    solve_em_dispersion(prefix, ssys, wguide, rcore, ncore, nclad)

    print(nbapp.final_report())


if __name__ == '__main__':
    do_main()
