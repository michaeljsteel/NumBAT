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

import matplotlib.pyplot as plt
import numpy as np


sys.path.append("../backend/")

import numbat
import materials
import plotting
from numbattools import launch_worker_threads_and_wait
from nbanalytic import ElasticRod

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

        arad_SI = arad*1e-9  # fix this messiness
        rod_ac = ElasticRod(rho, c11, c12, c44, arad_SI)
        v_Om_modes = rod_ac.find_nu_for_q(q, nu_hi, p, nmodes)[1]
        #v_Om_modes = findroots_elastic_rod_chareq(Vl, nu_hi, p, q, #nmodes, rho, c11, c12, c44, arad)
        return (p, iq, v_Om_modes)


    # Two simplified char eqs apply for p=0. We use p=-1 for one of them.
    plo = -1
    phi = 7

    num_cores = os.cpu_count()
    num_cores=1
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


    nbapp=numbat.NumBATApp(prefix)

    # Elastic problem of a single rod in vacuum
    rcore = 1000
    acore = 2*rcore
    inc_shape = 'circ_onion1'
    unitcell_x = rcore*2.5  # system size in nm
    unitcell_y = unitcell_x
    #mat_PMMA = materials.make_material('PMMA')  #

    mat_core = mat_As2S3
    mat_bkg  = mat_vac

    wguide = nbapp.make_structure(unitcell_x, acore, inc_shape=inc_shape,
                            unitcell_y=unitcell_y, inc_b_x  =rcore*.1,
                            material_bkg=mat_bkg,
                            material_a=mat_core,
                            lc_bkg=.1, lc_refine_1=3.0*refine_fac, lc_refine_2=3*refine_fac)

    prefix =pref0
    #wguide.plot_mesh(prefix)
    sim_EM= wguide.calc_EM_modes(40, 1550,1.5) # solve one EM step to prep the waveguide meshing

    solve_elastic_dispersion(prefix, ssys, wguide, sim_EM, rcore, mat_core)

    print(nbapp.final_report())




if __name__ == '__main__':
    do_main()
