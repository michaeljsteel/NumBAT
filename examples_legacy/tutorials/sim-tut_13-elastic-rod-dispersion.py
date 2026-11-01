""" Example showing how the 'onion' geometry template can be used
    to calculate the dispersion profile of standard SMF-28 optical fibre.
"""

import os
import time
import sys
import copy
import functools
import multiprocessing

import matplotlib.pyplot as plt
import numpy as np


from pathlib import Path

sys.path.append(str(Path("../../backend")))


import numbat
import materials
from nbthread import launch_worker_processes_and_wait
from nbanalytic.elasticcylinder import ElasticIsotropicCylinder
from nbtypes import twopi

import starter


# a helper function to make consistent pretty plots
def decorate_and_save_plot(
    fig, ax, fname, title, xlab, ylab, xlims, ylims, legend_loc=None, legend_ncol=1
):
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_xlim(xlims[0], xlims[1])
    ax.set_ylim(ylims[0], ylims[1])
    ax.set_title(title)
    if legend_loc is not None:
        ax.legend(loc=legend_loc, fontsize="xx-small", ncol=legend_ncol)
    else:
        ax.legend(fontsize="xx-small", ncol=legend_ncol)
    fig.savefig(fname)


def plot_and_label(ax, vx, mat, sty, lab, col=None):
    if col is None:
        ax.plot(vx, mat, sty, ms=5)
        ax.plot(vx, mat[:, 0], sty, label=lab, ms=5)
    else:
        ax.plot(vx, mat, sty, color=col, ms=5)
        ax.plot(vx, mat[:, 0], sty, label=lab, color=col, ms=5)


QLABEL = r"$q\, [\mathrm{μm}^{-1}]$ "
NU_LABEL = r"$\nu$ [GHz]"
NEFF_LABEL = r"$\bar{n}=v_s/v$"

# Two simplified char eqs apply for p=0. We use p=-1 for one of them.
P_LO = -1
P_HI = 7


def p_label(p):
    """Legend label for azimuthal order p."""
    if p == -1:
        return r"$p=0$, Torsional"
    if p == 0:
        return r"$p=0$, Poch."
    return rf"$p={p}$"


def calc_neff(qvec, nu, Vref):
    """Effective index Vref*q/Omega, with a small offset to avoid division by zero."""
    return Vref * qvec[:, np.newaxis] / (nu * twopi + 1e-5)


def add_sound_lines(ax, qvec, Vl, Vs, **kwargs):
    """Dotted longitudinal and shear sound lines on a (q, nu) plot."""
    ax.plot(qvec * 1e-6, qvec * Vl / (twopi * 1e9), ":", **kwargs)
    ax.plot(qvec * 1e-6, qvec * Vs / (twopi * 1e9), ":", **kwargs)


def add_neff_lines(ax, Vl, Vs):
    """Reference lines at the sound-line indices on a (q, neff) plot."""
    ax.axhline(1, color="gray", linestyle="dotted")
    ax.axhline(Vs / Vl, color="gray", linestyle="dotted")


def save_dispersion_txt(fname, qvec, nu, nmodes):
    """Write dispersion data with qvec as 0th column."""
    mdata = np.hstack((qvec[:, np.newaxis], nu))
    np.savetxt(fname, mdata, fmt="%.8e", header="q [/m]    " + "nu_j [Hz]   " * nmodes)


# Module-level (picklable) worker function to perform one analytic task.
def callee_elastic_rod_analytical(args, rho, c11, c12, c44, arad, nu_hi, nmodes):
    (p, iq, q) = args  # Matches queues passed to launch_worker...

    arad_SI = arad * 1e-9  # fix this messiness
    rod_ac = ElasticIsotropicCylinder(rho, c11, c12, c44, arad_SI)
    v_Om_modes = rod_ac.find_Omega_at_q(q, nu_hi * 2 * np.pi, p, nmodes)[1]
    return (p, iq, v_Om_modes)


def plot_chareq_debug(rod_ac, qvec):
    """Plot the char equation at q as func of Omega, to check it behaves nicely for root-finding."""
    ysc_by_p = {
        -1: 5e-9,
        0: 1e-5,
        1: 2e-3,
        2: 1e3,
        3: 0.1e-9,
        4: 1e-0,
        5: 0.1e-9,
        6: 1e-8,
        7: 1,
    }
    Omvec = np.linspace(1.0e5, twopi * 2e9, 200)

    for p in range(P_LO, P_HI + 1):
        figp, axs = plt.subplots(2, 1)
        for iq, q in enumerate(qvec):
            if iq % 5:
                continue

            drvec = np.array([rod_ac.chareq_elastic_cylinder(Om, p, q) for Om in Omvec])
            lab = rf"$q={q * 1e-6:.2f}/$ μm"
            axs[0].plot(Omvec / (twopi * 1e9), drvec, ".", label=lab)
            axs[1].plot(Omvec / (twopi * 1e9), np.log10(np.abs(drvec)), ".", label=lab)

        # Now the drvec=1 values representing no sols will be skipped
        axs[0].set_ylim(-ysc_by_p[p], ysc_by_p[p])
        for ax in axs:
            ax.set_xlabel(r"$\Omega/2\pi$  [GHz]")
            ax.axhline(0, color="gray", linestyle="dotted")
            ax.legend()
        figp.savefig(f"omdr_p{p}.png")


def compute_analytic_dispersion(qvec, nmodes, rho, c11, c12, c44, arad, nu_hi):
    """Find the analytic elastic modes for all p and q. Returns nu_an[p-P_LO, iq, mode]."""
    task_func = functools.partial(
        callee_elastic_rod_analytical,
        rho=rho, c11=c11, c12=c12, c44=c44, arad=arad, nu_hi=nu_hi, nmodes=nmodes,
    )

    num_cores = 1
    manager = multiprocessing.Manager()
    q_work = manager.Queue()  # for assigning the work
    q_result = manager.Queue()  # for returning the answers

    for p in range(P_LO, P_HI + 1):
        for iq, q in enumerate(qvec):
            q_work.put((p, iq, q))

    print(f"Assigned {q_work.qsize()} tasks to {num_cores} processes.")

    # No plotting is done in this calc so no plot queue is required.
    launch_worker_processes_and_wait(num_cores, task_func, q_result, q_work)

    nu_an = np.zeros([P_HI - P_LO + 1, len(qvec), nmodes], dtype=float)
    while not q_result.empty():  # workers have finished, so no blocking to worry about
        (p, iq, v_Om_modes) = q_result.get()
        nu_an[p - P_LO, iq, : len(v_Om_modes)] = v_Om_modes / twopi
    return nu_an


def plot_analytic_dispersion(prefix, qvec, nu_an, Vl, Vs, nu_hi):
    fig, ax = plt.subplots()
    fign, axn = plt.subplots()
    cmap = plt.get_cmap("tab10")

    neff_an = calc_neff(qvec, nu_an, Vs)

    for p in range(P_LO, P_HI + 1):
        lab = p_label(p)
        col = cmap(p - P_LO)
        plot_and_label(ax, qvec * 1e-6, nu_an[p - P_LO, :, :] * 1e-9, ".", lab, col)
        plot_and_label(axn, qvec * 1e-6, neff_an[p - P_LO, :, :], ".", lab, col)

    add_sound_lines(ax, qvec, Vl, Vs, color="gray", lw=0.5)
    decorate_and_save_plot(
        fig, ax, prefix + "-acdisp_qnu_exact.png",
        "Analytic acoustic dispersion", QLABEL, NU_LABEL,
        (0, qvec[-1] * 1e-6), (0, nu_hi * 1e-9 / 2),
        legend_ncol=2,
    )

    add_neff_lines(axn, Vl, Vs)
    decorate_and_save_plot(
        fign, axn, prefix + "-acdisp_qneff_exact.png",
        "Analytic acoustic dispersion", QLABEL, NEFF_LABEL,
        (0, qvec[-1] * 1e-6), (0.4, 1.2),
        legend_ncol=2,
    )


def solve_elastic_rod_analytical(prefix, qvec, nmodes, coremat, arad):
    cmat = coremat.stiffness_c_IJ
    c11 = cmat[1, 1].real
    c12 = cmat[1, 2].real
    c44 = cmat[4, 4].real
    rho = coremat.rho
    Vl = coremat.Vac_longitudinal()
    Vs = coremat.Vac_shear()
    nu_hi = 1.75 * Vs * qvec[-1] / twopi  # Gives a suitable number of modes

    print("Material properties:")
    print(
        f"c11: {c11/1e9} GPa, c12: {c12/1e9} GPa, c44: {c44/1e9} GPa, Vs {Vs:.2f} m/s, Vl {Vl:.2f} m/s"
    )

    plot_chareqs = False
    if plot_chareqs:
        rod_ac = ElasticIsotropicCylinder(rho, c11, c12, c44, arad * 1e-9)
        plot_chareq_debug(rod_ac, qvec)

    nu_an = compute_analytic_dispersion(qvec, nmodes, rho, c11, c12, c44, arad, nu_hi)
    plot_analytic_dispersion(prefix, qvec, nu_an, Vl, Vs, nu_hi)
    return nu_an


# Module-level (picklable) worker function to perform one FEM task.
def callee_elastic_rod_numerical(args, nbasis_AC, Vs, simres_EM, nmodes, nq, prefix):
    (iq, tq, doplot, wg) = args  # Matches queues passed to launch_worker...

    # look for modes not to far from the shear sound line
    shift_Hz = 0.8 * tq * Vs / twopi
    simres_AC = wg.calc_AC_modes(nbasis_AC, tq, shift_Hz=shift_Hz, EM_sim=simres_EM)
    v_nu_num = np.sort(np.real(simres_AC.nu_AC_all()))[:nmodes]

    if doplot:
        print(
            f"{multiprocessing.current_process().name} is plotting elastic modes "
            f"at iq = {iq:d} of [0..{nq - 1:d}]."
        )
        simres_AC.plot_modes(mode_indices=range(nmodes), prefix=f"{prefix}_{iq}")

    return (iq, tq, v_nu_num)


def compute_numerical_dispersion(prefix, qvec, nmodes, wguide, simres_EM, Vs):
    """Run the FEM elastic calculation at every q. Returns nu_num[iq, mode]."""
    nbasis_AC = min(80, nmodes * 2)  # 80 ish

    n_field_outs = 2  # how many ksteps to write full set of fields
    field_out_skip = max(int(len(qvec) / n_field_outs), 1)

    task_func = functools.partial(
        callee_elastic_rod_numerical,
        nbasis_AC=nbasis_AC, Vs=Vs, simres_EM=simres_EM, nmodes=nmodes,
        nq=len(qvec), prefix=prefix,
    )

    # Multiprocess-safe queues for communicating with workers
    manager = multiprocessing.Manager()
    q_work = manager.Queue()  # for assigning the work
    q_result = manager.Queue()  # for returning the answers

    for iq, tq in enumerate(qvec):
        doplot = iq % field_out_skip == 0  # Time for some field plots!
        # wguide and sim are not safe to share when we plot mode profiles
        q_work.put((iq, tq, doplot, copy.deepcopy(wguide)))

    # Each process has its own matplotlib state, so plotting in workers is safe
    num_cores = max(4, int(os.cpu_count() / 4))
    launch_worker_processes_and_wait(num_cores, task_func, q_result, q_work)

    nu_num = np.zeros([len(qvec), nmodes], dtype=float)
    while not q_result.empty():  # workers have finished, so no blocking to worry about
        (iq, tq, nu_modes) = q_result.get()
        nu_num[iq, :] = nu_modes
    return nu_num


def plot_numerical_dispersion(prefix, qvec, nu_num, Vl, Vs, nu_hi):
    fig, ax = plt.subplots()
    plot_and_label(ax, qvec * 1e-6, nu_num * 1e-9, "xk", "FEM")
    add_sound_lines(ax, qvec, Vl, Vs, color="gray")
    decorate_and_save_plot(
        fig, ax, prefix + "-acdisp_qnu_num.png",
        "Elastic dispersion of single rod", QLABEL, NU_LABEL,
        (0, qvec[-1] * 1e-6), (0, nu_hi * 1e-9),
    )

    fign, axn = plt.subplots()
    plot_and_label(axn, qvec * 1e-6, calc_neff(qvec, nu_num, Vs), "xk", "FEM")
    add_neff_lines(axn, Vl, Vs)
    decorate_and_save_plot(
        fign, axn, prefix + "-acdisp_qneff_num.png",
        "Elastic dispersion of single rod", QLABEL, NEFF_LABEL,
        (0, qvec[-1] * 1e-6), (0.4, 1.2),
        "upper right",
    )


def solve_elastic_rod_numerical(prefix, qvec, nmodes, wguide, simres_EM, cmat):
    Vl = cmat.Vac_longitudinal()
    Vs = cmat.Vac_shear()
    nu_hi = 1.75 * Vs * qvec[-1] / twopi  # Gives a suitable number of modes

    nu_num = compute_numerical_dispersion(prefix, qvec, nmodes, wguide, simres_EM, Vs)
    plot_numerical_dispersion(prefix, qvec, nu_num, Vl, Vs, nu_hi)
    return nu_num


def plot_joint_dispersion(prefix, qvec_an, nu_an, qvec_num, nu_num, Vl, Vs, nu_hi, nmodes):
    fig, ax = plt.subplots()
    fign, axn = plt.subplots()
    cmap = plt.get_cmap("tab10")

    neff_an = calc_neff(qvec_an, nu_an, Vs)
    neff_num = calc_neff(qvec_num, nu_num, Vs)

    # Row ip of nu_an is p = ip + P_LO because of the 2 p=0 cases
    for ip in range(nu_an.shape[0]):
        p = ip + P_LO
        lab = p_label(p)
        col = cmap(ip)
        plot_and_label(ax, qvec_an * 1e-6, nu_an[ip, :, :] * 1e-9, ".", lab, col)
        plot_and_label(axn, qvec_an * 1e-6, neff_an[ip, :, :], ".", lab, col)
        save_dispersion_txt(f"{prefix}-acdisp_qnu_p{p}.txt", qvec_an, nu_an[ip, :, :], nmodes)

    # Numerical part of joint (q,nu) plot
    plot_and_label(ax, qvec_num * 1e-6, nu_num * 1e-9, "xk", "FEM")
    add_sound_lines(ax, qvec_num, Vl, Vs, lw=0.5)
    save_dispersion_txt(f"{prefix}-acdisp_qnu_num.txt", qvec_num, nu_num, nmodes)

    decorate_and_save_plot(
        fig, ax, prefix + "-acdisp_qnu.png",
        "Elastic dispersion of single rod", QLABEL, NU_LABEL,
        (0, qvec_num[-1] * 1e-6), (0, nu_hi * 1e-9),
        legend_ncol=2,
    )

    # Joint (q, neff) plot
    plot_and_label(axn, qvec_num * 1e-6, neff_num, "xk", "FEM")
    add_neff_lines(axn, Vl, Vs)
    decorate_and_save_plot(
        fign, axn, prefix + "-acdisp_qneff.png",
        "Acoustic effective index of single rod", QLABEL, NEFF_LABEL,
        (0, qvec_num[-1] * 1e-6), (0.4, 1.2),
        legend_loc="upper left", legend_ncol=2,
    )


def solve_elastic_dispersion(prefix, ssys, wguide, simres_EM, rcore, mat_core):
    print("\n\nAcoustic dispersion problem")
    print("---------------------------")

    Vl = mat_core.Vac_longitudinal()
    Vs = mat_core.Vac_shear()

    qsteps_num = 51
    qsteps_analytic = 201

    qlo = 0.0001e6
    qhi = 5e6
    qvec_an = np.linspace(qlo, qhi, qsteps_analytic)
    qvec_num = np.linspace(qlo, qhi, qsteps_num)
    nu_hi = 1.75 * Vs * qhi / twopi  # Gives a suitable number of modes

    nmodes = 40

    print("Doing analytic problem")
    nu_an = solve_elastic_rod_analytical(prefix, qvec_an, nmodes, mat_core, rcore)

    print("Doing numerical problem")
    nu_num = solve_elastic_rod_numerical(prefix, qvec_num, nmodes, wguide, simres_EM, mat_core)

    plot_joint_dispersion(prefix, qvec_an, nu_an, qvec_num, nu_num, Vl, Vs, nu_hi, nmodes)


def do_main():

    start = time.time()

    pref0, refine_fac = starter.read_args(13, sys.argv, refine=4)

    # Geometric Parameters - all in nm.

    smf28 = 'smf28' in sys.argv  # which system to study: SMF-28 or chalc rod in silica


    mat_SiO2_GeO2 = materials.make_material("SiO2GeO2_smf28")
    mat_SiO2 = materials.make_material("SiO2_smf28")
    mat_vac = materials.make_material("Vacuum")
    # mat_Si= materials.make_material("Si_2021_Poulton")
    mat_As2S3 = materials.make_material("As2S3_2023_Steel")

    if smf28:
        ssys = "SMF-28 optical fibre"
        dcore = 8200  # 8 micron core diameter
        dclad = 125000  # 125 micron cladding diameter
        mat_core = mat_SiO2_GeO2
        mat_clad = mat_SiO2
        mat_bkg = mat_vac
        prefix = pref0 + "-smf28"
    else:
        ssys = "chalc rod in silica"
        dcore = 300  # 150 nm core diameter
        dclad = 800  # 500 nm cladding diameter
        mat_core = mat_As2S3
        mat_clad = mat_SiO2
        mat_bkg = mat_vac
        prefix = pref0 + "-chalcsil"

    nbapp = numbat.NumBATApp(prefix, prefix+'-out')

    # Elastic problem of a single rod in vacuum
    rcore = 1000
    acore = 2 * rcore
    inc_shape = "circ_onion1"
    domain_x = rcore * 2.5  # system size in nm
    domain_y = domain_x
    # mat_PMMA = materials.make_material('PMMA')  #

    mat_core = mat_As2S3
    mat_bkg = mat_vac

    wguide = nbapp.make_structure(
        inc_shape,
        domain_x,
        domain_y,
        acore,
        inc_b_x=rcore * 0.1,
        material_bkg=mat_bkg,
        material_a=mat_core,
        lc_bkg=0.1,
        lc_refine_1=3.0 * refine_fac,
        lc_refine_2=3 * refine_fac,
    )

    #prefix = pref0
    # wguide.plot_mesh(prefix)

    # move to solve_elastic_rod_numerical
    # solve one EM step to prep the waveguide meshing
    simres_EM = wguide.calc_EM_modes(40, 1550, 1.5)

    solve_elastic_dispersion(prefix, ssys, wguide, simres_EM, rcore, mat_core)

    print(nbapp.final_report())


if __name__ == "__main__":
    do_main()
