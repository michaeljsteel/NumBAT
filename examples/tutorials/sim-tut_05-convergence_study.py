""" Calculate the convergence as a function of FEM mesh for
    backward SBS gain spectra of a silicon waveguide
    surrounded in air.
"""


import sys
import time

import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
sys.path.append(str(Path('../../backend')))

import numbat
import materials

import integration

import starter

def start_plot():
    ax1 = plt.subplots()
    ax2 = ax1.twinx()

    return fig, ax1, ax2

def finish_plot(fig, ax1, ax2, ylabl, ylabr, suffix):
    xlabel = 'Mesh refinement factor'

    ax2.yaxis.tick_left()
    ax2.yaxis.set_label_position("left")
    ax1.yaxis.tick_right()
    ax1.spines['right'].set_color('red')
    ax1.yaxis.label.set_color('red')
    ax1.yaxis.set_label_position("right")
    ax1.tick_params(axis='y', colors='red')

    ax1.set_xscale('log')
    ax2.set_xscale('log')

    handles, labels = ax2.get_legend_handles_labels()

    ax2.legend(handles, labels)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabl)
    ax2.set_ylabel(ylabr)
    ax2.set_yscale('log')  # , nonposx='clip')

    fig.savefig(str(nbapp.outpath())+'-convergence-'+suffix+'.png', bbox_inches='tight')



# Geometric Parameters - all in nm.
lambda_nm = 1550
#domain_x = 2.5*lambda_nm
domain_x = 1000
domain_y = domain_x
inc_a_x = 300
inc_a_y = 280
inc_shape = 'rectangular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 20
EM_mode_index_pump = 0
EM_mode_index_Stokes = 0
AC_mode_index = 'All'

prefix, refine_fac = starter.read_args(5, sys.argv, refine=1.5)

nbapp = numbat.NumBATApp(prefix, prefix+'-out')


lc_list = [1, 2, 4, 8, 16, 32]
nu_lcs = len(lc_list)
lc_bkg_list = .25*np.ones(nu_lcs)
x_axis = lc_list
conv_list = []
time_list = []

for i_lc, lc_ref in enumerate(lc_list):
    start = time.time()
    print(f'\nRunning simulation {i_lc+1}/{nu_lcs} with mesh reduction factor {lc_ref} x {refine_fac} = {lc_ref*refine_fac}.')
    lc_bkg = lc_bkg_list[i_lc]
    wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                               material_bkg=materials.make_material("Vacuum"),
                               material_a=materials.make_material("Si_2016_Smith"),
                               #lc_bkg=lc_bkg, lc_refine_1=lc_ref*refine_fac,
                               #               lc_refine_2=(lc_ref/2)*refine_fac)
                               lc_bkg=lc_bkg/lc_ref, lc_refine_1=refine_fac, lc_refine_2=refine_fac/2)

    # Expected effective index of fundamental guided mode.
    n_eff = wguide.get_material('a').refindex_n-0.1

    # Calculate Electromagnetic modes.
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff)
    sim_EM_Stokes = sim_EM_pump.clone_as_backward_modes()

    # Calculate Acoustic modes.
    q_AC = np.real(sim_EM_pump.kz_EM(EM_mode_index_pump) -
                   sim_EM_Stokes.kz_EM(EM_mode_index_Stokes))
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump)

    # Calculate interaction integrals and SBS gain.
    gain_box= integration.get_gains_and_qs(
        sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
        EM_mode_index_pump=EM_mode_index_pump, EM_mode_index_Stokes=EM_mode_index_Stokes, AC_mode_index=AC_mode_index)

    conv_list.append([sim_EM_pump, sim_AC, gain_box.gain_total_all(),
                      gain_box.gain_PE_all(), gain_box.gain_MB_all()])
    end = time.time()
    time_list.append(end - start)

# It is crucial that you preselect modes with significant gain!
# Otherwise you will observe large relative errors similar to dividing by zero.
rel_modes = [3, 4, 8, 10]

# If you do not know the mode numbers of the significant AC modes you may wish to simply plot them all
# by uncommenting the line below and check if the modes with large gain have low relative errors.
# rel_modes = np.linspace(0,num_modes_AC-1,num_modes_AC)
rel_mode_kz_EM = np.zeros(nu_lcs, dtype=complex)
rel_mode_freq_AC = np.zeros((nu_lcs, len(rel_modes)), dtype=complex)
rel_mode_gain = np.zeros((nu_lcs, len(rel_modes)), dtype=complex)
rel_mode_gain_MB = np.zeros((nu_lcs, len(rel_modes)), dtype=complex)
rel_mode_gain_PE = np.zeros((nu_lcs, len(rel_modes)), dtype=complex)
for i_conv, conv_obj in enumerate(conv_list):
    rel_mode_kz_EM[i_conv] = conv_obj[0].kz_EM(0)
    for i_m, rel_mode in enumerate(rel_modes):
        rel_mode_freq_AC[i_conv, i_m] = conv_obj[1].nu_AC(rel_mode)
        rel_mode_gain[i_conv, i_m] = conv_obj[2][rel_mode]
        rel_mode_gain_PE[i_conv, i_m] = conv_obj[3][rel_mode]
        rel_mode_gain_MB[i_conv, i_m] = conv_obj[4][rel_mode]


fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
EM_plot_Mk = rel_mode_kz_EM*1e-6
error0 = np.abs((np.array(EM_plot_Mk[0:-1])-EM_plot_Mk[-1])/EM_plot_Mk[-1])
ax1.plot(x_axis, np.real(EM_plot_Mk), 'r-.o', label=r'EM k$_z$')
ax2.plot(x_axis[0:-1], error0, 'b-v', label='Mode #%i' % EM_mode_index_pump)
finish_plot(fig, ax1, ax2, r"EM k$_z$ ($\times 10^6$ 1/m)", r"Relative Error EM k$_z$", 'freq_EM')


fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
for i_m, rel_mode in enumerate(rel_modes):
    rel_mode_freq_AC_plot_GHz = rel_mode_freq_AC[:, i_m]*1e-9
    error0 = np.abs((np.array(
        rel_mode_freq_AC_plot_GHz[0:-1])-rel_mode_freq_AC_plot_GHz[-1])/rel_mode_freq_AC_plot_GHz[-1])
    ax1.plot(x_axis, np.real(rel_mode_freq_AC_plot_GHz), '-.o', label=r'AC Freq mode #%i' % rel_mode)
    ax2.plot(x_axis[0:-1], error0, '-v', label='Mode #%i' % rel_mode)
finish_plot(fig, ax1, ax2, r"AC Freq (GHz)", r"Relative Error AC Freq", 'freq_AC')


fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
for i_m, rel_mode in enumerate(rel_modes):
    rel_mode_gain_plot = rel_mode_gain[:, i_m]
    error0 = np.abs(
        (np.array(rel_mode_gain_plot[0:-1])-rel_mode_gain_plot[-1])/rel_mode_gain_plot[-1])
    ax1.plot(x_axis, -np.real(rel_mode_gain_plot), '-.o', label=r'Gain mode #%i' % rel_mode)
    ax2.plot(x_axis[0:-1], error0, '-v', label=r'Mode #%i' % rel_mode)
finish_plot(fig, ax1, ax2, r'-(Gain)  $\mathrm{mW}^{-1}$', r"Relative Error Gain", 'gain')

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
for i_m, rel_mode in enumerate(rel_modes):
    rel_mode_gain_PE_plot = rel_mode_gain_PE[:, i_m]
    error0 = np.abs((np.array(
        rel_mode_gain_PE_plot[0:-1])-rel_mode_gain_PE_plot[-1])/rel_mode_gain_PE_plot[-1])
    ax1.plot(x_axis, -np.real(rel_mode_gain_PE_plot), '-.o', label=r'Gain mode #%i' % rel_mode)
    ax2.plot(x_axis[0:-1], error0, '-v', label=r'Mode #%i' % rel_mode)
finish_plot(fig, ax1, ax2, r'-(PE Gain)  $\mathrm{mW}^{-1}$', r"Relative Error Gain (PE)", 'gain_PE')



fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
for i_m, rel_mode in enumerate(rel_modes):
    rel_mode_gain_MB_plot = rel_mode_gain_MB[:, i_m]
    error0 = np.abs((np.array(
        rel_mode_gain_MB_plot[0:-1])-rel_mode_gain_MB_plot[-1])/rel_mode_gain_MB_plot[-1])
    ax1.plot(x_axis, -np.real(rel_mode_gain_MB_plot), '-.o', label=r'Gain mode #%i' % rel_mode)
    ax2.plot(x_axis[0:-1], error0, '-v', label=r'Mode #%i' % rel_mode)
finish_plot(fig, ax1, ax2, r'-(MB Gain)  $\mathrm{mW}^{-1}$', r"Relative Error Gain (MB)", 'gain_MB')

print("Calculation times (secs.): ", ', '.join(
    map(lambda x: f'{x:.2f}', time_list)))

print(nbapp.final_report())
