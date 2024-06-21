""" Calculate the convergence as a function of FEM mesh for
    backward SBS gain spectra of a silicon waveguide
    surrounded in air.
"""


import sys
import time

import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
sys.path.append(str(Path('../backend')))

import numbat
import materials
import mode_calcs
import integration

import starter


# Geometric Parameters - all in nm.
lambda_nm = 1550
unitcell_x = 2.5*lambda_nm
unitcell_y = unitcell_x
inc_a_x = 300
inc_a_y = 280
inc_shape = 'rectangular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 20
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix, refine_fac = starter.read_args(5, sys.argv)

nbapp = numbat.NumBATApp(prefix, prefix+'-out')


# Warning: The fine grids in this list will take considerable time to run!
lc_list = [4, 10, 20, 50, 100, 200, 300, 400]
nu_lcs = len(lc_list)
lc_bkg_list = 1*np.ones(nu_lcs)
x_axis = lc_list
conv_list = []
time_list = []

for i_lc, lc_ref in enumerate(lc_list):
    start = time.time()
    print("\n Running simulation", i_lc+1, "/", nu_lcs)
    lc_refine_2 = lc_ref/2
    lc_bkg = lc_bkg_list[i_lc]
    wguide = nbapp.make_structure(inc_shape, unitcell_x, unitcell_y, inc_a_x, inc_a_y,
                               material_bkg=materials.make_material("Vacuum"),
                               material_a=materials.make_material(
                                   "Si_2016_Smith"),
                               lc_bkg=lc_bkg, lc_refine_1=lc_ref*refine_fac, lc_refine_2=lc_refine_2*refine_fac, force_mesh=True)

    # Expected effective index of fundamental guided mode.
    n_eff = wguide.get_material('a').refindex_n-0.1

    # Calculate Electromagnetic modes.
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff)
    sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

    # Calculate Acoustic modes.
    q_AC = np.real(sim_EM_pump.kz_EM(EM_ival_pump) -
                   sim_EM_Stokes.kz_EM(EM_ival_Stokes))
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump)

    # Calculate interaction integrals and SBS gain.
    gain_box= integration.get_gains_and_qs(
        sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
        EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)

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
                                                   


xlabel = "Mesh Refinement Factor"
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax2.yaxis.tick_left()
ax2.yaxis.set_label_position("left")
EM_plot_Mk = rel_mode_kz_EM*1e-6
error0 = np.abs((np.array(EM_plot_Mk[0:-1])-EM_plot_Mk[-1])/EM_plot_Mk[-1])
ax2.plot(x_axis[0:-1], error0, 'b-v', label='Mode #%i' % EM_ival_pump)
ax1.plot(x_axis, np.real(EM_plot_Mk), 'r-.o', label=r'EM k$_z$')
ax1.yaxis.tick_right()
ax1.spines['right'].set_color('red')
ax1.yaxis.label.set_color('red')
ax1.yaxis.set_label_position("right")
ax1.tick_params(axis='y', colors='red')
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels)
ax1.set_xlabel(xlabel)
ax1.set_ylabel(r"EM k$_z$ ($\times 10^6$ 1/m)")
ax2.set_ylabel(r"Relative Error EM k$_z$")
ax2.set_yscale('log')  # , nonposx='clip')
fig.savefig(nbapp.outpath()+'-convergence-freq_EM.png', bbox_inches='tight')

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax2.yaxis.tick_left()
ax2.yaxis.set_label_position("left")
for i_m, rel_mode in enumerate(rel_modes):
    rel_mode_freq_AC_plot_GHz = rel_mode_freq_AC[:, i_m]*1e-9
    error0 = np.abs((np.array(
        rel_mode_freq_AC_plot_GHz[0:-1])-rel_mode_freq_AC_plot_GHz[-1])/rel_mode_freq_AC_plot_GHz[-1])
    ax2.plot(x_axis[0:-1], error0, '-v', label='Mode #%i' % rel_mode)
    ax1.plot(x_axis, np.real(rel_mode_freq_AC_plot_GHz),
             '-.o', label=r'AC Freq mode #%i' % rel_mode)
ax1.yaxis.tick_right()
ax1.spines['right'].set_color('red')
ax1.yaxis.label.set_color('red')
ax1.yaxis.set_label_position("right")
ax1.tick_params(axis='y', colors='red')
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels)
ax1.set_xlabel(xlabel)
ax1.set_ylabel(r"AC Freq (GHz)")
ax2.set_ylabel(r"Relative Error AC Freq")
ax2.set_yscale('log')  # , nonposx='clip')
fig.savefig(nbapp.outpath()+'-convergence-freq_AC.png', bbox_inches='tight')

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax2.yaxis.tick_left()
ax2.yaxis.set_label_position("left")
for i_m, rel_mode in enumerate(rel_modes):
    rel_mode_gain_plot = rel_mode_gain[:, i_m]
    error0 = np.abs(
        (np.array(rel_mode_gain_plot[0:-1])-rel_mode_gain_plot[-1])/rel_mode_gain_plot[-1])
    ax2.plot(x_axis[0:-1], error0, '-v', label=r'Mode #%i' % rel_mode)
    ax1.plot(x_axis, -np.real(rel_mode_gain_plot),
             '-.o', label=r'Gain mode #%i' % rel_mode)
ax1.yaxis.tick_right()
ax1.spines['right'].set_color('red')
ax1.yaxis.label.set_color('red')
ax1.yaxis.set_label_position("right")
ax1.tick_params(axis='y', colors='red')
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels)
ax1.set_xlabel(xlabel)
ax1.set_ylabel(r'-Gain $\mathrm{mW}^{-1}$')
ax2.set_ylabel(r"Relative Error Gain")
ax2.set_yscale('log')  # , nonposx='clip')
fig.savefig(nbapp.outpath()+'-convergence-gain.png', bbox_inches='tight')

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax2.yaxis.tick_left()
ax2.yaxis.set_label_position("left")
for i_m, rel_mode in enumerate(rel_modes):
    rel_mode_gain_PE_plot = rel_mode_gain_PE[:, i_m]
    error0 = np.abs((np.array(
        rel_mode_gain_PE_plot[0:-1])-rel_mode_gain_PE_plot[-1])/rel_mode_gain_PE_plot[-1])
    ax2.plot(x_axis[0:-1], error0, '-v', label=r'Mode #%i' % rel_mode)
    ax1.plot(x_axis, -np.real(rel_mode_gain_PE_plot),
             '-.o', label=r'Gain mode #%i' % rel_mode)
ax1.yaxis.tick_right()
ax1.spines['right'].set_color('red')
ax1.yaxis.label.set_color('red')
ax1.yaxis.set_label_position("right")
ax1.tick_params(axis='y', colors='red')
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels)
ax1.set_xlabel(xlabel)
ax1.set_ylabel(r'-(PE Gain)  $\mathrm{mW}^{-1}$')
ax2.set_ylabel(r"Relative Error Gain (PE)")
ax2.set_yscale('log')  # , nonposx='clip')
fig.savefig(nbapp.outpath()+'-convergence-gain_PE.png', bbox_inches='tight')

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax2.yaxis.tick_left()
ax2.yaxis.set_label_position("left")
for i_m, rel_mode in enumerate(rel_modes):
    rel_mode_gain_MB_plot = rel_mode_gain_MB[:, i_m]
    error0 = np.abs((np.array(
        rel_mode_gain_MB_plot[0:-1])-rel_mode_gain_MB_plot[-1])/rel_mode_gain_MB_plot[-1])
    ax2.plot(x_axis[0:-1], error0, '-v', label=r'Mode #%i' % rel_mode)
    ax1.plot(x_axis, -np.real(rel_mode_gain_MB_plot),
             '-.o', label=r'Gain mode #%i' % rel_mode)
ax1.yaxis.tick_right()
ax1.spines['right'].set_color('red')
ax1.yaxis.label.set_color('red')
ax1.yaxis.set_label_position("right")
ax1.tick_params(axis='y', colors='red')
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels)
ax1.set_xlabel(xlabel)
ax1.set_ylabel(r'-(MB Gain)  $\mathrm{mW}^{-1}$')
ax2.set_ylabel(r"Relative Error Gain (MB)")
ax2.set_yscale('log')  # , nonposx='clip')
fig.savefig(nbapp.outpath()+'-convergence-gain_MB.png', bbox_inches='tight')

print("Calculation times (secs.): ", ', '.join(
    map(lambda x: f'{x:.2f}', time_list)))

print(nbapp.final_report())
