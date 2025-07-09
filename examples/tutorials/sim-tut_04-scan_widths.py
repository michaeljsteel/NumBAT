""" Calculate the backward SBS gain spectra as a function of
    waveguide width, for silicon waveguides surrounded in air.

    Also shows how to use python multiprocessing library.
"""

import os
import sys
from multiprocessing import Pool
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

from pathlib import Path
sys.path.append(str(Path('../backend')))


import numbat
from nbtypes import SI_GHz
import materials

import integration
import plotgain
import starter



# Select the number of CPUs to use in simulation.
num_cores = os.cpu_count()
num_cores = 4

# Geometric Parameters - all in nm.
lambda_nm = 1550
inc_shape = 'rectangular'

num_modes_EM_pump = 30
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 50
EM_mode_index_pump = 0
EM_mode_index_Stokes = 0
AC_mode_index = 'All'

prefix, refine_fac = starter.read_args(4, sys.argv, refine=2)

nbapp = numbat.NumBATApp(prefix, prefix+'-out')

use_multiproc = num_cores >1  and nbapp.is_linux()

# Width previous simo's done for, with known meshing params
known_geo = 315.


def modes_n_gain(wwguide):
    width, wguide = wwguide
    print('Process %d: commencing mode calculation for width a_x = %f' % (
        os.getpid(), width))
    # Expected effective index of fundamental guided mode.
    n_eff = (wguide.get_material('a').refindex_n-0.1) * width/known_geo

    # Calculate Electromagnetic modes.
    simres_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff)
    simres_EM_Stokes = simres_EM_pump.clone_as_backward_modes()
    q_AC = np.real(simres_EM_pump.kz_EM(EM_mode_index_pump) -
                   simres_EM_Stokes.kz_EM(EM_mode_index_Stokes))

    # Calculate Acoustic modes.
    simres_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=simres_EM_pump)

    # Calculate interaction integrals and SBS gain.
    gain_box = integration.get_gains_and_qs(simres_EM_pump, simres_EM_Stokes, simres_AC, q_AC,
        EM_mode_index_pump=EM_mode_index_pump, EM_mode_index_Stokes=EM_mode_index_Stokes, AC_mode_index=AC_mode_index)

    print('Process %d: completed mode calculation for width a_x = %.3f' % (
        os.getpid(), width))
    return [width, simres_EM_pump, simres_AC, gain_box]


n_widths = 6
waveguide_widths = np.linspace(300, 350, n_widths)
l_wguides = []
# Scale meshing to new structures.
for width in waveguide_widths:
    msh_ratio = width/known_geo
    domain_x = 1.5*lambda_nm*msh_ratio
    domain_y = domain_x
    inc_a_x = width
    inc_a_y = 0.9*inc_a_x

    wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                               material_bkg=materials.make_material("Vacuum"),
                               material_a=materials.make_material(
                                   "Si_2016_Smith"),
                               lc_bkg=.05, lc_refine_1=4.0*refine_fac, lc_refine_2=4.0*refine_fac)

    l_wguides.append((width, wguide))

    #wguide.plot_refractive_index_profile(prefix+'_%d'%int(width))

# Run widths in parallel across num_cores CPUs using multiprocessing package.
if use_multiproc:
    with Pool(num_cores) as pool:
        l_width_data = pool.map(modes_n_gain, l_wguides)
else:
    l_width_data = list(map(modes_n_gain, l_wguides))

v_width_data = np.array(l_width_data, dtype=object)



n_effs = []
freqs_gains = []
interp_grid_points = 10000
int_min = 10e9
int_max = 26e9
interp_grid = np.linspace(int_min, int_max, interp_grid_points)

for i_w, width_obj in enumerate(v_width_data):
    interp_values = np.zeros(interp_grid_points)
    width, simres_EM_pump, simres_AC, gain_box = width_obj

    # Calculate the EM effective index of the waveguide (q_AC = 2*k_EM).
    n_eff_sim = simres_EM_pump.neff(0)
    n_effs.append(n_eff_sim)

    # Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
    freq_min = np.real(simres_AC.nu_AC_all()[0]) - 5e9  # Hz
    freq_max = np.real(simres_AC.nu_AC_all()[-1]) + 5e9  # Hz
    decorator = plotgain.Decorator()
    decorator.set_title(f'Gain for width $w={width:.2f}$ nm')
    gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max,
                      suffix=f'_wscan_{i_w}' ,  # include scan step in file name
                      decorator=decorator, logy=True)

    # Repeat calc to collect data for waterfall plot.
    tune_steps = 50000
    tune_range = 10e9  # GHz
    detuning_range = np.append(np.linspace(-1*tune_range, 0, tune_steps),
                               np.linspace(0, tune_range, tune_steps)[1:])  # Hz

    # Linewidth of Lorentzian is half the FWHM style linewidth.
    linewidth = gain_box.linewidth_Hz_all()/2
    num_modes = len(linewidth)
    gain_box.set_EM_modes(EM_mode_index_pump, EM_mode_index_Stokes)
    for m in range(num_modes):
        gain_list = np.real(gain_box.gain_total(m) * linewidth[m]**2/(linewidth[m]**2 + detuning_range**2))
        freq_list = np.real(simres_AC.nu_AC(m) + detuning_range)
        interp_spectrum = np.interp(interp_grid, freq_list, gain_list)
        interp_values += interp_spectrum
    freqs_gains.append(list(zip(interp_grid/SI_GHz, abs(interp_values))))
fgs = np.array(freqs_gains)

fig, ax = plt.subplots(1,2, figsize=(14,6))
ax[0].set_xlabel('Frequency (GHz)', fontsize=14)
ax[0].set_ylabel('|Gain| (1/Wm)', fontsize=14)
ax[0].set_yscale('log')
for ifg,fg in enumerate(fgs):
    ax[0].plot(fg[:,0], fg[:,1], lw=1, label=r'$d=$'+f'{waveguide_widths[ifg]:.1f} nm')
ax[0].legend(loc='upper left', fontsize='x-small')

ax[1].set_xlabel('Width', fontsize=16)
ax[1].set_ylabel('Effective index', fontsize=16)
ax[1].plot(waveguide_widths, n_effs, '-o')


plt.savefig(str(nbapp.outpath())+'-gain_spectra-scan.png')

print(nbapp.final_report())
