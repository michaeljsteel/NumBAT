""" Calculate the backward SBS gain spectra of a Si
    slot waveguide containing As2S3 on a SiO2 slab.

    This time include a capping layer of SiO2 and
    investigate the effect of this layer's thickness.
"""

import os
import sys
from multiprocessing import Pool
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
unitcell_x = 4*lambda_nm
unitcell_y = 0.3*unitcell_x
inc_shape = 'slot_coated'
inc_a_x = 150
inc_a_y = 190
inc_b_x = 250
# Current mesh template assume inc_b_y = inc_a_y
slab_a_x = 1000
slab_a_y = 100

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix, refine_fac = starter.read_args(8, sys.argv)

nbapp = numbat.NumBATApp(prefix)

# Function to return ac freqs for given coating thickness
def ac_mode_freqs(coat_y):
    print(f'Commencing mode calculation for coat_y = {coat_y}')

    wguide = nbapp.make_structure(inc_shape, unitcell_x, unitcell_y, inc_a_x, inc_a_y,
                            slab_a_x=slab_a_x, slab_a_y=slab_a_y, inc_b_x=inc_b_x,
                            coat_y=coat_y,
                            material_bkg=materials.make_material("Vacuum"),            # background
                            material_a=materials.make_material("As2S3_2017_Morrison"), # slot
                            material_b=materials.make_material("SiO2_2013_Laude"),     # slab
                            material_c=materials.make_material("Si_2016_Smith"),       # walls of slot
                            material_d=materials.make_material("SiO2_2013_Laude"),     # coating
                            lc_bkg=1, lc_refine_1=100.0*refine_fac, lc_refine_2=50.0*refine_fac)

    # Expected effective index of fundamental guided mode.
    n_eff = wguide.get_material('a').refindex_n-0.1

    # Calculate Electromagnetic modes.
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff=n_eff)
    sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

    q_AC = np.real(sim_EM_pump.kz_EM(EM_ival_pump) - sim_EM_Stokes.kz_EM(EM_ival_Stokes))

    shift_Hz = 4e9

    # Calculate Acoustic modes.
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)

    if coat_y == 20.0: # Shouldn't really test equality on floats like this
        sim_AC.plot_modes(suffix='_%i' %int(coat_y))

    set_q_factor = 1000.

    gain_box= integration.get_gains_and_qs(
        sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
        EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)

    # Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
    freq_min = 4e9
    freq_max = 14e9

    gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max,
         suffix='_%i' %int(coat_y))

    # Convert to GHz
    mode_freqs = sim_AC.nu_AC_all()*1.e-9
    # Clear memory
    wguide = sim_EM_pump = sim_EM_Stokes = sim_AC = None

    print(f'Completed mode calculation for coating coat_y = {coat_y}')

    # Return the frequencies and simulated k_ac value in a list
    return mode_freqs


n_coats = 11
coat_min = 20
coat_max = 200
coat_y_list = np.linspace(coat_min, coat_max, n_coats)

num_cores = os.cpu_count()  # should be appropriate for individual machine/vm, and memory!

use_multiproc = num_cores >1 and not nbapp.is_macos()

if use_multiproc:
    pool = Pool(num_cores)
    pooled_mode_freqs = pool.map(ac_mode_freqs, coat_y_list)
else:
    pooled_mode_freqs = map(ac_mode_freqs, coat_y_list)

# We will pack the above values into a single array for plotting purposes, initialise first
freq_arr = np.empty((n_coats, num_modes_AC))
for i_w, sim_freqs in enumerate(pooled_mode_freqs):
    # Set the value to the values in the frequency array
    freq_arr[i_w] = np.real(sim_freqs)

# Also plot a figure for reference

fig, ax = plt.subplots()
for idx in range(min(15,num_modes_AC)):
    # slicing in the row direction for plotting purposes
    freq_slice = freq_arr[:, idx]
    ax.plot(coat_y_list, freq_slice, '.g')

# Set the limits and plot axis labels
ax.set_xlim(coat_min,coat_max)
ax.set_xlabel(r'Coating Thickness (nm)')
ax.set_ylabel(r'Frequency (GHz)')
ax.set_ylim(2.5, 9)
fig.savefig(prefix+'-acdisp_coating.png', bbox_inches='tight')


print(nbapp.final_report())
