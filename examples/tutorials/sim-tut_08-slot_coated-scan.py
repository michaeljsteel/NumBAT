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
sys.path.append(str(Path('../../backend')))

import numbat
import materials

import integration


import starter


# Geometric Parameters - all in nm.
lambda_nm = 1550
domain_x = 2500
domain_x = 5000
domain_y = 1500
inc_shape = 'slot_coated'

slot_w = 150
rib_h = 190
rib_w = 250

coat_t= 50

slab_a_x = 2000
slab_a_y = 300

num_modes_EM_pump = 40
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_mode_index_pump = 0
EM_mode_index_Stokes = 0
AC_mode_index = 'All'

prefix, refine_fac = starter.read_args(8, sys.argv, refine=3)

nbapp = numbat.NumBATApp(prefix, prefix+'-out')

mat_vac=materials.make_material("Vacuum")  # background
mat_slot=materials.make_material("As2S3_2017_Morrison")  # slot
mat_slab=materials.make_material("SiO2_2013_Laude")  # slab
mat_ribs=materials.make_material("Si_2016_Smith")  # walls of
mat_coat= mat_slab


# Function to return ac freqs for given coating thickness
def ac_mode_freqs(wid_x):

    print(f'Commencing mode calculation for wid_x = {wid_x}')

    nbapp.set_outprefix(f'wid_{wid_x}')

    wguide = nbapp.make_structure(inc_shape, domain_x, domain_y,
                     slot_w=wid_x, rib_w=rib_w, rib_h=rib_h, slab_w=slab_a_x, slab_h=slab_a_y, coat_t=coat_t,
                      material_bkg=mat_vac, material_a=mat_slot, material_b=mat_slab, material_c=mat_ribs, material_d=mat_coat,
                            #lc_bkg=.05/refine_fac, lc_refine_1=5*refine_fac, lc_refine_2=5*refine_fac, lc_refine_3=5*refine_fac)
                            lc_bkg=.05/refine_fac, lc_refine_1=5, lc_refine_2=5, lc_refine_3=5)

    pref = prefix+f'_{wid_x:.1f}'
    wguide.plot_refractive_index_profile(pref)
    #wguide.plot_mesh(pref)

    # Expected effective index of fundamental guided mode.
    n_eff = wguide.get_material('a').refindex_n-0.1

    # Calculate Electromagnetic modes.
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff=n_eff)
    sim_EM_Stokes = sim_EM_pump.clone_as_backward_modes()
    #sim_EM_pump.plot_modes(mode_indices=range(5))

    # Calculate Acoustic modes.
    q_AC = np.real(sim_EM_pump.kz_EM(EM_mode_index_pump) - sim_EM_Stokes.kz_EM(EM_mode_index_Stokes))
    shift_Hz = 4e9
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)

    if wid_x == 150.0:
        sim_AC.plot_modes(prefix=pref, comps=('x','y','z', 'a','t'))

    # Calculate gain
    set_q_factor = 1000.
    gain_box = integration.get_gains_and_qs(
        sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
        EM_mode_index_pump=EM_mode_index_pump, EM_mode_index_Stokes=EM_mode_index_Stokes, AC_mode_index=AC_mode_index, fixed_Q=set_q_factor)

    # Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
    freq_min = 4e9
    freq_max = 14e9

    gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, suffix='_%i' %int(wid_x))

    # Convert to GHz
    mode_freqs = sim_AC.nu_AC_all()*1.e-9
    # Clear memory
    wguide = sim_EM_pump = sim_EM_Stokes = sim_AC = None

    print(f'Completed mode calculation for coating wid_x = {wid_x}')

    # Return the frequencies and simulated k_ac value in a list
    return mode_freqs


n_widths = 21
wid_min = 150
wid_max = 300
wid_x_list = np.linspace(wid_min, wid_max, n_widths)

num_cores = os.cpu_count()  # should be appropriate for individual machine/vm, and memory!

num_cores = 1
use_multiproc = num_cores >1 and not nbapp.is_macos()

if use_multiproc:
    pool = Pool(num_cores)
    pooled_mode_freqs = pool.map(ac_mode_freqs, wid_x_list)
else:
    pooled_mode_freqs = map(ac_mode_freqs, wid_x_list)

# We will pack the above values into a single array for plotting purposes, initialise first
freq_arr = np.empty((n_widths, num_modes_AC))
for i_w, sim_freqs in enumerate(pooled_mode_freqs):
    # Set the value to the values in the frequency array
    freq_arr[i_w] = np.real(sim_freqs)


fig, ax = plt.subplots()
for idx in range(min(15,num_modes_AC)): # Plot at most the lowest 15 modes
    # slicing in the row direction for plotting purposes
    freq_slice = freq_arr[:, idx]
    ax.plot(wid_x_list, freq_slice, '.g')

# Set the limits and plot axis labels
ax.set_xlim(wid_min,wid_max)
ax.set_xlabel(r'Slot width (nm)')
ax.set_ylabel(r'Frequency (GHz)')
#ax.set_ylim(2.5, 9)
fig.savefig(prefix+'-acdisp_slotwidth.png', bbox_inches='tight')


print(nbapp.final_report())
