""" We've now covered most of the features of NumBAT.
    In the following tutorials we'll show how to study different
    geometries and materials.

    Calculate the backward SBS gain spectra of a silicon waveguide
    surrounded by vacuum (air).
"""

import sys
import numpy as np

from pathlib import Path
sys.path.append(str(Path('../../backend')))

import numbat
import materials

import integration

import starter


# Geometric Parameters - all in nm.
lambda_nm = 1550
domain_x = 4500
domain_y = 1.0*domain_x
inc_a_x = 550
inc_a_y = inc_a_x
inc_shape = 'circular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_mode_index_pump = 0
EM_mode_index_Stokes = 0
AC_mode_index = 'All'

prefix, refine_fac = starter.read_args(6, sys.argv, sub='a')

nbapp = numbat.NumBATApp(prefix)

wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=materials.make_material("SiO2_2016_Smith"),
                        lc_bkg=.05, lc_refine_1=1.5*refine_fac, lc_refine_2=1.5*refine_fac)

wguide.plot_refractive_index_profile(prefix)

# Expected effective index of fundamental guided mode.
n_eff = 1.4

use_old_fields=False     # run the calculation from scratch
#use_old_fields=True   # reuse saved fields from previous calculation

# Calculate Electromagnetic Modes
if use_old_fields:
    simres_EM_pump = numbat.load_simulation(prefix+'_pump')
    simres_EM_Stokes = numbat.load_simulation(prefix+'_stokes')
else:
    simres_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff=n_eff)
    simres_EM_Stokes = simres_EM_pump.clone_as_backward_modes()
    simres_EM_pump.save_simulation(prefix+'_pump')
    simres_EM_Stokes.save_simulation(prefix+'_stokes')

print('\nPlotting EM fields')
trim=0.0
#trim=0.0

# Display the wavevectors of EM modes.
v_kz=simres_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz): print(f'{i:3d}  {np.real(kz):.4e}')

simres_EM_pump.plot_modes(mode_indices=range(5), xlim_min=trim, xlim_max=trim, ylim_min=trim, ylim_max=trim)

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(simres_EM_pump.neff(0))
print(f"n_eff = {n_eff_sim:.4e}")

# Acoustic wavevector
q_AC = np.real(simres_EM_pump.kz_EM(EM_mode_index_pump) - simres_EM_Stokes.kz_EM(EM_mode_index_Stokes))

shift_Hz = 4e9

# Calculate Acoustic modes.
if use_old_fields:
    simres_AC = numbat.load_simulation(prefix+'_acoustic')
else:
    simres_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=simres_EM_pump, shift_Hz=shift_Hz)
    simres_AC.save_simulation(prefix+'_acoustic')

# Print the frequencies of AC modes.
v_nu=simres_AC.nu_AC_all()
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(v_nu): print(f'{i:3d}  {np.real(nu)*1e-9:.5f}')

simres_AC.plot_modes(quiver_points=20)

set_q_factor = 1000.
#set_q_factor = None

print('\nCalculating gains')
# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
gain_box = integration.get_gains_and_qs(
    simres_EM_pump, simres_EM_Stokes, simres_AC, q_AC,
    EM_mode_index_pump=EM_mode_index_pump, EM_mode_index_Stokes=EM_mode_index_Stokes,
    AC_mode_index=AC_mode_index, fixed_Q=set_q_factor)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5e9  # Hz
freq_max = 15e9  # Hz

# Generate gain spectra on linear and log vertical scales.
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True)

print(nbapp.final_report())
