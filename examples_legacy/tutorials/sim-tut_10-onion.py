""" Example showing how the 'onion' geometry template can be used
    to simulate a circular Si waveguide clad in SiO2.
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
domain_x = 3.5*lambda_nm
domain_y = domain_x
inc_shape = 'onion'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 20
EM_mode_index_pump = 0
EM_mode_index_Stokes = 0
AC_mode_index = 'All'

prefix, refine_fac = starter.read_args(10, sys.argv)

nbapp = numbat.NumBATApp(prefix)

d_hi = 200
d_lo = 400
mat_hi = materials.make_material("Si_2016_Smith")
mat_lo = materials.make_material("SiO2_2016_Smith")
# Use of a more refined mesh to produce field plots.
wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x=d_lo,
                           inc_b_x=d_hi,
                           inc_c_x=d_lo, inc_d_x=d_hi, inc_e_x=d_lo, inc_f_x=d_hi,
                           inc_g_x=d_lo, inc_h_x=d_hi, inc_i_x=d_lo, inc_j_x=d_hi,
                           inc_k_x=d_lo, inc_l_x=d_hi,
                           material_bkg=materials.make_material("Vacuum"),
                           material_a=mat_lo, material_b=mat_hi,
                           material_c=mat_lo, material_d=mat_hi,
                           material_e=mat_lo, material_f=mat_hi,
                           material_g=mat_lo, material_h=mat_hi,
                           material_i=mat_lo, material_j=mat_hi,
                           material_k=mat_lo, material_l=mat_hi,
                           lc_bkg=.1, lc_refine_1=5.0*refine_fac, lc_refine_2=5.0*refine_fac)

# wguide.check_mesh()
# Move origin to nominal centre of waveguide, y adjustment needs refinement based on mesh template.
# Note different shifts are allowed for EM and acoustic, because acoustic domain excludes vacuum regions
# Shifts are in nm
wguide.set_xyshift_em(-domain_x*.5, domain_y*.5)
wguide.set_xyshift_ac(-domain_x*.5, domain_y*.5)

# wguide.plot_mesh()

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

new_calcs = True

# Calculate Electromagnetic modes.
if new_calcs:
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff)
    sim_EM_Stokes = sim_EM_pump.clone_as_backward_modes()

    #sim_EM_pump.save_simulation('tut_10_pump')
    #sim_EM_Stokes.save_simulation('tut_10_stokes')
else:
    sim_EM_pump = numbat.load_simulation('tut_10_pump')
    sim_EM_Stokes = numbat.load_simulation('tut_10_stokes')

# Display the wavevectors of EM modes.
v_kz = sim_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz):
    print(f'{i:3d}  {np.real(kz):.4e}')

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff(0))
print("n_eff", np.round(n_eff_sim, 4))

# # Plot the E fields of the EM modes fields
# # Zoom in on the central region (of big unitcell) with xlim_, ylim_ args.
# # Only plot fields of fundamental (mode_index = 0) mode.
sim_EM_pump.plot_modes(xlim_min=0.3, xlim_max=0.3, ylim_min=0.3,
                          ylim_max=0.3, mode_indices=range(10), contours=True,
                           ticks=True, quiver_points=20)

sim_EM_pump.plot_modes(xlim_min=0.3, xlim_max=0.3, ylim_min=0.3,
                          ylim_max=0.3, mode_indices=range(10), contours=True, field_type='EM_H',
                           ticks=True, quiver_points=20)

# Acoustic wavevector
q_AC = np.real(sim_EM_pump.kz_EM(EM_mode_index_pump) -
               sim_EM_Stokes.kz_EM(EM_mode_index_Stokes))

# Calculate Acoustic modes.
if new_calcs:
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump)
    sim_AC.save_simulation('tut_10_acoustic')
else:
    sim_AC = numbat.load_simulation('tut_10_acoustic')

# Print the frequencies of AC modes.
v_nu = sim_AC.nu_AC_all()
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(v_nu):
    print(f'{i:3d}  {np.real(nu)*1e-9:.4e}')

sim_AC.plot_modes(contours=False, ticks=True, mode_indices=[10], quiver_points=20)

# Calculate the acoustic loss from our fields.
# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
gain_box = integration.get_gains_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC, EM_mode_index_pump=EM_mode_index_pump,
    EM_mode_index_Stokes=EM_mode_index_Stokes, AC_mode_index=AC_mode_index)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = np.real(sim_AC.nu_AC_all()[0]) - 2e9  # Hz
freq_max = np.real(sim_AC.nu_AC_all()[-1]) + 2e9  # Hz
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max)

print(nbapp.final_report())
