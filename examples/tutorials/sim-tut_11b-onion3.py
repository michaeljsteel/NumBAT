""" Example showing how the 'onion3' geometry template can be used
    to simulate a circular Si waveguide clad in SiO2 with a coating.
"""

import sys
import numpy as np

from pathlib import Path
sys.path.append(str(Path('../backend')))


import numbat
import materials

import integration


import starter


# Geometric Parameters - all in nm.
lambda_nm = 1550

domain_x = 2500
domain_y = domain_x
inc_a_x = 700
inc_b_x = 350
inc_c_x = 50

inc_shape = 'onion3'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 80
EM_mode_index_pump = 1
EM_mode_index_Stokes = 1
AC_mode_index = 'All'

prefix, refine_fac = starter.read_args(11, sys.argv, sub='b')

nbapp = numbat.NumBATApp(prefix)

# Use of a more refined mesh to produce field plots.
wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x,
                           inc_b_x=inc_b_x, inc_c_x=inc_c_x,
                           material_bkg=materials.make_material("Vacuum"),
                           material_a=materials.make_material("As2S3_2021_Poulton"),
                           material_b=materials.make_material("SiO2_2021_Poulton"),
                           material_c=materials.make_material("PMMA"),
                           lc_bkg=.1, lc_refine_1=1.0*refine_fac, lc_refine_2=1*refine_fac)

# wguide.check_mesh()

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

new_calcs = True

# Calculate Electromagnetic modes.
if new_calcs:
    simres_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff)
    simres_EM_Stokes = simres_EM_pump.clone_as_backward_modes()

    # simres_EM_pump.save_simulation('tut_11b_pump')
    # simres_EM_Stokes.save_simulation('tut_11b_stokes')
else:
    simres_EM_pump = numbat.load_simulation('tut_11b_pump')
    simres_EM_Stokes = numbat.load_simulation('tut_11b_stokes')

# Display the wavevectors of EM modes.
v_kz = simres_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz):
    print(f'{i:3d}  {np.real(kz):.4e}')

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(simres_EM_pump.neff(0))
print("n_eff", np.round(n_eff_sim, 4))

# # Plot the E fields of the EM modes fields
# # Zoom in on the central region (of big unitcell) with xlim_, ylim_ args.
# # Only plot fields of fundamental (mode_index = 0) mode.

for ft in ('EM_E', 'EM_H'):
    simres_EM_pump.plot_modes(mode_indices=range(10), field_type = ft)

# Acoustic wavevector
q_AC = np.real(simres_EM_pump.kz_EM(EM_mode_index_pump) -
               simres_EM_Stokes.kz_EM(EM_mode_index_Stokes))

# Calculate Acoustic modes.
if new_calcs:
  simres_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=simres_EM_pump)
#  simres_AC.save_simulation('tut_11b_acoustic')
else:
  simres_AC = numbat.load_simulation('tut_11b_acoustic')

# Print the frequencies of AC modes.
v_nu = simres_AC.nu_AC_all()
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(v_nu):
  print(f'{i:3d}  {np.real(nu)*1e-9:.4e}')

simres_AC.plot_modes(xlim_min=-.2, xlim_max=-.2, ylim_min=-.2, ylim_max=-.2,
                           mode_indices=range(num_modes_AC), quiver_points=20)

# Calculate the acoustic loss from our fields.
# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
gain_box = integration.get_gains_and_qs(
    simres_EM_pump, simres_EM_Stokes, simres_AC, q_AC, EM_mode_index_pump=EM_mode_index_pump,
    EM_mode_index_Stokes=EM_mode_index_Stokes, AC_mode_index=AC_mode_index)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 4.e9
freq_max = 10.e9
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True )

print(nbapp.final_report())
