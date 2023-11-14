""" Example showing how the 'onion3' geometry template can be used
    to simulate a circular Si waveguide clad in SiO2 with a coating.
"""

import sys
import numpy as np

sys.path.append("../backend/")

import numbat
import materials
import objects
import mode_calcs
import integration
import plotting


import starter


# Geometric Parameters - all in nm.
lambda_nm = 1550

unitcell_x = 2500
unitcell_y = unitcell_x
inc_a_x = 700
inc_b_x = 350
inc_c_x = 50

inc_shape = 'onion3'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 80
EM_ival_pump = 1
EM_ival_Stokes = 1
AC_ival = 'All'

prefix, refine_fac = starter.read_args(11, sys.argv, sub='b')

numbat = numbat.NumBAT()

# Use of a more refined mesh to produce field plots.
wguide = objects.Structure(unitcell_x, inc_a_x, inc_shape=inc_shape,
                           inc_b_x=inc_b_x, inc_c_x=inc_c_x,
                           unitcell_y=unitcell_y,
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
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff)
    sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

    # sim_EM_pump.save_simulation('tut_11b_pump')
    # sim_EM_Stokes.save_simulation('tut_11b_stokes')
else:
    sim_EM_pump = mode_calcs.load_simulation('tut_11b_pump')
    sim_EM_Stokes = mode_calcs.load_simulation('tut_11b_stokes')

# Display the wavevectors of EM modes.
v_kz = sim_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz):
    print('{0:3d}  {1:.4e}'.format(i, np.real(kz)))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff(0))
print("n_eff", np.round(n_eff_sim, 4))

# # Plot the E fields of the EM modes fields - specified with EM_AC='EM_E'.
# # Zoom in on the central region (of big unitcell) with xlim_, ylim_ args.
# # Only plot fields of fundamental (ival = 0) mode.

for em_ac in ('EM_E', 'EM_H'):
    plotting.plot_mode_fields(
        sim_EM_pump, ivals=range(10), EM_AC=em_ac, prefix=prefix)

# Acoustic wavevector
q_AC = np.real(sim_EM_pump.kz_EM(EM_ival_pump) -
               sim_EM_Stokes.kz_EM(EM_ival_Stokes))

# Calculate Acoustic modes.
if new_calcs:
  sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump)
#  sim_AC.save_simulation('tut_11b_acoustic')
else:
  sim_AC = mode_calcs.load_simulation('tut_11b_acoustic')

# Print the frequencies of AC modes.
v_nu = sim_AC.nu_AC_all()
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(v_nu):
  print(f'{i:3d}  {np.real(nu)*1e-9:.4e}')

plotting.plot_mode_fields(sim_AC, xlim_min=-.2, xlim_max=-.2, ylim_min=-.2, ylim_max=-.2,
                          prefix=prefix, ivals=range(num_modes_AC), quiver_points=20)

# Calculate the acoustic loss from our fields.
# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC, EM_ival_pump=EM_ival_pump,
    EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 4.e9
freq_max = 10.e9
plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
                           EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
                           prefix=prefix)

print(numbat.final_report())
