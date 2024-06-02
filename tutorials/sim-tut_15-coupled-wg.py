""" Example showing how the 'onion' geometry template can be used
    to simulate a circular Si waveguide clad in SiO2.
"""

import sys
import numpy as np

sys.path.append("../backend/")

import numbat
import materials
import mode_calcs

import plotting

import starter



# Geometric Parameters - all in nm.
lambda_nm = 1550
unitcell_x = 2.*lambda_nm
unitcell_y = .5*unitcell_x
inc_a_x = 800
inc_a_y = 500
inc_b_x = 800
inc_b_y = 500
sep = 200
inc_shape = 'twoincl'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 20
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix_str, refine_fac = starter.read_args(14, sys.argv)

nbapp=numbat.NumBATApp(prefix)

# Use of a more refined mesh to produce field plots.
wguide = nbapp.make_structure(unitcell_x, inc_a_x, inc_shape=inc_shape,
                           inc_a_y=inc_a_y,
                           inc_b_x=inc_b_x,
                           inc_b_y=inc_b_y,
                           two_inc_sep=sep,
                           unitcell_y=unitcell_y,
                           material_bkg=materials.make_material(
                               "SiO2_2021_Poulton"),
                           material_a=materials.make_material(
                               "As2S3_2021_Poulton"),
                           material_b=materials.make_material(
                               "As2S3_2021_Poulton"),
                           lc_bkg=.1, lc_refine_1=3.0*refine_fac, lc_refine_2=3*refine_fac)

wguide.plot_mesh(prefix_str)

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

new_calcs = True

# Calculate Electromagnetic modes.
if new_calcs:
  sim_EM_pump = wguide.calc_EM_modes(
      num_modes_EM_pump, lambda_nm, n_eff, calc_EM_mode_energy=True)
  sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

  sim_EM_pump.save_simulation('tut_14_pump')
  sim_EM_Stokes.save_simulation('tut_14_stokes')
else:
  sim_EM_pump = mode_calcs.load_simulation('tut_14_pump')
  sim_EM_Stokes = mode_calcs.load_simulation('tut_14_stokes')

print('EM modes:\n')
kz_EM_mu = np.real(sim_EM_pump.kz_EM_all())*1e-6
neff_EM = sim_EM_pump.neff_all()
ng_EM = sim_EM_pump.ngroup_EM_all()
print('m    |   k_z [1/micron]  | neff  | ng')
for m in range(num_modes_EM_pump):
  print(f'{m:4d}  {kz_EM_mu[m]:12.6f}  {neff_EM[m]:12.6f}  {ng_EM[m]:12.6f}')

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff(0))
print("n_eff", np.round(n_eff_sim, 4))

sim_EM_pump.plot_modes(ivals=range(8), contours=True, ticks=True, quiver_points=20)

sim_EM_pump.plot_modes(ivals=range(8), contours=True, field_type='EM_H', ticks=True, quiver_points=20)

# Acoustic wavevector
q_AC = np.real(sim_EM_pump.kz_EM(EM_ival_pump) -
               sim_EM_Stokes.kz_EM(EM_ival_Stokes))

# Calculate Acoustic modes.
if new_calcs:
  sim_AC = wguide.calc_AC_modes(
      num_modes_AC, q_AC, EM_sim=sim_EM_pump, calc_AC_mode_power=True)
  sim_AC.save_simulation('tut_14_acoustic')
else:
  sim_AC = mode_calcs.load_simulation('tut_14_acoustic')

print('AC mode properties (GHz) \n')
nu_AC = np.real(sim_AC.nu_AC_all())*1e-9
vp_AC = np.real(sim_AC.vp_AC_all())
vg_AC = np.real(sim_AC.vg_AC_all())
print('m    |   nu [GHz]  | vp [m/s]  | vg [m/s]')
for m in range(num_modes_AC):
  print(f'{m:4d}  {nu_AC[m]:12.6f}  {vp_AC[m]:12.2f}  {vg_AC[m]:12.2f}')

sim_AC.calc_acoustic_losses()

sim_AC.plot_modes(contours=False,
                          ticks=True, ivals=range(10), quiver_points=20)

print(nbapp.final_report())
