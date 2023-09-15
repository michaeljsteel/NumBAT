""" Example showing how the 'onion2' geometry template can be used
    to simulate a circular Si waveguide clad in SiO2 surrounded by air.

    Uses the two-layer onion structure which with vacuum surrounds makes a three material structure.

    With a silicon core and silica cladding, we do not expect this structure to guide elastic
    waves in the core.

    ***** Something is broken in assigning material properties at the central grid point.
"""

import time
import datetime
import numpy as np
import sys
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

sys.path.append("../backend/")
import materials
import objects
import mode_calcs
import integration
import plotting
from fortran import NumBAT

import starter



start = time.time()

# Geometric Parameters - all in nm.
lambda_nm = 1550
unitcell_x = 2500
unitcell_y = unitcell_x
inc_a_x = 700
inc_b_x = 350
inc_shape = 'onion2'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 80
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix, refine_fac = starter.read_args(11, sys.argv, sub='a')

# Use of a more refined mesh to produce field plots.
wguide = objects.Structure(unitcell_x,
                           inc_a_x, # inner diameter
                           inc_shape=inc_shape,
                        inc_b_x=inc_b_x, # first annulus width
                        unitcell_y=unitcell_y,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=materials.make_material("As2S3_2021_Poulton"),
                        material_b=materials.make_material("SiO2_2021_Poulton"),
                        lc_bkg=.1, lc_refine_1=2.0*refine_fac, lc_refine_2=2*refine_fac)

wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

new_calcs=True

# Calculate Electromagnetic modes.
if new_calcs:
  sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff, calc_EM_mode_energy=True)
  sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

  #sim_EM_pump.save_simulation('tut_11_pump')
  #sim_EM_Stokes.save_simulation('tut_11_stokes')
else:
  sim_EM_pump = mode_calcs.load_simulation('tut_11_pump')
  sim_EM_Stokes = mode_calcs.load_simulation('tut_11_stokes')

print('EM modes:\n')
kz_EM_mu =np.real(sim_EM_pump.kz_EM_all())*1e-6
neff_EM =sim_EM_pump.neff_all()
ng_EM =sim_EM_pump.ngroup_EM_all()
print('m    |   k_z [1/micron]  | neff  | ng')
for m in range(num_modes_EM_pump):
  print('{0:4d}  {1:12.6f}  {2:12.6f}  {3:12.6f}'.format(m, kz_EM_mu[m], neff_EM[m], ng_EM[m]))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff(0))
print("n_eff", np.round(n_eff_sim, 4))

# # Plot the E fields of the EM modes fields - specified with EM_AC='EM_E'.
# # Zoom in on the central region (of big unitcell) with xlim_, ylim_ args.
# # Only plot fields of fundamental (ival = 0) mode.
#npzfile = np.load('wguide_data.npz', allow_pickle=True)
#npzfile = np.load('wguide_data2.npz', allow_pickle=True)
#sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

for em_ac in ('EM_E', 'EM_H'):
    plotting.plot_mode_fields(sim_EM_pump, ivals=range(10), EM_AC=em_ac, prefix=prefix)

# Acoustic wavevector
q_AC = np.real(sim_EM_pump.kz_EM(EM_ival_pump) - sim_EM_Stokes.kz_EM(EM_ival_Stokes))

# Calculate Acoustic modes.
if new_calcs:
  sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, calc_AC_mode_power=True)
  #sim_AC.save_simulation('tut_11_acoustic')
else:
  sim_AC = mode_calcs.load_simulation('tut_11_acoustic')

print('AC mode properties (GHz) \n')
nu_AC = np.real(sim_AC.nu_AC_all())*1e-9
vp_AC = np.real(sim_AC.vp_AC_all())
vg_AC = np.real(sim_AC.vg_AC_all())
print('m    |   nu [GHz]  | vp [m/s]  | vg [m/s]')
for m in range(num_modes_AC):
  print('{0:4d}  {1:12.6f}  {2:12.2f}  {3:12.2f}'.format(m, nu_AC[m], vp_AC[m], vg_AC[m]))

sim_AC.calc_acoustic_losses()

plotting.plot_mode_fields(sim_AC, xlim_min=-.1, xlim_max=-.1, ylim_min=-.1, ylim_max=-.1,
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

end = time.time()
print("\nSimulation time: {0:10.3f} (secs.)".format(end - start))
print('\n\n')

