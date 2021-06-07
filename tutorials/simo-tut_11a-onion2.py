""" Example showing how the 'onion' geometry template can be used
    to simulate a circular Si waveguide clad in SiO2.
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


start = time.time()

# Geometric Parameters - all in nm.
wl_nm = 1550
unitcell_x = 3.5*wl_nm
unitcell_y = unitcell_x
inc_a_x = 800
inc_b_x = 500
inc_shape = 'onion2'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 20
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

if len(sys.argv)>1 and sys.argv[1]=='fast=1':  # choose between faster or more accurate calculation
  print('\n\nCommencing NumBAT tutorial 11a - fast mode')
  prefix_str = 'ftut_11a-'
  refine_fac=1
else:
  print('\n\nCommencing NumBAT tutorial 11a')
  prefix_str = 'tut_11a-'
  refine_fac=5

# Use of a more refined mesh to produce field plots.
wguide = objects.Struct(unitcell_x,inc_a_x,inc_shape=inc_shape,
                        inc_b_x=inc_b_x,
                        unitcell_y=unitcell_y,
                        material_bkg=materials.get_material("Vacuum"),
                        material_a=materials.get_material("Si_2016_Smith"),
                        material_b=materials.get_material("SiO2_2016_Smith"),
                        lc_bkg=.25, lc_refine_1=5.0*refine_fac, lc_refine_2=5*refine_fac, plt_mesh=True)


# Expected effective index of fundamental guided mode.
n_eff = wguide.material_a.n-0.1

new_calcs=True

# Calculate Electromagnetic modes.
if new_calcs:
  sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff, calc_EM_mode_energy=True)
  np.savez('wguide_data', sim_EM_pump=sim_EM_pump)

  sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
  np.savez('wguide_data2', sim_EM_Stokes=sim_EM_Stokes)
else:
  npzfile = np.load('wguide_data.npz', allow_pickle=True)
  sim_EM_pump = npzfile['sim_EM_pump'].tolist()
  npzfile = np.load('wguide_data2.npz', allow_pickle=True)
  sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

print('EM modes:\n')
kz_EM_mu =np.real(sim_EM_pump.kz_EM_all())*1e-6
neff_EM =sim_EM_pump.neff_all()
ng_EM =sim_EM_pump.ngroup_all()
print('m    |   k_z [1/micron]  | neff  | ng')
for m in range(num_modes_EM_pump):
  print('{0:4d}  {1:12.6f}  {2:12.6f}  {3:12.6f}'.format(m, kz_EM_mu[m], neff_EM[m], ng_EM[m]))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff(0))
print("n_eff", np.round(n_eff_sim, 4))

# # Plot the E fields of the EM modes fields - specified with EM_AC='EM_E'.
# # Zoom in on the central region (of big unitcell) with xlim_, ylim_ args.
# # Only plot fields of fundamental (ival = 0) mode.
plotting.plot_mode_fields(sim_EM_pump, xlim_min=0.3, xlim_max=0.3, ylim_min=0.3,
                         ylim_max=0.3, ivals=[0], contours=True, EM_AC='EM_E', 
                         prefix_str=prefix_str, ticks=True, quiver_points=20, comps=['Et'])

plotting.plot_mode_fields(sim_EM_pump, xlim_min=0.3, xlim_max=0.3, ylim_min=0.3,
                         ylim_max=0.3, ivals=[0], contours=True, EM_AC='EM_H', 
                         prefix_str=prefix_str, ticks=True, quiver_points=20, comps=['Ht'])

# Acoustic wavevector
k_AC = np.real(sim_EM_pump.kz_EM(EM_ival_pump) - sim_EM_Stokes.kz_EM(EM_ival_Stokes))

# Calculate Acoustic modes.
if new_calcs:
  sim_AC = wguide.calc_AC_modes(num_modes_AC, k_AC, EM_sim=sim_EM_pump, calc_AC_mode_power=True)
  np.savez('wguide_data_AC', sim_AC=sim_AC)
else:
  npzfile = np.load('wguide_data_AC.npz', allow_pickle=True)
  sim_AC = npzfile['sim_AC'].tolist()

print('AC mode properties (GHz) \n')
nu_AC = np.real(sim_AC.nu_AC_all())*1e-9
vp_AC = np.real(sim_AC.vp_AC_all())
vg_AC = np.real(sim_AC.vg_AC_all())
print('m    |   nu [GHz]  | vp [m/s]  | vg [m/s]')
for m in range(num_modes_AC):
  print('{0:4d}  {1:12.6f}  {2:12.2f}  {3:12.2f}'.format(m, nu_AC[m], vp_AC[m], vg_AC[m]))

sim_AC.calc_acoustic_losses()

plotting.plot_mode_fields(sim_AC, EM_AC='AC', pdf_png='png', contours=False, 
                         prefix_str=prefix_str, ticks=True, ivals=[0], quiver_points=20)

# Calculate the acoustic loss from our fields.
# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC, EM_ival_pump=EM_ival_pump, 
    EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = np.real(sim_AC.nu_AC_all()[0])*1e-9 - 2  # GHz
freq_max = np.real(sim_AC.nu_AC_all()[-1])*1e-9 + 2  # GHz
plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max, 
    prefix_str=prefix_str)

end = time.time()
print("\n Simulation time (sec.)", (end - start))

