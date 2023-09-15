""" Replicating the results of
    Generation of phonons from electrostriction in 
    small-core optical waveguides
    Laude et al.
    http://dx.doi.org/10.1063/1.4801936

    Replicating silica example for backwards SBS.
     -- Compare gain in log plot with Fig 4. in paper. 
         (Different though related quantities showing very similar resonant mode pattern.)
     -- Mode profiles m=5 (5.7 GHz) and
"""

import time
import datetime
import numpy as np
import sys
import matplotlib
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
wl_nm = 1550
unitcell_x = 7*wl_nm
unitcell_y = unitcell_x
inc_a_x = 1500
inc_a_y = 1000
inc_shape = 'rectangular'

# Optical Parameters
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 120
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix, refine_fac = starter.read_args(1, sys.argv)

# Use all specified parameters to create a waveguide object.
wguide = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=materials.make_material("SiO2_2013_Laude"),
                        lc_bkg=.1, lc_refine_1=5.0*refine_fac, lc_refine_2=5.0*refine_fac)

wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = 1.3

# Calculate Electromagnetic modes.
recalc=True
#recalc=False
if recalc:
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)
    sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
    sim_EM_pump.save_simulation(prefix+'_pump')
    sim_EM_Stokes.save_simulation(prefix+'_pump')

    plotting.plot_mode_fields(sim_EM_pump, xlim_min=0.4, xlim_max=0.4, ivals=[EM_ival_pump],
                         ylim_min=0.4, ylim_max=0.4, prefix=prefix)
else:
    sim_EM_pump = mode_calcs.load_simulation(prefix+'_pump')
    sim_EM_Stokes = mode_calcs.load_simulation(prefix+'_pump')


# Print the wavevectors of EM modes.
kzs = sim_EM_pump.kz_EM_all()
print('k_z of EM modes \n', np.round(np.real(kzs), 4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff_all()) 
print("n_eff = ", np.round(n_eff_sim, 4))

q_AC = np.real(sim_EM_pump.kz_EM_all()[EM_ival_pump] - sim_EM_Stokes.kz_EM_all()[EM_ival_Stokes])

shift_Hz = 8e9

# Calculate Acoustic modes.
if recalc:
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)
    sim_AC.save_simulation(prefix+'_ac')
else:
    sim_AC = mode_calcs.load_simulation(prefix+'_ac')


# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)

# Print the frequencies and gains of AC modes.
nus = sim_AC.nu_AC_all()
gains = SBS_gain[EM_ival_pump, EM_ival_Stokes, :]

print('Acoustic modes')
print('m      Freq (GHz)      Total gain (1/(Wm))')
for m in range(len(nus)):
    print('{0:3d} {1:11.6f}       {2:8.4f}'.format(m, np.real(nus[m])*1e-9, gains[m]))

# find indices selection of nplot highest gain modes to plot
nplot=min(20, len(nus))
high_g_indices = (np.abs(gains).argsort()[-nplot:])
high_g_indices.sort()

if recalc:
    plotting.plot_mode_fields(sim_AC, prefix=prefix, ivals=high_g_indices)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 4.e9  # Hz
freq_max = 13.e9  # Hz
plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    semilogy=True, prefix=prefix, mode_comps=True, dB=True)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5.3e9  # Hz
freq_max = 6.6e9  # Hz
plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    prefix=prefix, suffix='_zoom')

end = time.time()
print("\nSimulation time: {0:10.3f} secs.\n\n".format(end - start))
print("--------------------------------------------------------------------\n\n\n")

