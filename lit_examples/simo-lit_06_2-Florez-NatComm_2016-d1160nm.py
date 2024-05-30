""" Replicating the results of
    Brillouin scattering self-cancellation
    Florez et al.
    http://dx.doi.org/10.1038/ncomms11759
"""

import time
import numpy as np
import sys

sys.path.append("../backend/")
import numbat
import materials
import mode_calcs
import integration
import plotting

import starter



# Geometric Parameters - all in nm.
wl_nm = 1550
unitcell_x = 3*wl_nm
unitcell_y = unitcell_x
inc_a_x = 1160 # Diameter
inc_a_y = inc_a_x
inc_shape = 'circular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

if len(sys.argv)>1 and sys.argv[1]=='fast=1':  # choose between faster or more accurate calculation
  prefix = 'flit_06b-'
  refine_fac=1
  print('\n\nCommencing NumBAT literature example 6b - fast mode')
else:
  prefix = 'lit_06b-'
  refine_fac=5
  print('\n\nCommencing NumBAT literature example 6b')

prefix, refine_fac = starter.read_args(6, sys.argv, sub='b')

nbapp = numbat.NumBATApp(prefix)

# Use all specified parameters to create a waveguide object.
wguide = nbapp.make_structure(inc_shape, unitcell_x, unitcell_y, inc_a_x, inc_a_y,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=materials.make_material("SiO2_2013_Laude"),
                        #lc_bkg=1, lc_refine_1=600.0, lc_refine_2=300.0)
                        lc_bkg=.1, lc_refine_1=5*refine_fac, lc_refine_2=5.0*refine_fac)

wguide.plot_mesh(prefix)
# Expected effective index of fundamental guided mode.
n_eff = 1.4

# Calculate Electromagnetic Modes
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)
sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

plotting.plot_mode_fields(sim_EM_pump, xlim_min=0.2, xlim_max=0.2, ivals=range(5),
                         ylim_min=0.2, ylim_max=0.2, field_type='EM_E',
                         )

# Print the wavevectors of EM modes.
kzs = sim_EM_pump.kz_EM_all()
print('k_z of EM modes \n', np.round(np.real(kzs), 4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff_all())
print("n_eff = ", np.round(n_eff_sim, 4))

q_AC = np.real(sim_EM_pump.kz_EM_all()[EM_ival_pump] - sim_EM_Stokes.kz_EM_all()[EM_ival_Stokes])

shift_Hz = 4e9

# Calculate Acoustic Modes
sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)

plotting.plot_mode_fields(sim_AC,  )

# Print the frequencies of AC modes.
print('Freq of AC modes (GHz) \n', np.round(np.real(sim_AC.nu_AC_all())*1e-9, 4))

set_q_factor = 1000.

# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5e9  # Hz
freq_max = 12e9  # Hz
plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    )

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5.2e9  # GHz
freq_max = 5.7e9  # GHz
plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    semilogy=True, )
print(nbapp.final_report())