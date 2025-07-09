""" Replicating the results of
    Brillouin scattering self-cancellation
    Florez et al.
    http://dx.doi.org/10.1038/ncomms11759
"""


import sys

import numpy as np

from pathlib import Path
sys.path.append(str(Path('../backend')))

import numbat
from nbtypes import SI_GHz
import integration
import modecalcs
import materials

import starter



# Geometric Parameters - all in nm.
wl_nm = 1550
domain_x = 3*wl_nm
domain_y = domain_x
inc_a_x = 1160  # Diameter
inc_a_y = inc_a_x
inc_shape = 'circular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_mode_index_pump = 0
EM_mode_index_Stokes = 0
AC_mode_index = 'All'

# choose between faster or more accurate calculation
if len(sys.argv) > 1 and sys.argv[1] == 'fast=1':
    prefix = 'flit_06b-'
    refine_fac = 1
    print('\n\nCommencing NumBAT literature example 6b - fast mode')
else:
    prefix = 'lit_06b-'
    refine_fac = 5
    print('\n\nCommencing NumBAT literature example 6b')

prefix, refine_fac = starter.read_args(6, sys.argv, sub='b')

nbapp = numbat.NumBATApp(prefix)

# Use all specified parameters to create a waveguide object.
wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                              material_bkg=materials.make_material("Vacuum"),
                              material_a=materials.make_material(
                                  "SiO2_2013_Laude"),
                              # lc_bkg=1, lc_refine_1=600.0, lc_refine_2=300.0)
                              lc_bkg=.1, lc_refine_1=5*refine_fac, lc_refine_2=5.0*refine_fac)

wguide.plot_mesh(prefix)
# Expected effective index of fundamental guided mode.
n_eff = 1.4

# Calculate Electromagnetic Modes
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)
sim_EM_Stokes = sim_EM_pump.clone_as_backward_modes()

sim_EM_pump.plot_modes(xlim_min=0.2, xlim_max=0.2, mode_indices=range(5))

# Print the wavevectors of EM modes.
kzs = sim_EM_pump.kz_EM_all()
print('k_z of EM modes \n', np.round(np.real(kzs), 4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff_all())
print("n_eff = ", np.round(n_eff_sim, 4))

q_AC = np.real(sim_EM_pump.kz_EM_all()[
               EM_mode_index_pump] - sim_EM_Stokes.kz_EM_all()[EM_mode_index_Stokes])

shift_Hz = 4e9

# Calculate Acoustic Modes
sim_AC = wguide.calc_AC_modes(
    num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)

sim_AC.plot_modes()

# Print the frequencies of AC modes.
print('Freq of AC modes (GHz) \n', np.round(
    np.real(sim_AC.nu_AC_all())*1e-9, 4))

set_q_factor = 1000.

# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
gain_box = integration.get_gains_and_qs(sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_mode_index_pump=EM_mode_index_pump, EM_mode_index_Stokes=EM_mode_index_Stokes, AC_mode_index=AC_mode_index)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5 * SI_GHz
freq_max = 12* SI_GHz

gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5.2* SI_GHz
freq_max = 5.7* SI_GHz
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True, suffix='-5')

print(nbapp.final_report())
