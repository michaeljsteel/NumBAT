""" Replicating the results of
    Compact Brillouin devices through hybrid
    integration on Silicon
    Morrison et al.
    https://doi.org/10.1364/OPTICA.4.000847
"""


import sys

import numpy as np

from pathlib import Path
sys.path.append(str(Path('../backend')))
import numbat

import materials
import mode_calcs
import integration
from nbtypes import SI_GHz

import starter

# Naming conventions
# AC: acoustic
# EM: electromagnetic
# q_AC: acoustic wavenumber

# Geometric Parameters - all in nm.
wl_nm = 1550 # Wavelength of EM wave in vacuum.
# Unit cell must be large to ensure fields are zero at boundary.
domain_x = 6*wl_nm
domain_y = 0.75*domain_x
# Waveguide widths.
inc_a_x = 1900
inc_a_y = 680
# Shape of the waveguide.
inc_shape = 'rib_coated'

slab_a_x = 4000
slab_a_y = 1000

coat_x = 200
coat_y = 1000


# Number of electromagnetic modes to solve for.
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
# Number of acoustic modes to solve for.
num_modes_AC = 100
# The EM pump mode(s) for which to calculate interaction with AC modes.
# Can specify a mode number (zero has lowest propagation constant) or 'All'.
EM_ival_pump = 0
# The EM Stokes mode(s) for which to calculate interaction with AC modes.
EM_ival_Stokes = 0
# The AC mode(s) for which to calculate interaction with EM modes.
AC_ival = 'All'


prefix, refine_fac = starter.read_args(9, sys.argv)

nbapp = numbat.NumBATApp(prefix)

# Use specified parameters to create a waveguide object.
# Note use of rough mesh for demonstration purposes.

wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                        slab_a_x=slab_a_x, slab_a_y=slab_a_y, coat_x=coat_x, coat_y=coat_y,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=materials.make_material("As2S3_2017_Morrison"), # waveguide
                        material_b=materials.make_material("SiO2_2016_Smith"),     # slab
                        material_c=materials.make_material("SiO2_2016_Smith"),     # coating
                        lc_bkg=.1, lc_refine_1=10.0, lc_refine_2=10.0, lc_refine_3=10)
wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

# Calculate the Electromagnetic modes of the pump field.
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff)
# # np.savez('wguide_data', sim_EM_pump=sim_EM_pump)
# npzfile = np.load('wguide_data.npz')
# sim_EM_pump = npzfile['sim_EM_pump'].tolist()

sim_EM_pump.plot_modes(xlim_min=0.4, xlim_max=0.4, ivals=[EM_ival_pump],
                         ylim_min=0.3, ylim_max=0.3, field_type='EM_E', num_ticks=3,
                         )

# Calculate the Electromagnetic modes of the Stokes field.
sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
# np.savez('wguide_data2', sim_EM_Stokes=sim_EM_Stokes)
# npzfile = np.load('wguide_data2.npz')
# sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

# Print the wavevectors of EM modes.
print('\n k_z of EM modes \n', np.round(np.real(sim_EM_pump.kz_EM_all()),4))

# Calculate the EM effective index of the waveguide.
print("\n n_eff = ", np.round(sim_EM_pump.neff_all(), 4))

q_AC = np.real(sim_EM_pump.kz_EM_all()[EM_ival_pump] - sim_EM_Stokes.kz_EM_all()[EM_ival_Stokes])
print('\n AC wavenumber (1/m) = ', np.round(q_AC, 4))


q_AC= 2.*9173922.1698

# Calculate Acoustic modes.
shift_Hz = 7.5*1e9 # select the lowest frequency to start FEM search from.
sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)
# # np.savez('wguide_data_AC', sim_AC=sim_AC)
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC = npzfile['sim_AC'].tolist()

sim_AC.plot_modes( num_ticks=3, xlim_min=0.1, xlim_max=0.1)

# Print the frequencies of AC modes.
print('\n Freq of AC modes (GHz) \n', np.round(np.real(sim_AC.nu_AC_all())*1e-9, 4))

set_Q_factor = 190 # set the mechanic Q manually

# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB. Also calculate acoustic loss alpha.
gain_box = integration.get_gains_and_qs(sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_Q_factor)

print('Gains by acoustic mode:')
print('Ac. mode | Freq (GHz) | G_tot (1/mW) | G_PE (1/mW) | G_MB (1/mW)')
v_nu = sim_AC.nu_AC_all()
for (m, nu) in enumerate(v_nu):
    print(f'{m:7d}    {np.real(nu)*1e-9:9.4e} {gain_box.gain_total(m):13.3e} ',
          f'{gain_box.gain_PE(m):13.3e} {gain_box.gain_MB(m):13.3e} ')



# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 7.2*SI_GHz
freq_max = 8.1*SI_GHz

gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True)

print(nbapp.final_report())
