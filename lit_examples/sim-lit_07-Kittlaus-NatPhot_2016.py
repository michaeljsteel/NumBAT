""" Replicating the results of
    Large Brillouin amplification in silicon
    Kittlaus et al.
    http://dx.doi.org/10.1038/nphoton.2016.112
"""

import sys
import copy
import numpy as np

from pathlib import Path
sys.path.append(str(Path('../backend')))
import numbat

import materials
import modecalcs
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
domain_x = 3*wl_nm
domain_y = 0.4*domain_x
# Waveguide widths.
#inc_a_x = 1000
#inc_a_y = 80
rib_w = 1000
rib_h = 80
# Shape of the waveguide.
inc_shape = 'rib'

#slab_a_x = 3000
#slab_a_y = 130
slab_w = 3000
slab_h = 130

# Number of electromagnetic modes to solve for.
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
# Number of acoustic modes to solve for.
num_modes_AC = 80
# The EM pump mode(s) for which to calculate interaction with AC modes.
# Can specify a mode number (zero has lowest propagation constant) or 'All'.
EM_mode_index_pump = 0
# The EM Stokes mode(s) for which to calculate interaction with AC modes.
EM_mode_index_Stokes = 0
# The AC mode(s) for which to calculate interaction with EM modes.
AC_mode_index = 'All'

Si_110 = copy.deepcopy(materials.make_material("Si_2016_Smith"))
Si_110.rotate('y-axis', np.pi/4, save_rotated_tensors=True)

mat_vac = materials.make_material("Vacuum")
mat_rib = Si_110
mat_slab = Si_110

prefix, refine_fac = starter.read_args(7, sys.argv)

nbapp = numbat.NumBATApp(prefix)

# Use specified parameters to create a waveguide object.
# Note use of rough mesh for demonstration purposes.
wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, rib_w=rib_w, rib_h=rib_h, slab_w=slab_w, slab_h=slab_h,
                              material_bkg=mat_vac, material_a=mat_rib, material_b=mat_slab, symmetry_flag=False,
                              lc_bkg=.05, lc_refine_1=5.0, lc_refine_2=5.0)

wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

# Calculate Electromagnetic Modes
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)
# np.savez('wguide_data', sim_EM_pump=sim_EM_pump)
# npzfile = np.load('wguide_data.npz')
# sim_EM_pump = npzfile['sim_EM_pump'].tolist()
sim_EM_Stokes = sim_EM_pump.clone_as_forward_modes()

sim_EM_pump.plot_modes(xlim_min=0.3, xlim_max=0.3, mode_indices=range(5),
                         ylim_min=0.2, ylim_max=0.2, field_type='EM_E', )

# Print the wavevectors of EM modes.
print('k_z of EM modes \n', np.round(np.real(sim_EM_pump.kz_EM_all()), 4))

# Calculate the EM effective index of the waveguide.
print("n_eff = ", np.round(sim_EM_pump.neff_all(), 4))

q_AC = 5e6   #TODO: where is this from?

shift_Hz = 2e9

# Calculate Acoustic Modes
sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)
# np.savez('wguide_data_AC', sim_AC=sim_AC)
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC = npzfile['sim_AC'].tolist()

# Print the frequencies of AC modes.
print('Freq of AC modes (GHz) \n', np.round(np.real(sim_AC.nu_AC_all())*1e-9, 4))

sim_AC.plot_modes( mode_indices=range(40))

set_q_factor = 680.

# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
gain_box = integration.get_gains_and_qs(sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_mode_index_pump=EM_mode_index_pump, EM_mode_index_Stokes=EM_mode_index_Stokes, AC_mode_index=AC_mode_index,fixed_Q=set_q_factor)

print('Gains by acoustic mode:')
print('Ac. mode | Freq (GHz) | G_tot (1/mW) | G_PE (1/mW) | G_MB (1/mW)')
v_nu = sim_AC.nu_AC_all()
for (m, nu) in enumerate(v_nu):
    print(f'{m:7d}    {np.real(nu)*1e-9:9.4e} {gain_box.gain_total(m):13.3e} ',
          f'{gain_box.gain_PE(m):13.3e} {gain_box.gain_MB(m):13.3e} ')

freq_min = 2.0* SI_GHz
freq_max = 20.* SI_GHz
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True)


freq_min = 4.5* SI_GHz
freq_max = 5.5* SI_GHz
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True, suffix='zoom')

print(nbapp.final_report())
