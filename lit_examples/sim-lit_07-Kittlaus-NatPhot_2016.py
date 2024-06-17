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
import mode_calcs
import integration
import plotting

import starter


# Naming conventions
# AC: acoustic
# EM: electromagnetic
# q_AC: acoustic wavenumber


# Geometric Parameters - all in nm.
wl_nm = 1550 # Wavelength of EM wave in vacuum.
# Unit cell must be large to ensure fields are zero at boundary.
unitcell_x = 3*wl_nm
unitcell_y = 0.4*unitcell_x
# Waveguide widths.
inc_a_x = 1000
inc_a_y = 80
# Shape of the waveguide.
inc_shape = 'rib'

slab_a_x = 3000
slab_a_y = 130

# Number of electromagnetic modes to solve for.
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
# Number of acoustic modes to solve for.
num_modes_AC = 80
# The EM pump mode(s) for which to calculate interaction with AC modes.
# Can specify a mode number (zero has lowest propagation constant) or 'All'.
EM_ival_pump = 0
# The EM Stokes mode(s) for which to calculate interaction with AC modes.
EM_ival_Stokes = 0
# The AC mode(s) for which to calculate interaction with EM modes.
AC_ival = 'All'

Si_110 = copy.deepcopy(materials.make_material("Si_2016_Smith"))
Si_110.rotate_axis('y-axis', np.pi/4, save_rotated_tensors=True)


prefix, refine_fac = starter.read_args(7, sys.argv)

nbapp = numbat.NumBATApp(prefix)

# Use specified parameters to create a waveguide object.
# Note use of rough mesh for demonstration purposes.
wguide = nbapp.make_structure(inc_shape, unitcell_x, unitcell_y, inc_a_x, inc_a_y,
                        slab_a_x=slab_a_x, slab_a_y=slab_a_y,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=Si_110,
                        material_b=Si_110, symmetry_flag=False,
                        lc_bkg=.1, lc_refine_1=10.0, lc_refine_2=10.0)

wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

# Calculate Electromagnetic Modes
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)
# np.savez('wguide_data', sim_EM_pump=sim_EM_pump)
# npzfile = np.load('wguide_data.npz')
# sim_EM_pump = npzfile['sim_EM_pump'].tolist()

sim_EM_Stokes = mode_calcs.fwd_Stokes_modes(sim_EM_pump)

sim_EM_pump.plot_modes(xlim_min=0.3, xlim_max=0.3, ivals=range(5),
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

sim_AC.plot_modes( ivals=range(40))

set_q_factor = 680.

# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)

# Mask negligible gain values to improve clarity of print out.
threshold = 1e-3
masked_PE = np.ma.masked_inside(SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked_MB = np.ma.masked_inside(SBS_gain_MB[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked = np.ma.masked_inside(SBS_gain[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)

print("\n Displaying results with negligible components masked out")
print("SBS_gain [1/(Wm)] PE contribution \n", masked_PE)
print("SBS_gain [1/(Wm)] MB contribution \n", masked_MB)
print("SBS_gain [1/(Wm)] total \n", masked)

freq_min = 2.0e9  # GHz
freq_max = 20.e9  # GHz
plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
    EM_ival_pump, EM_ival_Stokes, AC_ival='All', freq_min=freq_min, freq_max=freq_max,
     suffix='')


freq_min = 4.5e9  # GHz
freq_max = 5.5e9  # GHz
plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
    EM_ival_pump, EM_ival_Stokes, AC_ival='All', freq_min=freq_min, freq_max=freq_max,
     suffix='zoom')

print(nbapp.final_report())
