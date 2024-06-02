""" Replicating the results of
    Germanium as a material for stimulated Brillouin scattering in the mid-infrared
    Wolff et al.
    https://doi.org/10.1364/OE.22.030735
"""

import time
import datetime
import numpy as np
import sys
import copy

sys.path.append("../backend/")
import materials
import mode_calcs
import integration
import plotting
from fortran import NumBAT
from nbtypes import PointGroup

import starter

# Naming conventions
# AC: acoustic
# EM: electromagnetic
# q_AC: acoustic wavenumber

start = time.time()

# Geometric Parameters - all in nm.
wl_nm = 4000 # Wavelength of EM wave in vacuum.
# Unit cell must be large to ensure fields are zero at boundary.
unitcell_x = 1.5*wl_nm
unitcell_y = 0.75*unitcell_x
# Waveguide widths.
inc_a_x = 1020
inc_a_y = 700
# Shape of the waveguide.
inc_shape = 'rectangular'


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

#reuse_fields=True   # use saved data
reuse_fields=False  # calculate from scratch

Ge_110 = copy.deepcopy(materials.make_material("Ge_cubic_2014_Wolff"))
print('Initial Ge_100:', Ge_110.full_str())
Ge_110.rotate_axis( 'y-axis', np.pi/4., save_rotated_tensors=True)
print('Rotated Ge_110:', Ge_110.full_str())

prefix, refine_fac = starter.read_args(10, sys.argv, sub='a')
nbapp = numbat.NumBATApp(prefix)

# Use specified parameters to create a waveguide object.
# Note use of rough mesh for demonstration purposes.
wguide = nbapp.make_structure(inc_shape, unitcell_x, unitcell_y, inc_a_x, inc_a_y,
                        material_bkg=materials.make_material("Si3N4_2014_Wolff"),
                        material_a=Ge_110,
                        lc_bkg=.1, lc_refine_1=5.0, lc_refine_2=5.0)
wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

# Calculate the Electromagnetic modes of the pump field.
if not reuse_fields:
  sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff)
  np.savez('wguide_data', sim_EM_pump=sim_EM_pump)
else:
  npzfile = np.load('wguide_data.npz', allow_pickle=True)
  sim_EM_pump = npzfile['sim_EM_pump'].tolist()

sim_EM_pump.analyse_symmetries(PointGroup.C2V)
sim_EM_pump.set_r0_offset(3.0e-6, -2.250e-6)

plotting.plot_modes(sim_EM_pump, xlim_min=0.2, xlim_max=0.2, ivals=[EM_ival_pump],
                         ylim_min=0.2, ylim_max=0.2, EM_AC='EM_E', num_ticks=3, ticks=True,
                         )

if not reuse_fields:
  sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
  np.savez('wguide_data2', sim_EM_Stokes=sim_EM_Stokes)
else:
  npzfile = np.load('wguide_data2.npz', allow_pickle=True)
  sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

# Print the wavevectors of EM modes.
v_kz=sim_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz): print('{0:3d}  {1:.4e}'.format(i, np.real(kz)))


# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff(0))
print("\n n_eff = ", np.round(n_eff_sim, 4))

q_AC = np.real(sim_EM_pump.kz_EM(EM_ival_pump) - sim_EM_Stokes.kz_EM(EM_ival_Stokes))
print('\n AC wavenumber (1/m) = ', np.round(q_AC, 4))


# Calculate Acoustic modes.
shift_Hz = 5.5*1e9 # select the lowest frequency to start FEM search from.
if not reuse_fields:
  sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)
  np.savez('wguide_data_AC', sim_AC=sim_AC)
else:
  npzfile = np.load('wguide_data_AC.npz', allow_pickle=True)
  sim_AC = npzfile['sim_AC'].tolist()

sim_AC.analyse_symmetries(PointGroup.C2V)
sim_AC.set_r0_offset(3.0e-6, -2.250e-6)

# Print the frequencies of AC modes.
v_nu=sim_AC.nu_AC_all()
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(v_nu): print('{0:3d}  {1:.4e}'.format(i, np.real(nu)*1e-9))


#set_Q_factor = 190 # set the mechanic Q manually

# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB. Also calculate acoustic loss alpha.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival #, fixed_Q=set_Q_factor
    )
# Mask negligible gain values to improve clarity of print out.
threshold = -1e-3
masked_PE = np.ma.masked_inside(SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked_MB = np.ma.masked_inside(SBS_gain_MB[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked = np.ma.masked_inside(SBS_gain[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
# Print the Backward SBS gain of the AC modes.
print("\n Displaying results with negligible components masked out")
print("SBS_gain [1/(Wm)] PE contribution \n", masked_PE)
print("SBS_gain [1/(Wm)] MB contribution \n", masked_MB)
print("SBS_gain [1/(Wm)] total \n", masked)

plotting.plot_modes(sim_AC,
     num_ticks=3, xlim_min=0.1, xlim_max=0.1,
     #modal_gains_PE=SBS_gain_PE[EM_ival_pump, EM_ival_Stokes,:],
     #modal_gains_MB=SBS_gain_MB[EM_ival_pump, EM_ival_Stokes,:],
     #modal_gains=SBS_gain[EM_ival_pump, EM_ival_Stokes,:]
                          )

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5.0e9  # Hz
freq_max = 11.0e9  # Hz

plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    )

print(nbapp.final_report())