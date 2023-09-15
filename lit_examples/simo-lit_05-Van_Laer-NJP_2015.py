""" Replicating the results of
    Net on-chip Brillouin gain based on suspended 
    silicon nanowires
    Van Laer et al.
    http://dx.doi.org/10.1088/1367-2630/17/11/115005
"""

import time
import datetime
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
import copy

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
unitcell_x = 5*wl_nm
unitcell_y = 0.5*unitcell_x
inc_a_x = 450
inc_a_y = 230
inc_shape = 'rectangular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 60
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'


prefix, refine_fac = starter.read_args(5, sys.argv)

# Rotate crystal axis of Si from <100> to <110>, starting with same Si_2016_Smith data.
Si_110 = copy.deepcopy(materials.make_material("Si_2016_Smith"))
Si_110.rotate_axis(np.pi/4,'y-axis', save_rotated_tensors=True)
# Use all specified parameters to create a waveguide object.
wguide = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=Si_110, symmetry_flag=False,
                        lc_bkg=.1, lc_refine_1=15.0*refine_fac, lc_refine_2=15.0*refine_fac)
wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

# Calculate Electromagnetic Modes
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)
# np.savez('wguide_data', sim_EM_pump=sim_EM_pump)
# npzfile = np.load('wguide_data.npz')
# sim_EM_pump = npzfile['sim_EM_pump'].tolist()

sim_EM_Stokes = mode_calcs.fwd_Stokes_modes(sim_EM_pump)
# np.savez('wguide_data2', sim_EM_Stokes=sim_EM_Stokes)
# npzfile = np.load('wguide_data2.npz')
# sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

plotting.plot_mode_fields(sim_EM_pump, xlim_min=0.45, xlim_max=0.45, 
                         ivals=[EM_ival_pump], ylim_min=0.45, ylim_max=0.45, 
                         EM_AC='EM_E', n_points=1500, prefix=prefix)

# Print the wavevectors of EM modes.
kzs = sim_EM_pump.kz_EM_all()
print('k_z of EM modes \n', np.round(np.real(kzs), 4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff_all()) 
print("n_eff = ", np.round(n_eff_sim, 4))

q_AC = 5 # close but not quite zero

# Calculate Acoustic Modes
sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump)
# np.savez('wguide_data_AC', sim_AC=sim_AC)
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC = npzfile['sim_AC'].tolist()

# Print the frequencies of AC modes.
print('Freq of AC modes (GHz) \n', np.round(np.real(sim_AC.nu_AC_all())*1e-9, 4))

plotting.plot_mode_fields(sim_AC, prefix=prefix)

set_q_factor = 230 # NJP

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

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5.0e9 # Hz
freq_max = 20.0e9 # Hz
plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max, prefix=prefix)

end = time.time()
print("\n Simulation time (sec.)", (end - start))
print("--------------------------------------------------------------------\n\n\n")
