# Simulation of SBS in LNOI waveguides 

import os
import sys
import threading
import numpy as np
import scipy.io as sc
import csv

import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
#from matplotlib.colors import colorConverter

#sys.path.append("../../NumBAT/backend/")
sys.path.append("../backend/")


import numbat
import materials
import mode_calcs
import integration
import plotting
#import starter

# Geometric Parameters - all in nm.
lambda_nm = 1550
unitcell_x = 5000
unitcell_y = 0.6*unitcell_x

slab_b_x = 3000    # width of buried part of rib
slab_a_y = 100     # height of buried part of rib

inc_a_x  = 1000     # width of top of rib
inc_a_y  = 500     # height of top of rib -> rib_height 
theta = 70
slab_a_x = 2*(inc_a_y/np.tan(theta*np.pi/180))+inc_a_x    # width at main surface


inc_shape = 'trapezoidal_rib' 


num_modes_EM_pump = 20                  
num_modes_EM_Stokes = num_modes_EM_pump   
num_modes_AC = 80                        
EM_ival_pump = 0                         
EM_ival_Stokes = 0
AC_ival = 'All'

prefix = 'tmplisa2'

nbapp = numbat.NumBATApp(prefix)

# inc_a_x inc_a_y slab_a_x slab_a_y 

mat_upper_surr = materials.make_material("Vacuum")
mat_rib=materials.make_material("LiNbO3_2023_Comsol")
mat_substrate=materials.make_material("SiO2_2023_Steel")


# Material tensors are by default  x->,   y up,   z -out (along waveguide)
# Rotate crystal so that x=right,   y= in (along waveguide), z=up
#mat_rib.rotate_axis(3*np.pi/2, 'x-axis')  
#mat_rib.rotate_axis(np.pi/2, [1,0,0])


#creating the waveguide here 
#wguide = nbapp.make_structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
#                   slab_a_x=slab_a_x, slab_a_y=slab_a_y, slab_b_x=slab_b_x, 
#                   material_bkg = mat_upper_surr, 
#                   material_a   = mat_rib, material_b   = mat_rib, material_c   = mat_substrate,
#                   lc_bkg=.1, lc_refine_1=10.0, lc_refine_2=20.0) 
wguide = nbapp.make_structure(inc_shape, unitcell_x,unitcell_y,inc_a_x, inc_a_y,
                   slab_a_x=slab_a_x, slab_a_y=slab_a_y, slab_b_x=slab_b_x, 
                   material_bkg = mat_upper_surr, 
                   material_a   = mat_rib, material_b   = mat_rib, material_c   = mat_substrate,
                   lc_bkg=.1, lc_refine_1=10.0, lc_refine_2=20.0) 

wguide.plot_mesh(prefix)
wguide.check_mesh()


# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1                            


recalc_fields=True     # run the calculation from scratch

# Calculate Electromagnetic modes.
if recalc_fields:
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff=n_eff)
    sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

    #sim_EM_pump.save_simulation('tut_07_pump')
    #sim_EM_Stokes.save_simulation('tut_07_stokes')
else:
    sim_EM_pump = mode_calcs.load_simulation('tut_07_pump')
    sim_EM_Stokes = mode_calcs.load_simulation('tut_07_stokes')

#plotting.plot_mode_fields(sim_EM_pump, quiver_points = 50, EM_AC='EM_E', ivals=range(6))
sim_EM_pump.plot_modes(quiver_points = 50, ivals=range(6))

# Display the wavevectors of EM modes.
v_kz=sim_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz): print('{0:3d}  {1:.4e}'.format(i, np.real(kz)))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff(0))
print("n_eff = ", np.round(n_eff_sim, 4))



q_AC = np.real(sim_EM_pump.kz_EM(EM_ival_pump) - sim_EM_Stokes.kz_EM(EM_ival_Stokes))
#q_AC = 14.37e6   # that is when we want to use a fixed wavenumber as in an example that we tried to reproduce, to get closer to the results that they had 

# Specify the expected acoustic frequency (chosen slightly lower than likely resonances).
shift_Hz = 6e9                   

# Calculate Acoustic modes.
if recalc_fields:
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)


# Print the frequencies of AC modes.
v_nu=sim_AC.nu_AC_all()
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(v_nu): print(f'{i:3d}  {np.real(nu)*1e-9:.5f}')


#plotting.plot_mode_fields(sim_AC, quiver_points=50, ivals=range(40),
#                          xlim_min=-.1, ylim_min=-.1, xlim_max=-.1, ylim_max=-.1)
sim_AC.plot_modes(quiver_points=50, ivals=range(40), xlim_min=-.1, ylim_min=-.1, xlim_max=-.1, ylim_max=-.1)



# set_q_factor = 1000.
set_q_factor = 355.

gainbox = integration.get_gains_and_qs( sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)

gainbox.set_allowed_EM_pumps(EM_ival_pump)
gainbox.set_allowed_EM_Stokes(EM_ival_Stokes)
gainbox.set_EM_modes(EM_ival_pump, EM_ival_Stokes)

SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)


# reshaping numpy arrays 
v_nu_reshape = v_nu.reshape(1, -1)
v_nu_reshape = v_nu_reshape.reshape(-1, 1)

SBS_gain_compressed = SBS_gain[0,0,:]
SBS_gain_PE_compressed = SBS_gain_PE[0,0,:]
SBS_gain_MB_compressed = SBS_gain_MB[0,0,:]

SBS_gain_compressed = SBS_gain_compressed.reshape(1, -1)
SBS_gain_compressed = SBS_gain_compressed.reshape(-1, 1)
SBS_gain_PE_compressed = SBS_gain_PE_compressed.reshape(1, -1)
SBS_gain_PE_compressed = SBS_gain_PE_compressed.reshape(-1, 1)
SBS_gain_MB_compressed = SBS_gain_MB_compressed.reshape(1, -1)
SBS_gain_MB_compressed = SBS_gain_MB_compressed.reshape(-1, 1)

# saving data as mat-files
results = {"header":prefix}
results["SBS_gain"] = SBS_gain_compressed
results["SBS_gain_PE"] = SBS_gain_PE_compressed
results["SBS_gain_MB"] = SBS_gain_MB_compressed
results["Acoustic_frequencies"] = v_nu_reshape
results["linewidth_Hz"] = linewidth_Hz
sc.savemat(prefix + ".mat", results)

# saving data as csv-files
results = np.concatenate((v_nu_reshape, SBS_gain_compressed, SBS_gain_PE_compressed, SBS_gain_MB_compressed), axis=1)
np.savetxt(prefix + ".csv", results, delimiter=',', newline='\n', header='Acoustic freq, SBS gain, SBS PE gain, SBS MB gain')


# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
#freq_min = np.real(sim_AC.nu_AC_all()[0]) - 2e9  # GHz
#freq_max = np.real(sim_AC.nu_AC_all()[-1]) + 2e9  # GHz
freq_min = 8.5e9
freq_max = 9e9


(nu_gain_tot, nu_gain_PE, nu_gain_MB) = gainbox.plot_spectra(
        freq_min, freq_max, semilogy=True,
        prefix = prefix)

#plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
#    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
#    )

print(nbapp.final_report())


