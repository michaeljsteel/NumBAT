""" Calculate the backward SBS gain spectra of a
    silicon waveguide surrounded in air.

    Show how to save simulation objects 
    (eg. EM mode calcs) to expedite the process 
    of altering later parts of simulations.

    Show how to implement integrals in python
    and how to load data from Comsol.
"""

import time
import datetime
import numpy as np
import sys

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
unitcell_x = 2.5*wl_nm
unitcell_y = unitcell_x
inc_a_x = 300
inc_a_y = 280
inc_shape = 'rectangular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 25
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

if len(sys.argv)>1 and sys.argv[1]=='fast=1':  # choose between faster or more accurate calculation
  prefix_str = 'ftut_02-'
  refine_fac=1
else:
  prefix_str = 'tut_02-'
  refine_fac=5

print('\n\nCommencing NumBAT tutorial 2')

# Use of a more refined mesh to produce field plots.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        material_bkg=materials.get_material("Vacuum"),
                        material_a=materials.get_material("Si_2016_Smith"),
                        lc_bkg=1, lc_refine_1=120.0*refine_fac, lc_refine_2=60.0*refine_fac)


# Estimate expected effective index of fundamental guided mode.
n_eff = wguide.material_a.n-0.1

recalc_fields=True     # run the calculation from scratch
recalc_fields=False   # reuse saved fields from previous calculation

if recalc_fields:
  # Calculate Electromagnetic modes.
  sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff)
  # Save calculated :Simmo: object for EM calculation.
  np.savez('wguide_data', sim_EM_pump=sim_EM_pump)
  sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
  np.savez('wguide_data2', sim_EM_Stokes=sim_EM_Stokes)
else:
  # Once npz files have been saved from one simulation run,
  # set recalc_fields=True to use the saved data
  #This provides precisely the same objects  for the remainder of the simulation.
  npzfile = np.load('wguide_data.npz', allow_pickle=True)
  sim_EM_pump = npzfile['sim_EM_pump'].tolist()
  npzfile = np.load('wguide_data2.npz', allow_pickle=True)
  sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

# Print the wavevectors of EM modes.
v_kz=sim_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz): print('{0:3d}  {1:.4e}'.format(i, np.real(kz)))

# Plot the E fields of the EM modes fields - specified with EM_AC='EM_E'.
# Zoom in on the central region (of big unitcell) with xlim_, ylim_ args,
# which specify the fraction of the axis to remove from the plot.
# For instance xlim_min=0.4 will remove 40% of the x axis from the left outer edge
# to the center. xlim_max=0.4 will remove 40% from the right outer edge towards the center.
# This leaves just the inner 20% of the unit cell displayed in the plot.
# The ylim variables perform the equivalent actions on the y axis.

# Let's plot fields for only the fundamental (ival = 0) mode.
#decorator=plotting.Decorator()
#decorator.set_multiplot_axes_property('subplots_wspace',.4)

plotting.plot_mode_fields(sim_EM_pump, xlim_min=0.4, xlim_max=0.4, ylim_min=0.4,
                         ylim_max=0.4, ivals=[EM_ival_pump], contours=True, EM_AC='EM_E', 
                         prefix_str=prefix_str, ticks=True, suppress_imimre=True)

#Repeat this plot in pdf output format
plotting.plot_mode_fields(sim_EM_pump, xlim_min=0.4, xlim_max=0.4, ylim_min=0.4,
                         ylim_max=0.4, ivals=[EM_ival_pump], contours=True, EM_AC='EM_E', 
                         pdf_png='pdf', prefix_str=prefix_str, ticks=True, suppress_imimre=True)

# Plot the H fields of the EM modes - specified with EM_AC='EM_H'.
plotting.plot_mode_fields(sim_EM_pump, xlim_min=0.4, xlim_max=0.4, ylim_min=0.4,
                         ylim_max=0.4, ivals=[EM_ival_pump], EM_AC='EM_H', 
                         prefix_str=prefix_str, ticks=True, suppress_imimre=True)

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff(0))
print("n_eff", np.round(n_eff_sim, 4))

# Acoustic wavevector
k_AC = np.real(sim_EM_pump.kz_EM(0) - sim_EM_Stokes.kz_EM(0))

if recalc_fields:
  # Calculate Acoustic modes.
  sim_AC = wguide.calc_AC_modes(num_modes_AC, k_AC, EM_sim=sim_EM_pump)
  # # Save calculated :Simmo: object for AC calculation.
  np.savez('wguide_data_AC', sim_AC=sim_AC)
else:
  npzfile = np.load('wguide_data_AC.npz', allow_pickle=True)
  sim_AC = npzfile['sim_AC'].tolist()

# Print the frequencies of AC modes.
v_nu=sim_AC.nu_AC_all()
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(v_nu): print('{0:3d}  {1:.4e}'.format(i, np.real(nu)*1e-9))

# Plot the AC modes fields, important to specify this with EM_AC='AC'.
# The AC modes are calculated on a subset of the full unitcell,
# which excludes vacuum regions, so there is usually no need to restrict the area plotted
# with xlim_min, xlim_max etc.
plotting.plot_mode_fields(sim_AC, EM_AC='AC', contours=False, prefix_str=prefix_str, 
    ticks=True, suppress_imimre=True, quiver_points=20, ivals=[0])

if recalc_fields:
  # Calculate the acoustic loss from our fields.
  # Calculate interaction integrals and SBS gain for PE and MB effects combined, 
  # as well as just for PE, and just for MB.
  SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC, EM_ival_pump=EM_ival_pump, 
    EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)
  # Save the gain calculation results
  np.savez('wguide_data_AC_gain', SBS_gain=SBS_gain, SBS_gain_PE=SBS_gain_PE, 
            SBS_gain_MB=SBS_gain_MB, linewidth_Hz=linewidth_Hz)
else:
  npzfile = np.load('wguide_data_AC_gain.npz', allow_pickle=True)
  SBS_gain = npzfile['SBS_gain']
  SBS_gain_PE = npzfile['SBS_gain_PE']
  SBS_gain_MB = npzfile['SBS_gain_MB']
  linewidth_Hz = npzfile['linewidth_Hz']

# The following function shows how integrals can be implemented purely in python,
# which may be of interest to users wanting to calculate expressions not currently
# included in NumBAT. Note that the Fortran routines are much faster!
# Also shows how field data can be imported (in this case from Comsol) and used.
comsol_ivals = 5 # Number of modes contained in data file.
SBS_gain_PE_py, alpha_py, SBS_gain_PE_comsol, alpha_comsol = integration.gain_python(
    sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC, 'Comsol_ac_modes_1-5.dat', 
    comsol_ivals=comsol_ivals)

# Print the PE contribution to gain SBS gain of the AC modes.
print("\n Displaying results with negligible components masked out")
# Mask negligible gain values to improve clarity of print out.
threshold = -1e-3
masked_PE = np.ma.masked_inside(SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,:comsol_ivals], 0, threshold)
print("SBS_gain [1/(Wm)] PE NumBAT default (Fortran)\n", masked_PE)
masked = np.ma.masked_inside(SBS_gain_PE_py[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
print("SBS_gain [1/(Wm)] python integration routines \n", masked)
masked = np.ma.masked_inside(SBS_gain_PE_comsol[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
print("SBS_gain [1/(Wm)] from loaded Comsol data \n", masked)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = np.real(sim_AC.Eig_values[0])*1e-9 - 2  # GHz
freq_max = np.real(sim_AC.Eig_values[-1])*1e-9 + 2  # GHz
plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max, 
    prefix_str=prefix_str)

# Repeat this plot focusing on one frequency range 
freq_min = 12  # GHz
freq_max = 14  # GHz
plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max, 
    prefix_str=prefix_str, suffix_str='_zoom')

end = time.time()
print("\n Simulation time (sec.)", (end - start))

