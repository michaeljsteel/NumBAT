""" Calculate the backward SBS gain for modes in a
    silicon waveguide surrounded in air.
"""

# Step 1

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

# Naming conventions
# AC: acoustic
# EM: electromagnetic
# k_AC: acoustic wavevector

start = time.time()
print('\n\nCommencing NumBAT tutorial 1')

# Step 2
# Geometric Parameters - all in nm.
wl_nm = 1550 # Wavelength of EM wave in vacuum.
# Unit cell must be large to ensure fields are zero at boundary.
unitcell_x = 2.5*wl_nm
unitcell_y = unitcell_x
# Waveguide widths.
inc_a_x = 300
inc_a_y = 280
# Shape of the waveguide.
inc_shape = 'rectangular'

# Step 3
# Number of electromagnetic modes to solve for.
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
# Number of acoustic modes to solve for.
num_modes_AC = 20
# The EM pump mode(s) for which to calculate interaction with AC modes.
# Can specify a mode number (zero has lowest propagation constant) or 'All'.
EM_ival_pump = 0
# The EM Stokes mode(s) for which to calculate interaction with AC modes.
EM_ival_Stokes = 0
# The AC mode(s) for which to calculate interaction with EM modes.
AC_ival = 'All'

# Step 4
# Use specified parameters to create a waveguide object.
# Note use of rough mesh for demonstration purposes, and use plt_mesh=True
# to save the geometry and mesh as png files in backend/fortran/msh/
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        material_bkg=materials.get_material("Vacuum"),
                        material_a=materials.get_material("Si_2016_Smith"),
                        lc_bkg=1, # in vacuum background
                        lc_refine_1=600.0, # on cylinder surfaces
                        lc_refine_2=300.0, # on cylinder center
                        plt_mesh=False)

# Explicitly remind ourselves what data we're using.
print('\nUsing material data: ', wguide.material_a)

# Step 5
# Estimate expected effective index of fundamental guided mode.
n_eff = wguide.material_a.n-0.1

# Calculate the Electromagnetic modes of the pump field.
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff)

# Display the wavevectors of EM modes.
v_kz=sim_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz): print('{0:3d}  {1:.4e}'.format(i, np.real(kz)))

# Calculate the Electromagnetic modes of the Stokes field.
# For an idealised backward SBS simulation the Stokes modes are identical 
# to the pump modes but travel in the opposite direction.
sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
# # Alt
# sim_EM_Stokes = wguide.calc_EM_modes(wl_nm, num_modes_EM_Stokes, n_eff, Stokes=True)

# Step 6
# Find the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff(0))
print("\n Fundamental optical mode ")
print(" n_eff = ", np.round(n_eff_sim, 4))

# Acoustic wavevector
k_AC = np.real(sim_EM_pump.kz_EM(0) - sim_EM_Stokes.kz_EM(0))

print('\n AC wavenumber (1/m) = ', np.round(k_AC, 4))

# Step 7
# Calculate Acoustic modes, using the mesh from the EM calculation.
sim_AC = wguide.calc_AC_modes(num_modes_AC, k_AC, EM_sim=sim_EM_pump)

# Print the frequencies of AC modes.
v_nu=sim_AC.nu_AC_all()
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(v_nu): print('{0:3d}  {1:.5f}'.format(i, np.real(nu)*1e-9))

# Do not calculate the acoustic loss from our fields, instead set a Q factor.
set_q_factor = 1000.

# Step 8
# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB. Also calculate acoustic loss alpha.
SBS_gain_tot, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC, EM_ival_pump=EM_ival_pump, 
    EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)

# SBS_gain_tot, SBS_gain_PE, SBS_gain_MB are 3D arrays indexed by pump, Stokes and acoustic mode
# Extract those of interest as a 1D array:
SBS_gain_PE_ij = SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,:]
SBS_gain_MB_ij = SBS_gain_MB[EM_ival_pump,EM_ival_Stokes,:]
SBS_gain_tot_ij = SBS_gain_tot[EM_ival_pump,EM_ival_Stokes,:]

# Print the Backward SBS gain of the AC modes.
print("\nContributions to SBS gain [1/(WM)]")
print("AC Mode number | Photoelastic (PE) | Moving boundary(MB) | Total")

for (m, gpe, gmb, gt) in zip(range(num_modes_AC), SBS_gain_PE_ij, SBS_gain_MB_ij, SBS_gain_tot_ij):
    print('{0:12d}  {1:19.6e}  {2:19.6e}  {3:16.6e}'.format(m, gpe, gmb, gt))

#with np.printoptions(formatter={ %    'float': '{: 12.4e}'.format, 'int': '{: 12d}'.format }):
#    print(np.array([ range(num_modes_AC), SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,:],
#                     SBS_gain_MB[EM_ival_pump,EM_ival_Stokes,:], SBS_gain[EM_ival_pump,EM_ival_Stokes,:]]
#                     ).T)


#print("\n SBS_gain [1/(Wm)] PE contribution \n", SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,:])
#print("SBS_gain [1/(Wm)] MB contribution \n", SBS_gain_MB[EM_ival_pump,EM_ival_Stokes,:])
#print("SBS_gain [1/(Wm)] total \n", SBS_gain[EM_ival_pump,EM_ival_Stokes,:])

# Mask negligible gain values to improve clarity of print out.
threshold = 1e-3
masked_PE = np.where(np.abs(SBS_gain_PE_ij)>threshold, SBS_gain_PE_ij, 0)
masked_MB = np.where(np.abs(SBS_gain_MB_ij)>threshold, SBS_gain_MB_ij, 0)
masked_tot = np.where(np.abs(SBS_gain_tot_ij)>threshold, SBS_gain_tot_ij, 0)

print("\n Displaying gain results with negligible components masked out:")

print("AC Mode number | Photoelastic (PE) | Moving boundary(MB) | Total")
for (m, gpe, gmb, gt) in zip( range(num_modes_AC), masked_PE, masked_MB, masked_tot):
    print('{0:12d}  {1:19.6e}  {2:19.6e}  {3:16.6e}'.format(m, gpe, gmb, gt))


end = time.time()
print("\n Simulation time (sec.)", (end - start))

