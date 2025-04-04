""" Sanity check implementation of fully anisotropic
    tensors by feeding in same parameters of simo_tut_06.
"""

import sys
import numpy as np

from pathlib import Path
sys.path.append(str(Path('../backend')))

import numbat
import materials

import integration

import starter


# Geometric Parameters - all in nm.
lambda_nm = 1550
domain_x = 2.0*lambda_nm
domain_y = domain_x
#inc_a_x = 314.7
#inc_a_y = 0.9*inc_a_x
inc_a_x = 300
inc_a_y = 280

inc_shape = 'rectangular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 20
EM_mode_index_pump = 0
EM_mode_index_Stokes = 0
AC_mode_index = 'All'

prefix, refine_fac = starter.read_args(9, sys.argv)

nbapp = numbat.NumBATApp(prefix)

# Use of a more refined mesh to produce field plots.
wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=materials.make_material("Si_test_anisotropic"),
                        lc_bkg=1, lc_refine_1=200.0*refine_fac, lc_refine_2=1.0*refine_fac)


# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

# Calculate the Electromagnetic modes of the pump field.
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff)

# Display the wavevectors of EM modes.
v_kz=sim_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz): print(f'{i:3d}  {np.real(kz):.4e}')

# Calculate the Electromagnetic modes of the Stokes field.
# For an idealised backward SBS simulation the Stokes modes are identical
# to the pump modes but travel in the opposite direction.
sim_EM_Stokes = sim_EM_pump.clone_as_backward_modes()
# # Alt
# sim_EM_Stokes = wguide.calc_EM_modes(lambda_nm, num_modes_EM_Stokes, n_eff, Stokes=True)

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff(0))
print("\n Fundamental optical mode ")
print(" n_eff = ", np.round(n_eff_sim, 4))

# Acoustic wavevector
q_AC = np.real(sim_EM_pump.kz_EM(EM_mode_index_pump) - sim_EM_Stokes.kz_EM(EM_mode_index_Stokes))
print('\n AC wavenumber (1/m) = ', np.round(q_AC, 4))

# Calculate Acoustic modes, using the mesh from the EM calculation.
sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump)

# Print the frequencies of AC modes.
v_nu=sim_AC.nu_AC_all()
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(v_nu): print(f'{i:3d}  {np.real(nu)*1e-9:.4e}')


# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB. Also calculate acoustic loss alpha.
gain_box = integration.get_gains_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_mode_index_pump=EM_mode_index_pump, EM_mode_index_Stokes=EM_mode_index_Stokes, AC_mode_index=AC_mode_index)

gain_box.plot_spectra(freq_min=0, freq_max=30e9, )

# SBS_gain_tot, SBS_gain_PE, SBS_gain_MB are 3D arrays indexed by pump, Stokes and acoustic mode
# Extract those of interest as a 1D array:
SBS_gain_PE_ij = gain_box.gain_PE_all()
SBS_gain_MB_ij = gain_box.gain_MB_all()
SBS_gain_tot_ij = gain_box.gain_total_all()


# Print the Backward SBS gain of the AC modes.
print("\nContributions to SBS gain [1/(WM)]")
print("AC Mode number | Photoelastic (PE) | Moving boundary(MB) | Total")

for (m, gpe, gmb, gt) in zip(range(num_modes_AC), SBS_gain_PE_ij, SBS_gain_MB_ij, SBS_gain_tot_ij):
    print(f'{m:12d}  {gpe:19.6e}  {gmb:19.6e}  {gt:16.6e}')

# Mask negligible gain values to improve clarity of print out.
threshold = 1e-3
masked_PE = np.where(np.abs(SBS_gain_PE_ij)>threshold, SBS_gain_PE_ij, 0)
masked_MB = np.where(np.abs(SBS_gain_MB_ij)>threshold, SBS_gain_MB_ij, 0)
masked_tot = np.where(np.abs(SBS_gain_tot_ij)>threshold, SBS_gain_tot_ij, 0)

print("\n Displaying gain results with negligible components masked out:")

print("AC Mode number | Photoelastic (PE) | Moving boundary(MB) | Total")
for (m, gpe, gmb, gt) in zip( range(num_modes_AC), masked_PE, masked_MB, masked_tot):
    print(f'{m:12d}  {gpe:19.6e}  {gmb:19.6e}  {gt:16.6e}')


print(nbapp.final_report())
