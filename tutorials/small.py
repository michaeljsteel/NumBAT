"""
    NumBAT Tutorial 2

    Calculate the backward SBS gain spectra of a silicon waveguide surrounded in air.

    Show how to save simulation objects (eg. EM mode calcs) to expedite the process
    of altering later parts of simulations.

    Show how to implement integrals in python and how to load data from Comsol.
"""

import sys
import numpy as np

from pathlib import Path
sys.path.append(str(Path('../backend')))

import numbat
import integration
import mode_calcs
import materials

# Geometric Parameters - all in nm.
lambda_nm = 1550.0  # Wavelength of EM wave in vacuum.

# Unit cell must be large to ensure fields are near-zero at boundary.
domain_x = 2000.0
domain_y = domain_x

# Waveguide widths.
inc_a_x = 500.0
inc_a_y = 600.0



inc_shape = 'rectangular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 25
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix = 'tt_small'
refine_fac = 1

print('\nCommencing NumBAT tutorial 2\n')

nbapp = numbat.NumBATApp(prefix)

# Use of a more refined mesh to produce field plots.
wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                           material_bkg=materials.make_material("Vacuum"),
                           material_a=materials.make_material("Si_2016_Smith"),
                           lc_bkg=.5, lc_refine_1=1.0*refine_fac, lc_refine_2=1.0*refine_fac)

#wguide.plot_mesh(prefix)
#wguide.plot_mail_mesh(prefix)


# Estimate expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

# Calculate Electromagnetic modes.
reuse_old_fields = False
if reuse_old_fields:
    print('\nLoading EM fields')
    simres_EM_pump = numbat.load_simulation('tut02_em_pump')
    simres_EM_Stokes = numbat.load_simulation('tut02_em_stokes')
else:
    simres_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff)
    simres_EM_Stokes = mode_calcs.bkwd_Stokes_modes(simres_EM_pump)
    print('\nSaving EM fields')
    simres_EM_pump.save_simulation('tut02_em_pump')
    simres_EM_Stokes.save_simulation('tut02_em_stokes')

# Print the wavevectors of EM modes.
v_kz = simres_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz):
    print(f'{i:3d}  {np.real(kz):.4e}')


# Find the EM effective index of the waveguide.
n_eff_sim = np.real(simres_EM_pump.neff(0))
print(f'\nThe fundamental optical mode has effective index n_eff = {n_eff_sim:.6f}')

# Zoom in on the central region (of big unitcell) with xlim_, ylim_ args,
# which specify the fraction of the axis to remove from the plot.
# For instance xlim_min=0.2 will remove 20% of the x axis from the left outer edge
# to the center. xlim_max=0.2 will remove 20% from the right outer edge towards the center.
# This leaves just the inner 60% of the unit cell displayed in the plot.
# The ylim variables perform the equivalent actions on the y axis.

# Let's plot fields for only the first few modes (ival=range(4)=0--3)


simres_EM_pump.get_mode(0).plot_mode_raw_fem(['x','y'])

print('\nPlotting EM fields')
# Plot the E field of the pump mode
#simres_EM_pump.plot_modes(xlim_min=0.2, xlim_max=0.2, ylim_min=0.2,
#                          ylim_max=0.2, ivals=range(4))

# Plot the H field of the pump mode
#simres_EM_pump.plot_modes(xlim_min=0.2, xlim_max=0.2, ylim_min=0.2,
#                          ylim_max=0.2, ivals=range(4),
#                          field_type='EM_H')

# Acoustic wavevector
q_AC = np.real(simres_EM_pump.kz_EM(0) - simres_EM_Stokes.kz_EM(0))


if reuse_old_fields:
    simres_AC = numbat.load_simulation('tut_02_ac')
else:
    # Calculate and save acoustic modes.
    simres_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=simres_EM_pump)
    print('Saving AC fields')
    simres_AC.save_simulation('tut_02_ac')

# Print the frequencies of AC modes.
v_nu = simres_AC.nu_AC_all()
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(v_nu):
    print(f'{i:3d}  {np.real(nu)*1e-9:.5f}')

# The AC modes are calculated on a subset of the full unitcell,
# which excludes vacuum regions, so there is usually no need to restrict the area plotted
# with xlim_min, xlim_max etc.

print('\nPlotting acoustic modes')
#simres_AC.plot_modes(contours=True, quiver_points=20, ivals=range(10))

#if reuse_old_fields:
# Calculate the acoustic loss from our fields.
# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
gain = integration.get_gains_and_qs(
    simres_EM_pump, simres_EM_Stokes, simres_AC, q_AC, EM_ival_pump=EM_ival_pump,
    EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)
# Save the gain calculation results
#np.savez('tut02_wguide_data_AC_gain', SBS_gain=gain)
#else:
 #   npzfile = np.load('wguide_data_AC_gain.npz', allow_pickle=True)
  #  gain = npzfile['SBS_gain']

# The following function shows how integrals can be implemented purely in python,
# which may be of interest to users wanting to calculate expressions not currently
# included in NumBAT. Note that the Fortran routines are much faster!
# Also shows how field data can be imported (in this case from Comsol) and used.
comsol_ivals = 5  # Number of modes contained in data file.
#SBS_gain_PE_py, alpha_py, SBS_gain_PE_comsol, alpha_comsol = integration.gain_python(
#    simres_EM_pump, simres_EM_Stokes, simres_AC, q_AC, 'comsol_ac_modes_1-5.dat',
#    comsol_ivals=comsol_ivals)

# Print the PE contribution to gain SBS gain of the AC modes.
print("\n Displaying results of first five modes with negligible components masked out")
SBS_gain_PE = gain.gain_PE_all()
SBS_gain_MB = gain.gain_MB_all()

# Mask negligible gain values to improve clarity of print out.
threshold = -1e-3
masked_PE = np.ma.masked_inside(SBS_gain_PE[:comsol_ivals], 0, threshold)
masked_MB = np.ma.masked_inside(SBS_gain_MB[:comsol_ivals], 0, threshold)
np.set_printoptions(precision=4)
print("SBS_gain [1/(Wm)] PE NumBAT default (Fortran)\n", masked_PE)
print("SBS_gain [1/(Wm)] MB NumBAT default (Fortran)\n", masked_MB)

#masked = np.ma.masked_inside(
#    SBS_gain_PE_py[EM_ival_pump, EM_ival_Stokes, :], 0, threshold)
#print("SBS_gain [1/(Wm)] python integration routines \n", masked)
#masked = np.ma.masked_inside(
#    SBS_gain_PE_comsol[EM_ival_pump, EM_ival_Stokes, :], 0, threshold)
#print("SBS_gain [1/(Wm)] from loaded Comsol data \n", masked)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = np.real(simres_AC.nu_AC_all()[0]) - 2e9  # Hz
freq_max = np.real(simres_AC.nu_AC_all()[-1]) + 2e9  # Hz
gain.plot_spectra(freq_min=freq_min, freq_max=freq_max, dB=True, logy=True)

# Repeat this plot focusing on one frequency range
freq_min = 11.5e9  # Hz
freq_max = 13.5e9  # Hz
gain.plot_spectra(freq_min=freq_min, freq_max=freq_max, suffix='_zoom')

print(nbapp.final_report())
