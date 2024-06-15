""" We've now covered most of the features of NumBAT.
    In the following tutorials we'll show how to study different
    geometries and materials.

    Calculate the backward SBS gain spectra of a silicon waveguide
    surrounded by vacuum (air).
"""

import sys
import numpy as np

from pathlib import Path
sys.path.append(str(Path('../backend')))

import numbat
import materials
import mode_calcs
import integration
import plotting

import starter


# Geometric Parameters - all in nm.
lambda_nm = 1550
unitcell_x = 5*lambda_nm
unitcell_y = unitcell_x
inc_a_x = 550
inc_a_y = inc_a_x
inc_shape = 'circular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix, refine_fac = starter.read_args(6, sys.argv)

nbapp = numbat.NumBATApp(prefix)

wguide = nbapp.make_structure(inc_shape, unitcell_x, unitcell_y, inc_a_x, inc_a_y,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=materials.make_material("SiO2_2016_Smith"),
                        #lc_bkg=.1, lc_refine_1=12.0*refine_fac, lc_refine_2=4.0*refine_fac)
                        #lc_bkg=.1, lc_refine_1=4.0*refine_fac, lc_refine_2=4.0*refine_fac)
                        lc_bkg=.1, lc_refine_1=3.0*refine_fac, lc_refine_2=3.0*refine_fac)

# Expected effective index of fundamental guided mode.
n_eff = 1.4

recalc_fields=True     # run the calculation from scratch
#recalc_fields=False   # reuse saved fields from previous calculation

# Calculate Electromagnetic Modes
if recalc_fields:
    simres_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff=n_eff)
    simres_EM_Stokes = mode_calcs.bkwd_Stokes_modes(simres_EM_pump)

    #simres_EM_pump.save_simulation('tut_06_pump')
    #simres_EM_Stokes.save_simulation('tut_06_stokes')
else:
    simres_EM_pump = mode_calcs.load_simulation('tut_06_pump')
    simres_EM_Stokes = mode_calcs.load_simulation('tut_06_stokes')

print('\nPlotting EM fields')
trim=0.4
#trim=0.0

# Display the wavevectors of EM modes.
v_kz=simres_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz): print(f'{i:3d}  {np.real(kz):.4e}')

simres_EM_pump.plot_modes(ivals=range(5), xlim_min=trim, xlim_max=trim, ylim_min=trim, ylim_max=trim)

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(simres_EM_pump.neff(0))
print(f"n_eff = {n_eff_sim:.4e}")

# Acoustic wavevector
q_AC = np.real(simres_EM_pump.kz_EM(EM_ival_pump) - simres_EM_Stokes.kz_EM(EM_ival_Stokes))

shift_Hz = 4e9

# Calculate Acoustic modes.
if recalc_fields:
    simres_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=simres_EM_pump, shift_Hz=shift_Hz)
    #simres_AC.save_simulation('tut_06_acoustic')
else:
    simres_AC = mode_calcs.load_simulation('tut_06_acoustic')


# Print the frequencies of AC modes.
v_nu=simres_AC.nu_AC_all()
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(v_nu): print(f'{i:3d}  {np.real(nu)*1e-9:.5f}')

simres_AC.set_r0_offset(0, -0.5e-9*unitcell_y)  # ensure plots identify centre as (0,0)
simres_AC.plot_modes()

set_q_factor = 1000.

print('\nCalculating gains')
# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    simres_EM_pump, simres_EM_Stokes, simres_AC, q_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)

# np.savez('wguide_data_AC_gain', SBS_gain=SBS_gain, SBS_gain_PE=SBS_gain_PE, SBS_gain_MB=SBS_gain_MB, alpha=alpha)
# npzfile = np.load('wguide_data_AC_gain.npz')
# SBS_gain = npzfile['SBS_gain']
# SBS_gain_PE = npzfile['SBS_gain_PE']
# SBS_gain_MB = npzfile['SBS_gain_MB']
# alpha = npzfile['alpha']

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5e9  # Hz
freq_max = 12e9  # Hz

plotting.plot_gain_spectra(simres_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    )

print(nbapp.final_report())
