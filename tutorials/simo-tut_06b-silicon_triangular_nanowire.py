""" We've now covered most of the features of NumBAT.
    In the following tutorials we'll show how to study different
    geometries and materials.

    Calculate the backward SBS gain spectra of a silicon waveguide
    surrounded by vacuum (air).
"""

import sys
import numpy as np

sys.path.append("../backend/")
import numbat
import materials
import mode_calcs
import integration
import plotting

import starter


# Geometric Parameters - all in nm.
#lambda_nm = 1550
lambda_nm = 1000
unitcell_x = lambda_nm
unitcell_y = unitcell_x
inc_a_x = 600      # base length (always horizontal)
inc_a_y = inc_a_x  # not used
inc_b_x = 500      # displacement of peak from left end of base
inc_b_y = 500      # height of peak from base

inc_shape = 'triangular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix, refine_fac = starter.read_args(6, sys.argv, sub='b')

nbapp = numbat.NumBATApp(prefix)

refine_fac = 1;
lc_bkg = .1 * refine_fac
lc_norm = 3
lc_corner = 6

lc_norm = 1
lc_corner = 1

wguide = nbapp.make_structure(inc_shape, unitcell_x, unitcell_y, inc_a_x, inc_a_y,
                              inc_b_x = inc_b_x, inc_b_y = inc_b_y,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=materials.make_material("Si_2021_Poulton"),
                        lc_bkg=lc_bkg, lc_refine_1=lc_norm, lc_refine_2=lc_corner)

# Expected effective index of fundamental guided mode.
n_eff = 2.5

recalc_fields=True     # run the calculation from scratch
#recalc_fields=False   # reuse saved fields from previous calculation

# Calculate Electromagnetic Modes
if recalc_fields:
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff=n_eff)
    sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

    sim_EM_pump.save_simulation('tut_06_pump')
    sim_EM_Stokes.save_simulation('tut_06_stokes')
else:
    sim_EM_pump = mode_calcs.load_simulation('tut_06_pump')
    sim_EM_Stokes = mode_calcs.load_simulation('tut_06_stokes')

sim_EM_pump.set_r0_offset(0, -0.5e-9*unitcell_y)  # ensure plots identify centre as (0,0)
sim_EM_Stokes.set_r0_offset(0, -0.5e-9*unitcell_y)  # ensure plots identify centre as (0,0)

print('\nPlotting EM fields')
plotting.plot_modes(sim_EM_pump, EM_AC='EM_E', ivals=[0])

# Display the wavevectors of EM modes.
v_kz=sim_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz): print('{0:3d}  {1:.4e}'.format(i, np.real(kz)))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff(0))
print("n_eff = {0:.4e}".format(n_eff_sim))

# Acoustic wavevector
q_AC = np.real(sim_EM_pump.kz_EM(EM_ival_pump) - sim_EM_Stokes.kz_EM(EM_ival_Stokes))

shift_Hz = 4e9

# Calculate Acoustic modes.
if recalc_fields:
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)
    sim_AC.save_simulation('tut_06_acoustic')
else:
    sim_AC = mode_calcs.load_simulation('tut_06_acoustic')

sim_AC.set_r0_offset(0, -0.5e-9*unitcell_y)  # ensure plots identify centre as (0,0)
plotting.plot_modes(sim_AC, EM_AC='AC', )

# Print the frequencies of AC modes.
v_nu=sim_AC.nu_AC_all()
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(v_nu): print('{0:3d}  {1:.5f}'.format(i, np.real(nu)*1e-9))

set_q_factor = 1000.

print('\nCalculating gains')
# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
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

plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    )

print(nbapp.final_report())
