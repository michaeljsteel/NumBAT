""" Calculate the backward SBS gain spectra of a Si
    slot waveguide containing As2S3 on a SiO2 slab.
"""

import sys
import numpy as np

sys.path.append("../backend/")
import numbat
import materials
import objects
import mode_calcs
import integration
import plotting

import starter



# Geometric Parameters - all in nm.
lambda_nm = 1550
unitcell_x = 4*lambda_nm
unitcell_y = 0.3*unitcell_x
inc_shape = 'slot'
inc_a_x = 150
inc_a_y = 190
inc_b_x = 250
# Current mesh template assume inc_b_y = inc_a_y
slab_a_x = 2000
slab_a_y = 100

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix, refine_fac = starter.read_args(7, sys.argv)

nbapp = numbat.NumBATApp()

wguide = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        slab_a_x=slab_a_x, slab_a_y=slab_a_y, inc_b_x=inc_b_x,
                        material_bkg=materials.make_material("Vacuum"),            # background
                        material_a=materials.make_material("As2S3_2017_Morrison"), # slot
                        material_b=materials.make_material("SiO2_2013_Laude"),     # slab
                        material_c=materials.make_material("Si_2016_Smith"),       # walls of slot
                        lc_bkg=.1, lc_refine_1=15.0*refine_fac, lc_refine_2=15.0*refine_fac)

wguide.plot_mesh(prefix)
# Move origin to nominal centre of waveguide, y adjustment needs refinement based on mesh template.
# Note different shifts are allowed for EM and acoustic, because acoustic domain excludes vacuum regions
# Shifts are in nm
wguide.set_xyshift_em(-unitcell_x*.5, unitcell_y*.5)
wguide.set_xyshift_ac(-unitcell_x*.5, unitcell_y*.5)

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

recalc_fields=True     # run the calculation from scratch
#recalc_fields=False   # reuse saved fields from previous calculation

# Calculate Electromagnetic modes.
if recalc_fields:
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff=n_eff)
    sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

    sim_EM_pump.save_simulation('tut_07_pump')
    sim_EM_Stokes.save_simulation('tut_07_stokes')
else:
    sim_EM_pump = mode_calcs.load_simulation('tut_07_pump')
    sim_EM_Stokes = mode_calcs.load_simulation('tut_07_stokes')

plotting.plot_mode_fields(sim_EM_pump, quiver_points = 20, xlim_min=0.2, xlim_max=0.2,
                           ylim_min=0.0, ylim_max=0.0, EM_AC='EM_E',
                           prefix=prefix)

# Display the wavevectors of EM modes.
v_kz=sim_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz): print('{0:3d}  {1:.4e}'.format(i, np.real(kz)))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff(0))
print("n_eff = ", np.round(n_eff_sim, 4))

q_AC = np.real(sim_EM_pump.kz_EM(EM_ival_pump) - sim_EM_Stokes.kz_EM(EM_ival_Stokes))

# Specify the expected acoustic frequency (chosen slightly lower than likely resonances).
shift_Hz = 4e9

# Calculate Acoustic modes.
if recalc_fields:
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)
    sim_AC.save_simulation('tut_07_acoustic')
else:
    sim_AC = mode_calcs.load_simulation('tut_07_acoustic')

plotting.plot_mode_fields(sim_AC, quiver_points=20, prefix=prefix)

# Print the frequencies of AC modes.
v_nu=sim_AC.nu_AC_all()
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(v_nu): print('f{i:3d}  {np.real(nu)*1e-9:.5f}')

set_q_factor = 1000.

SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)
# np.savez('wguide_data_AC_gain', SBS_gain=SBS_gain, SBS_gain_PE=SBS_gain_PE, SBS_gain_MB=SBS_gain_MB, alpha=alpha)
# npzfile = np.load('wguide_data_AC_gain.npz', allow_pickle=True)
# SBS_gain = npzfile['SBS_gain']
# SBS_gain_PE = npzfile['SBS_gain_PE']
# SBS_gain_MB = npzfile['SBS_gain_MB']
# alpha = npzfile['alpha']

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = np.real(sim_AC.nu_AC_all()[0]) - 2e9  # GHz
freq_max = np.real(sim_AC.nu_AC_all()[-1]) + 2e9  # GHz

plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    prefix=prefix)

print(nbapp.final_report())