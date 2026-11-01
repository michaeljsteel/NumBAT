""" We've now covered most of the features of NumBAT.
    In the following tutorials we'll show how to study different
    geometries and materials.

    Calculate the backward SBS gain spectra of a silicon waveguide
    surrounded by vacuum (air).
"""

import sys
import numpy as np

from pathlib import Path
sys.path.append(str(Path('../../backend')))

import numbat
import materials

import integration

import starter


# Geometric Parameters - all in nm.
lambda_nm = 1550
#lambda_nm = 1000
domain_x = lambda_nm
domain_y = domain_x
basewid = 1000      # base length (always horizontal)
inc_a_y = basewid  # not used
peak_xoff = 500      # displacement of peak from left end of base
peak_ht = 500      # height of peak from base

inc_shape = 'triangular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_mode_index_pump = 0
EM_mode_index_Stokes = 0
AC_mode_index = 'All'

prefix, refine_fac = starter.read_args(6, sys.argv, sub='b')

nbapp = numbat.NumBATApp(prefix)

refine_fac = 1;
lc_bkg = .05 * refine_fac

lc_norm = 1
lc_corner = 1


wguide = nbapp.make_structure(inc_shape, domain_x, domain_y,
                              base_width=basewid, peak_height=peak_ht,
                              peak_xoff = peak_xoff,
                              material_bkg=materials.make_material("Vacuum"),
                              material_a=materials.make_material("Si_2021_Poulton"),
                              lc_bkg=lc_bkg, lc_refine_1=lc_norm, lc_refine_2=lc_corner)


wguide.plot_refractive_index_profile(prefix)

# Expected effective index of fundamental guided mode.
n_eff = 2.5

#use_old_fields=True   # reuse saved fields from previous calculation
use_old_fields=False     # run the calculation from scratch

# Calculate Electromagnetic Modes
if use_old_fields:
    simres_EM_pump = numbat.load_simulation('tut_06bpump')
    simres_EM_Stokes = numbat.load_simulation('tut_06bstokes')
else:
    simres_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff=n_eff)
    simres_EM_Stokes = simres_EM_pump.clone_as_backward_modes()

    simres_EM_pump.save_simulation('tut_06bpump')
    simres_EM_Stokes.save_simulation('tut_06bstokes')

#simres_EM_pump.set_r0_offset(0, -0.5e-9*domain_y)  # ensure plots identify centre as (0,0)
#simres_EM_Stokes.set_r0_offset(0, -0.5e-9*domain_y)  # ensure plots identify centre as (0,0)

print('\nPlotting EM fields')
simres_EM_pump.plot_modes(mode_indices=[0])

# Display the wavevectors of EM modes.
v_kz=simres_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz): print(f'{i:3d}  {np.real(kz):.4e}')

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(simres_EM_pump.neff(0))
print(f"n_eff = {n_eff_sim:.4e}")

# Acoustic wavevector
q_AC = np.real(simres_EM_pump.kz_EM(EM_mode_index_pump) - simres_EM_Stokes.kz_EM(EM_mode_index_Stokes))

shift_Hz = 4e9

# Calculate Acoustic modes.
if use_old_fields:
    simres_AC = numbat.load_simulation('tut_06bacoustic')
else:
    simres_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=simres_EM_pump, shift_Hz=shift_Hz)
    simres_AC.save_simulation('tut_06bacoustic')

#simres_AC.set_r0_offset(0, -0.5e-9*domain_y)  # ensure plots identify centre as (0,0)

# Print the frequencies of AC modes.
v_nu=simres_AC.nu_AC_all()
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(v_nu): print(f'{i:3d}  {np.real(nu) * 1e-09:.5f}')

#simres_AC.plot_modes()
v_widy = np.zeros(num_modes_AC)
v_w2 = np.zeros(num_modes_AC)
v_x0 = np.zeros(num_modes_AC)
v_y0 = np.zeros(num_modes_AC)

nm_to_um = 0.001
m_likely_wedge = 0
peak_ht_um = peak_ht*nm_to_um
m_apex_dist=5*peak_ht_um # comfortably large to be displaced. Convert peak_ht to microns.
for m in range(num_modes_AC):
    md = simres_AC.get_mode(m)
    md = simres_AC.mode_set[m]
    md.set_width_r0_reference(0, peak_ht_um/2) # set width reference to top of triangle # this part done in m!
    md.analyse_mode()
    v_widy[m] = md.wy()
    v_w2[m] = md.w0()
    v_x0[m] = md.center_of_mass_x()
    v_y0[m] = md.center_of_mass_y()

    # look for mode with energy most concentrated near the apex
    if abs(peak_ht_um-v_y0[m])<m_apex_dist:
        m_apex_dist = abs(peak_ht_um-v_y0[m])
        m_likely_wedge = m


    print(f'mode {m:2d}: r0=({md.center_of_mass_x(): .4f}, {md.center_of_mass_y(): .4f}), '
    + f'wid= ({md.wx():.4f},{md.wy():.4f},{md.w0():.4f})')

print(f'The apex wedge mode is most likely  mode {m_likely_wedge}')

set_q_factor = 1000.


print('\nCalculating gains')
# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
gain_box = integration.get_gains_and_qs(
    simres_EM_pump, simres_EM_Stokes, simres_AC, q_AC,
    EM_mode_index_pump=EM_mode_index_pump, EM_mode_index_Stokes=EM_mode_index_Stokes, AC_mode_index=AC_mode_index, fixed_Q=set_q_factor)

# np.savez('wguide_data_AC_gain', SBS_gain=SBS_gain, SBS_gain_PE=SBS_gain_PE, SBS_gain_MB=SBS_gain_MB, alpha=alpha)
# npzfile = np.load('wguide_data_AC_gain.npz')
# SBS_gain = npzfile['SBS_gain']
# SBS_gain_PE = npzfile['SBS_gain_PE']
# SBS_gain_MB = npzfile['SBS_gain_MB']
# alpha = npzfile['alpha']

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5e9  # Hz
freq_max = 42e9  # Hz

gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True)

print(nbapp.final_report())
