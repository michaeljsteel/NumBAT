""" Example showing how the 'onion2' geometry template can be used
    to simulate a circular Si waveguide clad in SiO2 surrounded by air.

    Uses the two-layer onion structure which with vacuum surrounds makes a three material structure.

    With a silicon core and silica cladding, we do not expect this structure to guide elastic
    waves in the core.

    ***** Something is broken in assigning material properties at the central grid point.
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
domain_x = 2500
domain_y = domain_x
inc_a_x = 700
inc_b_x = 350
inc_shape = 'onion2'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 80
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix, refine_fac = starter.read_args(11, sys.argv, sub='a')

nbapp = numbat.NumBATApp(prefix)

mat_bkg=materials.make_material("Vacuum")
mat_a=materials.make_material("As2S3_2021_Poulton")
mat_b=materials.make_material("SiO2_2021_Poulton")

wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x,
                           inc_b_x=inc_b_x,  # first annulus width
                           material_bkg=mat_bkg, material_a=mat_a, material_b=mat_b,
                           lc_bkg=.1, lc_refine_1=2.0*refine_fac, lc_refine_2=2*refine_fac)

wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

new_calcs = True

# Calculate Electromagnetic modes.
if new_calcs:
    simres_EM_pump = wguide.calc_EM_modes(
        num_modes_EM_pump, lambda_nm, n_eff, calc_EM_mode_energy=True)
    simres_EM_Stokes = simres_EM_pump.bkwd_Stokes_modes()

  # simres_EM_pump.save_simulation('tut_11_pump')
  # simres_EM_Stokes.save_simulation('tut_11_stokes')
else:
    simres_EM_pump = numbat.load_simulation('tut_11_pump')
    simres_EM_Stokes = numbat.load_simulation('tut_11_stokes')

print('EM modes:\n')
kz_EM_mu = np.real(simres_EM_pump.kz_EM_all())*1e-6
neff_EM = simres_EM_pump.neff_all()
ng_EM = simres_EM_pump.ngroup_EM_all()
print('m    |   k_z [1/micron]  | neff  | ng')
for m in range(num_modes_EM_pump):
    print(f'{m:4d}  {kz_EM_mu[m]:12.6f}  {neff_EM[m]:12.6f}  {ng_EM[m]:12.6f}')

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(simres_EM_pump.neff(0))
print("n_eff", np.round(n_eff_sim, 4))

# # Plot the E fields of the EM modes fields - specified with field_type='EM_E'.
# # Zoom in on the central region (of big unitcell) with xlim_, ylim_ args.
# # Only plot fields of fundamental (ival = 0) mode.
# npzfile = np.load('wguide_data.npz', allow_pickle=True)
# npzfile = np.load('wguide_data2.npz', allow_pickle=True)
# simres_EM_Stokes = npzfile['simres_EM_Stokes'].tolist()

simres_EM_pump.plot_modes(ivals=range(10), field_type='EM_E')
simres_EM_pump.plot_modes(ivals=range(10), field_type='EM_H')

# Acoustic wavevector
q_AC = np.real(simres_EM_pump.kz_EM(EM_ival_pump) -
               simres_EM_Stokes.kz_EM(EM_ival_Stokes))

# Calculate Acoustic modes.
if new_calcs:
    simres_AC = wguide.calc_AC_modes(
        num_modes_AC, q_AC, EM_sim=simres_EM_pump, calc_AC_mode_power=True)
  # simres_AC.save_simulation('tut_11_acoustic')
else:
    simres_AC = numbat.load_simulation('tut_11_acoustic')

print('AC mode properties (GHz) \n')
nu_AC = np.real(simres_AC.nu_AC_all())*1e-9
vp_AC = np.real(simres_AC.vp_AC_all())
vg_AC = np.real(simres_AC.vg_AC_all())
print('m    |   nu [GHz]  | vp [m/s]  | vg [m/s]')
for m in range(num_modes_AC):
    print(f'{m:4d}  {nu_AC[m]:12.6f}  {vp_AC[m]:12.2f}  {vg_AC[m]:12.2f}')

#simres_AC.calc_acoustic_losses()

simres_AC.plot_modes(xlim_min=-.1, xlim_max=-.1, ylim_min=-.1, ylim_max=-.1,
                           ivals=range(num_modes_AC), quiver_points=20)

# Calculate the acoustic loss from our fields.
# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
gain_box = integration.get_gains_and_qs(
    simres_EM_pump, simres_EM_Stokes, simres_AC, q_AC, EM_ival_pump=EM_ival_pump,
    EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 4.e9
freq_max = 10.e9
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max )

print(nbapp.final_report())
