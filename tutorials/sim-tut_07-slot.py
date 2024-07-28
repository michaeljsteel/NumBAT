"""Calculate the backward SBS gain spectra of a Si
slot waveguide containing As2S3 on a SiO2 slab.
"""

import sys
import numpy as np

from pathlib import Path

sys.path.append(str(Path("../backend")))

import numbat
import materials
import mode_calcs
import integration
from nbtypes import SI_GHz

import starter


# Geometric Parameters - all in nm.
lambda_nm = 1550
domain_x = 4 * lambda_nm
domain_y = 0.3 * domain_x
inc_shape = "slot"
inc_a_x = 150
inc_a_y = 190
inc_b_x = 250
# Current mesh template assume inc_b_y = inc_a_y
slab_a_x = 2000
slab_a_y = 100

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_ival_pump = 1
EM_ival_Stokes = 1
AC_ival = "All"

prefix, refine_fac = starter.read_args(7, sys.argv)

nbapp = numbat.NumBATApp(prefix)

mat_vac=materials.make_material("Vacuum")  # background
mat_slot=materials.make_material("As2S3_2017_Morrison")  # slot
mat_slab=materials.make_material("SiO2_2013_Laude")  # slab
mat_walls=materials.make_material("Si_2016_Smith")  # walls of slot


wguide = nbapp.make_structure(inc_shape, domain_x, domain_y,
    inc_a_x, inc_a_y, slab_a_x=slab_a_x, slab_a_y=slab_a_y, inc_b_x=inc_b_x,
    material_bkg=mat_vac, material_a=mat_slot, material_b=mat_slab, material_c=mat_walls,
    lc_bkg=0.05, lc_refine_1=7.0 * refine_fac, lc_refine_2=7.0 * refine_fac,
)

# wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material("a").refindex_n - 0.1

recalc_fields = True  # run the calculation from scratch
# recalc_fields=False   # reuse saved fields from previous calculation

# Calculate Electromagnetic modes.
if recalc_fields:
    simres_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff=n_eff)
    simres_EM_Stokes = mode_calcs.bkwd_Stokes_modes(simres_EM_pump)

    simres_EM_pump.save_simulation("tut_07_pump")
    simres_EM_Stokes.save_simulation("tut_07_stokes")
else:
    simres_EM_pump = mode_calcs.load_simulation("tut_07_pump")
    simres_EM_Stokes = mode_calcs.load_simulation("tut_07_stokes")

simres_EM_pump.plot_modes(
    quiver_points=20,
    xlim_min=0.2,
    xlim_max=0.2,
    ylim_min=0.0,
    ylim_max=0.0,
    aspect=1.5,
    hide_vector_field=True,
)

# Display the wavevectors of EM modes.
v_kz = simres_EM_pump.kz_EM_all()
print("\n k_z of EM modes [1/m]:")
for i, kz in enumerate(v_kz):
    print(f"{i:3d}  {np.real(kz):.4e}")

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(simres_EM_pump.neff(0))
print("n_eff = ", np.round(n_eff_sim, 4))

q_AC = np.real(
    simres_EM_pump.kz_EM(EM_ival_pump) - simres_EM_Stokes.kz_EM(EM_ival_Stokes)
)

# Specify the expected acoustic frequency (chosen slightly lower than likely resonances).
shift_Hz = 4e9

# Calculate Acoustic modes.
if recalc_fields:
    simres_AC = wguide.calc_AC_modes(
        num_modes_AC, q_AC, EM_sim=simres_EM_pump, shift_Hz=shift_Hz
    )
    simres_AC.save_simulation("tut_07_acoustic")
else:
    simres_AC = mode_calcs.load_simulation("tut_07_acoustic")

simres_AC.plot_modes(quiver_points=20, aspect=4)

# Print the frequencies of AC modes.
v_nu = simres_AC.nu_AC_all()
print("\n Freq of AC modes (GHz):")
for m, nu in enumerate(v_nu):
    print(f"{m:3d}  {np.real(nu)*1e-9:.5f}")

set_q_factor = 1000.0

gain_box = integration.get_gains_and_qs(
    simres_EM_pump,
    simres_EM_Stokes,
    simres_AC,
    q_AC,
    EM_ival_pump=EM_ival_pump,
    EM_ival_Stokes=EM_ival_Stokes,
    AC_ival=AC_ival,
    fixed_Q=set_q_factor,
)

print("\nContributions to SBS gain [1/(WM)]")
print("AC Mode number | Photoelastic (PE) | Moving boundary(MB) | Total")
for m, nu in enumerate(v_nu):
    print(
        f"{m:7d}    {np.real(nu)/SI_GHz:9.4e} {gain_box.gain_total(m):13.3e} {gain_box.gain_PE(m):13.3e} {gain_box.gain_MB(m):13.3e} "
    )


# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = np.real(simres_AC.nu_AC_all()[0]) - 2e9  # GHz
freq_max = np.real(simres_AC.nu_AC_all()[-1]) + 2e9  # GHz

gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max)
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True)

print(nbapp.final_report())
