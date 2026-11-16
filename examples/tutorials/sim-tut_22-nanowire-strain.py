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
from numbat import materials

from numbat.solvers import sbsgain

import starter


# Geometric Parameters - all in nm.
lambda_nm = 1550
domain_x = 4500
domain_y = 1.0*domain_x
inc_a_x = 550
inc_a_y = inc_a_x
inc_shape = 'circular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_mode_index_pump = 0
EM_mode_index_Stokes = 0
AC_mode_index = 'All'

prefix, refine_fac = starter.read_args(22, sys.argv)

nbapp = numbat.NumBATApp(prefix)

wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=materials.make_material("SiO2_2016_Smith"),
                        lc_bkg=.05, lc_refine_1=1.5*refine_fac, lc_refine_2=1.5*refine_fac)

wguide.plot_refractive_index_profile(prefix)

# Expected effective index of fundamental guided mode.
n_eff = 1.4

use_old_fields=False     # run the calculation from scratch
#use_old_fields=True   # reuse saved fields from previous calculation

# Calculate Electromagnetic Modes
simres_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff=n_eff)
simres_EM_Stokes = simres_EM_pump.clone_as_backward_modes()

print('\nPlotting EM fields')
trim=0.0

# Display the wavevectors of EM modes.
v_kz=simres_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz): print(f'{i:3d}  {np.real(kz):.4e}')

#simres_EM_pump.plot_modes(mode_indices=range(5), trim=(trim, trim, trim, trim), comps="exyt")

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(simres_EM_pump.neff(0))
print(f"n_eff = {n_eff_sim:.4e}")

# Acoustic wavevector
q_AC = np.real(simres_EM_pump.kz_EM(EM_mode_index_pump) - simres_EM_Stokes.kz_EM(EM_mode_index_Stokes))

shift_Hz = 4e9

# Calculate Acoustic modes.
simres_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=simres_EM_pump, shift_Hz=shift_Hz)

# Print the frequencies of AC modes.
v_nu=simres_AC.nu_AC_all()
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(v_nu): print(f'{i:3d}  {np.real(nu)*1e-9:.5f}')

md = simres_AC.get_mode(0)
md.plot_strain_1D("x", 0)
md.plot_strain_1D("y", 0)

md.write_strain()

print(nbapp.final_report())
