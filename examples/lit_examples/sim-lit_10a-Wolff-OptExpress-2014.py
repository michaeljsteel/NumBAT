""" Replicating the results of
    Germanium as a material for stimulated Brillouin scattering in the mid-infrared
    Wolff et al.
    https://doi.org/10.1364/OE.22.030735
"""

import sys
import copy
import time
import numpy as np
from pathlib import Path

sys.path.append(str(Path('../backend')))

import materials
import integration
from nbtypes import SI_GHz
import numbat

from nbtypes import PointGroup

import starter

# Naming conventions
# AC: acoustic
# EM: electromagnetic
# q_AC: acoustic wavenumber

start = time.time()

# Geometric Parameters - all in nm.
wl_nm = 4000  # Wavelength of EM wave in vacuum.
# Unit cell must be large to ensure fields are zero at boundary.
domain_x = 1.5 * wl_nm
domain_y = 0.75 * domain_x
# Waveguide widths.
inc_a_x = 1020
inc_a_y = 700
# Shape of the waveguide.
inc_shape = "rectangular"


# Number of electromagnetic modes to solve for.
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
# Number of acoustic modes to solve for.
num_modes_AC = 100
# The EM pump mode(s) for which to calculate interaction with AC modes.
# Can specify a mode number (zero has lowest propagation constant) or 'All'.
EM_mode_index_pump = 0
# The EM Stokes mode(s) for which to calculate interaction with AC modes.
EM_mode_index_Stokes = 0
# The AC mode(s) for which to calculate interaction with EM modes.
AC_mode_index = "All"

# reuse_fields=True   # use saved data
reuse_fields = False  # calculate from scratch

Ge_110 = copy.deepcopy(materials.make_material("Ge_cubic_2014_Wolff"))
print("Initial Ge_100:", Ge_110.full_str())
Ge_110.rotate("y-axis", np.pi / 4.0, save_rotated_tensors=True)
print("Rotated Ge_110:", Ge_110.full_str())

prefix, refine_fac = starter.read_args(10, sys.argv, sub="a")
nbapp = numbat.NumBATApp(prefix)

# Use specified parameters to create a waveguide object.
# Note use of rough mesh for demonstration purposes.
wguide = nbapp.make_structure(
    inc_shape,
    domain_x,
    domain_y,
    inc_a_x,
    inc_a_y,
    material_bkg=materials.make_material("Si3N4_2014_Wolff"),
    material_a=Ge_110,
    lc_bkg=0.1,
    lc_refine_1=5.0,
    lc_refine_2=5.0,
)
wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material("a").refindex_n - 0.1

# Calculate the Electromagnetic modes of the pump field.
if not reuse_fields:
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff)
    np.savez("wguide_data", sim_EM_pump=sim_EM_pump)
else:
    npzfile = np.load("wguide_data.npz", allow_pickle=True)
    sim_EM_pump = npzfile["sim_EM_pump"].tolist()

sim_EM_pump.analyse_symmetries(PointGroup.C2V)
#sim_EM_pump.set_r0_offset(3.0e-6, -2.250e-6)

sim_EM_pump.plot_modes( xlim_min=0.2, xlim_max=0.2, mode_indices=[EM_mode_index_pump], ylim_min=0.2, ylim_max=0.2, num_ticks=3, ticks=True,)

if not reuse_fields:
    sim_EM_Stokes = sim_EM_pump.clone_as_backward_modes()
    np.savez("wguide_data2", sim_EM_Stokes=sim_EM_Stokes)
else:
    npzfile = np.load("wguide_data2.npz", allow_pickle=True)
    sim_EM_Stokes = npzfile["sim_EM_Stokes"].tolist()

# Print the wavevectors of EM modes.
v_kz = sim_EM_pump.kz_EM_all()
print("\n k_z of EM modes [1/m]:")
for i, kz in enumerate(v_kz):
    print(f"{i:3d}  {np.real(kz):.4e}")


# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff(0))
print("\n n_eff = ", np.round(n_eff_sim, 4))

q_AC = np.real(sim_EM_pump.kz_EM(EM_mode_index_pump) - sim_EM_Stokes.kz_EM(EM_mode_index_Stokes))
print("\n AC wavenumber (1/m) = ", np.round(q_AC, 4))


# Calculate Acoustic modes.
shift_Hz = 5.5 * 1e9  # select the lowest frequency to start FEM search from.
if not reuse_fields:
    sim_AC = wguide.calc_AC_modes(
        num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz
    )
    np.savez("wguide_data_AC", sim_AC=sim_AC)
else:
    npzfile = np.load("wguide_data_AC.npz", allow_pickle=True)
    sim_AC = npzfile["sim_AC"].tolist()

sim_AC.analyse_symmetries(PointGroup.C2V)
#sim_AC.set_r0_offset(3.0e-6, -2.250e-6)
sim_AC.plot_modes(num_ticks=3, xlim_min=0.1, xlim_max=0.1)

# Print the frequencies of AC modes.
v_nu = sim_AC.nu_AC_all()
print("\n Freq of AC modes (GHz):")
for i, nu in enumerate(v_nu):
    print(f"{i:3d}  {np.real(nu) * 1e-09:.4e}")


# set_Q_factor = 190 # set the mechanic Q manually

# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB. Also calculate acoustic loss alpha.
gain_box = integration.get_gains_and_qs(sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_mode_index_pump=EM_mode_index_pump, EM_mode_index_Stokes=EM_mode_index_Stokes, AC_mode_index=AC_mode_index)


print('Gains by acoustic mode:')
print('Ac. mode | Freq (GHz) | G_tot (1/mW) | G_PE (1/mW) | G_MB (1/mW)')
v_nu = sim_AC.nu_AC_all()
for (m, nu) in enumerate(v_nu):
    print(f'{m:7d}    {np.real(nu)*1e-9:9.4e} {gain_box.gain_total(m):13.3e} ',
          f'{gain_box.gain_PE(m):13.3e} {gain_box.gain_MB(m):13.3e} ')



# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5.0 * SI_GHz
freq_max = 11.0 * SI_GHz

gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True)

print(nbapp.final_report())
