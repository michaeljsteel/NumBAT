"""
Script to evaluate intermodal forward Brillouin scattering in a cylindrical SiO2 waveguide
"""

# Import the necessary packages


import sys
import math
import numpy as np

from pathlib import Path
sys.path.append(str(Path('../../backend')))
import numbat
import materials
import modecalcs
import integration
import plotgain


import starter

# Naming conventions
# AC: acoustic
# EM: electromagnetic
# q_AC: acoustic wavenumber


# Specify Geometric Parameters - all in [nm].
wl_nm = 1550  # Wavelength of EM wave in vacuum.
# Unit cell dimensions must be sufficiently large to ensure fields are zero at outermost boundary.
domain_x = 4.01 * wl_nm  # be careful to ensure not whole integer(8) multiples
domain_y = domain_x
inc_a_x = 1000  # Waveguide width.
inc_a_y = inc_a_x
inc_shape = "circular"  # Shape of the waveguide.

# Specify number of electromagnetic modes, acoustic modes, and which EM indices
# are involved in the calculation for intermodal FSBS
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 100  # Number of acoustic modes to solve for.
# The EM pump mode(s) for which to calculate interaction with AC modes.
# Can specify a mode number (zero has lowest propagation constant) or 'All'.
EM_mode_index_pump = 1
# The EM Stokes mode(s) for which to calculate interaction with AC modes.
EM_mode_index_Stokes = 0
# The AC mode(s) for which to calculate interaction with EM modes.
AC_mode_index = "All"

prefix, refine_fac = starter.read_args(5, sys.argv)

nbapp = numbat.NumBATApp(prefix)

# Use all specified parameters to create a waveguide object
wguide = nbapp.make_structure(
    inc_shape,
    domain_x,
    domain_y,
    inc_a_x,
    inc_a_y,
    material_bkg=materials.make_material("Vacuum"),
    material_a=materials.make_material("SiO2_2021_Poulton"),
    lc_bkg=0.05,  # mesh coarseness in background, larger lc_bkg = coarser along horizontal outer edge
    lc_refine_1=4
    * refine_fac,  # mesh refinement factor near the interface of waveguide, larger lc2 = finer along horizontal interface
    lc_refine_2=5 * refine_fac,
)  # mesh refinement factor near the origin/centre of waveguide

# Initial guess for the EM effective index of the waveguide
n_eff = wguide.get_material("a").refindex_n - 0.1

# Calculate Electromagnetic Modes
print("Starting EM pump modes")
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff, debug=False)

print("Starting EM Stokes modes")
sim_EM_Stokes = sim_EM_pump.clone_as_forward_modes()

# Generate images for the EM modes involved in the calculation
print("Starting EM field plotting ")
sim_EM_pump.plot_modes(
    mode_indices=[EM_mode_index_pump, EM_mode_index_Stokes],
    field_type="EM_E",
    num_ticks=3,
    xlim_min=0.2,
    xlim_max=0.2,
    ylim_min=0.2,
    ylim_max=0.2,
    quiver_points=40,
    n_points=1000,
    colorbar=True,
)

# A computation interruption if needed
# sys.exit("We interrupt your regularly scheduled computation to bring you something completely different... for now")

# Print the wavevectors of EM modes.
v_kz = sim_EM_pump.kz_EM_all()
print("\n k_z of EM modes [1/m]:")
for i, kz in enumerate(v_kz):
    print(f"{i:3d}  {np.real(kz):.4e}")

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff(0))

print("n_eff = ", np.round(n_eff_sim, 4))

# Calculate and print the acoustic wave vector
q_AC = np.real(sim_EM_pump.kz_EM(EM_mode_index_pump) - sim_EM_Stokes.kz_EM(EM_mode_index_Stokes))
print("Intermode q_AC (Hz) \n", q_AC)

# Calculate Acoustic Modes
print("Starting acoustic modes")
sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, debug=False)

# Print the frequencies of AC modes.
AC_freqs_GHz = sim_AC.nu_AC_all() * 1e-9
print("\n Freq of AC modes (GHz):")
for i, nu in enumerate(AC_freqs_GHz):
    print(f"{i:3d}  {np.real(nu):.4e}")

# Calculate total SBS gain, photoelastic and moving boundary contributions, as
# well as other important quantities
gain_box= integration.get_gains_and_qs(
        sim_EM_pump,
        sim_EM_Stokes,
        sim_AC,
        q_AC,
        EM_mode_index_pump=EM_mode_index_pump,
        EM_mode_index_Stokes=EM_mode_index_Stokes,
        AC_mode_index=AC_mode_index,
    )

freq_min = 2.5e9
freq_max = 3.0e9
gain_box.plot_spectra(freq_min=freq_min,freq_max=freq_max)

# Mask negligible gain values to improve clarity of print out.
threshold = 1e-3
masked_PE = np.ma.masked_inside(gain_box.gain_PE_all(), 0, threshold)
masked_MB = np.ma.masked_inside(gain_box.gain_MB_all(), 0, threshold)
masked = np.ma.masked_inside(gain_box.gain_total_all(), 0, threshold)

# Display these in terminal
print("\n Displaying results with negligible components masked out")
print("SBS_gain [1/(Wm)] PE contribution \n", masked_PE)
print("SBS_gain [1/(Wm)] MB contribution \n", masked_MB)
print("SBS_gain [1/(Wm)] total \n", masked)
# determining the location of the maximum gain
maxGainloc = 6
# note sometimes its necessary to manually specify as certain values are NOT possible by symmetry arguments

print("Plotting acoustic mode corresponding to maximum")
sim_AC.plot_modes(mode_indices=range(15), num_ticks=3, quiver_points=40, colorbar=True)

# Displaying results for the maximum found in the selection
print("\n\n-----------------")
print('Displaying results for maximum (physically realisable) "gain" value found:')
print(f"Greatest total SBS_gain  [1/(Wm)] = {masked.data[maxGainloc]:4e}) for  acoustic mode number {maxGainloc}.")
print(f"EM Pump Power [Watts]:       {np.real(sim_EM_pump.EM_mode_power[EM_mode_index_pump]):.5e}")
print(f"EM Stokes Power [Watts]:     {np.real(sim_EM_Stokes.EM_mode_power[EM_mode_index_Stokes]):.5e}")
print(f"EM angular frequency [THz]:  {sim_EM_pump.omega_EM / 1e12:.5e}")
print(f"AC Energy Density [J/m]:     {sim_AC.AC_mode_energy[maxGainloc]:.5e}")
print(f"AC loss alpha [1/s]:         {gain_box.alpha[maxGainloc]:.5e}")
print(f"AC frequency [GHz]:          {np.real(sim_AC.Omega_AC[maxGainloc]) / (1e9 * 2 * math.pi):.5e}")
print(f"AC linewidth [MHz]:          {gain_box.linewidth_Hz[maxGainloc] / 1e6:.5e}")

# since the overlap is not returned directly we'll have to deduce it
absQtot2 = (
    gain_box.alpha[maxGainloc]
    * sim_EM_pump.EM_mode_power[EM_mode_index_pump]
    * sim_EM_Stokes.EM_mode_power[EM_mode_index_Stokes]
    * sim_AC.AC_mode_energy[maxGainloc]
    * masked.data[maxGainloc]
) / (2 * sim_EM_pump.omega_EM * sim_AC.Omega_AC[maxGainloc])
absQtot = pow(absQtot2, 1 / 2)
print(f"Total coupling |Qtot| [J/m]: {np.real(absQtot):.4e}")

print(nbapp.final_report())
