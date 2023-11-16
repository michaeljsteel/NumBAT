"""
Script to evaluate intermodal forward Brillouin scattering in a cylindrical SiO2 waveguide
"""

# Import the necessary packages


import numpy as np
import sys
import math
sys.path.append("../backend/")
import numbat
import materials
import mode_calcs
import integration
import plotting


import starter

# Naming conventions
# AC: acoustic
# EM: electromagnetic
# q_AC: acoustic wavenumber



# Specify Geometric Parameters - all in [nm].
wl_nm = 1550 # Wavelength of EM wave in vacuum.
# Unit cell dimensions must be sufficiently large to ensure fields are zero at outermost boundary.
unitcell_x = 4.01*wl_nm #be careful to ensure not whole integer multiples
unitcell_y = unitcell_x
inc_a_x = 1000 # Waveguide width.
inc_a_y = inc_a_x
inc_shape = 'circular' # Shape of the waveguide.

# Specify number of electromagnetic modes, acoustic modes, and which EM indices
# are involved in the calculation for intermodal FSBS
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 100 # Number of acoustic modes to solve for.
# The EM pump mode(s) for which to calculate interaction with AC modes.
# Can specify a mode number (zero has lowest propagation constant) or 'All'.
EM_ival_pump = 1
# The EM Stokes mode(s) for which to calculate interaction with AC modes.
EM_ival_Stokes = 0
# The AC mode(s) for which to calculate interaction with EM modes.
AC_ival = 'All'

prefix, refine_fac = starter.read_args(5, sys.argv)

nbapp = numbat.NumBATApp(prefix)

# Use all specified parameters to create a waveguide object
wguide = nbapp.make_structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=materials.make_material("SiO2_2021_Poulton"),
                        lc_bkg=0.05, # mesh coarseness in background, larger lc_bkg = coarser along horizontal outer edge
                        lc_refine_1=4*refine_fac, # mesh refinement factor near the interface of waveguide, larger lc2 = finer along horizontal interface
                        lc_refine_2=5*refine_fac) # mesh refinement factor near the origin/centre of waveguide

# Initial guess for the EM effective index of the waveguide
n_eff = wguide.get_material('a').refindex_n-0.1

# Calculate Electromagnetic Modes
print("Starting EM pump modes")
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff, debug=False)

print("Starting EM Stokes modes")
sim_EM_Stokes = mode_calcs.fwd_Stokes_modes(sim_EM_pump)

# Generate images for the EM modes involved in the calculation
print("Starting EM field plotting ")
plotting.plot_mode_fields(sim_EM_pump,
                         ivals=[EM_ival_pump,EM_ival_Stokes],
                         EM_AC='EM_E', num_ticks=3,xlim_min=0.2, xlim_max=0.2, ylim_min=0.2, ylim_max=0.2,
                          quiver_points=40,
                         n_points=1000, colorbar=True)

# A computation interruption if needed
# sys.exit("We interrupt your regularly scheduled computation to bring you something completely different... for now")

# Print the wavevectors of EM modes.
v_kz=sim_EM_pump.kz_EM_all()
print('\n k_z of EM modes [1/m]:')
for (i, kz) in enumerate(v_kz): print('{0:3d}  {1:.4e}'.format(i, np.real(kz)))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff(0))

print("n_eff = ", np.round(n_eff_sim, 4))

# Calculate and print the acoustic wave vector
q_AC = np.real(sim_EM_pump.kz_EM(EM_ival_pump) - sim_EM_Stokes.kz_EM(EM_ival_Stokes))
print('Intermode q_AC (Hz) \n', q_AC)

# Calculate Acoustic Modes
print("Starting acoustic modes")
sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, debug=False)

# Print the frequencies of AC modes.
AC_freqs_GHz=sim_AC.nu_AC_all()*1e-9
print('\n Freq of AC modes (GHz):')
for (i, nu) in enumerate(AC_freqs_GHz): print('{0:3d}  {1:.4e}'.format(i, np.real(nu)))

# Calculate total SBS gain, photoelastic and moving boundary contributions, as
# well as other important quantities
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)

freq_min = 2.5e9
freq_max = 3.0e9
plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max, )

# Mask negligible gain values to improve clarity of print out.
threshold = 1e-3
masked_PE = np.ma.masked_inside(SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked_MB = np.ma.masked_inside(SBS_gain_MB[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked = np.ma.masked_inside(SBS_gain[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)

# Display these in terminal
print("\n Displaying results with negligible components masked out")
print("SBS_gain [1/(Wm)] PE contribution \n", masked_PE)
print("SBS_gain [1/(Wm)] MB contribution \n", masked_MB)
print("SBS_gain [1/(Wm)] total \n", masked)
# determining the location of the maximum gain
maxGainloc=6; #note sometimes its necessary to manually specify as certain values are NOT possible by symmetry arguments

print("Plotting acoustic mode corresponding to maximum")
plotting.plot_mode_fields(sim_AC,  ivals=range(15),
                          num_ticks=3, quiver_points=40, colorbar=True)

# Displaying results for the maximum found in the selection
print("-----------------")
print("Displaying results for maximum gain value found:")
print("Greatest SBS_gain  [1/(Wm)] total \n", masked.data[maxGainloc])
print("displaying corresponding acoustic mode number (i.e., AC_field_#) for reference \n",maxGainloc )
print("EM Pump Power [Watts] \n", sim_EM_pump.EM_mode_power[EM_ival_pump] )
print("EM Stokes Power [Watts] \n", sim_EM_Stokes.EM_mode_power[EM_ival_Stokes] )
print("EM angular frequency [THz] \n", sim_EM_pump.omega_EM/1e12 )
print("AC Energy Density [J*m^{-1}] \n", sim_AC.AC_mode_energy[maxGainloc] )
print("AC loss alpha [1/s] \n", alpha[maxGainloc] )
print("AC frequency [GHz] \n", sim_AC.Omega_AC[maxGainloc]/(1e9*2*math.pi) )
print("AC linewidth [MHz] \n", linewidth_Hz[maxGainloc]/1e6)

#since the overlap is not returned directly we'll have to deduce it
absQtot2 = (alpha[maxGainloc]*sim_EM_pump.EM_mode_power[EM_ival_pump]*sim_EM_Stokes.EM_mode_power[EM_ival_Stokes]*sim_AC.AC_mode_energy[maxGainloc]*masked.data[maxGainloc])/(2*sim_EM_pump.omega_EM*sim_AC.Omega_AC[maxGainloc]);
absQtot = pow(absQtot2,1/2)
print("Total coupling |Qtot| [W*m^{-1}*s] \n", absQtot )

print(nbapp.final_report())
