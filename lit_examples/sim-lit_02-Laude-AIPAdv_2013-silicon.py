""" Replicating the results of
    Generation of phonons from electrostriction in
    small-core optical waveguides
    Laude et al.
    http://dx.doi.org/10.1063/1.4801936

    Replicating silicon example.
    Note requirement for lots of modes and therefore lots of memory.
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
wl_nm = 1550
domain_x = 4*wl_nm
domain_y = domain_x*2/3
inc_a_x = 1500
inc_a_y = 1000
inc_shape = 'rectangular'

# Optical Parameters
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
#num_modes_AC = 800
num_modes_AC = 300
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'


prefix, refine_fac = starter.read_args(2, sys.argv)

nbapp = numbat.NumBATApp(prefix)

# Use all specified parameters to create a waveguide object.
wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=materials.make_material("Si_2013_Laude"),
                        lc_bkg=.1, lc_refine_1=8*refine_fac, lc_refine_2=8.0*refine_fac)
wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = 3.4

# Calculate Electromagnetic modes.
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)
sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

sim_EM_pump.plot_modes(xlim_min=0.2, xlim_max=0.2, ivals=[EM_ival_pump],
                         ylim_min=0.2, ylim_max=0.2, )

# Print the wavevectors of EM modes.
kzs = sim_EM_pump.kz_EM_all()
print('k_z of EM modes \n', np.round(np.real(kzs), 4))

# Calculate the EM effective index of the waveguide.
print("n_eff = ", np.round(sim_EM_pump.neff_all(), 4))

q_AC = np.real(sim_EM_pump.kz_EM_all()[EM_ival_pump] - sim_EM_Stokes.kz_EM_all()[EM_ival_Stokes])
print(q_AC)

shift_Hz = 31e9

# Calculate Acoustic modes.
sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)

sim_AC.plot_modes()

# Print the frequencies of AC modes.
print('Freq of AC modes (GHz) \n', np.round(np.real(sim_AC.nu_AC_all())*1e-9, 4))

# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 20e9
freq_max = 45e9
plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    logy=True, )

print(nbapp.final_report())
