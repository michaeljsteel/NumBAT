""" Replicating the results of
    Generation of phonons from electrostriction in
    small-core optical waveguides
    Laude et al.
    http://dx.doi.org/10.1063/1.4801936

    Replicating silica example for backwards SBS.
     -- Compare gain in log plot with Fig 4. in paper.
         (Different though related quantities showing very similar resonant mode pattern.)
     -- Mode profiles m=5 (5.7 GHz) and
"""

import sys
import numpy as np
from pathlib import Path

sys.path.append(str(Path('../backend')))

import numbat
import materials
import modecalcs
import integration

import starter


# Geometric Parameters - all in nm.
wl_nm = 1550
domain_x = 7*wl_nm
domain_y = domain_x
inc_a_x = 1500
inc_a_y = 1000
inc_shape = 'rectangular'

# Optical Parameters
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 120
EM_mode_index_pump = 0
EM_mode_index_Stokes = 0
AC_mode_index = 'All'

prefix, refine_fac = starter.read_args(1, sys.argv)

nbapp=numbat.NumBATApp(prefix)

wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=materials.make_material("SiO2_2013_Laude"),
                        lc_bkg=.05, lc_refine_1=2.5*refine_fac, lc_refine_2=2.5*refine_fac)

#wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = 1.3

# Calculate Electromagnetic modes.
recalc=True
#recalc=False
if recalc:
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)
    sim_EM_Stokes = sim_EM_pump.clone_as_backward_modes()
    sim_EM_pump.save_simulation(prefix+'_pump')
    sim_EM_Stokes.save_simulation(prefix+'_pump')

    sim_EM_pump.plot_modes(xlim_min=0.2, xlim_max=0.2, mode_indices=[EM_mode_index_pump],
                         ylim_min=0.2, ylim_max=0.2, )
else:
    sim_EM_pump = numbat.load_simulation(prefix+'_pump')
    sim_EM_Stokes = numbat.load_simulation(prefix+'_pump')


# Print the wavevectors of EM modes.
kzs = sim_EM_pump.kz_EM_all()
print('k_z of EM modes \n', np.round(np.real(kzs), 4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff_all())
print("n_eff = ", np.round(n_eff_sim, 4))

q_AC = np.real(sim_EM_pump.kz_EM_all()[EM_mode_index_pump] - sim_EM_Stokes.kz_EM_all()[EM_mode_index_Stokes])

shift_Hz = 8e9

# Calculate Acoustic modes.
if recalc:
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)
    sim_AC.save_simulation(prefix+'_ac')
else:
    sim_AC = numbat.load_simulation(prefix+'_ac')

print('maintest',
    np.max(np.abs(sim_AC.fem_evecs[0,:6,0,:])),
    np.max(np.abs(sim_AC.fem_evecs[1,:6,0,:])),
    np.max(np.abs(sim_AC.fem_evecs[2,:6,0,:])),
      '\n',
    np.max(np.abs(sim_AC.fem_evecs[0,:6,1,:])),
    np.max(np.abs(sim_AC.fem_evecs[1,:6,1,:])),
    np.max(np.abs(sim_AC.fem_evecs[2,:6,1,:])),
    )



gain_box = integration.get_gains_and_qs(sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_mode_index_pump=EM_mode_index_pump, EM_mode_index_Stokes=EM_mode_index_Stokes, AC_mode_index=AC_mode_index)

print('maintest2',
    np.max(np.abs(sim_AC.fem_evecs[0,:6,0,:])),
    np.max(np.abs(sim_AC.fem_evecs[1,:6,0,:])),
    np.max(np.abs(sim_AC.fem_evecs[2,:6,0,:])),
      '\n',
    np.max(np.abs(sim_AC.fem_evecs[0,:6,1,:])),
    np.max(np.abs(sim_AC.fem_evecs[1,:6,1,:])),
    np.max(np.abs(sim_AC.fem_evecs[2,:6,1,:])),
    )


print('Gains by acoustic mode:')
print('Ac. mode | Freq (GHz) | G_tot (1/mW) | G_PE (1/mW) | G_MB (1/mW)')
v_nu = sim_AC.nu_AC_all()
for (m, nu) in enumerate(v_nu):
    print(f'{m:7d}    {np.real(nu)*1e-9:9.4e} {gain_box.gain_total(m):13.3e} ',
          f'{gain_box.gain_PE(m):13.3e} {gain_box.gain_MB(m):13.3e} ')



# find indices selection of nplot highest gain modes to plot
#nplot=min(20, len(v_nu))
nplot=min(100, len(v_nu))
high_g_indices = np.abs(gain_box.gain_total_all()).argsort()[-nplot:]
high_g_indices.sort()

if recalc:
    sim_AC.plot_modes(mode_indices=high_g_indices)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 4.e9  # Hz
freq_max = 13.e9  # Hz

gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True,  mode_comps=True, dB=True)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5.3e9  # Hz
freq_max = 6.6e9  # Hz
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, suffix='_zoom')

print(nbapp.final_report())
