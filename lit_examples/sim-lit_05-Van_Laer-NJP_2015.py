""" Replicating the results of
    Net on-chip Brillouin gain based on suspended
    silicon nanowires
    Van Laer et al.
    http://dx.doi.org/10.1088/1367-2630/17/11/115005
"""

import sys
import copy
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
domain_x = 5*wl_nm
domain_y = 0.5*domain_x
inc_a_x = 450
inc_a_y = 230
inc_shape = 'rectangular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 60
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'


prefix, refine_fac = starter.read_args(5, sys.argv)

nbapp = numbat.NumBATApp(prefix)

# Rotate crystal axis of Si from <100> to <110>, starting with same Si_2016_Smith data.
Si_110 = copy.deepcopy(materials.make_material("Si_2016_Smith"))
Si_110.rotate_axis('y-axis', np.pi/4, save_rotated_tensors=True)
# Use all specified parameters to create a waveguide object.
wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                        material_bkg=materials.make_material("Vacuum"),
                        material_a=Si_110, symmetry_flag=False,
                        lc_bkg=.1, lc_refine_1=15.0*refine_fac, lc_refine_2=15.0*refine_fac)
wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

# Calculate Electromagnetic Modes
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)
# np.savez('wguide_data', sim_EM_pump=sim_EM_pump)
# npzfile = np.load('wguide_data.npz')
# sim_EM_pump = npzfile['sim_EM_pump'].tolist()

sim_EM_Stokes = modecalcs.fwd_Stokes_modes(sim_EM_pump)
# np.savez('wguide_data2', sim_EM_Stokes=sim_EM_Stokes)
# npzfile = np.load('wguide_data2.npz')
# sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

sim_EM_pump.plot_modes(xlim_min=0.45, xlim_max=0.45,
                         ivals=[EM_ival_pump], ylim_min=0.45, ylim_max=0.45, n_points=1500, )

# Print the wavevectors of EM modes.
kzs = sim_EM_pump.kz_EM_all()
print('k_z of EM modes \n', np.round(np.real(kzs), 4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.neff_all())
print("n_eff = ", np.round(n_eff_sim, 4))

q_AC = 5 # close but not quite zero

# Calculate Acoustic Modes
sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump)
# np.savez('wguide_data_AC', sim_AC=sim_AC)
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC = npzfile['sim_AC'].tolist()

# Print the frequencies of AC modes.
print('Freq of AC modes (GHz) \n', np.round(np.real(sim_AC.nu_AC_all())*1e-9, 4))

sim_AC.plot_modes()

set_q_factor = 230 # NJP

# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
gain_box = integration.get_gains_and_qs(sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)

print('Gains by acoustic mode:')
print('Ac. mode | Freq (GHz) | G_tot (1/mW) | G_PE (1/mW) | G_MB (1/mW)')
v_nu = sim_AC.nu_AC_all()
for (m, nu) in enumerate(v_nu):
    print(f'{m:7d}    {np.real(nu)*1e-9:9.4e} {gain_box.gain_total(m):13.3e} ',
          f'{gain_box.gain_PE(m):13.3e} {gain_box.gain_MB(m):13.3e} ')

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5.0e9 # Hz
freq_max = 20.0e9 # Hz
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True,  mode_comps=True)

print(nbapp.final_report())
