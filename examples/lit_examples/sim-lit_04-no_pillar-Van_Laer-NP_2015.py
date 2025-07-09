""" Replicating the results of
    Interaction between light and highly confined
    hypersound in a silicon photonic nanowire
    Van Laer et al.
    http://dx.doi.org/10.1038/nphoton.2015.11

    Making simplification of ignoring the pedestal.
"""

import sys
import copy
import numpy as np

from pathlib import Path
sys.path.append(str(Path('../../backend')))

from plotgain import Decorator
import integration
import modecalcs
import materials
import numbat

import starter


# use this class to add or alter features to the final plots
class EMDecorator(Decorator):
    def __init__(self):
        super().__init__()

    def extra_axes_commands(self, ax):
        ax.tick_params(axis='x', color='gray', which='both')
        ax.tick_params(axis='y', color='gray', which='both')
        ax.tick_params(axis='x', length=20)
        ax.tick_params(axis='y', length=20)
        ax.tick_params(axis='x', width=2)
        ax.tick_params(axis='y', width=2)

        ax.tick_params(axis='x', length=10, which='minor')
        ax.tick_params(axis='y', length=10, which='minor')
        ax.tick_params(axis='x', width=2, which='minor')
        ax.tick_params(axis='y', width=2, which='minor')


emdecorate = EMDecorator()
# acdecorate=ACDecorator()


# Geometric Parameters - all in nm.
wl_nm = 1550
domain_x = 5*wl_nm
domain_y = 0.5*domain_x
inc_a_x = 485
inc_a_y = 230
inc_shape = 'rectangular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 60
EM_mode_index_pump = 0
EM_mode_index_Stokes = 0
AC_mode_index = 'All'

prefix, refine_fac = starter.read_args(4, sys.argv, sub='a')

nbapp = numbat.NumBATApp(prefix)
# Rotate crystal axis of Si from <100> to <110>, starting with same Si_2016_Smith data.
Si_110 = copy.deepcopy(materials.make_material("Si_2016_Smith"))
# Si_110 = copy.deepcopy(materials.materials_dict["Si_2015_Van_Laer"])
Si_110.rotate('z-axis', np.pi/4)
# Use all specified parameters to create a waveguide object.
wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                              material_bkg=materials.make_material("Vacuum"),
                              material_a=Si_110, symmetry_flag=False,
                              lc_bkg=.1, lc_refine_1=20.0*refine_fac, lc_refine_2=20.0*refine_fac)

wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1


doem = True
doac = True
new_calcs = True

if doem:
    # Calculate Electromagnetic Modes
    if new_calcs:
        sim_EM_pump = wguide.calc_EM_modes(
            num_modes_EM_pump, wl_nm, n_eff=n_eff)
        np.savez(prefix+'-wguide_data', sim_EM_pump=sim_EM_pump)
    else:
        npzfile = np.load(prefix+'-wguide_data.npz', allow_pickle=True)
        sim_EM_pump = npzfile['sim_EM_pump'].tolist()

    sim_EM_Stokes = sim_EM_pump.clone_as_forward_modes()
    np.savez(prefix+'-wguide_data2', sim_EM_Stokes=sim_EM_Stokes)
    # npzfile = np.load(prefix+'-wguide_data2.npz', allow_pickle=True)
    # sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

    sim_EM_pump.plot_modes(xlim_min=0.43, xlim_max=0.43, mode_indices=[EM_mode_index_pump],
                        n_points=2000, quiver_points=10,  decorator=emdecorate)

    # Print the wavevectors of EM modes.
    print('k_z of EM modes \n', np.round(np.real(sim_EM_pump.kz_EM_all()), 4))

    # Calculate the EM effective index of the waveguide.
    print("n_eff = ", np.round(sim_EM_pump.neff_all(), 4))


if doac:
    q_AC = 5  # close but not quite zero

    # Calculate Acoustic Modes
    if new_calcs:
        sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump)
        # np.savez(prefix+'-wguide_data_AC', sim_AC=sim_AC)  # Pickle problem
    else:
        npzfile = np.load(prefix+'-wguide_data_AC.npz', allow_pickle=True)
        sim_AC = npzfile['sim_AC'].tolist()

    sim_AC.plot_modes(mode_indices=range(20))

set_q_factor = 306

# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
gain_box = integration.get_gains_and_qs(sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_mode_index_pump=EM_mode_index_pump, EM_mode_index_Stokes=EM_mode_index_Stokes, AC_mode_index=AC_mode_index)

print('Gains by acoustic mode:')
print('Ac. mode | Freq (GHz) | G_tot (1/mW) | G_PE (1/mW) | G_MB (1/mW)')
v_nu = sim_AC.nu_AC_all()
for (m, nu) in enumerate(v_nu):
    print(f'{m:7d}    {np.real(nu)*1e-9:9.4e} {gain_box.gain_total(m):13.3e} ',
          f'{gain_box.gain_PE(m):13.3e} {gain_box.gain_MB(m):13.3e} ')


# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5e9  # 9.1e9 # Hz
freq_max = 50e9  # 9.3e9 # Hz
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True,  mode_comps=True, dB=True)

print(nbapp.final_report())
