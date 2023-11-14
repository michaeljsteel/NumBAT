""" Replicating the results of
    Interaction between light and highly confined
    hypersound in a silicon photonic nanowire
    Van Laer et al.
    http://dx.doi.org/10.1038/nphoton.2015.11
"""

import sys
import time
import copy
import numpy as np

sys.path.append("../backend/")
import numbat
import materials
import objects
import mode_calcs
import integration
import plotting

import starter


# Geometric Parameters - all in nm.
wl_nm = 1550
unitcell_x = 4*wl_nm
unitcell_y = 0.5*unitcell_x
inc_a_x = 450
inc_a_y = 230
inc_shape = 'pedestal'
pillar_x = 15
pillar_y = 300
slab_a_x = 2000
slab_a_y = 800

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 60
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix, refine_fac = starter.read_args(4, sys.argv, sub='b')

nbapp = numbat.NumBAT()

# Rotate crystal axis of Si from <100> to <110>, starting with same Si_2016_Smith data.
Si_110 = copy.deepcopy(materials.make_material("Si_2015_Van_Laer"))
Si_110.rotate_axis(np.pi/4,'y-axis', save_rotated_tensors=True)

# Use all specified parameters to create a waveguide object.
wguide = objects.Structure(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        slab_a_x=slab_a_x, slab_a_y=slab_a_y,
                        pillar_x=pillar_x, pillar_y=pillar_y,
                        material_bkg=materials.make_material("Vacuum"),            # background
                        material_a=Si_110,                        # rib
                        material_b=materials.make_material("SiO2_2015_Van_Laer"),  # slab
                        material_c=materials.make_material("SiO2_2015_Van_Laer"),  # pillar
                        lc_bkg=.1, lc_refine_1=10.0*refine_fac, lc_refine_2=10.0*refine_fac)

wguide.plot_mesh(prefix)

# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

# Calculate Electromagnetic Modes
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff)
sim_EM_Stokes = mode_calcs.fwd_Stokes_modes(sim_EM_pump)

plotting.plot_mode_fields(sim_EM_pump, ivals=[EM_ival_pump],
                         xlim_min=0.4, xlim_max=0.4, ylim_min=0.4, ylim_max=0.2,
                         EM_AC='EM_E', prefix=prefix)

# Print the wavevectors of EM modes.
print('k_z of EM modes \n', np.round(np.real(sim_EM_pump.kz_EM_all()), 4))

q_AC = 5
shift_Hz = 8e9

# Calculate Acoustic Modes
sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)

plotting.plot_mode_fields(sim_AC, prefix=prefix)

set_q_factor = 306

# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = np.real(sim_AC.nu_AC_all()[0]) - 2e9  # GHz
freq_max = np.real(sim_AC.nu_AC_all()[-1]) + 2e9  # GHz
plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    prefix=prefix)

print(nbapp.final_report())
