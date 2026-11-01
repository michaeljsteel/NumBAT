""" Replicating the results of
    On-chip inter-modal Brillouin scattering
    Kittlaus et al.
    http://dx.doi.org/10.1038/ncomms15819
"""

import time
#import datetime
import sys
import copy
import numpy as np

#from matplotlib.ticker import AutoMinorLocator


from pathlib import Path
sys.path.append(str(Path('../../backend')))

import numbat
from nbtypes import SI_GHz
from plotgain import Decorator
import plotgain
import integration
import modecalcs
import materials

import starter


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

        ax.set_xticks(np.arange(-1.00, 1.01, .1), minor=True)
        ax.set_yticks([-.3, 0, .3])
        ax.set_yticks([-.5, -.4, -.2, -.1, .1, .2, .4, .5], minor=True)

        ax.set_aspect('equal')


emdecorate = EMDecorator()


class ACDecorator(plotgain.Decorator):
    def __init__(self):
        super().__init__()

        self.set_property('ax_label_fs', 30)
        self.set_property('subplot_title_fs', 30)
        self.set_property('cbar_tick', 20)
        self.set_property('ax_tick', 30)

        # compensate for call to set_aspect() below
        self.set_property('cbar_pad', '-30%')
        # compensate for call to set_aspect() below
        self.set_property('cbar_pad', '-30%')
        # compensate for call to set_aspect() below
        self.set_property('cbar_size', '2%')
        # compensate for call to set_aspect() below
        self.set_property('cbar_size', '2%')

    def extra_axes_commands(self, ax):
        ax.set_aspect(3)


acdecorate = ACDecorator()
# Naming conventions
# AC: acoustic
# EM: electromagnetic
# q_AC: acoustic wavenumber

start = time.time()

# Geometric Parameters - all in nm.
wl_nm = 1550  # Wavelength of EM wave in vacuum.

# Shape of the waveguide.
# Use double coated geometry to control meshing around rib waveguide.
inc_shape = 'rib_double_coated'

# Waveguide widths.
rib_w = 1500
rib_h = 80

slab_w = 2850
slab_h = 135

# Unit cell must be large to ensure fields are zero at boundary.
domain_x = 5000
domain_y = 0.7*domain_x

# areas included purely
slab_b_h = 500
coat_x = 50
coat_y = 100
coat2_x = 100
coat2_y = 200

lc_bkg = .05  # background
lc_refine_1 = 5  # edge of rib
lc_refine_2 = 5   # edge of slab_a
lc_refine_3 = 5    # edge of coat
lc_refine_4 = 5    # edge of slab_b
lc_refine_5 = 5     # edge of coat2


# ##working set
# lc_bkg = .5  # background
# lc_refine_1 = 1000  # edge of rib
# lc_refine_2 = 400   # edge of slab_a
# lc_refine_3 = 10    # edge of coat
# lc_refine_4 = 5   # edge of slab_b
# lc_refine_5 = 1   # edge of coat2

# #scaled working set: doesn't work
# lc_bkg = 1  # background
# lc_refine_1 = 2000  # edge of rib
# lc_refine_2 = 600   # edge of slab_a
# lc_refine_3 = 10    # edge of coat
# lc_refine_4 = 4   # edge of slab_b
# lc_refine_5 = 1  # edge of coat2

# scaled working set: doesn't work
# lc_bkg = .1  # background
# lc_refine_1 = 200  # edge of rib
# lc_refine_2 = 75   # edge of slab_a
# lc_refine_3 = 50    # edge of coat
# c5  = 50  # edge of slab_b
# lc_refine_5 = 50  # edge of coat2


# Number of electromagnetic modes to solve for.
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
# Number of acoustic modes to solve for.
num_modes_AC = 35
# The EM pump mode(s) for which to calculate interaction with AC modes.
# Can specify a mode number (zero has lowest propagation constant) or 'All'.
EM_mode_index_pump = 0
# The EM Stokes mode(s) for which to calculate interaction with AC modes.
EM_mode_index_Stokes = 1  # INTERMODE SBS TE0 to TE1
# The AC mode(s) for which to calculate interaction with EM modes.
AC_mode_index = 'All'

# Si_110 = copy.deepcopy(materials.make_material("Si_2015_Van_Laer")
Si_110 = copy.deepcopy(materials.make_material("Si_2016_Smith"))
Si_110.rotate('z-axis', np.pi/4, save_rotated_tensors=True)

prefix, refine_fac = starter.read_args(8, sys.argv)

nbapp = numbat.NumBATApp(prefix)

vac = materials.make_material("Vacuum")
# Use specified parameters to create a waveguide object.
wguide = nbapp.make_structure(inc_shape, domain_x, domain_y,
                              rib_w=rib_w, rib_h=rib_h,
                              slab_w=slab_w, slab_h=slab_h, slab_b_h=slab_b_h,
                              coat_w=coat_x, coat_h=coat_y, coat2_w=coat2_x, coat2_h=coat2_y,
                              material_bkg=vac,
                              material_a=Si_110, material_b=Si_110, material_c=vac,
                              material_d=vac, material_e=vac,
                              lc_bkg=lc_bkg, lc_refine_1=lc_refine_1, lc_refine_2=lc_refine_2,
                              lc_refine_3=lc_refine_3, lc_refine_4=lc_refine_4, lc_refine_5=lc_refine_5)
wguide.plot_mesh(prefix)
# Expected effective index of fundamental guided mode.
n_eff = wguide.get_material('a').refindex_n-0.1

# Calculate Electromagnetic Modes
print("starting EM pump modes")
sim_EM_pump = wguide.calc_EM_modes(
    num_modes_EM_pump, wl_nm, n_eff=n_eff, debug=True)
# np.savez('wguide_data', sim_EM_pump=sim_EM_pump)
# npzfile = np.load('wguide_data.npz', allow_pickle=True)
# sim_EM_pump = npzfile['sim_EM_pump'].tolist()

print("starting EM Stokes modes")
sim_EM_Stokes = sim_EM_pump.clone_as_forward_modes()
# np.savez('wguide_data2', sim_EM_Stokes=sim_EM_Stokes)
# npzfile = np.load('wguide_data2.npz', allow_pickle=True)
# sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

print("starting EM field plotting ")
sim_EM_pump.plot_modes(xlim_min=0.4, xlim_max=0.4,
                       mode_indices=[EM_mode_index_pump, EM_mode_index_Stokes],
                       ylim_min=0.435, ylim_max=0.435, field_type='EM_E', num_ticks=3,

                       decorator=emdecorate, quiver_points=20,
                       comps=('x', 'y', 'z', 'a', 't'), n_points=2000, colorbar=True)

sim_EM_pump.plot_modes(xlim_min=0.4, xlim_max=0.4,
                       mode_indices=[EM_mode_index_pump, EM_mode_index_Stokes],
                       ylim_min=0.435, ylim_max=0.435, field_type='EM_H', num_ticks=3,

                       decorator=emdecorate, quiver_points=20,
                       comps=('x', 'y', 'z', 'a', 't'), n_points=2000, colorbar=True)

# Print the wavevectors of EM modes.
print('k_z of EM modes \n', np.round(np.real(sim_EM_pump.kz_EM_all()), 4))

# Calculate the EM effective index of the waveguide.
print("n_eff = ", np.round(sim_EM_pump.neff_all(), 4))

q_AC = np.real(sim_EM_pump.kz_EM_all()[
               EM_mode_index_pump] - sim_EM_Stokes.kz_EM_all()[EM_mode_index_Stokes])
print('Intermode q_AC (Hz) \n', q_AC)

shift_Hz = 2e9

# Calculate Acoustic Modes
print("starting acoustic modes")
sim_AC = wguide.calc_AC_modes(
    num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz, debug=True)
# np.savez('wguide_data_AC', sim_AC=sim_AC)
# npzfile = np.load('wguide_data_AC.npz', allow_pickle=True)
# sim_AC = npzfile['sim_AC'].tolist()



selected_AC_modes = [7, 13, 23]
print("AC modes selected for field plotting", selected_AC_modes)
print("plotting acoustic modes")
sim_AC.plot_modes(mode_indices=selected_AC_modes,
                  num_ticks=3, xlim_min=-.05, xlim_max=-0.05, ylim_min=-.1, ylim_max=-0.1,
                  quiver_points=20, decorator=acdecorate, colorbar=True)


set_q_factor = 460.

# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
gain_box = integration.get_gains_and_qs(sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_mode_index_pump=EM_mode_index_pump, EM_mode_index_Stokes=EM_mode_index_Stokes, AC_mode_index=AC_mode_index,fixed_Q=set_q_factor)

print('Gains by acoustic mode:')
print('Ac. mode | Freq (GHz) | G_tot (1/mW) | G_PE (1/mW) | G_MB (1/mW)')
v_nu = sim_AC.nu_AC_all()
for (m, nu) in enumerate(v_nu):
    print(f'{m:7d}    {np.real(nu)*1e-9:9.4e} {gain_box.gain_total(m):13.3e} ',
          f'{gain_box.gain_PE(m):13.3e} {gain_box.gain_MB(m):13.3e} ')


# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 0.5 * SI_GHz
freq_max = 9.5 * SI_GHz
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True)

print(nbapp.final_report())
