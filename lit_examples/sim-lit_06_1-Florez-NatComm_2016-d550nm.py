""" Replicating the results of
    Brillouin scattering self-cancellation
    Florez et al.
    http://dx.doi.org/10.1038/ncomms11759
"""

import sys

import time
import numpy as np

import matplotlib.pyplot as plt

from pathlib import Path
sys.path.append(str(Path('../backend')))

import numbat
from plotgain import Decorator
import integration
import modecalcs
import materials
from nbtypes import SI_GHz


import starter


# use this class to add or alter features to the final plots

class EMDecorator(Decorator):
    def __init__(self):
        super().__init__()

        # title_font=24
        # self._multi_sizes= {'title':title_font-2, 'subplot_title':title_font-5, 'cbar_tick':title_font-10, 'ax_tick':title_font-10, 'ax_label':title_font-10 }
        # self._single_sizes= {'ax_label':80, 'subplot_title':80, 'cbar_tick':60, 'ax_tick':70}
        # self._single_sizes= {'ax_label':60, 'subplot_title':60, 'cbar_tick':40, 'ax_tick':40}
        # self._is_single=True

    def extra_axes_commands(self, ax):
        circle1 = plt.Circle((0, 0), 0.275, color='black', fill=False)
        ax.add_artist(circle1)


class ACDecorator(Decorator):
    def __init__(self):
        super().__init__()
        # title_font=24
        # self._multi_sizes= {'title':title_font-2, 'subplot_title':title_font-5, 'cbar_tick':title_font-10, 'ax_tick':title_font-10, 'ax_label':title_font-10 }
#
#    self._single_sizes= {'ax_label':30, 'subplot_title':40, 'cbar_tick':20, 'ax_tick':30, 'title_pad':20}
#    self._single_sizes= {'ax_label':80, 'subplot_title':80, 'cbar_tick':60, 'ax_tick':70, 'title_pad':25}
#    self._is_single=True

    def extra_axes_commands(self, ax):
        circle1 = plt.Circle((0, 0), 0.275, color='black', fill=False)
        ax.add_artist(circle1)


emdecorate = EMDecorator()
acdecorate = ACDecorator()


start = time.time()

# Geometric Parameters - all in nm.
wl_nm = 1550
domain_x = 2*wl_nm
domain_y = domain_x
inc_a_x = 550  # Diameter
inc_a_y = inc_a_x
inc_shape = 'circular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix, refine_fac = starter.read_args(6, sys.argv, sub='a')

nbapp = numbat.NumBATApp(prefix)

# Use all specified parameters to create a waveguide object.
wguide = nbapp.make_structure(inc_shape, domain_x, domain_y, inc_a_x, inc_a_y,
                              material_bkg=materials.make_material("Vacuum"),
                              material_a=materials.make_material(
                                  "SiO2_2013_Laude"),
                              lc_bkg=.1, lc_refine_1=5*refine_fac, lc_refine_2=5.0*refine_fac)

wguide.plot_mesh(prefix)
# Expected effective index of fundamental guided mode.
n_eff = 1.4

doem = True
doac = True
new_calcs = True

if doem:
    # Calculate Electromagnetic Modes
    if new_calcs:
        sim_EM_pump = wguide.calc_EM_modes(
            num_modes_EM_pump, wl_nm, n_eff=n_eff)
        np.savez(prefix+'-wguide_data_florez', sim_EM_pump=sim_EM_pump)
    else:
        npzfile = np.load(prefix+'-wguide_data_florez.npz', allow_pickle=True)
        sim_EM_pump = npzfile['sim_EM_pump'].tolist()
    sim_EM_Stokes = sim_EM_pump.bkwd_Stokes_modes()

    sim_EM_pump.plot_modes(xlim_min=0.2, xlim_max=0.2, ivals=range(6),
                        ylim_min=0.2, ylim_max=0.2, decorator=emdecorate,
                        )

    # Print the wavevectors of EM modes.
    kzs = sim_EM_pump.kz_EM_all()
    print('k_z of EM modes \n', np.round(np.real(kzs), 4))

    # Calculate the EM effective index of the waveguide.
    n_eff_sim = np.real(sim_EM_pump.neff_all())
    print("n_eff = ", np.round(n_eff_sim, 4))

if not doac:
    sys.exit(0)

q_AC = np.real(sim_EM_pump.kz_EM_all()[
               EM_ival_pump] - sim_EM_Stokes.kz_EM_all()[EM_ival_Stokes])

shift_Hz = 4*SI_GHz

# Calculate Acoustic Modes
if new_calcs:
    sim_AC = wguide.calc_AC_modes(
        num_modes_AC, q_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)
    # np.savez('wguide_data_florez_AC', sim_AC=sim_AC)  # triangulation pickle
else:
    npzfile = np.load('wguide_data_florez_AC.npz', allow_pickle=True)
    sim_AC = npzfile['sim_AC'].tolist()

sim_AC.plot_modes(ivals=range(10),
                    xlim_min=-.1, ylim_min=-.1, xlim_max=-.1, ylim_max=-.1,
                    decorator=acdecorate)

# Print the frequencies of AC modes.
print('Freq of AC modes (GHz) \n', np.round(
    np.real(sim_AC.nu_AC_all())*1e-9, 4))

set_q_factor = 1000.

# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB.
gain_box = integration.get_gains_and_qs(sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5*SI_GHz
freq_max = 12*SI_GHz
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True)

freq_min = 5.06*SI_GHz
freq_max = 6.0*SI_GHz
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True, suffix='-5')

freq_min = 6.28*SI_GHz
freq_max = 6.32*SI_GHz
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True, suffix='-6')


freq_min = 8.09*SI_GHz
freq_max = 8.13*SI_GHz
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True, suffix='-8')

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 11.65*SI_GHz
freq_max = 11.69*SI_GHz
gain_box.plot_spectra(freq_min=freq_min, freq_max=freq_max, logy=True, suffix='-11')

print(nbapp.final_report())
