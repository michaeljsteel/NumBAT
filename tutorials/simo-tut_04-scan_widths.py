""" Calculate the backward SBS gain spectra as a function of
    waveguide width, for silicon waveguides surrounded in air.

    Also shows how to use python multiprocessing library.
"""

import os
import sys
from multiprocessing import Pool
import threading
import numpy as np

import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
#from matplotlib.colors import colorConverter

sys.path.append("../backend/")

import numbat
import materials
import mode_calcs
import integration
import plotting
import starter



# Select the number of CPUs to use in simulation.
num_cores = os.cpu_count()
num_cores = 2

# Geometric Parameters - all in nm.
lambda_nm = 1550
inc_shape = 'rectangular'

num_modes_EM_pump = 30
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 50
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix, refine_fac = starter.read_args(4, sys.argv, refine=3)

nbapp = numbat.NumBATApp(prefix, prefix+'-out')

use_multiproc = num_cores >1  and not nbapp.is_macos()

# Width previous simo's done for, with known meshing params
known_geo = 315.


def modes_n_gain(wwguide):
    width, wguide = wwguide
    thread_nm = threading.current_thread().name
    print('Process %d, thread %s: commencing mode calculation for width a_x = %f' % (
        os.getpid(), thread_nm, wguide.inc_a_x))
    # Expected effective index of fundamental guided mode.
    n_eff = (wguide.get_material('a').refindex_n-0.1) * \
        wguide.inc_a_x/known_geo

    # Calculate Electromagnetic modes.
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, lambda_nm, n_eff)
    sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
    q_AC = np.real(sim_EM_pump.kz_EM(EM_ival_pump) -
                   sim_EM_Stokes.kz_EM(EM_ival_Stokes))

    # Calculate Acoustic modes.
    sim_AC = wguide.calc_AC_modes(num_modes_AC, q_AC, EM_sim=sim_EM_pump)

    # Calculate interaction integrals and SBS gain.
    SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
        sim_EM_pump, sim_EM_Stokes, sim_AC, q_AC,
        EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)

    print('Process %d, thread %s: completed mode calculation for width a_x = %.3f' % (
        os.getpid(), thread_nm, wguide.inc_a_x))
    return [width, sim_EM_pump, sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz]


n_widths = 6
waveguide_widths = np.linspace(300, 350, n_widths)
l_wguides = []
# Scale meshing to new structures.
for width in waveguide_widths:
    msh_ratio = width/known_geo
    unitcell_x = 1.5*lambda_nm*msh_ratio
    unitcell_y = unitcell_x
    inc_a_x = width
    inc_a_y = 0.9*inc_a_x

    wguide = nbapp.make_structure(inc_shape, unitcell_x, unitcell_y, inc_a_x, inc_a_y,
                               material_bkg=materials.make_material("Vacuum"),
                               material_a=materials.make_material(
                                   "Si_2016_Smith"),
                               lc_bkg=.1, lc_refine_1=8.0*refine_fac, lc_refine_2=8.0*refine_fac)

    l_wguides.append((width, wguide))

    # wguide.plot_mesh(prefix+'_%d'%int(width))

new_calcs = True
if new_calcs:
  # Run widths in parallel across num_cores CPUs using multiprocessing package.
    if use_multiproc:
        with Pool(num_cores) as pool:
            l_width_data = pool.map(modes_n_gain, l_wguides)
    else:
        l_width_data = list(map(modes_n_gain, l_wguides))

    v_width_data = np.array(l_width_data, dtype=object)

    # This generates a warning abut ragged nested sequences. Is there an option to pool.map that would clean this up?
    np.savez('%s_simo_results' % prefix, width_objs=v_width_data)

else:
    npzfile = np.load(f'{prefix}_simo_results.npz', allow_pickle=True)
    v_width_data = npzfile['width_objs'].tolist()

n_effs = []
freqs_gains = []
interp_grid_points = 10000
int_min = 10e9
int_max = 26e9
interp_grid = np.linspace(int_min, int_max, interp_grid_points)

for i_w, width_obj in enumerate(v_width_data):
    interp_values = np.zeros(interp_grid_points)
    width, sim_EM_pump, sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz = width_obj

    sim_EM_pump.plot_modes(suffix='_wid_%d' % i_w, ivals=range(5))

    sim_AC.plot_modes(suffix='_wid_%d' % i_w, ivals=range(20))

    # Calculate the EM effective index of the waveguide (q_AC = 2*k_EM).
    # np.round(np.real((q_AC/2.)*((lambda_nm*1e-9)/(2.*np.pi))), 4)
    n_eff_sim = sim_EM_pump.neff(0)
    n_effs.append(n_eff_sim)

    # Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
    freq_min = np.real(sim_AC.nu_AC_all()[0]) - 5e9  # Hz
    freq_max = np.real(sim_AC.nu_AC_all()[-1]) + 5e9  # Hz
    decorator = plotting.Decorator()
    decorator.set_title(f'Gain for width $w={width:.2f}.2f$ nm')
    plotting.plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
                               EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
                                suffix=f'_wscan_{i_w}' ,  # include scan step in file name
                               decorator=decorator)

    # Repeat calc to collect data for waterfall plot.
    tune_steps = 50000
    tune_range = 10e9  # GHz
    detuning_range = np.append(np.linspace(-1*tune_range, 0, tune_steps),
                               np.linspace(0, tune_range, tune_steps)[1:])  # Hz

    # Linewidth of Lorentzian is half the FWHM style linewidth.
    linewidth = linewidth_Hz/2
    num_modes = len(linewidth)
    for m in range(num_modes):
        gain_list = np.real(SBS_gain[EM_ival_Stokes, EM_ival_pump, m]
                            * linewidth[m]**2/(linewidth[m]**2 + detuning_range**2))
        freq_list = np.real(sim_AC.nu_AC(m) + detuning_range)
        interp_spectrum = np.interp(interp_grid, freq_list, gain_list)
        interp_values += interp_spectrum
    freqs_gains.append(list(zip(interp_grid*1e-9, abs(interp_values))))

print('Widths', waveguide_widths)
print('n_effs', n_effs)

# Plot a 'waterfall' plot.
ax = plt.figure().add_subplot(projection='3d')
poly = PolyCollection(freqs_gains)
poly.set_alpha(0.7)
ax.add_collection3d(poly, zs=waveguide_widths, zdir='y')
ax.set_xlabel('Frequency (GHz)', fontsize=14)
ax.set_xlim3d(int_min*1e-9, int_max*1e-9)
ax.set_ylabel('Width (nm)', fontsize=14)
ax.set_ylim3d(waveguide_widths[0], waveguide_widths[-1])
ax.set_zlabel('|Gain| (1/Wm)', fontsize=14)
ax.set_zlim3d(0, 1500)
plt.tick_params(axis='both', which='major', labelsize=12, pad=-2)
plt.savefig(nbapp.outpath()+'-gain_spectra-waterfall.png')
plt.close()


print(nbapp.final_report())
