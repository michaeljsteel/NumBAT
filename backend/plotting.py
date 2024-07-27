# plotting.py is a subroutine of NumBAT that contains numerous plotting
# routines.

# Copyright (C) 2017  Bjorn Sturmberg.

# NumBAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import math

from pathlib import Path

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#import matplotlib.colors as mplcolors

#from mpl_toolkits.axes_grid1 import make_axes_locatable
#from collections.abc import Iterable

import numbat
from nbtypes import SI_GHz
import reporting

from plottools import save_and_close_figure
from plotmodes import Decorator


# setting this true is useful for use in jupyter style notebooks. TODO: Make a nicer interface
keep_plots_open = False

# TODO: make this some kind of global numbat setting and get it out of the startup space
try:
    # Location search path is at matplotlib.get_configdir()
    # Will load from ~/.config/matplotlib/NumBATstyle.mplstyle if found
    plt.style.use('NumBATstyle')
except (ValueError, IOError, AttributeError) as err:
    print("Preferred NumBAT matplotlib style file not found. Using matplotlib defaults.")

mycolors = [color['color'] for color in list(plt.rcParams['axes.prop_cycle'])]




def gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, q_AC,
                 EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=0., freq_max=50e9, num_interp_pts=3000,
                 dB=False, dB_peak_amp=10, mode_comps=False, logy=False,
                 pdf_png='png', save_txt=False, prefix='', suffix='', decorator=None,
                 show_gains='All'):
    reporting.report_and_exit('The function plotting.gain_spectra() is deprecated.\n'
                              + 'Please now use plotting.plot_gain_spectra() and observe the changes in calling convention described in the Release notes.')


def plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
                      EM_ival_pump, EM_ival_Stokes, AC_ival='All',
                      freq_min=0., freq_max=50e9, num_interp_pts=3000,
                      dB=False, dB_peak_amp=10, mode_comps=False, logy=False,
                      pdf_png='png', save_txt=False, prefix='', suffix='', decorator=None,
                      show_gains='All', mark_modes_threshold=0.02):
    r""" Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.


        Args:
            sim_AC : An AC ``Struct`` instance that has had calc_modes calculated

            SBS_gain  (array): Totlat SBS gain of modes.

            SBS_gain_PE  (array): Moving Bountary gain of modes.

            SBS_gain_MB  (array): Photoelastic gain of modes.

            linewidth_Hz  (array): Linewidth of each mode [Hz].

            EM_ival_pump  (int or 'All'): Which EM pump mode(s) to consider.

            EM_ival_Stokes  (int or 'All'): Which EM Stokes mode(s) to consider.

            AC_ival  (int or 'All'):  Which AC mode(s) to consider.

            freq_min  (float): Minimum of frequency range.

            freq_max  (float): Maximum of frequency range.

        Keyword Args:
            num_interp_pts  (int): Number of frequency points to interpolate to.

            dB  (bool): Save a set of spectra in dB units.

            dB_peak_amp  (float): Set the peak amplitude of highest gain mode in dB.

            mode_comps  (bool): Plot decomposition of spectra into individual modes.

            logy  (bool): PLot y-axis on log scale.

            pdf_png  (str): Save figures as 'png' or 'pdf'.

            save_txt  (bool): Save spectra data to txt file.

            prefix  (str): String to be appended to start of file name.

            suffix  (str): String to be appended to end of file name.
    """

    # TODO: there are a lot of calls to ax.plot() in here to get every resonance on its
    # own grid and matplotlib leaks badly. Do that now by just plotting the top value

    # TODO: give a no plot option 'calc_gain_spectra'

    # process = psutil.Process()

    if not prefix: prefix=numbat.NumBATApp().outprefix()

    pref = f'{prefix}-gain_spectra'
    pathpref = str(Path(numbat.NumBATApp().outdir(), pref))

    if decorator is None:
        decorator = Decorator()
    decorator.set_for_single()

    # This is expensive but helps to get the peaks in most cases.
    tune_steps = 50000
    # Need to do smarter interpolation to make this number lower
    tune_range = 10e9  # Hz

    # Construct an odd range of freqs guaranteed to include central resonance frequency.
    if tune_steps % 2 == 0:
        tune_steps += 1
    detuning_range = np.linspace(-tune_range, tune_range, tune_steps)
    detran2 = detuning_range**2

    nu_grid = np.linspace(freq_min, freq_max, num_interp_pts)

    nu_grid_GHz = nu_grid / SI_GHz
    nu_min_GHz = nu_grid_GHz[0]
    nu_max_GHz = nu_grid_GHz[-1]
    if freq_min:
        nu_min_GHz = freq_min * 1e-9
    if freq_max:
        nu_max_GHz = freq_max * 1e-9

    # Linewidth of Lorentzian is half the FWHM style linewidth.
    linewidth = linewidth_Hz/2
    num_modes = len(linewidth)

    # Plot decomposition of spectra into individual modes.
    v_gain_global = np.zeros(num_interp_pts)

    fig, ax = plt.subplots()

    lw = decorator.get_property('linewidth')
    fs = decorator.get_property('ax_label_fs')
    ts = decorator.get_property('ax_tick_fs')

    ivals = []
    if AC_ival is None or AC_ival == 'All':
        ivals = range(num_modes)
    else:
        try:
            ivals = [i for i in AC_ival if 0 <= i < num_modes]
        except Exception:
            reporting.report_and_exit(
                'AC_ival in gain_spectra() must be "All" or a list of mode numbers.')

    # Total gain via sum over all modes in the vicinity of their peak
    for m in ivals:  # TODO: this is same as top of next big loop? DELETE?
        # build lorentzian centered on mode m
        v_nu_loc = np.real(sim_AC.nu_AC(m) + detuning_range)
        v_Lorentz = linewidth[m]**2/(linewidth[m]**2 + detran2)
        v_gain_loc = np.real(SBS_gain[EM_ival_pump, EM_ival_Stokes, m]) * v_Lorentz

        if mode_comps:
            ax.plot(v_nu_loc, np.abs(v_gain_loc), linewidth=lw)
            if save_txt:
                save_array = np.array([v_nu_loc, v_gain_loc]).T
                np.savetxt(f'{pathpref}-mode_comps{suffix}-{m}.csv', save_array, delimiter=',')
        # set up an interpolation for summing all the gain peaks
        v_gain_global += np.interp(nu_grid, v_nu_loc, v_gain_loc)

    #return_interp_values = v_gain_global

    if mode_comps:  # TODO: delete this as no longer important?
        ax.plot(nu_grid_GHz, np.abs(v_gain_global),
                'b', linewidth=lw, label="Total")
        ax.legend(loc=0)
        ax.set_xlim(nu_min_GHz, nu_max_GHz)
        ax.set_xlabel('Frequency (GHz)', size=fs)
        ax.set_ylabel('|Gain| (1/Wm)', size=fs)
        ax.tick_params(labelsize=ts)
        decorator.add_title(ax)
        decorator.extra_axes_commands(ax)

        save_and_close_figure(fig, f'{pathpref}-mode_comps{suffix}.{pdf_png}')

        if save_txt:
            save_array = np.array([nu_grid, v_gain_global]).T
            np.savetxt(f'{pathpref}-mode_comps{suffix}-Total.csv' , save_array, delimiter=',')

    v_gain_global_tot = np.zeros(num_interp_pts)
    v_gain_global_PE = np.zeros(num_interp_pts)
    v_gain_global_MB = np.zeros(num_interp_pts)

    show_mode_indices = mark_modes_threshold > 0.0

    if show_gains not in ('All', 'PE', 'MB', 'Total'):
        show_gains = 'All'
    do_tot = show_gains in ('Total', 'All')
    do_PE = show_gains in ('PE', 'All')
    do_MB = show_gains in ('MB', 'All')

    maxG = 0  # this should be largest of the plotted gains
    # only maximise over the freqs that will be plotted
    mask = np.where((sim_AC.nu_AC_all() > freq_min) *
                    (sim_AC.nu_AC_all() < freq_max), 1, 0)

    gain = np.max(np.abs(np.real(mask*SBS_gain[EM_ival_pump, EM_ival_Stokes, :])))
    gain_PE = np.max(np.abs(np.real(mask*SBS_gain_PE[EM_ival_pump, EM_ival_Stokes, :])))
    gain_MB = np.max(np.abs(np.real(mask*SBS_gain_MB[EM_ival_pump, EM_ival_Stokes, :])))

    if do_tot:
        maxG = np.max(np.array([gain, gain_PE, gain_MB]))
    elif do_PE:
        maxG = gain_PE
    else:
        maxG = gain_MB

    fig, ax = plt.subplots()

    modelabs = []
    modelabs_logy = []

    for m in ivals:
        nu0_m = np.real(sim_AC.nu_AC(m))

        if not freq_min < nu0_m < freq_max:
            continue
        v_nu_loc = nu0_m + detuning_range

        v_Lorentz = linewidth[m]**2/(linewidth[m]**2 + detran2)

        if do_tot:
            gain_m = abs(np.real(SBS_gain[EM_ival_pump, EM_ival_Stokes, m]))
            v_gain_loc = gain_m * v_Lorentz
            interp_spectrum = np.interp(nu_grid, v_nu_loc, v_gain_loc)
            v_gain_global_tot += interp_spectrum
            if gain_m > mark_modes_threshold*maxG:
                ax.plot(nu0_m*1e-9, gain_m, 'ob')

        if do_PE:
            gain_PE_m = abs(np.real(SBS_gain_PE[EM_ival_pump, EM_ival_Stokes, m]))
            v_gain_loc = gain_PE_m * v_Lorentz
            v_gain_global_PE += np.interp(nu_grid, v_nu_loc, v_gain_loc)
            if gain_PE_m > mark_modes_threshold*maxG:
                ax.plot(nu0_m*1e-9, gain_PE_m, 'or')

        if do_MB:
            gain_MB_m = abs(np.real(SBS_gain_MB[EM_ival_pump, EM_ival_Stokes, m]))
            v_gain_loc = gain_MB_m * v_Lorentz
            v_gain_global_MB += np.interp(nu_grid, v_nu_loc, v_gain_loc)
            if gain_MB_m > mark_modes_threshold*maxG:
                ax.plot(nu0_m*1e-9, gain_MB_m, 'og')

        if show_mode_indices:  # mark modes with gains larger than 5% of the maximum found
            Gm = {'All': max(gain_m, gain_PE_m, gain_MB_m),
                  'Total': gain_m, 'PE': gain_PE_m, 'MB': gain_MB_m}[show_gains]
            xloc = nu0_m*1e-9 + .01*(nu_max_GHz-nu_min_GHz)
            if Gm > mark_modes_threshold*maxG and nu_min_GHz < xloc < nu_max_GHz:
                ax.text(xloc, abs(Gm), m, fontsize=fs, horizontalalignment='left',
                        verticalalignment='top')
                modelabs.append((xloc, Gm, m))
            # keep extra small ones for the log display
            if Gm > mark_modes_threshold*maxG/1000 and nu_min_GHz < xloc < nu_max_GHz:
                modelabs_logy.append((xloc, Gm, m))

    if do_PE:
        ax.plot(nu_grid_GHz, np.abs(v_gain_global_PE), 'r', linewidth=lw, label='PE')
    if do_MB:
        ax.plot(nu_grid_GHz, np.abs(v_gain_global_MB), 'g', linewidth=lw, label='MB')
    if do_tot:
        ax.plot(nu_grid_GHz, np.abs(v_gain_global_tot), 'b', linewidth=lw, label='Total gain')

    ax.legend(loc=0)
    ax.set_xlim(nu_min_GHz, nu_max_GHz)
    ax.set_xlabel('Frequency (GHz)', size=fs)
    ax.set_ylabel('|Gain| (1/Wm)', size=fs)
    ax.tick_params(labelsize=ts)

    decorator.extra_axes_commands(ax)
    decorator.add_title(ax)

    save_and_close_figure(fig, f'{pathpref}{suffix}.{pdf_png}')


    if save_txt:
        save_array = np.array([nu_grid, v_gain_global_tot]).T
        np.savetxt(f'{pathpref}-MB_PE_comps{suffix}-Total.csv', save_array, delimiter=',')
        save_array = np.array([nu_grid, v_gain_global_PE]).T
        np.savetxt(f'{pathpref}-MB_PE_comps{suffix}-PE.csv',  save_array, delimiter=',')
        save_array = np.array([nu_grid, v_gain_global_MB]).T
        np.savetxt(f'{pathpref}-MB_PE_comps{suffix}-MB.csv',  save_array, delimiter=',')


    if dB:
        fig, ax = plt.subplots()

        max_G = np.max(abs(v_gain_global_tot))
        Leff = math.log(10**(dB_peak_amp*.1))/max_G
        # dB_const = dB_peak_amp/(4.34*max_G)
        # ax.plot(nu_grid_GHz, np.abs(10*np.log10(np.exp(abs(interp_values)*dB_const))), 'b', linewidth=3, label="Total")
        v_scale = dB_peak_amp * Leff*math.log10(math.exp(1.0))
        #v_amp = v_scale * abs(v_gain_global_tot)
        if do_PE:
            ax.plot(nu_grid_GHz, v_scale*abs(v_gain_global_PE), 'r', linewidth=lw, label="PE")
        if do_MB:
            ax.plot(nu_grid_GHz, v_scale*abs(v_gain_global_MB), 'g', linewidth=lw, label="MB")
        if do_tot:
            ax.plot(nu_grid_GHz, v_scale*abs(v_gain_global_tot), 'b', linewidth=lw, label="Total")
        ax.legend(loc=0)
        ax.set_xlim(nu_min_GHz, nu_max_GHz)
        ax.set_xlabel('Frequency (GHz)', size=fs)
        ax.set_ylabel('Amplification (dB)', size=fs)
        ax.tick_params(labelsize=ts)
        for (nuloc, Gm, m) in modelabs:
            ax.text(nuloc, abs(Gm)*v_scale, m, fontsize=fs, horizontalalignment='left',
                    verticalalignment='top')

        #save_and_close_figure(fig, f'{pathpref}-gain_spectra-dB{suffix}.{pdf_png}')
        save_and_close_figure(fig, f'{pathpref}-dB{suffix}.{pdf_png}')

        if save_txt:
            dB_const=1.0
            save_array = (
                nu_grid, 10*np.log10(np.exp(abs(v_gain_global_tot)*dB_const)))
            np.savetxt(f'{pathpref}-gain_spectra-dB{suffix}.csv', save_array, delimiter=',')

    if logy:
        fig, ax = plt.subplots()
        ax.semilogy(nu_grid_GHz, abs(v_gain_global_PE), 'r', linewidth=lw, label="PE")
        ax.semilogy(nu_grid_GHz, abs(v_gain_global_MB), 'g', linewidth=lw, label="MB")
        ax.semilogy(nu_grid_GHz, abs(v_gain_global_tot), 'b', linewidth=lw, label="Total")
        ax.legend(loc=0)
        ax.set_xlim(nu_min_GHz, nu_max_GHz)
        ax.set_xlabel('Frequency (GHz)', size=fs)
        ax.set_ylabel('|Gain| (1/Wm)', size=fs)
        ax.tick_params(labelsize=ts)
        for (nuloc, Gm, m) in modelabs_logy:
            ax.text(nuloc, abs(Gm), m, fontsize=fs, horizontalalignment='left',
                    verticalalignment='top')

        save_and_close_figure(fig, f'{pathpref}-logy{suffix}.{pdf_png}')

    return v_gain_global, v_gain_global_PE, v_gain_global_MB


