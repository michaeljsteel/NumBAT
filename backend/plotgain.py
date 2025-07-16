# plotgain.py is a subroutine of NumBAT that contains numerous plotting
# routines.

# Copyright (C) 2017-2025  Michael Steel, Bjorn Sturmberg, Kokou Dossou.

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
import matplotlib.pyplot as plt


import numbat
from nbtypes import SI_GHz
import reporting

from plotting.plottools import save_and_close_figure
from plotmodes import Decorator, TidyAxes


def gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, v_linewidth_Hz, q_AC,
                 EM_mode_index_pump, EM_mode_index_Stokes, AC_mode_index, freq_min=0., freq_max=50e9,
                 num_interp_pts=0, dB=False, dB_peak_amp=10, mode_comps=False, logy=False,
                 save_txt=False, prefix='', suffix='', decorator=None,
                 show_gains='All'):
    reporting.deprecated_function('plotgain.gain_spectra',
                                  'Gain.plot_spectra')
    # reporting.report_and_exit('The function plotgain.gain_spectra() is deprecated.\n'
    #                           + 'Please now use Gain.plot_spectra() and observe the changes in calling convention described in the Release notes.')


def plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, v_linewidth_Hz,
                      EM_mode_index_pump, EM_mode_index_Stokes, AC_mode_index='All',
                      freq_min=0., freq_max=50e9, num_interp_pts=0,
                      dB=False, dB_peak_amp=10, mode_comps=False, logy=False,
                      save_txt=False, prefix='', suffix='', decorator=None,
                      show_gains='All', mark_modes_threshold=0.02):
    r""" Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.


        Args:
            sim_AC : An AC ``Struct`` instance that has had calc_modes calculated

            SBS_gain  (array): Totlat SBS gain of modes.

            SBS_gain_PE  (array): Moving Bountary gain of modes.

            SBS_gain_MB  (array): Photoelastic gain of modes.

            v_linewidth_Hz  (array): Linewidth of each mode [Hz].

            EM_mode_index_pump  (int or 'All'): Which EM pump mode(s) to consider.

            EM_mode_index_Stokes  (int or 'All'): Which EM Stokes mode(s) to consider.

            AC_mode_index  (int or 'All'):  Which AC mode(s) to consider.

            freq_min  (float): Minimum of frequency range.

            freq_max  (float): Maximum of frequency range.

        Keyword Args:
            num_interp_pts  (int): Number of frequency points to interpolate to.

            dB  (bool): Save a set of spectra in dB units.

            dB_peak_amp  (float): Set the peak amplitude of highest gain mode in dB.

            mode_comps  (bool): Plot decomposition of spectra into individual modes.

            logy  (bool): PLot y-axis on log scale.

            save_txt  (bool): Save spectra data to txt file.

            prefix  (str): String to be appended to start of file name.

            suffix  (str): String to be appended to end of file name.
    """

    nbapp = numbat.NumBATApp()

    plot_prefs = numbat.NumBATPlotPrefs()
    plot_fileext = plot_prefs._plot_extension

    pathpref = nbapp.outpath(prefix) + '-gain_spectra'

    if decorator is None:
        decorator = Decorator()
    decorator.set_for_single()

    # This is expensive but helps to get the peaks in most cases.
    tune_steps = 50000
    # Need to do smarter interpolation to make this number lower
    tune_range = 10 * SI_GHz

    # Construct an odd range of freqs guaranteed to include central resonance frequency.
    if tune_steps % 2 == 0:
        tune_steps += 1
    detuning_range = np.linspace(-tune_range, tune_range, tune_steps)
    detran2 = detuning_range**2


    # Linewidth of Lorentzian is half the FWHM style linewidth.
    v_linewidth = v_linewidth_Hz/2
    num_modes = len(v_linewidth)

    # use Linewdith to pick good number of interp points.
    if num_interp_pts == 0:
        num_interp_pts = 5000
        # min_lw = np.min(v_linewidth)
        # if min_lw == 0.0:
        #     num_interp_pts = 5000
        # else: # put 5 points across the narrowest peak
        #     d_nu = min_lw/5
        #     print('got dnu', d_nu,freq_max-freq_min)
        #     num_interp_pts = int((freq_max-freq_min)/d_nu)
        #     num_interp_pts = max(num_interp_pts, 1000)

    nu_grid = np.linspace(freq_min, freq_max, num_interp_pts)

    nu_grid_GHz = nu_grid / SI_GHz
    nu_min_GHz = freq_min / SI_GHz if freq_min else nu_grid_GHz[0]
    nu_max_GHz = freq_max / SI_GHz if freq_max else nu_grid_GHz[-1]


    # Plot decomposition of spectra into individual modes.
    v_gain_global = np.zeros(num_interp_pts)

    fig, ax = plt.subplots()

    tidy = TidyAxes(nax=1)
    mode_fs = tidy._props['mode_index_label_fs']

    lw = decorator.get_property('ax_linewidth')
    fs = decorator.get_property('ax_label_fs')
    ts = decorator.get_property('ax_tick_fs')

    mode_indices = []
    if AC_mode_index is None or AC_mode_index == 'All':
        mode_indices = range(num_modes)
    else:
        try:
            mode_indices = [i for i in AC_mode_index if 0 <= i < num_modes]
        except Exception:
            reporting.report_and_exit(
                'AC_mode_index in gain_spectra() must be "All" or a list of mode numbers.')

    v_gain_global = {'Total': np.zeros(num_interp_pts),
                     'PE': np.zeros(num_interp_pts),
                     'MB': np.zeros(num_interp_pts)}

    # v_gain_global_tot = np.zeros(num_interp_pts)
    # v_gain_global_PE = np.zeros(num_interp_pts)
    # v_gain_global_MB = np.zeros(num_interp_pts)

    show_mode_indices = mark_modes_threshold > 0.0

    if show_gains not in ('All', 'PE', 'MB', 'Total'):
        show_gains = 'All'
    do_tot = show_gains in ('Total', 'All')
    do_PE = show_gains in ('PE', 'All')
    do_MB = show_gains in ('MB', 'All')

    # Find max gain over the freq range that will be plotted
    mask = np.where((sim_AC.nu_AC_all() > freq_min) *
                    (sim_AC.nu_AC_all() < freq_max), 1, 0)

    gain = np.max(np.abs(np.real(mask*SBS_gain[EM_mode_index_pump, EM_mode_index_Stokes, :])))
    gain_PE = np.max(np.abs(np.real(mask*SBS_gain_PE[EM_mode_index_pump, EM_mode_index_Stokes, :])))
    gain_MB = np.max(np.abs(np.real(mask*SBS_gain_MB[EM_mode_index_pump, EM_mode_index_Stokes, :])))

    if do_tot:
        maxG = np.max(np.array([gain, gain_PE, gain_MB]))
    elif do_PE:
        maxG = gain_PE
    else:
        maxG = gain_MB

    fig, ax = plt.subplots()

    modelabs = []
    modelabs_logy = []

    for m in mode_indices:
        nu0_m = np.real(sim_AC.nu_AC(m))
        nu0_m_GHz = nu0_m / SI_GHz

        if not freq_min < nu0_m < freq_max:
            continue
        v_nu_loc = nu0_m + detuning_range

        v_Lorentz = v_linewidth[m]**2/(v_linewidth[m]**2 + detran2)

        gsets = (
            (do_tot, SBS_gain, 'ob', 'Total'),
            (do_PE, SBS_gain_PE, 'or', 'PE' ),
            (do_MB, SBS_gain_MB, 'og', 'MB'),
            )

        gain_peaks={}
        for (do_g, m_gain, style, gcont) in gsets:
            t_gain = abs(np.real(m_gain[EM_mode_index_pump, EM_mode_index_Stokes, m]))
            v_gain_loc = t_gain * v_Lorentz
            interp_spectrum = np.interp(nu_grid, v_nu_loc, v_gain_loc)
            v_gain_global[gcont] += interp_spectrum
            gain_peaks[gcont] = t_gain

            # mark peak gain with a dot
            if t_gain > mark_modes_threshold*maxG:
                ax.plot(nu0_m_GHz, t_gain, style)

        gain_peaks['All'] = max(gain_peaks.values())

        # mark modes with gains larger than 5% of the maximum found
        # Marks the largest of the dots for a given mode
        if show_mode_indices:
            Gm = gain_peaks[show_gains]  # extract the desired gain
            xloc = nu0_m_GHz + .01*(nu_max_GHz-nu_min_GHz)

            if nu_min_GHz < xloc < nu_max_GHz:
                if Gm > mark_modes_threshold*maxG:
                    ax.text(xloc, abs(Gm), m, fontsize=mode_fs, horizontalalignment='left',
                        verticalalignment='top')
                    modelabs.append((xloc, Gm, m))

                # keep extra small ones for the log display
                if Gm > mark_modes_threshold*maxG/1000:
                    modelabs_logy.append((xloc, Gm, m))


    gpsets = ( (do_PE, 'PE', 'r', 'PE'),
               (do_MB, 'MB', 'g', 'MB'),
               (do_tot,'Total', 'b', 'Total gain') )

    # plain graph
    for do_g, gcont, col, lab in gpsets:
        if do_g:
            ax.plot(nu_grid_GHz, np.abs(v_gain_global[gcont]), col, linewidth=lw, label=lab)

    ax.legend(loc=0)
    ax.set_xlim(nu_min_GHz, nu_max_GHz)
    ax.set_xlabel('Frequency (GHz)', size=fs)
    ax.set_ylabel('|Gain| (1/Wm)', size=fs)
    ax.tick_params(labelsize=ts)

    decorator.extra_axes_commands(ax)
    decorator.add_title(ax)

    tidy.apply_to_axes(ax)

    save_and_close_figure(fig, f'{pathpref}{suffix}{plot_fileext}')


    if save_txt:
        save_array = np.array([nu_grid, v_gain_global['Total']]).T
        np.savetxt(f'{pathpref}-MB_PE_comps{suffix}-Total.csv', save_array, delimiter=',')
        save_array = np.array([nu_grid, v_gain_global['PE']]).T
        np.savetxt(f'{pathpref}-MB_PE_comps{suffix}-PE.csv',  save_array, delimiter=',')
        save_array = np.array([nu_grid, v_gain_global['MB']]).T
        np.savetxt(f'{pathpref}-MB_PE_comps{suffix}-MB.csv',  save_array, delimiter=',')

    # dB graph
    if dB:
        fig, ax = plt.subplots()

        max_G = np.max(abs(v_gain_global['Total']))
        Leff = math.log(10**(dB_peak_amp*.1))/max_G
        v_scale = dB_peak_amp * Leff*math.log10(math.exp(1.0))

        for do_g, gcont, col, lab in gpsets:
            if do_g:
                ax.plot(nu_grid_GHz, v_scale*np.abs(v_gain_global[gcont]), col, linewidth=lw, label=lab)

        ax.legend(loc=0)
        ax.set_xlim(nu_min_GHz, nu_max_GHz)
        ax.set_xlabel('Frequency (GHz)', size=fs)
        ax.set_ylabel('Amplification (dB)', size=fs)
        ax.tick_params(labelsize=ts)
        for (nuloc, Gm, m) in modelabs:
            ax.text(nuloc, abs(Gm)*v_scale, m, fontsize=mode_fs, horizontalalignment='left',
                    verticalalignment='top')

        tidy.apply_to_axes(ax)
        save_and_close_figure(fig, f'{pathpref}-dB{suffix}{plot_fileext}')

        if save_txt:
            dB_const=1.0
            save_array = (
                nu_grid, 10*np.log10(np.exp(abs(v_gain_global['Total'])*dB_const)))
            np.savetxt(f'{pathpref}-gain_spectra-dB{suffix}.csv', save_array, delimiter=',')

    # log scale graph
    if logy:
        fig, ax = plt.subplots()
        for do_g, gcont, col, lab in gpsets:
            if do_g:
                ax.semilogy(nu_grid_GHz, np.abs(v_gain_global[gcont]), col, linewidth=lw, label=lab)

        ax.legend(loc=0)
        ax.set_xlim(nu_min_GHz, nu_max_GHz)
        ax.set_xlabel('Frequency (GHz)', size=fs)
        ax.set_ylabel('|Gain| (1/Wm)', size=fs)
        ax.tick_params(labelsize=ts)
        for (nuloc, Gm, m) in modelabs_logy:
            ax.text(nuloc, abs(Gm), m, fontsize=mode_fs, horizontalalignment='left',
                    verticalalignment='top')

        tidy.apply_to_axes(ax)
        save_and_close_figure(fig, f'{pathpref}-logy{suffix}{plot_fileext}')

    return v_gain_global['Total'], v_gain_global['PE'], v_gain_global['MB']


