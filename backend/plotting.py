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
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors

from mpl_toolkits.axes_grid1 import make_axes_locatable


import numbat
from nbtypes import QAcMethod, FieldType, component_t
from fortran import NumBAT
import reporting

# setting this true is useful for use in jupyter style notebooks. TODO: Make a nicer interface
keep_plots_open = False

# TODO: make this some kind of global numbat setting and get it out of the startup space
try:
    # Will load from ~/.config/matplotlib/NumBATstyle.mplstyle if found
    plt.style.use('NumBATstyle')
except (ValueError, IOError, AttributeError):
    print("Preferred NumBAT matplotlib style file not found. Using matplotlib defaults.")
mycolors = [color['color'] for color in list(plt.rcParams['axes.prop_cycle'])]


def savefig(fig, fname):
    fig.savefig(fname)


def plot_filename(plps, ival, label=None):
    filestart = '%(pre)s-fields/%(s)s_field_%(i)02i%(add)s' % {'pre': plps['prefix'],
                                                               's': plps['EM_AC'].name, 'i': ival, 'add': plps['suffix']}
    if label is not None:
        filestart += '_'+label

    if plps['pdf_png'] == 'png':
        fig_filename = filestart+'.png'
    elif plps['pdf_png'] == 'pdf':
        fig_filename = filestart+'.pdf'
    else:
        raise ValueError("pdf_png must be either 'png' or 'pdf'.")

    return fig_filename


class Decorator(object):

    def __init__(self):
        mp_base_font = 18

        self._multiplot_fontsizes = {'title': mp_base_font-2, 'subplot_title': mp_base_font-5,
                                     'ax_label': mp_base_font-10, 'data_label': mp_base_font-8,
                                     'ax_tick': mp_base_font-10, 'cbar_tick': mp_base_font-12}
        self._multiplot_axesprops = {'figsize': (10, 8),
                                     'subplots_hspace': .2, 'subplots_wspace': .2,
                                     'linewidth': '.75', 'edgecolor': 'gray', 'title_pad': 5, 'cbar_size': '5%', 'cbar_pad': '4%'}

        sp_base_font = 24
        self._singleplot_fontsizes = {
            'title': sp_base_font-2, 'ax_label': 20, 'subplot_title': 25, 'cbar_tick': 20, 'ax_tick': 20}
        self._singleplot_axesprops = {
            'linewidth': '.75', 'edgecolor': 'gray', 'title_pad': 20, 'cbar_size': '5%', 'cbar_pad': '2%'}

        self._is_single = True
        self._cmap_limits = {}
        self._title = ''
        self._frame_drawer = None

    def set_frame_drawer(self, fd):
        self._frame_drawer = fd

    def add_frame(self, ax):
        if self._frame_drawer is not None:
            self._frame_drawer.draw_mpl_frame(ax)

    def _fontsizes(self):
        if self._is_single:
            return self._singleplot_fontsizes
        else:
            return self._multiplot_fontsizes

    def _get_axes_prop(self):
        if self._is_single:
            return self._singleplot_axesprops
        else:
            return self._multiplot_axesprops

    def _set_for_single(self):
        self._is_single = True

    def _set_for_multi(self):
        self._is_single = False

    def set_cmap_limits(self, d):
        '''Specify the lower and upper contour plot limits for a field component plot.

          :param dict d: Dictionary mapping component ('x','y','z') to (zlo, zhi) tuple for contour plots.
          '''
        self._cmap_limits.update(d)

    def get_cmap_limits(self, comp):
        return self._cmap_limits.get(comp)

    def get_font_size(self, lab):
        fs = 10
        try:
            fs = self._fontsizes()[lab]
        except Exception:
            print(f'Warning: unknown fontsize label "{lab}" in Decorator::get_font_size()')
        return fs

    def get_axes_property(self, lab):
        return self._get_axes_prop().get(lab, None)
        # prop=''
        # try:
        #  prop=self._get_axes_prop()[lab]
        # except:
        #  print('Warning: unknown axes property label "{0}" in Decorator::get_axes_property()'.format(lab))
        #  print (self._get_axes_prop())
        #  sys.exit(1)
        # return prop

    def is_single_plot(self):
        '''Returns True if this Decorator is for a single axes plot such as a spectrum or spatial map of a single field component.
           '''
        return self._is_single

    def set_singleplot_axes_property(self, label, prop):
        '''Add or override an axes property for a single plot corresponding to the given label.'''
        self._singleplot_axesprops[label] = prop

    def set_multiplot_axes_property(self, label, prop):
        '''Add or override an axes property for a multiple axes plot corresponding to the given label.'''
        self._multiplot_axesprops[label] = prop

    def set_singleplot_fontsize(self, label, sz):
        '''Override a font size for a single plot corresponding to the given label.'''
        self._singleplot_fontsizes[label] = sz

    def set_multiplot_fontsize(self, label, sz):
        '''Override a font size for a mutiple axes plot corresponding to the given label.'''
        self._multiplot_fontsizes[label] = sz

    def set_title(self, t): self._title = t

    def add_title(self, ax):
        if self._title:
            ax.set_title(self._title)

    def extra_axes_commands(self, ax):
        '''Function to make additions to a standard plot.

          This function is called after all other plotting operations are performed.
          The default implementation does nothing. By subclassing :Decorator: and implementing this function,
          users may add extra features to a plot.
          '''
        pass


def gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, q_AC,
                 EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=0., freq_max=50e9, num_interp_pts=3000,
                 dB=False, dB_peak_amp=10, mode_comps=False, semilogy=False,
                 pdf_png='png', save_txt=False, prefix='', suffix='', decorator=None,
                 show_gains='All'):
    reporting.report_and_exit('The function plotting.gain_spectra() is deprecated.\n'
                              + 'Please now use plotting.plot_gain_spectra() and observe the changes in calling convention described in the Release notes.')


def plot_gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz,
                      EM_ival_pump, EM_ival_Stokes, AC_ival='All',
                      freq_min=0., freq_max=50e9, num_interp_pts=3000,
                      dB=False, dB_peak_amp=10, mode_comps=False, semilogy=False,
                      pdf_png='png', save_txt=False, prefix='', suffix='', decorator=None,
                      show_gains='All', mark_modes_thresh=0.02):
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

            semilogy  (bool): PLot y-axis on log scale.

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

    if decorator is None:
        decorator = Decorator()
    decorator._set_for_single()

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

    nu_grid_GHz = nu_grid * 1e-9
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

    lw = decorator.get_axes_property('linewidth')
    fs = decorator.get_font_size('ax_label')
    ts = decorator.get_font_size('ax_tick')

    ivals = []
    if AC_ival is None or AC_ival == 'All':
        ivals = range(num_modes)
    else:
        try:
            ivals = [i for i in AC_ival if (0 <= i < num_modes)]
        except Exception:
            reporting.report_and_exit(
                'AC_ival in gain_spectra() must be "All" or a list of mode numbers.')

    # Total gain via sum over all modes in the vicinity of their peak
    for m in ivals:  # TODO: this is same as top of next big loop? DELETE?
        # build lorentzian centered on mode m
        v_nu_loc = np.real(sim_AC.nu_AC(m) + detuning_range)
        v_Lorentz = linewidth[m]**2/(linewidth[m]**2 + detran2)
        v_gain_loc = np.real(
            SBS_gain[EM_ival_pump, EM_ival_Stokes, m]) * v_Lorentz

        if mode_comps:
            ax.plot(v_nu_loc, np.abs(v_gain_loc), linewidth=lw)
            if save_txt:
                save_array = np.array([v_nu_loc, v_gain_loc]).T
                np.savetxt('%(pre)s-mode_comps%(add)s-%(mode)i.csv'
                           % {'pre': pref, 'add': suffix, 'mode': m},
                           save_array, delimiter=',')
        # set up an interpolation for summing all the gain peaks
        v_gain_global += np.interp(nu_grid, v_nu_loc, v_gain_loc)

    return_interp_values = v_gain_global

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

        fig.savefig('%(pre)s-mode_comps%(add)s.%(png)s' % {
            'pre': pref, 'add': suffix, 'png': pdf_png})

        if save_txt:
            save_array = np.array([nu_grid, v_gain_global]).T
            np.savetxt('%(pre)s-mode_comps%(add)s-Total.csv'
                       % {'pre': pref, 'add': suffix}, save_array, delimiter=',')

    v_gain_global_tot = np.zeros(num_interp_pts)
    v_gain_global_PE = np.zeros(num_interp_pts)
    v_gain_global_MB = np.zeros(num_interp_pts)

    show_mode_indices = mark_modes_thresh > 0.0

    if show_gains not in ('All', 'PE', 'MB', 'Total'):
        show_gains = 'All'
    do_tot = show_gains in ('Total', 'All')
    do_PE = show_gains in ('PE', 'All')
    do_MB = show_gains in ('MB', 'All')

    maxG = 0  # this should be largest of the plotted gains
    # only maximise over the freqs that will be plotted
    mask = np.where((sim_AC.nu_AC_all() > freq_min) *
                    (sim_AC.nu_AC_all() < freq_max), 1, 0)

    gain = np.max(
        np.abs(np.real(mask*SBS_gain[EM_ival_pump, EM_ival_Stokes, :])))
    gain_PE = np.max(
        np.abs(np.real(mask*SBS_gain_PE[EM_ival_pump, EM_ival_Stokes, :])))
    gain_MB = np.max(
        np.abs(np.real(mask*SBS_gain_MB[EM_ival_pump, EM_ival_Stokes, :])))

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
            if gain_m > mark_modes_thresh*maxG:
                ax.plot(nu0_m*1e-9, gain_m, 'ob')

        if do_PE:
            gain_PE_m = abs(
                np.real(SBS_gain_PE[EM_ival_pump, EM_ival_Stokes, m]))
            v_gain_loc = gain_PE_m * v_Lorentz
            v_gain_global_PE += np.interp(nu_grid, v_nu_loc, v_gain_loc)
            if gain_PE_m > mark_modes_thresh*maxG:
                ax.plot(nu0_m*1e-9, gain_PE_m, 'or')

        if do_MB:
            gain_MB_m = abs(
                np.real(SBS_gain_MB[EM_ival_pump, EM_ival_Stokes, m]))
            v_gain_loc = gain_MB_m * v_Lorentz
            v_gain_global_MB += np.interp(nu_grid, v_nu_loc, v_gain_loc)
            if gain_MB_m > mark_modes_thresh*maxG:
                ax.plot(nu0_m*1e-9, gain_MB_m, 'og')

        if show_mode_indices:  # mark modes with gains larger than 5% of the maximum found
            Gm = {'All': max(gain_m, gain_PE_m, gain_MB_m),
                  'Total': gain_m, 'PE': gain_PE_m, 'MB': gain_MB_m}[show_gains]
            xloc = nu0_m*1e-9 + .01*(nu_max_GHz-nu_min_GHz)
            if Gm > mark_modes_thresh*maxG and nu_min_GHz < xloc < nu_max_GHz:
                ax.text(xloc, abs(Gm), m, fontsize=fs, horizontalalignment='left',
                        verticalalignment='top')
                modelabs.append((xloc, Gm, m))
            # keep extra small ones for the log display
            if Gm > mark_modes_thresh*maxG/1000 and nu_min_GHz < xloc < nu_max_GHz:
                modelabs_logy.append((xloc, Gm, m))

    if do_PE:
        ax.plot(nu_grid_GHz, np.abs(v_gain_global_PE),
                'r', linewidth=lw, label='PE')
    if do_MB:
        ax.plot(nu_grid_GHz, np.abs(v_gain_global_MB),
                'g', linewidth=lw, label='MB')
    if do_tot:
        ax.plot(nu_grid_GHz, np.abs(v_gain_global_tot),
                'b', linewidth=lw, label='Total gain')

    ax.legend(loc=0)
    ax.set_xlim(nu_min_GHz, nu_max_GHz)
    ax.set_xlabel('Frequency (GHz)', size=fs)
    ax.set_ylabel('|Gain| (1/Wm)', size=fs)
    ax.tick_params(labelsize=ts)

    decorator.extra_axes_commands(ax)
    decorator.add_title(ax)

    # fig.tight_layout()

    savefig(fig, '%(pre)s%(add)s.%(png)s' %
            {'pre': pref, 'add': suffix, 'png': pdf_png})
    plt.close(fig)

    if save_txt:
        save_array = np.array([nu_grid, v_gain_global_tot]).T
        np.savetxt('%(pre)s-MB_PE_comps%(add)s-Total.csv' % {'pre': pref, 'add': suffix},
                   save_array, delimiter=',')
        save_array = np.array([nu_grid, v_gain_global_PE]).T
        np.savetxt('%(pre)s-MB_PE_comps%(add)s-PE.csv' % {'pre': pref, 'add': suffix},
                   save_array, delimiter=',')
        save_array = np.array([nu_grid, v_gain_global_MB]).T
        np.savetxt('%(pre)s-MB_PE_comps%(add)s-MB.csv' % {'pre': pref, 'add': suffix},
                   save_array, delimiter=',')

    if dB:
        fig, ax = plt.subplots()

        max_G = np.max(abs(v_gain_global_tot))
        Leff = math.log(10**(dB_peak_amp*.1))/max_G
        # dB_const = dB_peak_amp/(4.34*max_G)
        # ax.plot(nu_grid_GHz, np.abs(10*np.log10(np.exp(abs(interp_values)*dB_const))), 'b', linewidth=3, label="Total")
        v_scale = dB_peak_amp * Leff*math.log10(math.exp(1.0))
        v_amp = v_scale * abs(v_gain_global_tot)
        if do_PE:
            ax.plot(nu_grid_GHz, v_scale*abs(v_gain_global_PE),
                    'r', linewidth=lw, label="PE")
        if do_MB:
            ax.plot(nu_grid_GHz, v_scale*abs(v_gain_global_MB),
                    'g', linewidth=lw, label="MB")
        if do_tot:
            ax.plot(nu_grid_GHz, v_scale*abs(v_gain_global_tot),
                    'b', linewidth=lw, label="Total")
        ax.legend(loc=0)
        ax.set_xlim(nu_min_GHz, nu_max_GHz)
        ax.set_xlabel('Frequency (GHz)', size=fs)
        ax.set_ylabel('Amplification (dB)', size=fs)
        ax.tick_params(labelsize=ts)
        for (nuloc, Gm, m) in modelabs:
            ax.text(nuloc, abs(Gm)*v_scale, m, fontsize=fs, horizontalalignment='left',
                    verticalalignment='top')

        fig.savefig('%(pre)s-gain_spectra-dB%(add)s.%(png)s' % {
            'pre': prefix, 'add': suffix, 'png': pdf_png})
        plt.close(fig)

        if save_txt:
            save_array = (
                nu_grid, 10*np.log10(np.exp(abs(v_gain_global_tot)*dB_const)))
            np.savetxt('%(pre)s-gain_spectra-dB%(add)s.csv'
                       % {'pre': prefix, 'add': suffix},
                       save_array, delimiter=',')

    if semilogy:
        fig, ax = plt.subplots()
        ax.semilogy(nu_grid_GHz, abs(v_gain_global_PE),
                    'r', linewidth=lw, label="PE")
        ax.semilogy(nu_grid_GHz, abs(v_gain_global_MB),
                    'g', linewidth=lw, label="MB")
        ax.semilogy(nu_grid_GHz, abs(v_gain_global_tot),
                    'b', linewidth=lw, label="Total")
        ax.legend(loc=0)
        ax.set_xlim(nu_min_GHz, nu_max_GHz)
        ax.set_xlabel('Frequency (GHz)', size=fs)
        ax.set_ylabel('|Gain| (1/Wm)', size=fs)
        ax.tick_params(labelsize=ts)
        for (nuloc, Gm, m) in modelabs_logy:
            ax.text(nuloc, abs(Gm), m, fontsize=fs, horizontalalignment='left',
                    verticalalignment='top')

        fig.savefig('%(pre)s-gain_spectra-logy%(add)s.%(png)s' % {
            'pre': prefix, 'add': suffix, 'png': pdf_png})
        plt.close(fig)

    return v_gain_global, v_gain_global_PE, v_gain_global_MB


def plot_set_ticks(ax, plps, decorator):
    if plps['ticks']:
        ax.tick_params(labelsize=decorator.get_font_size('ax_tick'))
        ax.xaxis.set_tick_params(width=1.0)
        ax.yaxis.set_tick_params(width=1.0)
        ax.set_xlabel('$x$ [μm]', size=decorator.get_font_size('ax_label'))
        ax.set_ylabel('$y$ [μm]', size=decorator.get_font_size('ax_label'))
        ax.grid(False)
    else:
        ax.set_xticks([])
        ax.set_yticks([])


def plot_set_axes_style(ax, plps, decorator):
    if decorator.get_axes_property('axes.linewidth') is not None:
        plt.setp(ax.spines.values(),
                 linewidth=decorator.get_axes_property('axes.linewidth'))


def plot_set_title(ax, comp_label, plps, decorator):
    if plps.get('add_title', True):
        ax.set_title(comp_label, fontsize=decorator.get_font_size('subplot_title'),
                     pad=decorator.get_axes_property('title_pad'))


def get_quiver_skip_range(npts, skip):
    if npts % 2:
        #mid = int((npts-1)/2)
        j0 = (int(round((npts-1)/2))) % skip
    else:
        if skip % 2 == 0:
            skip += 1
        jk = int(round(npts/2 - (skip+1)/2))
        j0 = jk % skip
    return np.array(range(j0, npts, skip))


def plot_one_component_axes_contour_and_quiver(ax, m_X, m_Y, l_fields, plps,
                                               cc_cont=None, cc_quiver=None):

    do_cont = not cc_cont is None
    do_quiv = not cc_quiver is None

    decorator = plps['decorator']
    cmap_signed = 'seismic'   # This should be a plot_param
    cmap_unsigned = 'OrRd'

    # Adjustments to the visible plot domain
    xlmi = plps.get('xlim_min', 0)
    xlma = plps.get('xlim_max', 0)
    ylmi = plps.get('ylim_min', 0)
    ylma = plps.get('ylim_max', 0)

    if do_cont:
        comp_label = cc_cont.get_label()
        if do_quiv:
            comp_label += ',  '+cc_quiver.get_label()
    else:
        comp_label = cc_quiver.get_label()

    if do_cont:
        cont_signed = cc_cont.is_signed_field()
        cmap = cmap_unsigned
        if cont_signed:
            cmap = cmap_signed

        # Extract and tidy the field
        # transpose because arrays
        cont_field = (l_fields[cc_cont._f_code]).copy().T

        if cc_cont.is_abs():
            # TODO: cleanup: plot |u| as |u|^2
            cont_field = np.abs(cont_field)**2

        # if the data is all noise, just plot zeros
        plot_threshold = 1e-8
        if np.max(np.abs(cont_field[~np.isnan(cont_field)])) < plot_threshold:
            cont_field = np.zeros(np.shape(cont_field))

        # set imshow plot (and tick) range to match the input x and y domain
        if True or plps['ticks']:
            x0 = m_X[0, 0]
            x1 = m_X[0, -1]
            y0 = m_Y[0, 0]
            y1 = m_Y[-1, 0]

            xm = x0+x1
            ym = y0+y1
            # Convert to length units in microns
            extents = np.array([x0-xm/2, x1-xm/2, y0-ym/2, y1-ym/2])*1e6
        else:
            extents = None

        interp = None
        # interp='bilinear';

        # User requested specific limits for each component x, y or z
        if decorator.get_cmap_limits(cc_cont._xyz) != None:
            (act_zlo, act_zhi) = decorator.get_cmap_limits(cc_cont._xyz)
            tsnorm = mplcolors.TwoSlopeNorm(
                vmin=act_zlo, vmax=act_zhi, vcenter=(act_zlo+act_zhi)/2)

            im = ax.imshow(cont_field, origin='lower', extent=extents,
                           interpolation=interp, cmap=cmap, norm=tsnorm)
        else:
            act_zlo = np.nanmin(cont_field)
            act_zhi = np.nanmax(cont_field)
            vma = max(abs(act_zlo), abs(act_zhi))

            if cont_signed:
                vmi = -vma
            else:
                vmi = 0
            im = ax.imshow(cont_field, origin='lower', extent=extents,
                           interpolation=interp, cmap=cmap, vmin=vmi, vmax=vma)

    if do_quiv:
        quiver_points = plps.get('quiver_points', 20)

        # shapes:
        # m_X and other matrices are stored with y points as rows but upside down
        # mesh_grid_objects:
        # m_X =  [[x0, x1, x2, ...]
        #         [x0, x1, x2, ...]
        #         [x0, x1, x2, ...]
        # m_Y =  [[y0, y0, y0, ...]
        #         [y1, y1, y1, ...]
        #         [y2, y2, y2, ...]

        n_pts_y, n_pts_x = m_X.shape
        # delprint('\n\nshape 1', n_pts_x, n_pts_y)
        # delprint('shape 2', m_X[0,:5])
        # delprint('shape 2', m_X[1,:5])
        # delprint('shape 2', m_Y[0,:5])
        # delprint('shape 2', m_Y[1,:5])

        # grid points to skip for each arrow
        # this could probably be chosen nicer
        # quiver_skip_x=int(round(min(n_pts_x,n_pts_y)/quiver_points *(1-xlmi-xlma) ))
        # quiver_skip_y=int(round(min(n_pts_x,n_pts_y)/quiver_points *(1-ylmi-ylma) ))
        quiver_skip_x = int(round(n_pts_x/quiver_points * (1-xlmi-xlma)))
        quiver_skip_y = int(round(n_pts_y/quiver_points * (1-ylmi-ylma)))

        # delprint('shape 3', quiver_points, quiver_skip_x, quiver_skip_y)
        # getting a nice symmetric pattern of points to do quivers centred around the middle
        qrange_x = get_quiver_skip_range(n_pts_x, quiver_skip_x)
        qrange_y = get_quiver_skip_range(n_pts_y, quiver_skip_y)

        # TODO: ensure the points are symmetric about the centre of the domain
        # is quiver_skip_x and _y around the right way given the .T operation?
        v_x_q = m_X.T[qrange_x[:, np.newaxis], qrange_y]
        v_y_q = m_Y.T[qrange_x[:, np.newaxis], qrange_y]
        # v_x_q = m_X.T[0::quiver_skip_x, 0::quiver_skip_y]
        # v_y_q = m_Y.T[0::quiver_skip_x, 0::quiver_skip_y]

        # TODO: why no transpose on these fields?
        # m_ReEx_q = l_fields['Fxr'][0::quiver_skip_x, 0::quiver_skip_y]
        # m_ReEy_q = l_fields['Fyr'][0::quiver_skip_x, 0::quiver_skip_y]
        # m_ImEx_q = l_fields['Fxi'][0::quiver_skip_x, 0::quiver_skip_y]
        # m_ImEy_q = l_fields['Fyi'][0::quiver_skip_x, 0::quiver_skip_y]
        m_ReEx_q = l_fields['Fxr'][qrange_x[:, np.newaxis], qrange_y]
        m_ReEy_q = l_fields['Fyr'][qrange_x[:, np.newaxis], qrange_y]
        m_ImEx_q = l_fields['Fxi'][qrange_x[:, np.newaxis], qrange_y]
        m_ImEy_q = l_fields['Fyi'][qrange_x[:, np.newaxis], qrange_y]

        # convert to microns
        v_x_q_um = v_x_q*1e6
        v_y_q_um = v_y_q*1e6

        # centre at zero
        xm = v_x_q_um[-1, 0]+v_x_q_um[0, 0]
        ym = v_y_q_um[0, -1]+v_y_q_um[0, 0]
        v_x_q_um -= xm/2
        v_y_q_um -= ym/2

    # Ignore all imaginary values. If there are significant imag values,
    # then instaneous vector plots don't make much sense anyway

        if do_cont:  # no colours in the quiver
            ax.quiver(v_x_q_um, v_y_q_um, m_ReEx_q, m_ReEy_q,
                      linewidths=(0.2,), color='gray', edgecolors=('gray'), pivot='mid', headlength=5)
        else:
            m_arrcolour = np.sqrt(m_ReEx_q*m_ReEx_q + m_ReEy_q*m_ReEy_q)
            ax.quiver(v_x_q_um, v_y_q_um, m_ReEx_q, m_ReEy_q, m_arrcolour,
                      linewidths=(0.2,), edgecolors=('gray'), pivot='mid', headlength=5)

        if not do_cont:
            ax.set_aspect('equal')
            # this step is needed because quiver doesn't seem
            # to use its input x and y vectors to set range limits
            # clean this up so as to avoid seemingly circular calls following
            ax.set_xlim(v_x_q_um[0, 0], v_x_q_um[-1, 0])
            ax.set_ylim(v_y_q_um[0, 0], v_y_q_um[0, -1])

    if xlmi > 0 or xlma > 0:
        xmin, xmax = ax.get_xlim()
        width_x = xmax-xmin
        ax.set_xlim(xmin+xlmi*width_x, xmax-xlma*width_x)

    if ylmi > 0 or ylma > 0:
        ymin, ymax = ax.get_ylim()
        width_y = ymax-ymin
        ax.set_ylim(ymin+ylmi*width_y, ymax-ylma*width_y)

    plot_set_ticks(ax, plps, decorator)

    plot_set_title(ax, comp_label, plps, decorator)

    plot_set_axes_style(ax, plps, decorator)

    decorator.add_frame(ax)
    # plot_set_contours(ax, plops, decorator)

    if do_cont:
        # contours and colorbars
        # colorbar
        do_cbar = plps['colorbar']
        do_contours = plps['contours']

        if do_cbar or do_contours:
            if plps['contour_lst']:
                cbarticks = plps['contour_lst']
            else:
                nt = plps.get('num_ticks', 7)
                if nt is None:
                    cbarticks = None
                else:
                    cbarticks = np.linspace(act_zlo, act_zhi, nt)

        if do_cbar:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size=decorator.get_axes_property('cbar_size'),
                                      pad=decorator.get_axes_property('cbar_pad'))
            cbar = plt.colorbar(im, cax=cax)
            if not cbarticks is None:
                cbar.set_ticks(cbarticks)
                cbarlabels = ['%.2f' % t for t in cbarticks]
                cbar.set_ticklabels(cbarlabels)
            cbar.ax.tick_params(labelsize=decorator.get_font_size('cbar_tick'))

        if do_contours:
            if np.max(np.abs(cont_field[~np.isnan(cont_field)])) > plot_threshold:
                CS2 = ax.contour(m_X, m_Y, cont_field, levels=cbarticks,
                                 colors=mycolors[::-1], linewidths=(1.5,))
                if do_cbar:
                    cbar.add_lines(CS2)

    if decorator is not None:
        decorator.extra_axes_commands(ax)


# contour plots a single component specified by cc
def delme_plot_one_component_axes(ax, m_X, m_Y, l_fields, plps, cc):
    decorator = plps['decorator']
    cmap_signed = 'seismic'
    cmap_unsigned = 'OrRd'

    # Adjustments to the visible plot domain
    xlmi = plps.get('xlim_min', 0)
    xlma = plps.get('xlim_max', 0)
    ylmi = plps.get('ylim_min', 0)
    ylma = plps.get('ylim_max', 0)

    comp_label = cc.get_label()  # for plotting
    signed = cc.is_signed_field()

    cmap = cmap_unsigned
    if signed:
        cmap = cmap_signed

    # Extract and tidy the field
    field = (l_fields[cc._f_code]).copy().T  # transpose because arrays

    if cc.is_abs():
        field = np.abs(field)**2  # TODO: cleanup: plot |u| as |u|^2

    # if the data is all noise, just plot zeros
    plot_threshold = 1e-8
    if np.max(np.abs(field[~np.isnan(field)])) < plot_threshold:
        field = np.zeros(np.shape(field))

    # set imshow plot (and tick) range to match the input x and y domain
    if True or plps['ticks']:
        x0 = m_X[0, 0]
        x1 = m_X[0, -1]
        y0 = m_Y[0, 0]
        y1 = m_Y[-1, 0]

        xm = x0+x1
        ym = y0+y1
        # Convert to length units in microns
        extents = np.array([x0-xm/2, x1-xm/2, y0-ym/2, y1-ym/2])*1e6
    else:
        extents = None

    interp = None
    # interp='bilinear';

    # User requested specific limits for each component x, y or z
    if decorator.get_cmap_limits(cc._xyz) != None:
        (act_zlo, act_zhi) = decorator.get_cmap_limits(cc._xyz)
        tsnorm = mplcolors.TwoSlopeNorm(
            vmin=act_zlo, vmax=act_zhi, vcenter=(act_zlo+act_zhi)/2)

        im = ax.imshow(field, origin='lower', extent=extents,
                       interpolation=interp, cmap=cmap, norm=tsnorm)
    else:
        act_zlo = np.nanmin(field)
        act_zhi = np.nanmax(field)
        vma = max(abs(act_zlo), abs(act_zhi))

        if signed:
            vmi = -vma
        else:
            vmi = 0
        im = ax.imshow(field, origin='lower', extent=extents,
                       interpolation=interp, cmap=cmap, vmin=vmi, vmax=vma)

    if xlmi > 0 or xlma > 0:
        xmin, xmax = ax.get_xlim()
        width_x = xmax-xmin
        ax.set_xlim(xmin+xlmi*width_x, xmax-xlma*width_x)

    if ylmi > 0 or ylma > 0:
        ymin, ymax = ax.get_ylim()
        width_y = ymax-ymin
        ax.set_ylim(ymin+ylmi*width_y, ymax-ylma*width_y)

    plot_set_ticks(ax, plps, decorator)

    plot_set_title(ax, comp_label, plps, decorator)

    plot_set_axes_style(ax, plps, decorator)

    # plot_set_contours(ax, plops, decorator)

    # contours and colorbars
    # colorbar
    do_cbar = plps['colorbar']
    do_contours = plps['contours']

    if do_cbar or do_contours:
        if plps['contour_lst']:
            cbarticks = plps['contour_lst']
        else:
            nt = plps.get('num_ticks', 7)
            cbarticks = np.linspace(act_zlo, act_zhi, nt)

    if do_cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=decorator.get_axes_property('cbar_size'),
                                  pad=decorator.get_axes_property('cbar_pad'))
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_ticks(cbarticks)
        cbarlabels = [f'{t:.2f}' for t in cbarticks]
        cbar.set_ticklabels(cbarlabels)
        cbar.ax.tick_params(labelsize=decorator.get_font_size('cbar_tick'))

    if do_contours:
        if np.max(np.abs(field[~np.isnan(field)])) > plot_threshold:
            CS2 = ax.contour(m_X, m_Y, field.T, levels=cbarticks,
                             colors=mycolors[::-1], linewidths=(1.5,))
            if do_cbar:
                cbar.add_lines(CS2)

    if decorator is not None:
        decorator.extra_axes_commands(ax)


def plot_one_component_quiver(ax, m_X, m_Y, v_fields, plps, cc):
    decorator = plps['decorator']

    quiver_points = plps.get('quiver_points', 20)

    # Adjustments to the visible plot domain
    xlmi = plps.get('xlim_min', 0)
    xlma = plps.get('xlim_max', 0)
    ylmi = plps.get('ylim_min', 0)
    ylma = plps.get('ylim_max', 0)

    n_pts_y, n_pts_x = m_X.shape

    # this could probably be chosen nicer
    quiver_skip_x = int(round(min(n_pts_x, n_pts_y) /
                        quiver_points * (1-xlmi-xlma)))
    # this could probably be chosen nicer
    quiver_skip_y = int(round(min(n_pts_x, n_pts_y) /
                        quiver_points * (1-ylmi-ylma)))

    # TODO: ensure the points are symmetric about the centre of the domain
    # is quiver_skip_x and _y around the right way given the .T operation?

    v_x_q = m_X.T[0::quiver_skip_x, 0::quiver_skip_y]
    v_y_q = m_Y.T[0::quiver_skip_x, 0::quiver_skip_y]

    # TODO: why no transpose on these fields?
    m_ReEx_q = v_fields['Fxr'][0::quiver_skip_x, 0::quiver_skip_y]
    m_ReEy_q = v_fields['Fyr'][0::quiver_skip_x, 0::quiver_skip_y]
    m_ImEx_q = v_fields['Fxi'][0::quiver_skip_x, 0::quiver_skip_y]
    m_ImEy_q = v_fields['Fyi'][0::quiver_skip_x, 0::quiver_skip_y]

    # convert to microns
    v_x_q_um = v_x_q*1e6
    v_y_q_um = v_y_q*1e6

    # centre at zero
    xm = v_x_q_um[-1, 0]+v_x_q_um[0, 0]
    ym = v_y_q_um[0, -1]+v_y_q_um[0, 0]
    v_x_q_um -= xm/2
    v_y_q_um -= ym/2


# Ignore all imaginary values. If there are significant imag values,
# then instaneous vector plots don't make much sense anyway
    m_arrcolour = np.sqrt(m_ReEx_q*m_ReEx_q + m_ReEy_q*m_ReEy_q)
    ax.quiver(v_x_q_um, v_y_q_um, m_ReEx_q, m_ReEy_q, m_arrcolour,
              linewidths=(0.2,), edgecolors=('k'), pivot='mid', headlength=5)  # length of the arrows

    ax.set_aspect('equal')
    # this step is needed because quiver doesn't seem
    ax.set_xlim(v_x_q_um[0, 0], v_x_q_um[-1, 0])
    # to use its input x and y vectors to set range limits
    ax.set_ylim(v_y_q_um[0, 0], v_y_q_um[0, -1])
    # clean this up so as to avoid seemingly circular calls following

    plot_set_ticks(ax, plps, decorator)

    comp_label = cc.get_label()
    plot_set_title(ax, comp_label, plps, decorator)

    plot_set_axes_style(ax, plps, decorator)

    # if plps['xlim_min'] != None:
    #    ax.set_xlim(xmin+plps['xlim_min']*width_x,xmax-plps['xlim_max']*width_x)
    # if plps['ylim_min'] != None:
    #    ax.set_ylim(ymin+plps['ylim_min']*width_y,ymax-plps['ylim_max']*width_y)

    if True or xlmi > 0 or xlma > 0:
        xmin, xmax = ax.get_xlim()
        width_x = xmax-xmin
        ax.set_xlim(xmin+xlmi*width_x, xmax-xlma*width_x)

    if True or ylmi > 0 or ylma > 0:
        ymin, ymax = ax.get_ylim()
        width_y = ymax-ymin
        ax.set_ylim(ymin+ylmi*width_y, ymax-ylma*width_y)

    if decorator is not None:
        decorator.extra_axes_commands(ax)


def plot_mode_data(ax, v_fields, plps, sim_wguide, ival, v_x, v_y):  # mode data summary

    # TODO: Move this analysis to a better place
    [m_ReEx, m_ReEy, m_ReEz, m_ImEx, m_ImEy, m_ImEz, m_AbsE] = [v_fields['Fxr'], v_fields['Fyr'],
                                                                v_fields['Fzr'], v_fields['Fxi'], v_fields['Fyi'], v_fields['Fzi'], v_fields['Fabs']]
    tms = sim_wguide.get_modes()
    modeobj = sim_wguide.get_modes()[ival]
    modeobj._analyse_mode(v_x, v_y, m_ReEx, m_ReEy, m_ReEz,
                          m_ImEx, m_ImEy, m_ImEz, m_AbsE)

    yhi = .99  # try and make these ranges match those of the physical domain?
    dy = .10
    decorator = plps['decorator']
    ax.set_xlim(0, 1)
    ax.set_ylim(0, yhi)
    ax.set_aspect('equal')
    ax.axis('off')
    fs = decorator.get_font_size('data_label')

    #  Beward that this is a very coarse measure of component fraction and is not consistent with energy fraction

    [f_x, f_y, f_t, f_z] = modeobj.field_fracs()
    (r0x, r0y) = modeobj.center_of_mass()*1e6
    (wx, wy, w0) = modeobj.second_moment_widths()*1e6

    r = yhi-dy
    x0 = .05
    ax.text(x0-.05, r, r'Mode properties: m={0}'.format(
        ival), transform=ax.transAxes, fontsize=fs+2)
    r -= dy
    if sim_wguide.is_EM():
        ax.text(x0, r, r'$\omega/(2\pi)$: {0:.5f} THz'.format(
            sim_wguide.omega_EM/(2*np.pi*1.e12)), transform=ax.transAxes, fontsize=fs)
        r -= dy
        ax.text(x0, r, r'$k$: {0:.5f} $\mu\mathrm{{m}}^{{-1}}$'.format(
            sim_wguide.kz_EM(ival)/1.e6), transform=ax.transAxes, fontsize=fs)
        r -= dy
        if sim_wguide.ngroup_EM_available():
            ax.text(x0, r, r'$\bar{{n}}$: {0:.6f}, $n_g$: {1:.6f}'.format(sim_wguide.neff(
                ival), sim_wguide.ngroup_EM(ival)), transform=ax.transAxes, fontsize=fs)
            r -= dy  # TODO:
        else:
            ax.text(x0, r, r'$\bar{{n}}$: {0:.6f}'.format(
                sim_wguide.neff(ival)), transform=ax.transAxes, fontsize=fs)
            r -= dy
    else:
        ax.text(x0, r, r'$q$: {0:.5f} $\mu\mathrm{{m}}^{{-1}}$, $\lambda:$ {1:.5f} $\mu$m'.format(
            sim_wguide.q_AC/1.e6, 2*np.pi/sim_wguide.q_AC*1e6), transform=ax.transAxes, fontsize=fs)
        r -= dy
        ax.text(x0, r, r'$q/2\pi$: {0:.5f} $\mu\mathrm{{m}}^{{-1}}$'.format(
            sim_wguide.q_AC/(2.e6*np.pi)), transform=ax.transAxes, fontsize=fs)
        r -= dy
        ax.text(x0, r, r'$\Omega/(2\pi)$: {0:.5f} GHz'.format(
            np.real(sim_wguide.nu_AC(ival))/1.e9), transform=ax.transAxes, fontsize=fs)
        r -= dy
        if sim_wguide.vgroup_AC_available():
            ax.text(x0, r, r'$v_p$: {0:.2f} m/s, $v_g$: {1:.2f} m/s'.format(sim_wguide.vp_AC(
                ival), sim_wguide.vg_AC(ival)), transform=ax.transAxes, fontsize=fs)
            r -= dy
        else:
            ax.text(x0, r, r'$v_p$: {0:.2f} m/s'.format(
                sim_wguide.vp_AC(ival)), transform=ax.transAxes, fontsize=fs)
            r -= dy
    ax.text(x0, r, r'$f_x:$ {0:.3f}, $f_y$: {1:.3f}'.format(
        f_x, f_y), transform=ax.transAxes, fontsize=fs)
    r -= dy
    ax.text(x0, r, r'$f_t:$ {0:.3f}, $f_z$: {1:.3f}'.format(
        f_t, f_z), transform=ax.transAxes, fontsize=fs)
    r -= dy
    ax.text(x0, r, r'$\mathbf{{r}}_0:$ ({0:.2f}, {1:.3f}) $\mu$m'.format(
        r0x, r0y), transform=ax.transAxes, fontsize=fs)
    r -= dy
    ax.text(x0, r, r'$(w_x, w_y):$ ({0:.2f}, {1:.2f}) $\mu$m'.format(
        wx, wy), transform=ax.transAxes, fontsize=fs)
    r -= dy

    if modeobj.field_type == FieldType.EM_H:
        ax.text(x0, r, r'$H$ field multiplied by $Z_0=376.7\, \Omega$',
                transform=ax.transAxes, fontsize=fs)
        r -= dy

    sc = sim_wguide.symmetry_classification(ival)
    if len(sc):
        ax.text(x0, r, r'Sym: {0}'.format(sc),
                transform=ax.transAxes, fontsize=fs)
        r -= dy

    if sim_wguide.is_EM() and sim_wguide.Q_method != QAcMethod.NotSet:
        if sim_wguide.Q_method == QAcMethod.Intrinsic:
            ax.text(x0, r, 'Losses (intrinsic):',
                    transform=ax.transAxes, fontsize=fs)
            r -= dy
        else:
            ax.text(x0, r, 'Losses (fixed Q):',
                    transform=ax.transAxes, fontsize=fs)
            r -= dy

        ax.text(x0+.1, r, r'$\alpha$: {0:.3e} s$^{{-1}}$, {1:.2f} cm$^{{-1}}$, {2:.2f} dB/cm'.format(
            sim_wguide.alpha_t_AC(ival),
            sim_wguide.alpha_s_AC(ival)/100.,
            sim_wguide.alpha_s_AC(ival)/100./(np.log(10.0)/10.0)
        ), transform=ax.transAxes, fontsize=fs)
        r -= dy

        ax.text(x0+.1, r, r'$Q_m$: {0:.2f}'.format(
            sim_wguide.Qmech_AC(ival)), transform=ax.transAxes, fontsize=fs)
        r -= dy
        ax.text(x0+.1, r, r'$\Delta\Omega/(2\pi)$: {0:.4f} MHz'.format(
            1.e-6*sim_wguide.linewidth_AC(ival)), transform=ax.transAxes, fontsize=fs)
        r -= dy
        mg_pe = plps['modal_gain'].get('PE', 0)
        mg_mb = plps['modal_gain'].get('MB', 0)
        mg_tot = plps['modal_gain'].get('Tot', 0)
        ax.text(x0+.1, r, r'Gain (PE, MB, Tot): ({0:.4f}, {1:.4f}, {2:.4f} ) $\mathrm{{(mW)}}^{{-1}}$'.format(
            mg_pe, mg_mb, mg_tot), transform=ax.transAxes, fontsize=fs)
        r -= dy
    d_ext = modeobj.get_mode_data()
    if len(d_ext):
        ax.text(x0, r, r'Extra data:', transform=ax.transAxes, fontsize=fs)
        r -= dy
        for (k, v) in d_ext.items():
            ax.text(x0, r, r'   {0}: {1}'.format(k, v),
                    transform=ax.transAxes, fontsize=fs)
            r -= dy


def plot_all_components(v_x, v_y, m_X, m_Y, v_plots, plps, sim_wguide, ival):
    decorator = plps['decorator']
    figsz = decorator.get_axes_property('figsize')
    ws = decorator.get_axes_property('subplots_wspace')
    hs = decorator.get_axes_property('subplots_hspace')
    lw = decorator.get_axes_property('linewidth')
    ec = decorator.get_axes_property('edgecolor')

    decorator.set_frame_drawer(sim_wguide.structure.wg_geom)

    EM_AC = plps['EM_AC']
    plt.clf()
    fig = plt.figure(figsize=figsz)
    fig.subplots_adjust(wspace=ws, hspace=hs)
    plt.rc('axes', linewidth=lw)
    plt.rc('axes', edgecolor=ec)

    rows = 3
    if plps['suppress_imimre']:
        rows = 2

    axi = 1
    ax = plt.subplot(rows, 3, axi)
    plot_mode_data(ax, v_plots, plps, sim_wguide,
                   ival, v_x, v_y)  # mode data summary
    axi += 1

    cc_cont = {FieldType.EM_E: component_t('Eabs'), FieldType.EM_H: component_t(
        'Habs'), FieldType.AC: component_t('uabs')}[EM_AC]
    cc_quiv = {FieldType.EM_E: component_t('Et'), FieldType.EM_H: component_t(
        'Ht'), FieldType.AC: component_t('ut')}[EM_AC]

    ax = plt.subplot(rows, 3, axi)
    # plot_one_component_axes(ax, m_X, m_Y, v_plots, plps, cc)  # the intensity plot
    plot_one_component_axes_contour_and_quiver(ax, m_X, m_Y, v_plots, plps,
                                               cc_cont=cc_cont, cc_quiver=cc_quiv)  # the intensity plot
    axi += 1

    ax = plt.subplot(rows, 3, axi)
    # plot_one_component_quiver(ax, m_X, m_Y, v_plots, plps, cc)  # the transverse vector plot
    plot_one_component_axes_contour_and_quiver(ax, m_X, m_Y, v_plots, plps,
                                               cc_quiver=cc_quiv)  # the intensity plot
    axi += 1

    for (flab, field) in v_plots.items():
        cc = component_t.make_comp(EM_AC, flab)

        if plps['suppress_imimre'] and not cc.is_dominant():
            continue
        if not cc.is_signed_field():
            continue
        ax = plt.subplot(rows, 3, axi)
        plot_one_component_axes_contour_and_quiver(
            ax, m_X, m_Y, v_plots, plps, cc_cont=cc)  # the scalar plots
        axi += 1

    figfile = plot_filename(plps, ival)
    save_figure(plt, figfile)

    if not keep_plots_open:
        plt.close()

    if False and plps['EM_AC'] == FieldType.AC and plps['stress_fields']:
        # Interpolate onto rectangular Cartesian grid
        # TODO: connect mgrid with meshgrid picture above
        xy = list(zip(v_x6p, v_y6p))
        grid_x, grid_y = np.mgrid[x_min:x_max:n_pts_x *
                                  1j, y_min:y_max:n_pts_y*1j]
        m_ReEx = interpolate.griddata(
            xy, v_Ex6p.real, (grid_x, grid_y), method='linear')
        m_ReEy = interpolate.griddata(
            xy, v_Ey6p.real, (grid_x, grid_y), method='linear')
        m_ReEz = interpolate.griddata(
            xy, v_Ez6p.real, (grid_x, grid_y), method='linear')
        m_ImEx = interpolate.griddata(
            xy, v_Ex6p.imag, (grid_x, grid_y), method='linear')
        m_ImEy = interpolate.griddata(
            xy, v_Ey6p.imag, (grid_x, grid_y), method='linear')
        m_ImEz = interpolate.griddata(
            xy, v_Ez6p.imag, (grid_x, grid_y), method='linear')
        m_AbsE = interpolate.griddata(
            xy, v_E6p.real, (grid_x, grid_y), method='linear')
        dx = grid_x[-1, 0] - grid_x[-2, 0]
        dy = grid_y[0, -1] - grid_y[0, -2]
        m_Ex = m_ReEx + 1j*m_ImEx
        m_Ey = m_ReEy + 1j*m_ImEy
        m_Ez = m_ReEz + 1j*m_ImEz
        m_Ex = m_Ex.reshape(n_pts_x, n_pts_y)
        m_Ey = m_Ey.reshape(n_pts_x, n_pts_y)
        m_Ez = m_Ez.reshape(n_pts_x, n_pts_y)
        m_AbsE = m_AbsE.reshape(n_pts_x, n_pts_y)

        m_ReEx = np.real(m_Ex)
        m_ReEy = np.real(m_Ey)
        m_ReEz = np.real(m_Ez)
        m_ImEx = np.imag(m_Ex)
        m_ImEy = np.imag(m_Ey)
        m_ImEz = np.imag(m_Ez)

        del_x_Ex = np.gradient(m_Ex, dx, axis=0)
        del_y_Ex = np.gradient(m_Ex, dy, axis=1)
        del_x_Ey = np.gradient(m_Ey, dx, axis=0)
        del_y_Ey = np.gradient(m_Ey, dy, axis=1)
        del_x_Ez = np.gradient(m_Ez, dx, axis=0)
        del_y_Ez = np.gradient(m_Ez, dy, axis=1)
        del_z_Ex = 1j*sim_wguide.q_AC*m_Ex
        del_z_Ey = 1j*sim_wguide.q_AC*m_Ey
        del_z_Ez = 1j*sim_wguide.q_AC*m_Ez

        # Flip y order as imshow has origin at top left
        del_mat = np.array([del_x_Ex[:, ::-1].real, del_x_Ey[:, ::-1].real, del_x_Ez[:, ::-1].real, del_x_Ex[:, ::-1].imag, del_x_Ey[:, ::-1].imag, del_x_Ez[:, ::-1].imag, del_y_Ex[:, ::-1].real, del_y_Ey[:, ::-1].real, del_y_Ez[:, ::-1].real,
                           del_y_Ex[:, ::-1].imag, del_y_Ey[:, ::-1].imag, del_y_Ez[:, ::-1].imag, del_z_Ex[:, ::-1].real, del_z_Ey[:, ::-1].real, del_z_Ez[:, ::-1].real, del_z_Ex[:, ::-1].imag, del_z_Ey[:, ::-1].imag, del_z_Ez[:, ::-1].imag])
        v_labels = ["Re($S_{xx}$)", "Re($S_{xy}$)", "Re($S_{xz}$)", "Im($S_{xx}$)", "Im($S_{xy}$)", "Im($S_{xz}$)", "Re($S_{yx}$)", "Re($S_{yy}$)", "Re($S_{yz}$)",
                    "Im($S_{yx}$)", "Im($S_{yy}$)", "Im($S_{yz}$)", "Re($S_{zx}$)", "Re($S_{zy}$)", "Re($S_{zz}$)", "Im($S_{zx}$)", "Im($S_{zy}$)", "Im($S_{zz}$)"]

        # stress field plots
        plt.clf()
        fig = plt.figure(figsize=(15, 30))
        for i_p, plot in enumerate(del_mat):
            ax = plt.subplot(6, 3, i_p+1)
            im = plt.imshow(plot.T)
            # no ticks
            plt.xticks([])
            plt.yticks([])
            # limits
            if xlim_min != None:
                ax.set_xlim(xlim_min*n_points, (1-xlim_max)*n_points)
            if ylim_min != None:
                ax.set_ylim((1-ylim_min)*n_points, ylim_max*n_points)
            # titles
            plt.title(v_labels[i_p], fontsize=decorator.get_font_size(
                'subplot_title'))
            # colorbar
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.1)
            cbar = plt.colorbar(im, cax=cax, format='%.2e')
            if num_ticks:
                cbarticks = np.linspace(
                    np.min(plot), np.max(plot), num=num_ticks)
            elif ylim_min != None:
                if xlim_min/ylim_min > 3:
                    cbarlabels = np.linspace(np.min(plot), np.max(plot), num=3)
                if xlim_min/ylim_min > 1.5:
                    cbarlabels = np.linspace(np.min(plot), np.max(plot), num=5)
                else:
                    cbarlabels = np.linspace(np.min(plot), np.max(plot), num=7)
            else:
                cbarlabels = np.linspace(np.min(plot), np.max(plot), num=7)
            cbar.set_ticks(cbarlabels)
            cbarlabels = ['%.2f' % t for t in cbarlabels]
            cbar.set_ticklabels(cbarlabels)
            if contours:
                if contour_lst:
                    cbarticks = contour_lst
                if np.max(np.abs(plot[~np.isnan(plot)])) > plot_threshold:
                    CS2 = ax.contour(
                        m_X, m_Y, plot.T, levels=cbarticks, colors=colors[::-1], linewidths=(1.5,))
                cbar.add_lines(CS2)
            cbar.ax.tick_params(labelsize=decorator.get_font_size('cbar_tick'))
        fig.set_tight_layout(True)
        n_str = ''
        if np.imag(sim_wguide.Eig_values[ival]) < 0:
            k_str = r'$\Omega/2\pi = %(re_k)f %(im_k)f i$ GHz' % \
                {'re_k': np.real(sim_wguide.Eig_values[ival]*1e-9),
                 'im_k': np.imag(sim_wguide.Eig_values[ival]*1e-9)}
        else:
            k_str = r'$\Omega/2\pi = %(re_k)f + %(im_k)f i$ GHz' % \
                {'re_k': np.real(sim_wguide.Eig_values[ival]*1e-9),
                 'im_k': np.imag(sim_wguide.Eig_values[ival]*1e-9)}
        plt.suptitle('Mode #' + str(ival) + '   ' + k_str + '   ' +
                     n_str, fontsize=decorator.get_font_size('title'))

        if pdf_png == 'png':
            plt.savefig('%(pre)sfields/%(s)s_S_field_%(i)i%(add)s.png' %
                        {'pre': prefix, 's': EM_AC, 'i': ival, 'add': suffix})
        elif pdf_png == 'pdf':
            plt.savefig('%(pre)sfields/%(s)s_S_field_%(i)i%(add)s.pdf' %
                        {'pre': prefix, 's': EM_AC, 'i': ival, 'add': suffix}, bbox_inches='tight')
        if not keep_plots_open:
            plt.close()


def save_figure(plt, figfile):
    if figfile[-3:-1] == 'png':
        plt.savefig(figfile)
    else:
        plt.savefig(figfile, bbox_inches='tight')


def plot_one_component(m_X, m_Y, v_fields, plps, sim_wguide, ival, cc, axis=None):

    decorator = plps['decorator']

    if axis is None:
        plt.clf()  # TODO: surely this is not needed

        plt.rc('axes', linewidth=decorator.get_axes_property('linewidth'))
        plt.rc('axes', edgecolor=decorator.get_axes_property('edgecolor'))

        fig = plt.figure(figsize=(12, 10))
        ax = plt.subplot(111)
    else:
        ax = axis

    EM_AC = plps['EM_AC']
    cc_cont = {FieldType.EM_E: component_t('Eabs'), FieldType.EM_H: component_t(
        'Habs'), FieldType.AC: component_t('uabs')}[EM_AC]
    cc_quiv = {FieldType.EM_E: component_t('Et'), FieldType.EM_H: component_t(
        'Ht'), FieldType.AC: component_t('ut')}[EM_AC]

    # if cc.is_transverse():
    #  plot_one_component_quiver(ax, m_X, m_Y, v_fields, plps, cc)
    # else:
    #  plot_one_component_axes(ax, m_X, m_Y, v_fields, plps, cc)

    plot_one_component_axes_contour_and_quiver(ax, m_X, m_Y, v_fields, plps,
                                               cc_cont=cc_cont, cc_quiver=cc_quiv)

    if axis is None:  # If user passed in the axis, they can look after saving.
        figfile = plot_filename(plps, ival, cc._user_code)
        save_figure(plt, figfile)
        if not keep_plots_open:
            plt.close()


# deprecated spelling.
def plt_mode_fields(sim_wguide, ivals=None, n_points=501, quiver_points=50,
                    xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None,
                    EM_AC='EM_E', num_ticks=None, colorbar=True, contours=False, contour_lst=None,
                    stress_fields=False, pdf_png='png',
                    prefix='', suffix='', ticks=True, comps=[], decorator=None,
                    suppress_imimre=True,
                    modal_gains_PE=None,
                    modal_gains_MB=None,
                    modal_gains=None):

    print('Warning: "plt_mode_fields" is deprecated, use "plot_mode_fields"')
    plot_mode_fields(sim_wguide, ivals, n_points, quiver_points,
                     xlim_min, xlim_max, ylim_min, ylim_max,
                     EM_AC, num_ticks, colorbar, contours, contour_lst, stress_fields, pdf_png,
                     prefix, suffix, ticks, comps, decorator,
                     suppress_imimre, modal_gains_PE, modal_gains_MB, modal_gains)

#### Standard plotting of spectra #############################################


def plot_mode_fields(sim_wguide, ivals=None, n_points=501, quiver_points=30,
                     xlim_min=0, xlim_max=0, ylim_min=0, ylim_max=0,
                     EM_AC='EM_E',
                     num_ticks=None, colorbar=True, contours=False, contour_lst=None,
                     stress_fields=False, pdf_png='png',
                     prefix='', suffix='', ticks=True, comps=[], decorator=None,
                     suppress_imimre=True,
                     modal_gains_PE=None,
                     modal_gains_MB=None,
                     modal_gains=None):
    """ Plot E or H fields of EM mode, or the AC modes displacement fields.

        Args:
            sim_wguide : A ``Struct`` instance that has had calc_modes calculated

        Keyword Args:
            ivals  (list): mode numbers of modes you wish to plot

            n_points  (int): The number of points across unitcell to
                interpolate the field onto

            xlim_min  (float): Limit plotted xrange to xlim_min:(1-xlim_max) of unitcell

            xlim_max  (float): Limit plotted xrange to xlim_min:(1-xlim_max) of unitcell

            ylim_min  (float): Limit plotted yrange to ylim_min:(1-ylim_max) of unitcell

            ylim_max  (float): Limit plotted yrange to ylim_min:(1-ylim_max) of unitcell

            EM_AC  (str): Either 'EM' or 'AC' modes

            num_ticks  (int): Number of tick marks

            contours  (bool): Controls contours being overlaid on fields

            contour_lst  (list): Specify contour values

            stress_fields  (bool): Calculate acoustic stress fields

            pdf_png  (str): File type to save, either 'png' or 'pdf'

            prefix  (str): Add a string to start of file name

            suffix  (str): Add a string to end of file name.

            modal_gains (float array): Pre-calculated gain for each acoustic mode given chosen optical fields.
    """

    # TODO: ft_emac  ,a new var to distinguish EM_AC string and enum
    if sim_wguide.is_AC():
        EM_AC = FieldType.AC
    else:
        try:
            # TODO:ugly that this changes from string to enum
            EM_AC = FieldType.from_str(EM_AC)
        except:
            raise ValueError("EM_AC must be either 'AC', 'EM_E' or 'EM_H'.")

    # Do this elsewhere by function on Mode

    # Calculate the magnetic field from the electric field
    # Move this into plot_mode_H
    if EM_AC == FieldType.EM_H:
        nnodes = 6
        sim_wguide.sol1_H = NumBAT.h_mode_field_ez(sim_wguide.k_0, sim_wguide.n_modes,
                                                   sim_wguide.n_msh_el, sim_wguide.n_msh_pts, nnodes, sim_wguide.table_nod,
                                                   sim_wguide.mesh_xy, sim_wguide.Eig_values, sim_wguide.sol1)

    # assemble desired list of eigenmodes to plot
    if not ivals is None:
        ival_range = ivals
    else:
        ival_range = range(len(sim_wguide.Eig_values))

    mode_helper = sim_wguide.get_mode_helper()
    mode_helper.setup_plot_grid(n_points=n_points)

    if decorator is None:
        decorator = Decorator()

    # TODO: what is this for?
    modal_gain = {}
    if modal_gains is not None:
        modal_gain['Tot'] = modal_gains[ival]
    if modal_gains_PE is not None:
        modal_gain['PE'] = modal_gains_PE[ival]
    if modal_gains_MB is not None:
        modal_gain['MB'] = modal_gains_MB[ival]


#    plot_params={'xlim_min': xlim_min, 'xlim_max': xlim_max, 'ylim_min': ylim_min,
#                 'ylim_max': ylim_max, 'ticks': ticks, 'num_ticks':num_ticks,
#                  'colorbar':colorbar, 'contours':contours, 'contour_lst':contour_lst, 'EM_AC':EM_AC,
#                  'prefix': prefix, 'suffix': suffix, 'pdf_png': pdf_png,
#                  'stress_fields':stress_fields, 'modal_gain':modal_gain, 'decorator': decorator,
#                  'suppress_imimre':suppress_imimre,
#              # 'n_pts_x': n_pts_x, 'n_pts_y': n_pts_y,
#               'quiver_points': quiver_points }
#

    if not prefix: prefix=numbat.NumBATApp().outprefix()

    mode_helper.set_plot_params(xlim_min=xlim_min, xlim_max=xlim_max, ylim_min=ylim_min, ylim_max=ylim_max,
                                EM_AC='EM_E',
                                quiver_points=quiver_points,
                                num_ticks=num_ticks, colorbar=colorbar, contours=contours, contour_lst=contour_lst,
                                pdf_png=pdf_png,
                                prefix=prefix, suffix=suffix, ticks=ticks,  # comps=comps,
                                decorator=decorator,
                                suppress_imimre=suppress_imimre
                                )
    # modal_gains_PE=modal_gains_PE, #purpose, can these go elsewhere?
    # modal_gains_MB=modal_gains_MB,
    # modal_gains=modal_gains)

    v_modes = sim_wguide.get_modes()

    modetype = 'em'
    if EM_AC == FieldType.AC:
        modetype = 'acoustic'

    if len(ival_range) > 1:
        print('Plotting {} modes m={} to {}.'.format(
            modetype, ival_range[0], ival_range[-1]))
    else:
        print('Plotting {} mode m={}.'.format(modetype, ival_range[0]))

    # TODO: mode drawing and especially saving is very slow. could do this in threads?
    for ival in ival_range:
        v_modes[ival].plot_mode(comps, EM_AC)


#### Plot mesh #############################################
def plot_msh(mesh_xy, prefix='', suffix=''):
    """ Plot EM mode fields.

        Args:
            sim_wguide : A ``Struct`` instance that has had calc_modes calculated

        Keyword Args:
            n_points  (int): The number of points across unitcell to \
                interpolate the field onto.
    """

    # plt.clf()
    # plt.figure(figsize=(13,13))
    # ax = plt.subplot(1,1,1)
    fig, ax = plt.subplots()
    for node in range(np.shape(mesh_xy)[1]):
        plt.plot(mesh_xy[0, node], mesh_xy[1, node], 'og')
    ax.set_aspect('equal')
    # plt.savefig('%(pre)smsh_%(add)s.pdf' %
    #   {'pre' : prefix, 'add' : suffix}, bbox_inches='tight')
    fig.savefig('%(pre)smsh_%(add)s.pdf' %
                {'pre': prefix, 'add': suffix}, bbox_inches='tight')
    if not keep_plots_open:
        plt.close()


# Plot nodal arrangement on mesh triangle.
# plt.figure(figsize=(13,13))
# el = 1
# plt.clf()
# for i in range(0,6):
#     print table_nod[i][el] - 1
#     x = mesh_xy[0,table_nod[i][el] - 1]
#     y = mesh_xy[1,table_nod[i][el] - 1]
#     print 'x1 = ', mesh_xy[0,table_nod[i][el] - 1]
#     print 'y1 = ', mesh_xy[1,table_nod[i][el] - 1]
#     plt.plot(x, y, 'o')
#     plt.text(x+0.001, y+0.001, str(i))
# plt.savefig('triangle_%i.png' %el)
if not keep_plots_open:
    plt.close()
