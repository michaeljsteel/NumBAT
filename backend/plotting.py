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


import os
import sys
import math
import numpy as np
from scipy import sqrt
import subprocess
from scipy import interpolate
import matplotlib
#matplotlib.use('pdf')  # TODO: remove if ok
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mplcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker
import copy

from nbtypes import *
from fortran import NumBAT

keep_plots_open=False  # setting this true is useful for use in jupyter style notebooks. TODO: Make a nicer interface

try: 
    plt.style.use('NumBATstyle')
except (ValueError, IOError, AttributeError): print("Preferred matplotlib style file not found.")
mycolors = [color['color'] for color in list(plt.rcParams['axes.prop_cycle'])]


# font = {'family' : 'normal',
#         'weight' : 'bold',
#         'size'   : 18}
# matplotlib.rc('font', **font)
linesstrength = 2.5

class Decorator(object):

  def __init__(self):
    mp_base_font=18

    self._multiplot_fontsizes= {'title':mp_base_font-2, 'subplot_title':mp_base_font-5, 
      'ax_label':mp_base_font-10, 'data_label':mp_base_font-8,
      'ax_tick':mp_base_font-10, 'cbar_tick':mp_base_font-12}
    self._multiplot_axesprops={'figsize': (10,8), 
       'subplots_hspace': .2, 'subplots_wspace': .2, 
       'linewidth':'.75', 'edgecolor':'gray', 'title_pad':5 , 'cbar_size':'5%', 'cbar_pad': '4%'}

    sp_base_font=24
    self._singleplot_fontsizes= {'title':sp_base_font-2, 'ax_label':20, 'subplot_title':25,'cbar_tick':20, 'ax_tick':20 }
    self._singleplot_axesprops={'linewidth':'.75', 'edgecolor':'gray', 'title_pad':20, 'cbar_size':'5%', 'cbar_pad': '2%'}

    self._is_single=True
    self._cmap_limits={}

  def _fontsizes(self):
    if self._is_single: return self._singleplot_fontsizes
    else: return self._multiplot_fontsizes

  def _get_axes_prop(self):
    if self._is_single: return self._singleplot_axesprops
    else: return self._multiplot_axesprops

  def _set_for_single(self): 
    self._is_single=True

  def _set_for_multi(self):  
    self._is_single=False

  def set_cmap_limits(self,d):
    '''Specify the lower and upper contour plot limits for a field component plot.

      :param dict d: Dictionary mapping component ('x','y','z') to (zlo, zhi) tuple for contour plots.
      '''
    self._cmap_limits.update(d)

  def get_cmap_limits(self, comp):
    return self._cmap_limits.get(comp)

  def get_font_size(self, lab):
    fs=10
    try:
      fs=self._fontsizes()[lab]
    except:
      print('Warning: unknown fontsize label "{0}" in Decorator::get_font_size()'.format(lab))
    return fs

  def get_axes_property(self, lab):
    prop=''
    try:
      prop=self._get_axes_prop()[lab]
    except:
      print('Warning: unknown axes property label "{0}" in Decorator::get_axes_property()'.format(lab))
      print (self._get_axes_prop())
      sys.exit(1)
    return prop

  def is_single_plot(self): 
    '''Returns True if this Decorator is for a single axes plot such as a spectrum or spatial map of a single field component.
       '''
    return self._is_single

  def set_singleplot_axes_property(self, label, prop): 
    '''Add or override an axes property for a single plot corresponding to the given label.'''
    self._singleplot_axesprops[label]=prop

  def set_multiplot_axes_property(self, label, prop): 
    '''Add or override an axes property for a multiple axes plot corresponding to the given label.'''
    self._multiplot_axesprops[label]=prop

  def set_singleplot_fontsize(self, label, sz): 
    '''Override a font size for a single plot corresponding to the given label.'''
    self._singleplot_fontsizes[label]=sz

  def set_multiplot_fontsize(self, label, sz): 
    '''Override a font size for a mutiple axes plot corresponding to the given label.'''
    self._multiplot_fontsizes[label]=sz

  def extra_axes_commands(self, ax):
    '''Function to make additions to a standard plot.
      
      This function is called after all other plotting operations are performed.
      The default implementation does nothing. By subclassing :Decorator: and implementing this function,
      users may add extra features to a plot.
      '''
    pass



#### Natural constants ########################################################
ASTM15_tot_I   = 900.084            # Integral ASTM 1.5 solar irradiance W/m**2
Plancks_h      = 6.62607015e-34     # Planck's constant in Js (exact)
speed_c        = 299792458          # Speed of light in vacuum in m/s (exact)
charge_e       = 1.602176634e-19    # Charge of an electron in C (exact)
###############################################################################

class component_t(object):
  @staticmethod
  def Ecomp(c): return component_t('E'+c)

  @staticmethod
  def Hcomp(c): return component_t('H'+c)

  @staticmethod
  def ucomp(c): return component_t('u'+c)

  @staticmethod
  def make_comp(emac, fc): 
    if emac=='EM_E': 
      uc = { 'Fxr':'Ex', 'Fxi':'Ex', 'Fyr':'Ey', 'Fyi':'Ey', 'Fzr':'Ez', 'Fzi':'Ez', 'Fabs':'Eabs', 'Ft':'Et'}[fc]
    elif emac=='EM_H': 
      uc = { 'Fxr':'Hx', 'Fxi':'Hx', 'Fyr':'Hy', 'Fyi':'Hy', 'Fzr':'Hz', 'Fzi':'Hz', 'Fabs':'Habs', 'Ft':'Ht'}[fc]
    else:
      uc = { 'Fxr':'ux', 'Fxi':'ux', 'Fyr':'uy', 'Fyi':'uy', 'Fzr':'uz', 'Fzi':'uz', 'Fabs':'uabs', 'Ft':'ut'}[fc]
    cc= component_t(uc)
    cc._f_code=fc
    return cc

  def __init__(self, uc): #must be an Ex, Ey style, not an EM_AC,Fxr style
    self._user_code=uc
    self._E=uc[0] # E, H, or u
    self._Ei=uc[:2] # Ex, Ey, Ez, Ea, Et, Hx, Hy, etc
    self._xyz=uc[1] # x, y, z, a, t
    self._Eimaj=self.reim_major(self._Ei)
    self._Eimin=self.reim_minor(self._Ei)

    #default real/imag given the _E value                          
    self._f_code= {'Ex':'Fxr', 'Ey':'Fyr', 'Ez':'Fzi', 'Eabs':'Fabs', 'Et':'Ft',
                'Hx':'Fxr', 'Hy':'Fyr', 'Hz':'Fzi', 'Habs':'Fabs', 'Ht':'Ft',
                'ux':'Fxr', 'uy':'Fyr', 'uz':'Fzi', 'uabs':'Fabs','ut':'Ft',
                   }[self._user_code]


  def get_label(self):
    c=self._E
    lab= { 'Fx':r'Re($E_x$)', 'Fy':r'Re($E_y$)', 'Fz':r'Im($E_z$)', 'Fxr':r'Re($E_x$)',
              'Fyr':r'Re($E_y$)', 'Fzi':r'Im($E_z$)', 'Fxi':r'Im($E_x$)', 'Fyi':r'Im($E_y$)', 'Fzr':r'Re($E_z$)',
              'Fabs':r'$|E|$', 'Ft':r'$(E_x, E_y)$'}[self._f_code]
    return lab.replace('E', c)

  def is_signed_field(self): return self._f_code not in ('Ft', 'Fabs')
  def is_transverse(self): return self._user_code.endswith('t')
  def is_dominant(self): return self._f_code in ('Fxr', 'Fyr', 'Fzi')


  def reim_major(self, fi): 
    try:
      return { 'Ex':'Exr', 'Ey':'Eyr', 'Ez':'Ezi', 'Ea':'Ea',
             'Hx':'Hxr', 'Hy':'Hyr', 'Hz':'Hzi', 'Ha':'Ha',
             'ux':'uxr', 'uy':'uyr', 'uz':'uzi', 'ua':'ua' }[fi]
    except KeyError:
      return fi

  def reim_minor(self, fi): 
    try:
      return { 'Ex':'Exi', 'Ey':'Eyi', 'Ez':'Ezr', 'Ea':None,
             'Hx':'Hxi', 'Hy':'Hyi', 'Hz':'Hyr', 'Ha':None,
             'ux':'uxi', 'uy':'uyi', 'uz':'uyr', 'ua':None }[fi]
    except KeyError:
      return None



#### Short utility functions ##################################################
def zeros_int_str(zero_int):
    """ Convert integer into string with '0' in place of ' '. """
    # if zero_int == 0:
    #     fmt_string = '0000'
    # else:
    #     string = '%4.0f' % zero_int
    #     fmt_string = string.replace(' ','0')
    string = '%4.0f' % zero_int
    fmt_string = string.replace(' ','0')
    return fmt_string

###############################################################################


def gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
                EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min, freq_max, num_interp_pts=3000,
                save_fig=True, dB=False, dB_peak_amp=10, mode_comps=False, semilogy=False,
                pdf_png='png', save_txt=False, prefix_str='', suffix_str='', decorator=None,
                show_gains='All'):
    r""" Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
            

        Args:
            sim_AC : An AC ``Struct`` instance that has had calc_modes calculated

            SBS_gain  (array): Totlat SBS gain of modes.

            SBS_gain_PE  (array): Moving Bountary gain of modes.

            SBS_gain_MB  (array): Photoelastic gain of modes.

            linewidth_Hz  (array): Linewidth of each mode [Hz].

            k_AC  (float): Acoustic wavevector.

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

            save_fig  (bool): Save figure at all.

            pdf_png  (str): Save figures as 'png' or 'pdf'.

            save_txt  (bool): Save spectra data to txt file.

            prefix_str  (str): String to be appended to start of file name.

            suffix_str  (str): String to be appended to end of file name.
    """

    if decorator is None: decorator = Decorator()
    decorator._set_for_single()

    tune_steps = 50000
    tune_range = 10 # GHz
    # Construct an odd range of freqs guaranteed to include central resonance frequency.
    detuning_range = np.append(np.linspace(-1*tune_range, 0, tune_steps),
                       np.linspace(0, tune_range, tune_steps)[1:])*1e9 # GHz
    interp_grid = np.linspace(freq_min, freq_max, num_interp_pts)

    # Linewidth of Lorentzian is half the FWHM style linewidth.
    linewidth = linewidth_Hz/2

    # Plot decomposition of spectra into individual modes.
    interp_values = np.zeros(num_interp_pts)
    if mode_comps and save_fig:
        plt.figure()
        plt.clf()  #Why is this necessary?
    if AC_ival == 'All':
        for AC_i in range(len(linewidth)):
            gain_list = np.real(SBS_gain[EM_ival_pump,EM_ival_Stokes,AC_i]
                         *linewidth[AC_i]**2/(linewidth[AC_i]**2 + detuning_range**2))
            freq_list_GHz = np.real(sim_AC.Eig_values[AC_i] + detuning_range)*1e-9
            if mode_comps:
                plt.plot(freq_list_GHz, np.abs(gain_list), linewidth=2)
                if save_txt:
                    save_array = (freq_list_GHz, gain_list)
                    np.savetxt('gain_spectra-mode_comps%(add)s-%(mode)i.csv' 
                                % {'add' : suffix_str, 'mode' : AC_i}, 
                                save_array, delimiter=',')
            # set up an interpolation for summing all the gain peaks
            interp_spectrum = np.interp(interp_grid, freq_list_GHz, gain_list)
            interp_values += interp_spectrum
    else: 
      raise NotImplementedError("Spectrum plotting for limited AC modes not implemented.")

    return_interp_values = interp_values
    if mode_comps:
        plt.plot(interp_grid, np.abs(interp_values), 'b', linewidth=3, label="Total")
        plt.legend(loc=0)
        if freq_min and freq_max:
            plt.xlim(freq_min,freq_max)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('|Gain| (1/Wm)')
        if save_txt:
            save_array = (interp_grid, interp_values)
            np.savetxt('gain_spectra-mode_comps%(add)s-Total.csv' 
                        % {'add' : suffix_str}, save_array, delimiter=',')

    if mode_comps and save_fig:
        if pdf_png=='png':
            plt.savefig('%(pre)sgain_spectra-mode_comps%(add)s.png' % {'pre' : prefix_str, 'add' : suffix_str})
        elif pdf_png=='pdf':
            plt.savefig('%(pre)sgain_spectra-mode_comps%(add)s.pdf' % {'pre' : prefix_str, 'add' : suffix_str})
        if not keep_plots_open: plt.close()


    interp_values = np.zeros(num_interp_pts)
    interp_values_PE = np.zeros(num_interp_pts)
    interp_values_MB = np.zeros(num_interp_pts)
    plt.figure()
    plt.clf()
    ax=plt.gca()
    if AC_ival == 'All':
        for AC_i in range(len(linewidth)):
            gain_list = np.real(SBS_gain[EM_ival_pump,EM_ival_Stokes,AC_i]
                         *linewidth[AC_i]**2/(linewidth[AC_i]**2 + detuning_range**2))
            freq_list_GHz = np.real(sim_AC.Eig_values[AC_i] + detuning_range)*1e-9
            interp_spectrum = np.interp(interp_grid, freq_list_GHz, gain_list)
            interp_values += interp_spectrum

            gain_list_PE = np.real(SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,AC_i]
                         *linewidth[AC_i]**2/(linewidth[AC_i]**2 + detuning_range**2))
            interp_spectrum_PE = np.interp(interp_grid, freq_list_GHz, gain_list_PE)
            interp_values_PE += interp_spectrum_PE

            gain_list_MB = np.real(SBS_gain_MB[EM_ival_pump,EM_ival_Stokes,AC_i]
                         *linewidth[AC_i]**2/(linewidth[AC_i]**2 + detuning_range**2))
            interp_spectrum_MB = np.interp(interp_grid, freq_list_GHz, gain_list_MB)
            interp_values_MB += interp_spectrum_MB
    else: raise NotImplementedError("Spectrum plotting for limited AC modes not implemented.")
    if save_fig:
        if not show_gains in ('All', 'PE', 'MB', 'Total'): show_gains='All'

        lw=decorator.get_axes_property('linewidth')
        fs=decorator.get_font_size('ax_label')
        ts=decorator.get_font_size('ax_tick')

        if show_gains in ('PE', 'All'):
          ax.plot(interp_grid, np.abs(interp_values_PE), 'r', linewidth=lw, label='PE')
        if show_gains in ('MB', 'All'):
          ax.plot(interp_grid, np.abs(interp_values_MB), 'g', linewidth=lw, label='MB')
        if show_gains in ('Total', 'All'):
#ax.plot(interp_grid, np.abs(interp_values), 'b', linewidth=3, label='Total')
          ax.plot(interp_grid, np.abs(interp_values), 'b', linewidth=lw)
        ax.legend(loc=0)
        if freq_min and freq_max:
            ax.set_xlim(freq_min,freq_max)
        ax.set_xlabel('Frequency (GHz)',size=decorator.get_font_size('ax_label'))
        ax.set_ylabel('|Gain| (1/Wm)', size=decorator.get_font_size('ax_label'))

        ax.tick_params(labelsize=decorator.get_font_size('ax_tick'))
#ax.xaxis.set_tick_params(width=1.0)
#        ax.yaxis.set_tick_params(width=1.0)

        if not decorator is None:
#ax=plt.gca()
          decorator.extra_axes_commands(ax)


        if pdf_png=='png':
            plt.savefig('%(pre)sgain_spectra-MB_PE_comps%(add)s.png' % {'pre' : prefix_str, 'add' : suffix_str})
        elif pdf_png=='pdf':
            plt.savefig('%(pre)sgain_spectra-MB_PE_comps%(add)s.pdf' % {'pre' : prefix_str, 'add' : suffix_str})
        if not keep_plots_open: plt.close()

    if save_txt:
        save_array = (interp_grid, interp_values)
        np.savetxt('%(pre)sgain_spectra-MB_PE_comps%(add)s-Total.csv' 
                    % {'pre' : prefix_str, 'add' : suffix_str}, 
                    save_array, delimiter=',')
        save_array = (interp_grid, interp_values_PE)
        np.savetxt('%(pre)sgain_spectra-MB_PE_comps%(add)s-PE.csv' 
                    % {'pre' : prefix_str, 'add' : suffix_str}, 
                    save_array, delimiter=',')
        save_array = (interp_grid, interp_values_MB)
        np.savetxt('%(pre)sgain_spectra-MB_PE_comps%(add)s-MB.csv' 
                    % {'pre' : prefix_str, 'add' : suffix_str}, 
                    save_array, delimiter=',')

    if dB:
        plt.figure()
        plt.clf()

        max_G = np.max(interp_values)
        dB_const = dB_peak_amp/(4.34*max_G)
        plt.plot(interp_grid, np.abs(10*np.log10(np.exp(abs(interp_values)*dB_const))), 'b', linewidth=3, label="Total")
        plt.legend(loc=0)
        if freq_min and freq_max:
            plt.xlim(freq_min,freq_max)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Amplification (dB)')

        if pdf_png=='png':
            plt.savefig('%(pre)sgain_spectra-dB%(add)s.png' % {'pre' : prefix_str, 'add' : suffix_str})
        elif pdf_png=='pdf':
            plt.savefig('%(pre)sgain_spectra-dB%(add)s.pdf' % {'pre' : prefix_str, 'add' : suffix_str})
        if not keep_plots_open: plt.close()

        if save_txt:
            save_array = (interp_grid, 10*np.log10(np.exp(abs(interp_values)*dB_const)))
            np.savetxt('%(pre)sgain_spectra-dB%(add)s.csv' 
                        % {'pre' : prefix_str, 'add' : suffix_str}, 
                        save_array, delimiter=',')

    if semilogy and save_fig:
        plt.figure()
        plt.clf()
        plt.semilogy(interp_grid, abs(interp_values_PE), 'r', linewidth=3, label="PE")
        plt.semilogy(interp_grid, abs(interp_values_MB), 'g', linewidth=3, label="MB")
        plt.semilogy(interp_grid, abs(interp_values), 'b', linewidth=2, label="Total")
        plt.legend(loc=0)
        if freq_min and freq_max:
            plt.xlim(freq_min,freq_max)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('|Gain| (1/Wm)')

        if pdf_png=='png':
            plt.savefig('%(pre)sgain_spectra-MB_PE_comps-logy%(add)s.png' % {'pre' : prefix_str, 'add' : suffix_str})
        elif pdf_png=='pdf':
            plt.savefig('%(pre)sgain_spectra-MB_PE_comps-logy%(add)s.pdf' % {'pre' : prefix_str, 'add' : suffix_str})
        if not keep_plots_open: plt.close()

    return return_interp_values

def plot_one_component_axes(ax, m_X, m_Y, fields, plps, cc):
  decorator = plps['decorator']
  cmap_signed='seismic'
  cmap_unsigned='OrRd'

  plot_threshold = 1e-8 # set negligible components to explicitly zero

  comp_label=cc.get_label()  # for plotting
  signed=cc.is_signed_field()

  cmap=cmap_unsigned
  if signed: cmap=cmap_signed

  field=fields[cc._f_code]


  if plps['ticks']: #set tick range
#xm=v_x[-1]+v_x[0]
#ym=v_y[-1]+v_y[0]
#extents=np.array([v_x[0]-xm/2, v_x[-1]-xm/2, v_y[0]-ym/2, v_y[-1]-ym/2])*1e6  # Convert to length units in microns
    x0=m_X[0,0]
    x1=m_X[0,-1]
    y0=m_Y[0,0]
    y1=m_Y[-1,0]

    xm = x0+x1
    ym = y0+y1
    extents=np.array([x0-xm/2, x1-xm/2, y0-ym/2, y1-ym/2])*1e6  # Convert to length units in microns
  else:
    extents=None

  req_lims=False

  if decorator.get_cmap_limits(cc._xyz)!=None: #User requested specific limits for each component x, y or z
    req_lims=True
    (req_zlo, req_zhi)= decorator.get_cmap_limits(cc._xyz)
    req_tsnorm=mplcolors.TwoSlopeNorm(vmin=req_zlo, vmax=req_zhi, vcenter=(req_zlo+req_zhi)/2)

  act_zlo=0
  act_zhi=0

  field_final=field.T

  # if the data is all noise, just plot zeros
  if np.max(np.abs(field[~np.isnan(field)])) < plot_threshold: 
    field_final=np.zeros(np.shape(field.T))

  interp=None   #interp='bilinear';
  tsnorm=None

  if req_lims:
    tsnorm=req_tsnorm
    act_zlo=req_zlo
    act_zhi=req_zhi
  else: 
    zlo=np.nanmin(field)
    zhi=np.nanmax(field)
    act_zlo=zlo
    act_zhi=zhi
    if signed: # ensure that zero values get mapped to the right part of the colormap
      if zlo<0 and zhi >0:
        vma=max(abs(zlo),abs(zhi))
#        tsnorm=mplcolors.TwoSlopeNorm(vmin=-vma, vmax=vma, vcenter=0.0)
#tsnorm=mplcolors.CenteredNorm()
      else: # do this anyway, even if zlo and zhi are the same sign
        vma=max(abs(zlo),abs(zhi))
#        tsnorm=mplcolors.TwoSlopeNorm(vmin=-vma, vmax=vma, vcenter=0.0)
#tsnorm=mplcolors.CenteredNorm() # Requires mpl > 3.4.2
    else:
      vma=max(abs(zlo),abs(zhi))


#TODO: Clean this up!
#im = plt.imshow(field_final, origin='lower', extent=extents, interpolation=interp, cmap=cmap, norm=tsnorm)
    if req_lims:
      im = plt.imshow(field_final, origin='lower', extent=extents, interpolation=interp, cmap=cmap, norm=tsnorm)
    else:
      if signed:
        im = plt.imshow(field_final, origin='lower', extent=extents, interpolation=interp, cmap=cmap, vmin=-vma, vmax=vma)
      else:
        im = plt.imshow(field_final, origin='lower', extent=extents, interpolation=interp, cmap=cmap, vmin=0, vmax=vma)

  axes = plt.gca()
  xmin, xmax = axes.get_xlim()
  ymin, ymax = axes.get_ylim()
  width_x = xmax-xmin
  width_y = ymax-ymin
  
  if plps['xlim_min'] != None:
      ax.set_xlim(xmin+plps['xlim_min']*width_x,xmax-plps['xlim_max']*width_x)
  if plps['ylim_min'] != None:
      ax.set_ylim(ymin+plps['ylim_min']*width_y,ymax-plps['ylim_max']*width_y)

  if plps['ticks']:
    plt.xticks()
    plt.yticks()
    ax.tick_params(labelsize=decorator.get_font_size('ax_tick'))
    ax.xaxis.set_tick_params(width=1.0)
    ax.yaxis.set_tick_params(width=1.0)
    ax.set_xlabel('$x$ [$\mathrm{\mu}$m]', size=decorator.get_font_size('ax_label'))
    ax.set_ylabel('$y$ [$\mathrm{\mu}$m]', size=decorator.get_font_size('ax_label'))
    ax.grid(False)
  else:
    plt.xticks([])
    plt.yticks([])

  plt.rc('axes',titlepad=decorator.get_axes_property('title_pad'))
  plt.title(comp_label,fontsize=decorator.get_font_size('subplot_title'))


  # colorbar
  docbar=plps['colorbar']
  if docbar:
    divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="2%", pad=-6)
    cax = divider.append_axes("right", size=decorator.get_axes_property('cbar_size'), 
        pad=decorator.get_axes_property('cbar_pad'))
    cbar = plt.colorbar(im, cax=cax)
    nt=plps.get('num_ticks')
    if nt==None: nt=7
    cbarticks = np.linspace(act_zlo, act_zhi, nt)
    cbar.set_ticks(cbarticks)
    cbarlabels = ['%.2f' %t for t in cbarticks]
    cbar.set_ticklabels(cbarlabels)
    cbar.ax.tick_params(labelsize=decorator.get_font_size('cbar_tick'))

  if plps['contours']:
    if plps['contour_lst']:
        if docbar: cbarticks = plps['contour_lst']
    if np.max(np.abs(field[~np.isnan(field)])) > plot_threshold:
        CS2 = ax.contour(m_X, m_Y, field.T, levels=cbarticks, colors=mycolors[::-1], linewidths=(1.5,))
        if docbar: cbar.add_lines(CS2)

  if decorator!=None: decorator.extra_axes_commands(ax)

def plot_filename(plps, ival, label=None):
  filestart='%(pre)sfields/%(s)s_field_%(i)i%(add)s' % {'pre' : plps['prefix_str'], 
          's' : plps['EM_AC'], 'i' : ival, 'add' : plps['suffix_str']}
  if label!=None: filestart+='_'+label
 
  if plps['pdf_png']=='png': fig_filename=filestart+'.png'
  elif plps['pdf_png']=='pdf': fig_filename=filestart+'.pdf'
  else: raise ValueError("pdf_png must be either 'png' or 'pdf'.")

  return fig_filename

def plot_one_component_quiver(ax, m_X, m_Y, v_fields, plps, cc):
  decorator = plps['decorator']

  quiver_points=20
  xlmi, xlma, ylmi, ylma= plps['xlim_min'], plps['xlim_max'], plps['ylim_min'], plps['ylim_max']
  n_pts_x=plps['n_pts_x']
  n_pts_y=plps['n_pts_y']
  quiver_points=plps['quiver_points']

  if xlmi is None: xlmi=0
  if xlma is None: xlma=0
  if ylmi is None: ylmi=0
  if ylma is None: ylma=0
  quiver_steps_x=int(round(min(n_pts_x,n_pts_y)/quiver_points *(1-xlmi-xlma) )) # this could probably be chosen nicer
  quiver_steps_y=int(round(min(n_pts_x,n_pts_y)/quiver_points *(1-ylmi-ylma) )) # this could probably be chosen nicer

  v_x_q = m_X.T[0::quiver_steps_x, 0::quiver_steps_y] # is quiver_steps_x and _y around the right way given the .T operation?
  v_y_q = m_Y.T[0::quiver_steps_x, 0::quiver_steps_y]

  m_ReEx_q = v_fields['Fxr'][0::quiver_steps_x, 0::quiver_steps_y]
  m_ReEy_q = v_fields['Fyr'][0::quiver_steps_x, 0::quiver_steps_y]
  m_ImEx_q = v_fields['Fxi'][0::quiver_steps_x, 0::quiver_steps_y]
  m_ImEy_q = v_fields['Fyi'][0::quiver_steps_x, 0::quiver_steps_y]

  # convert to microns
  v_x_q_um=v_x_q*1e6
  v_y_q_um=v_y_q*1e6

  # centre at zero
  xm=v_x_q_um[-1,0]+v_x_q_um[0,0]
  ym=v_y_q_um[0,-1]+v_y_q_um[0,0]
  v_x_q_um-=xm/2
  v_y_q_um-=ym/2


# Ignore all imaginary values. If there are significant imag values, 
# then instaneous vector plots don't make much sense anyway
  m_arrcolour= np.sqrt(m_ReEx_q*m_ReEx_q +m_ReEy_q*m_ReEy_q)
  plt.quiver(v_x_q_um, v_y_q_um, m_ReEx_q, m_ReEy_q, m_arrcolour,
      linewidths=(0.2,), edgecolors=('k'), pivot='mid', headlength=5) # length of the arrows

  ax.set_aspect('equal')
  ax.set_xlim(v_x_q_um[0,0],v_x_q_um[-1,0]) # this step is needed because quiver doesn't seem
  ax.set_ylim(v_y_q_um[0,0],v_y_q_um[0,-1]) # to use its input x and y vectors to set range limits
                                            # clean this up so as to avoid seemingly circular calls following
  axes = plt.gca()
  xmin, xmax = axes.get_xlim()
  ymin, ymax = axes.get_ylim()
  width_x = xmax-xmin
  width_y = ymax-ymin

  if plps['ticks']:
    plt.xticks()
    plt.yticks()
    ax.tick_params(labelsize=decorator.get_font_size('ax_tick'))
    ax.xaxis.set_tick_params(width=1.0)
    ax.yaxis.set_tick_params(width=1.0)
    ax.set_xlabel('$x$ [$\mathrm{\mu}$m]', size=decorator.get_font_size('ax_label'))
    ax.set_ylabel('$y$ [$\mathrm{\mu}$m]', size=decorator.get_font_size('ax_label'))
    ax.grid(False)
  else:
    plt.xticks([])
    plt.yticks([])

#return
  if plps['xlim_min'] != None:
      ax.set_xlim(xmin+plps['xlim_min']*width_x,xmax-plps['xlim_max']*width_x)
  if plps['ylim_min'] != None:
      ax.set_ylim(ymin+plps['ylim_min']*width_y,ymax-plps['ylim_max']*width_y)

  plt.rc('axes',titlepad=decorator.get_axes_property('title_pad'))


  plt.title(cc.get_label(), fontsize=decorator.get_font_size('subplot_title'))

  decorator = plps['decorator']
  if decorator!=None: decorator.extra_axes_commands(ax)

def plot_mode_data(ax, v_fields, plps, sim_wguide, ival, v_x, v_y):  # mode data summary

  #TODO: Move this analysis to a better place                                                                                                 
  [m_ReEx, m_ReEy, m_ReEz, m_ImEx, m_ImEy, m_ImEz, m_AbsE]=[ v_fields['Fxr'], v_fields['Fyr'], 
     v_fields['Fzr'], v_fields['Fxi'], v_fields['Fyi'], v_fields['Fzi'], v_fields['Fabs'] ]
  tms=sim_wguide.get_modes()
  modeobj=sim_wguide.get_modes()[ival]
  modeobj.analyse_mode(v_x, v_y, m_ReEx, m_ReEy, m_ReEz, m_ImEx, m_ImEy, m_ImEz, m_AbsE)

  yhi=.99  # try and make these ranges match those of the physical domain?
  dy=.10
  decorator = plps['decorator']
  ax.set_xlim(0,1)
  ax.set_ylim(0,yhi)
  ax.set_aspect('equal')
  ax.axis('off')
  fs=decorator.get_font_size('data_label')
#plt.title('Mode {0}'.format(ival), fontsize=decorator.get_font_size('subplot_title'))

  #  Beward that this is a very coarse measure of component fraction and is not consistent with energy fraction

  [f_x, f_y, f_t, f_z] = modeobj.field_fracs()
  (r0x, r0y)=modeobj.center_of_mass()*1e6
  (wx, wy, w0)=modeobj.second_moment_widths()*1e6

  r=yhi-dy
  x0=.05
  ax.text(x0-.05,r,r'Mode properties: m={0}'.format(ival), transform=ax.transAxes, fontsize=fs+2); r-=dy
  if sim_wguide.EM_AC=='EM':
    ax.text(x0,r,r'$\omega/(2\pi)$: {0:.4f} THz'.format(sim_wguide.omega_EM/(2*np.pi*1e12)), transform=ax.transAxes, fontsize=fs); r-=dy
    ax.text(x0,r,r'$k$: {0:.2f} m$^{{-1}}$'.format(sim_wguide.kz_EM(ival)), transform=ax.transAxes, fontsize=fs); r-=dy
    if sim_wguide.ngroup_EM_available():
      ax.text(x0,r,r'$\bar{{n}}$: {0:.6f}, $n_g$: {1:.6f}'.format(sim_wguide.neff(ival), sim_wguide.ngroup_EM(ival)), transform=ax.transAxes, fontsize=fs); r-=dy  # TODO: 
    else:
      ax.text(x0,r,r'$\bar{{n}}$: {0:.6f}'.format(sim_wguide.neff(ival)), transform=ax.transAxes, fontsize=fs); r-=dy  
  else:
    ax.text(x0,r,r'$q$: {0:.4e} m$^{{-1}}$, $\lambda:$ {1:.4f} $\mu$m'.format(sim_wguide.k_AC, 2*np.pi/sim_wguide.k_AC*1e6), transform=ax.transAxes, fontsize=fs); r-=dy
    ax.text(x0,r,r'$\Omega/(2\pi)$: {0:.4f} GHz'.format(sim_wguide.nu_AC(ival)/1.e9), transform=ax.transAxes, fontsize=fs); r-=dy
    if sim_wguide.vgroup_AC_available():
      ax.text(x0,r,r'$v_p$: {0:.2f} m/s, $v_g$: {1:.2f} m/s'.format( sim_wguide.vp_AC(ival), sim_wguide.vg_AC(ival)), transform=ax.transAxes, fontsize=fs); r-=dy
    else:
      ax.text(x0,r,r'$v_p$: {0:.2f} m/s'.format( sim_wguide.vp_AC(ival)), transform=ax.transAxes, fontsize=fs); r-=dy
  ax.text(x0,r,r'$f_x:$ {0:.3f}, $f_y$: {1:.3f}'.format(f_x, f_y), transform=ax.transAxes, fontsize=fs); r-=dy
  ax.text(x0,r,r'$f_t:$ {0:.3f}, $f_z$: {1:.3f}'.format(f_t, f_z), transform=ax.transAxes, fontsize=fs); r-=dy
  ax.text(x0,r,r'$\mathbf{{r}}_0:$ ({0:.2f}, {1:.3f}) $\mu$m'.format(r0x, r0y), transform=ax.transAxes, fontsize=fs); r-=dy
  ax.text(x0,r,r'$(w_x, w_y):$ ({0:.2f}, {1:.2f}) $\mu$m'.format(wx, wy), transform=ax.transAxes, fontsize=fs); r-=dy
  sc = sim_wguide.symmetry_classification(ival)
  if len(sc): ax.text(x0,r,r'Sym: {0}'.format(sc), transform=ax.transAxes, fontsize=fs); r-=dy

  if sim_wguide.EM_AC=='AC' and sim_wguide.Q_method != QAcMethod.NotSet:
      if sim_wguide.Q_method==QAcMethod.Intrinsic:
        ax.text(x0,r,'Losses (intrinsic):', transform=ax.transAxes, fontsize=fs); r-=dy
      else:
        ax.text(x0,r,'Losses (fixed Q):', transform=ax.transAxes, fontsize=fs); r-=dy

      ax.text(x0+.1,r,r'$\alpha$: {0:.3e} s$^{{-1}}$, {1:.2f} cm$^{{-1}}$, {2:.2f} dB/cm'.format(
           sim_wguide.alpha_t_AC(ival),
           sim_wguide.alpha_s_AC(ival)/100.,
           sim_wguide.alpha_s_AC(ival)/100./(np.log(10.0)/10.0)
           ), transform=ax.transAxes, fontsize=fs); r-=dy
      ax.text(x0+.1,r,r'$Q_m$: {0:.2f}'.format(sim_wguide.Qmech_AC(ival)), transform=ax.transAxes, fontsize=fs); r-=dy
      ax.text(x0+.1,r,r'$\Delta\Omega/(2\pi)$: {0:.4f} MHz'.format(1.e-6*sim_wguide.linewidth_AC(ival)), transform=ax.transAxes, fontsize=fs); r-=dy
      mg_pe=plps['modal_gain'].get('PE',0)
      mg_mb=plps['modal_gain'].get('MB',0)
      mg_tot=plps['modal_gain'].get('Tot',0)
      ax.text(x0+.1,r,r'Gain (PE, MB, Tot): ({0:.4f}, {1:.4f}, {2:.4f} ) $\mathrm{{(mW)}}^{{-1}}$'.format(mg_pe, mg_mb, mg_tot), transform=ax.transAxes, fontsize=fs); r-=dy
  d_ext=modeobj.get_mode_data()
  if len(d_ext):
    ax.text(x0,r,r'Extra data:', transform=ax.transAxes, fontsize=fs); r-=dy
    for (k, v) in d_ext.items():
      ax.text(x0,r,r'   {0}: {1}'.format(k,v), transform=ax.transAxes, fontsize=fs); r-=dy





def plot_all_components(v_x, v_y, m_X, m_Y, v_plots, plps, sim_wguide, ival, suppress_imimre, EM_AC):
  decorator = plps['decorator']
  figsz=decorator.get_axes_property('figsize')
  ws=decorator.get_axes_property('subplots_wspace')
  hs=decorator.get_axes_property('subplots_hspace')
  lw=decorator.get_axes_property('linewidth')
  ec=decorator.get_axes_property('edgecolor')

  plt.clf()
  fig = plt.figure(figsize=figsz)
  fig.subplots_adjust(wspace=ws, hspace=hs)
  plt.rc('axes', linewidth=lw)
  plt.rc('axes', edgecolor=ec)

  rows=3
  if suppress_imimre: rows=2

  axi=1
  ax = plt.subplot(rows,3,axi)
  plot_mode_data(ax, v_plots, plps, sim_wguide, ival, v_x, v_y)  # mode data summary
  axi+=1

  cc={ 'EM_E':component_t('Eabs'), 'EM_H':component_t('Habs'), 'AC':component_t('uabs')}[EM_AC]

  ax = plt.subplot(rows,3,axi)
  plot_one_component_axes(ax, m_X, m_Y, v_plots, plps, cc)  # the intensity plot
  axi+=1

  cc={ 'EM_E':component_t('Et'), 'EM_H':component_t('Ht'), 'AC':component_t('ut')}[EM_AC]
  ax = plt.subplot(rows,3,axi)
  plot_one_component_quiver(ax, m_X, m_Y, v_plots, plps, cc)  # the transverse vector plot
  axi+=1

  for (flab, field) in v_plots.items():
     cc=component_t.make_comp(EM_AC, flab)

     if suppress_imimre and not cc.is_dominant() :  continue 
     if not cc.is_signed_field():continue
     ax = plt.subplot(rows,3,axi)
     plot_one_component_axes(ax, m_X, m_Y, v_plots, plps, cc)  # the scalar plots
     axi+=1

   
  figfile=plot_filename(plps, ival)
  save_figure(plt, figfile)

  if not keep_plots_open: plt.close()


  if plps['EM_AC']=='AC' and plps['stress_fields']:
     ### Interpolate onto rectangular Cartesian grid
     xy = list(zip(v_x6p, v_y6p))
     grid_x, grid_y = np.mgrid[x_min:x_max:n_pts_x*1j, y_min:y_max:n_pts_y*1j]
     m_ReEx = interpolate.griddata(xy, v_Ex6p.real, (grid_x, grid_y), method='linear')
     m_ReEy = interpolate.griddata(xy, v_Ey6p.real, (grid_x, grid_y), method='linear')
     m_ReEz = interpolate.griddata(xy, v_Ez6p.real, (grid_x, grid_y), method='linear')
     m_ImEx = interpolate.griddata(xy, v_Ex6p.imag, (grid_x, grid_y), method='linear')
     m_ImEy = interpolate.griddata(xy, v_Ey6p.imag, (grid_x, grid_y), method='linear')
     m_ImEz = interpolate.griddata(xy, v_Ez6p.imag, (grid_x, grid_y), method='linear')
     m_AbsE = interpolate.griddata(xy, v_E6p.real, (grid_x, grid_y), method='linear')
     dx = grid_x[-1,0] - grid_x[-2,0]
     dy = grid_y[0,-1] - grid_y[0,-2]
     m_Ex = m_ReEx + 1j*m_ImEx
     m_Ey = m_ReEy + 1j*m_ImEy
     m_Ez = m_ReEz + 1j*m_ImEz
     m_Ex = m_Ex.reshape(n_pts_x,n_pts_y)
     m_Ey = m_Ey.reshape(n_pts_x,n_pts_y)
     m_Ez = m_Ez.reshape(n_pts_x,n_pts_y)
     m_AbsE = m_AbsE.reshape(n_pts_x,n_pts_y)

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
     del_z_Ex = 1j*sim_wguide.k_AC*m_Ex
     del_z_Ey = 1j*sim_wguide.k_AC*m_Ey
     del_z_Ez = 1j*sim_wguide.k_AC*m_Ez

     # Flip y order as imshow has origin at top left
     del_mat = np.array([del_x_Ex[:,::-1].real, del_x_Ey[:,::-1].real, del_x_Ez[:,::-1].real, del_x_Ex[:,::-1].imag, del_x_Ey[:,::-1].imag, del_x_Ez[:,::-1].imag, del_y_Ex[:,::-1].real, del_y_Ey[:,::-1].real, del_y_Ez[:,::-1].real, del_y_Ex[:,::-1].imag, del_y_Ey[:,::-1].imag, del_y_Ez[:,::-1].imag, del_z_Ex[:,::-1].real, del_z_Ey[:,::-1].real, del_z_Ez[:,::-1].real, del_z_Ex[:,::-1].imag, del_z_Ey[:,::-1].imag, del_z_Ez[:,::-1].imag])
     v_labels = ["Re($S_{xx}$)","Re($S_{xy}$)","Re($S_{xz}$)", "Im($S_{xx}$)","Im($S_{xy}$)","Im($S_{xz}$)","Re($S_{yx}$)","Re($S_{yy}$)","Re($S_{yz}$)","Im($S_{yx}$)","Im($S_{yy}$)","Im($S_{yz}$)","Re($S_{zx}$)","Re($S_{zy}$)","Re($S_{zz}$)","Im($S_{zx}$)","Im($S_{zy}$)","Im($S_{zz}$)"]

     # stress field plots
     plt.clf()
     fig = plt.figure(figsize=(15,30))
     for i_p,plot in enumerate(del_mat):
         ax = plt.subplot(6,3,i_p+1)
         im = plt.imshow(plot.T);
         # no ticks
         plt.xticks([])
         plt.yticks([])
         # limits
         if xlim_min != None:
             ax.set_xlim(xlim_min*n_points,(1-xlim_max)*n_points)
         if ylim_min != None:
             ax.set_ylim((1-ylim_min)*n_points,ylim_max*n_points)
         # titles
         plt.title(v_labels[i_p],fontsize=decorator.get_font_size('subplot_title'))
         # colorbar
         divider = make_axes_locatable(ax)
         cax = divider.append_axes("right", size="5%", pad=0.1)
         cbar = plt.colorbar(im, cax=cax, format='%.2e')
         if num_ticks:
             cbarticks = np.linspace(np.min(plot), np.max(plot), num=num_ticks)                
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
         cbarlabels = ['%.2f' %t for t in cbarlabels]
         cbar.set_ticklabels(cbarlabels)
         if contours:
             if contour_lst:
                 cbarticks = contour_lst
             if np.max(np.abs(plot[~np.isnan(plot)])) > plot_threshold:
                 CS2 = ax.contour(m_X, m_Y, plot.T, levels=cbarticks, colors=colors[::-1], linewidths=(1.5,))
             cbar.add_lines(CS2)
         cbar.ax.tick_params(labelsize=decorator.get_font_size('cbar_tick'))
     fig.set_tight_layout(True)
     n_str = ''
     if np.imag(sim_wguide.Eig_values[ival]) < 0:
         k_str = r'$\Omega/2\pi = %(re_k)f %(im_k)f i$ GHz'% \
             {'re_k' : np.real(sim_wguide.Eig_values[ival]*1e-9),
             'im_k' : np.imag(sim_wguide.Eig_values[ival]*1e-9)}
     else:
         k_str = r'$\Omega/2\pi = %(re_k)f + %(im_k)f i$ GHz'% \
             {'re_k' : np.real(sim_wguide.Eig_values[ival]*1e-9),
             'im_k' : np.imag(sim_wguide.Eig_values[ival]*1e-9)}
     plt.suptitle('Mode #' + str(ival) + '   ' + k_str + '   ' + n_str, fontsize=decorator.get_font_size('title'))

     if pdf_png=='png':
         plt.savefig('%(pre)sfields/%(s)s_S_field_%(i)i%(add)s.png' %
             {'pre' : prefix_str, 's' : EM_AC, 'i' : ival, 'add' : suffix_str})
     elif pdf_png=='pdf':
         plt.savefig('%(pre)sfields/%(s)s_S_field_%(i)i%(add)s.pdf' %
             {'pre' : prefix_str, 's' : EM_AC, 'i' : ival, 'add' : suffix_str}, bbox_inches='tight')
     if not keep_plots_open: plt.close()


def save_figure(plt, figfile):
  if figfile[-3:-1]=='png': plt.savefig(figfile)
  else: plt.savefig(figfile, bbox_inches='tight')


def plot_one_component(m_X, m_Y, v_fields, plps, sim_wguide, ival, cc):
  plt.clf()
  fig = plt.figure(figsize=(12,10))

  decorator = plps['decorator']
  plt.rc('axes', linewidth=decorator.get_axes_property('linewidth'))
  plt.rc('axes', edgecolor=decorator.get_axes_property('edgecolor'))

  ax = plt.subplot(111)
  if cc.is_transverse():
    plot_one_component_quiver(ax, m_X, m_Y, v_fields, plps, cc)
  else: 
    plot_one_component_axes(ax, m_X, m_Y, v_fields, plps, cc)

  figfile=plot_filename(plps, ival, cc._user_code)
  save_figure(plt, figfile)
  if not keep_plots_open: plt.close()


#deprecated spelling.
def plt_mode_fields(sim_wguide, ivals=None, n_points=501, quiver_points=50, 
                  xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None,
                  EM_AC='EM_E', num_ticks=None, colorbar=True, contours=False, contour_lst=None,
                  stress_fields=False, pdf_png='png', 
                  prefix_str='', suffix_str='', ticks=False, comps=[], decorator=None, 
                  suppress_imimre=True, 
                  modal_gains_PE=None,
                  modal_gains_MB=None,
                  modal_gains=None):
  plot_mode_fields(sim_wguide, ivals, n_points, quiver_points, 
                  xlim_min, xlim_max, ylim_min, ylim_max,
                  EM_AC, num_ticks, colorbar, contours, contour_lst, stress_fields, pdf_png, 
                  prefix_str, suffix_str, ticks, comps, decorator, 
                  suppress_imimre, modal_gains_PE, modal_gains_MB, modal_gains)

#### Standard plotting of spectra #############################################
def plot_mode_fields(sim_wguide, ivals=None, n_points=501, quiver_points=30, 
                  xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None,
                  EM_AC='EM_E', num_ticks=None, colorbar=True, contours=False, contour_lst=None,
                  stress_fields=False, pdf_png='png', 
                  prefix_str='', suffix_str='', ticks=False, comps=[], decorator=None, 
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

            prefix_str  (str): Add a string to start of file name
            
            suffix_str  (str): Add a string to end of file name.

            modal_gains (float array): Pre-calculated gain for each acoustic mode given chosen optical fields.
    """


    if EM_AC != 'EM_E' and EM_AC != 'EM_H' and EM_AC != 'AC':
        raise ValueError("EM_AC must be either 'AC', 'EM_E' or 'EM_H'.")

#if EM_AC == 'AC': sim_wguide.calc_acoustic_losses()

    # Calculate the magnetic field from the electric field
    if EM_AC == 'EM_H':
        nnodes = 6
        sim_wguide.sol1_H = NumBAT.h_mode_field_ez(sim_wguide.k_0, sim_wguide.num_modes, 
            sim_wguide.n_msh_el, sim_wguide.n_msh_pts, nnodes, sim_wguide.table_nod, 
            sim_wguide.x_arr, sim_wguide.Eig_values, sim_wguide.sol1)


    plt.clf()

    # field mapping
    x_tmp = []
    y_tmp = []
    for i in np.arange(sim_wguide.n_msh_pts):
        x_tmp.append(sim_wguide.x_arr[0,i])
        y_tmp.append(sim_wguide.x_arr[1,i])
    x_min = np.min(x_tmp); x_max=np.max(x_tmp)
    y_min = np.min(y_tmp); y_max=np.max(y_tmp)
    area = abs((x_max-x_min)*(y_max-y_min))
    n_pts_x = int(n_points*abs(x_max-x_min)/np.sqrt(area))
    n_pts_y = int(n_points*abs(y_max-y_min)/np.sqrt(area))

    v_x=np.linspace(x_min,x_max,n_pts_x)
    v_y=np.linspace(y_min,y_max,n_pts_y)

    m_X, m_Y = np.meshgrid(v_x, v_y)

    # unrolling data for the interpolators
    table_nod = sim_wguide.table_nod.T
    x_arr = sim_wguide.x_arr.T

    if ivals:
        ival_range = ivals
    else:
        ival_range = range(len(sim_wguide.Eig_values))

    for ival in ival_range:
        # dense triangulation with multiple points
        v_x6p = np.zeros(6*sim_wguide.n_msh_el)
        v_y6p = np.zeros(6*sim_wguide.n_msh_el)
        v_Ex6p = np.zeros(6*sim_wguide.n_msh_el, dtype=np.complex128)
        v_Ey6p = np.zeros(6*sim_wguide.n_msh_el, dtype=np.complex128)
        v_Ez6p = np.zeros(6*sim_wguide.n_msh_el, dtype=np.complex128)
        v_triang6p = []

        i = 0
        for i_el in np.arange(sim_wguide.n_msh_el):
            # triangles
            idx = np.arange(6*i_el, 6*(i_el+1))
            triangles = [[idx[0], idx[3], idx[5]],
                         [idx[1], idx[4], idx[3]],
                         [idx[2], idx[5], idx[4]],
                         [idx[3], idx[4], idx[5]]]
            v_triang6p.extend(triangles)

            for i_node in np.arange(6):
                # index for the coordinates
                i_ex = table_nod[i_el, i_node]-1
                # values
                v_x6p[i] = x_arr[i_ex, 0]
                v_y6p[i] = x_arr[i_ex, 1]

                if EM_AC == 'EM_E' or EM_AC == 'AC' :
                    v_Ex6p[i] = sim_wguide.sol1[0,i_node,ival,i_el]
                    v_Ey6p[i] = sim_wguide.sol1[1,i_node,ival,i_el]
                    v_Ez6p[i] = sim_wguide.sol1[2,i_node,ival,i_el]
                if EM_AC == 'EM_H':
                    v_Ex6p[i] = sim_wguide.sol1_H[0,i_node,ival,i_el]
                    v_Ey6p[i] = sim_wguide.sol1_H[1,i_node,ival,i_el]
                    v_Ez6p[i] = sim_wguide.sol1_H[2,i_node,ival,i_el]
                i += 1

        v_E6p = np.sqrt(np.abs(v_Ex6p)**2 + np.abs(v_Ey6p)**2 + np.abs(v_Ez6p)**2)

        ### Interpolate onto triangular grid - honest to FEM elements
        # dense triangulation with unique points
        v_triang1p = []
        for i_el in np.arange(sim_wguide.n_msh_el):
            # triangles
            triangles = [[table_nod[i_el,0]-1,table_nod[i_el,3]-1,table_nod[i_el,5]-1],
                         [table_nod[i_el,1]-1,table_nod[i_el,4]-1,table_nod[i_el,3]-1],
                         [table_nod[i_el,2]-1,table_nod[i_el,5]-1,table_nod[i_el,4]-1],
                         [table_nod[i_el,3]-1,table_nod[i_el,4]-1,table_nod[i_el,5]-1]]
            v_triang1p.extend(triangles)

        # triangulations
        triang6p = matplotlib.tri.Triangulation(v_x6p,v_y6p,v_triang6p)
        triang1p = matplotlib.tri.Triangulation(x_arr[:,0],x_arr[:,1],v_triang1p)
        
        #fields are called Ex, Ey, etc regardless of whether we are plotting E, H, or u/S


        # building interpolators: triang1p for the finder, triang6p for the values
        #TODO: could be more efficient only interpolating the fields which are ultimately to be used?
        v_x_flat=m_X.flatten('F')  # There might be a cleaner way of doing this
        v_y_flat=m_Y.flatten('F')
        finder = matplotlib.tri.TrapezoidMapTriFinder(triang1p)
        interper_f = lambda x: matplotlib.tri.LinearTriInterpolator(
             triang6p, x, trifinder=finder)(v_x_flat, v_y_flat).reshape(n_pts_x,n_pts_y)

        # interpolated fields.  Always need these ones.
        m_ReEx=interper_f(v_Ex6p.real)
        m_ReEy=interper_f(v_Ey6p.real)
        m_ImEz=interper_f(v_Ez6p.imag)
        m_AbsE=interper_f(v_E6p)

        # often not needed for plotting, but are used for measuring fractions. (Could fix taht?)
        m_ImEx=interper_f(v_Ex6p.imag)
        m_ImEy=interper_f(v_Ey6p.imag)
        m_ReEz=interper_f(v_Ez6p.real)

        v_plots = {'Fxr':m_ReEx, 'Fyr':m_ReEy, 'Fzi':m_ImEz, 'Fxi':m_ImEx, 'Fyi':m_ImEy, 'Fzr':m_ReEz, 'Fabs':m_AbsE}

        if decorator == None: decorator = Decorator()

        modal_gain={}
        if modal_gains is not None: modal_gain['Tot']=modal_gains[ival]
        if modal_gains_PE is not None: modal_gain['PE']=modal_gains_PE[ival]
        if modal_gains_MB is not None: modal_gain['MB']=modal_gains_MB[ival]


        plot_params={'xlim_min': xlim_min, 'xlim_max': xlim_max, 'ylim_min': ylim_min, 
                     'ylim_max': ylim_max, 'ticks': ticks, 'num_ticks':num_ticks,
                      'colorbar':colorbar, 'contours':contours, 'contour_lst':contour_lst, 'EM_AC':EM_AC,
                      'prefix_str': prefix_str, 'suffix_str': suffix_str, 'pdf_png': pdf_png, 
                      'stress_fields':stress_fields, 'modal_gain':modal_gain, 'decorator': decorator,
       'n_pts_x': n_pts_x, 'n_pts_y': n_pts_y, 'quiver_points': quiver_points }

        #TODO: should test that xlim_min + xlim_max < 1 so as some picture survies (and there is no flip)
        #TODO: should test that ylim_min + ylim_max < 1 so as some picture survies (and there is no flip)

        if not os.path.exists("%sfields" % prefix_str): os.mkdir("%sfields" % prefix_str)

        decorator._set_for_multi()
        plot_all_components(v_x, v_y, m_X, m_Y, v_plots, 
            plot_params, sim_wguide, ival, suppress_imimre, EM_AC)


        decorator._set_for_single()
        for comp in comps: #options are ['Ex', 'Hx', 'ux', 'Ey', 'Hy', 'uy', 'Ez', 'Hz', 'uz','Eabs', 'Habs', 'uabs', 'Et', 'Ht', 'ut'] 
          cc=component_t(comp)
          plot_one_component(m_X, m_Y, v_plots, plot_params, sim_wguide, ival, cc)



#### Plot mesh #############################################
def plot_msh(x_arr, prefix_str='', suffix_str=''):
    """ Plot EM mode fields.

        Args:
            sim_wguide : A ``Struct`` instance that has had calc_modes calculated

        Keyword Args:
            n_points  (int): The number of points across unitcell to \
                interpolate the field onto.
    """

    plt.clf()
    plt.figure(figsize=(13,13))
    ax = plt.subplot(1,1,1)
    for node in range(np.shape(x_arr)[1]):
        plt.plot(x_arr[0,node], x_arr[1,node], 'og')
    ax.set_aspect('equal')
    plt.savefig('%(pre)smsh_%(add)s.pdf' %
        {'pre' : prefix_str, 'add' : suffix_str}, bbox_inches='tight')
    if not keep_plots_open: plt.close()



### Plot nodal arrangement on mesh triangle.
# plt.figure(figsize=(13,13))
# el = 1
# plt.clf()
# for i in range(0,6):
#     print table_nod[i][el] - 1
#     x = x_arr[0,table_nod[i][el] - 1]
#     y = x_arr[1,table_nod[i][el] - 1]
#     print 'x1 = ', x_arr[0,table_nod[i][el] - 1]
#     print 'y1 = ', x_arr[1,table_nod[i][el] - 1]
#     plt.plot(x, y, 'o')
#     plt.text(x+0.001, y+0.001, str(i))
# plt.savefig('triangle_%i.png' %el)
if not keep_plots_open: plt.close()
