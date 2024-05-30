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

#import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
from pathlib import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable
from collections.abc import Iterable

import numbat
from numbattools import *
from nbtypes import QAcMethod, FieldType, component_t, twopi
#from fortran import NumBAT
#mport reporting





def plot_filename(plps, ival, label=None):
    #fullpref = str(Path(numbat.NumBATApp().outdir(), plps['prefix']))
    fullpref = str(numbat.NumBATApp().path_fields())

    comp = plps['EM_AC'].name
    suf = plps['suffix']
    filestart = f'{fullpref}/{comp}_field_{ival:02d}{suf}'

    if label is not None:
        filestart += '_'+label

    if plps['pdf_png'] == 'png':
        fig_filename = filestart+'.png'
    elif plps['pdf_png'] == 'pdf':
        fig_filename = filestart+'.pdf'
    else:
        raise ValueError("pdf_png must be either 'png' or 'pdf'.")

    return fig_filename


class TidyAxes:

    def __init__(self, nax=1, props={}):

        self._nax = nax
        self._set_defaults_for_num_axes(nax)
        self._props.update(props)
        
    def _set_defaults_for_num_axes(self, nax):
        props = {}
        
        if nax == 4:
            props['ax_label_fs'] = 8
            props['ax_ticklabel_fs'] = 8
            props['ax_tickwidth'] = .25
            props['ax_linewidth'] = 1
        elif nax == 2:
            props['ax_label_fs'] = 12
            props['ax_ticklabel_fs'] =  10
            props['ax_tickwidth'] = 1
            props['ax_linewthidth'] = 1
        else:
            props['ax_label_fs'] =  12
            props['ax_ticklabel_fs'] =  10
            props['ax_tickwidth'] =  1
            props['ax_linewidth'] =  1

        props['aspect'] = 0.0

        props['cb_linewidth'] = 1 
        props['cb_label_fs'] = 10 
        props['cb_ticklabel_fs'] = 8 
        props['cb_tickwidth'] = 0.25

        props['cb_shrink'] = 0  # possible?
        props['cb_pad'] = 0.    # possible?

                                     
        self._props = props
        
    def apply_to_axes(self, axs):

        if not isinstance(axs, Iterable):
            axs = (axs,)
        
        pr = self._props
        
        #for ax in axs:  #TODO: for some reason, this is not working and have to index explicitly?!
        for i in range(len(axs)):
            ax = axs[i]
            
            if pr['aspect'] >0 : ax.set_aspect(pr['aspect'])

            ax.tick_params(labelsize=pr['ax_ticklabel_fs'], 
                           width=pr['ax_tickwidth'])
            
            ax.xaxis.label.set_size(pr['ax_label_fs'])
            ax.yaxis.label.set_size(pr['ax_label_fs'])

            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(pr['ax_linewidth'])

    def apply_to_cbars(self, cbs):
        if not isinstance(cbs, Iterable):
            cbs = (cbs,)

        pr = self._props
        for cb in cbs:

            lab = cb.ax.get_ylabel()  # don't seem to be able to set label size except by setting label again
            cb.set_label(lab, size=pr['cb_label_fs'])
            
            cb.outline.set_linewidth(pr['cb_linewidth'])
            
            cb.ax.tick_params(labelsize=pr['cb_ticklabel_fs'], 
                           width=pr['cb_tickwidth'])        

            # Not sure how to do these
            
            #if pr['cb_shrink'] >0 :   
            #    cb.set_shrink(pr['cb_shrink'])
            #    cb.set_pad(pr['cb_pad'])

                
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


def plot_mode_data(ax, plps, sim_result, ival, v_x, v_y):  # mode data summary

    modeobj = sim_result.get_all_modes()[ival]
    modeobj.analyse_mode()

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
    if sim_result.is_EM():
        ax.text(x0, r, r'$\omega/(2\pi)$: {0:.5f} THz'.format(
            sim_result.omega_EM/(twopi*1.e12)), transform=ax.transAxes, fontsize=fs)
        r -= dy
        ax.text(x0, r, r'$k$: {0:.5f} $\mu\mathrm{{m}}^{{-1}}$'.format(
            sim_result.kz_EM(ival)/1.e6), transform=ax.transAxes, fontsize=fs)
        r -= dy
        if sim_result.ngroup_EM_available():
            ax.text(x0, r, r'$\bar{{n}}$: {0:.6f}, $n_g$: {1:.6f}'.format(sim_result.neff(
                ival), sim_result.ngroup_EM(ival)), transform=ax.transAxes, fontsize=fs)
            r -= dy  # TODO:
        else:
            ax.text(x0, r, r'$\bar{{n}}$: {0:.6f}'.format(
                sim_result.neff(ival)), transform=ax.transAxes, fontsize=fs)
            r -= dy
    else:
        ax.text(x0, r, r'$q$: {0:.5f} $\mu\mathrm{{m}}^{{-1}}$, $\lambda:$ {1:.5f} $\mu$m'.format(
            sim_result.q_AC/1.e6, twopi/sim_result.q_AC*1e6), transform=ax.transAxes, fontsize=fs)
        r -= dy
        ax.text(x0, r, r'$q/2\pi$: {0:.5f} $\mu\mathrm{{m}}^{{-1}}$'.format(
            sim_result.q_AC/(2.e6*np.pi)), transform=ax.transAxes, fontsize=fs)
        r -= dy
        ax.text(x0, r, r'$\Omega/(2\pi)$: {0:.5f} GHz'.format(
            np.real(sim_result.nu_AC(ival))/1.e9), transform=ax.transAxes, fontsize=fs)
        r -= dy
        if sim_result.vgroup_AC_available():
            ax.text(x0, r, r'$v_p$: {0:.2f} m/s, $v_g$: {1:.2f} m/s'.format(sim_result.vp_AC(
                ival), sim_result.vg_AC(ival)), transform=ax.transAxes, fontsize=fs)
            r -= dy
        else:
            ax.text(x0, r, r'$v_p$: {0:.2f} m/s'.format(
                sim_result.vp_AC(ival)), transform=ax.transAxes, fontsize=fs)
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

    sc = sim_result.symmetry_classification(ival)
    if len(sc):
        ax.text(x0, r, r'Sym: {0}'.format(sc),
                transform=ax.transAxes, fontsize=fs)
        r -= dy

    # TODO: reactivate
    # if sim_result.is_AC() and sim_result.Q_method != QAcMethod.NotSet:
    #     if sim_result.Q_method == QAcMethod.Intrinsic:
    #         ax.text(x0, r, 'Losses (intrinsic):',
    #                 transform=ax.transAxes, fontsize=fs)
    #         r -= dy
    #     else:
    #         ax.text(x0, r, 'Losses (fixed Q):',
    #                 transform=ax.transAxes, fontsize=fs)
    #         r -= dy

    #     ax.text(x0+.1, r, r'$\alpha$: {0:.3e} s$^{{-1}}$, {1:.2f} cm$^{{-1}}$, {2:.2f} dB/cm'.format(
    #         sim_result.alpha_t_AC(ival),
    #         sim_result.alpha_s_AC(ival)/100.,
    #         sim_result.alpha_s_AC(ival)/100./(np.log(10.0)/10.0)
    #     ), transform=ax.transAxes, fontsize=fs)
    #     r -= dy

    #     ax.text(x0+.1, r, r'$Q_m$: {0:.2f}'.format(
    #         sim_result.Qmech_AC(ival)), transform=ax.transAxes, fontsize=fs)
    #     r -= dy
    #     ax.text(x0+.1, r, r'$\Delta\Omega/(2\pi)$: {0:.4f} MHz'.format(
    #         1.e-6*sim_result.linewidth_AC(ival)), transform=ax.transAxes, fontsize=fs)
    #     r -= dy
    #     mg_pe = plps['modal_gain'].get('PE', 0)
    #     mg_mb = plps['modal_gain'].get('MB', 0)
    #     mg_tot = plps['modal_gain'].get('Tot', 0)
    #     ax.text(x0+.1, r, r'Gain (PE, MB, Tot): ({0:.4f}, {1:.4f}, {2:.4f} ) $\mathrm{{(mW)}}^{{-1}}$'.format(
    #         mg_pe, mg_mb, mg_tot), transform=ax.transAxes, fontsize=fs)
    #     r -= dy
    
    d_ext = modeobj.get_mode_data()
    if len(d_ext):
        ax.text(x0, r, r'Extra data:', transform=ax.transAxes, fontsize=fs)
        r -= dy
        for (k, v) in d_ext.items():
            ax.text(x0, r, r'   {0}: {1}'.format(k, v),
                    transform=ax.transAxes, fontsize=fs)
            r -= dy


def plot_all_components(v_x, v_y, m_X, m_Y, v_plots, plps, sim_result, ival):
    decorator = plps['decorator']
    figsz = decorator.get_axes_property('figsize')
    ws = decorator.get_axes_property('subplots_wspace')
    hs = decorator.get_axes_property('subplots_hspace')
    lw = decorator.get_axes_property('linewidth')
    ec = decorator.get_axes_property('edgecolor')

    decorator.set_frame_drawer(sim_result._structure.wg_geom)

    field_type = plps['EM_AC']
    plt.rc('axes', linewidth=lw)
    plt.rc('axes', edgecolor=ec)


    rows = 3
    if plps['suppress_imimre']:
        rows = 2

    fig, axs = plt.subplots(rows, 3, figsize=figsz)
    axi = 0
    ax = axs[axi//3, axi%3]


    plot_mode_data(ax, plps, sim_result, ival, v_x, v_y)  # mode data summary
    axi += 1

    cc_cont = {FieldType.EM_E: component_t('Eabs'), FieldType.EM_H: component_t(
        'Habs'), FieldType.AC: component_t('uabs')}[field_type]
    cc_quiv = {FieldType.EM_E: component_t('Et'), FieldType.EM_H: component_t(
        'Ht'), FieldType.AC: component_t('ut')}[field_type]

    ax = axs[axi//3, axi%3]
    # plot_one_component_axes(ax, m_X, m_Y, v_plots, plps, cc)  # the intensity plot
    plot_one_component_axes_contour_and_quiver(ax, m_X, m_Y, v_plots, plps,
                                               cc_cont=cc_cont, cc_quiver=cc_quiv)  # the intensity plot
    axi += 1

    ax = axs[axi//3, axi%3]
    # plot_one_component_quiver(ax, m_X, m_Y, v_plots, plps, cc)  # the transverse vector plot
    plot_one_component_axes_contour_and_quiver(ax, m_X, m_Y, v_plots, plps,
                                               cc_quiver=cc_quiv)  # the intensity plot
    axi += 1

    for (flab, field) in v_plots.items():
        cc = component_t.make_comp(field_type, flab)

        if plps['suppress_imimre'] and not cc.is_dominant():
            continue
        if not cc.is_signed_field():
            continue
        ax = axs[axi//3, axi%3]
        plot_one_component_axes_contour_and_quiver(
            ax, m_X, m_Y, v_plots, plps, cc_cont=cc)  # the scalar plots
        axi += 1

    fig_fname = plot_filename(plps, ival)
    save_and_close_figure(fig,  fig_fname)




def plot_one_component(m_X, m_Y, v_fields, plps, sim_result, ival, cc, axis=None):

    decorator = plps['decorator']

    if axis is None:
        plt.clf()  # TODO: surely this is not needed

        plt.rc('axes', linewidth=decorator.get_axes_property('linewidth'))
        plt.rc('axes', edgecolor=decorator.get_axes_property('edgecolor'))

        fig, ax = plt.subplots(figsize=(12, 10))

    else:
        ax = axis

    field_type = plps['EM_AC']
    cc_cont = {FieldType.EM_E: component_t('Eabs'), FieldType.EM_H: component_t(
        'Habs'), FieldType.AC: component_t('uabs')}[field_type]
    cc_quiv = {FieldType.EM_E: component_t('Et'), FieldType.EM_H: component_t(
        'Ht'), FieldType.AC: component_t('ut')}[field_type]

    # if cc.is_transverse():
    #  plot_one_component_quiver(ax, m_X, m_Y, v_fields, plps, cc)
    # else:
    #  plot_one_component_axes(ax, m_X, m_Y, v_fields, plps, cc)

    plot_one_component_axes_contour_and_quiver(ax, m_X, m_Y, v_fields, plps,
                                               cc_cont=cc_cont, cc_quiver=cc_quiv)

    if axis is None:  # If user passed in the axis, they can look after saving.
        fig_fname = plot_filename(plps, ival, cc._user_code)
        save_and_close_figure(fig, fig_fname)


# deprecated spelling.
def plt_mode_fields(sim_result, ivals=None, n_points=501, quiver_points=50,
                    xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None,
                    field_type='EM_E', num_ticks=None, colorbar=True, contours=False, contour_lst=None,
                    stress_fields=False, pdf_png='png',
                    prefix='', suffix='', ticks=True, comps=[], decorator=None,
                    suppress_imimre=True,
                    modal_gains_PE=None,
                    modal_gains_MB=None,
                    modal_gains=None):

    print('Warning: "plt_mode_fields" is deprecated, use "plot_mode_fields"')
    plot_mode_fields(sim_result, ivals, n_points, quiver_points,
                     xlim_min, xlim_max, ylim_min, ylim_max,
                     field_type, num_ticks, colorbar, contours, contour_lst, stress_fields, pdf_png,
                     prefix, suffix, ticks, comps, decorator,
                     suppress_imimre, modal_gains_PE, modal_gains_MB, modal_gains)

#### Standard plotting of spectra #############################################


def plot_mode_fields(sim_result, ivals=None, n_points=501, quiver_points=30,
                     xlim_min=0, xlim_max=0, ylim_min=0, ylim_max=0,
                     field_type='EM_E',
                     num_ticks=None, colorbar=True, contours=False, contour_lst=None,
                     stress_fields=False, pdf_png='png',
                     prefix='', suffix='', ticks=True, comps=[], decorator=None,
                     suppress_imimre=True,
                     modal_gains_PE=None,
                     modal_gains_MB=None,
                     modal_gains=None):
    """ Plot E or H fields of EM mode, or the AC modes displacement fields.

        Args:
            sim_result : A ``Struct`` instance that has had calc_modes calculated

        Keyword Args:
            ivals  (list): mode numbers of modes you wish to plot

            n_points  (int): The number of points across unitcell to
                interpolate the field onto

            xlim_min  (float): Limit plotted xrange to xlim_min:(1-xlim_max) of unitcell

            xlim_max  (float): Limit plotted xrange to xlim_min:(1-xlim_max) of unitcell

            ylim_min  (float): Limit plotted yrange to ylim_min:(1-ylim_max) of unitcell

            ylim_max  (float): Limit plotted yrange to ylim_min:(1-ylim_max) of unitcell

            field_type  (str): Either 'EM' or 'AC' modes

            num_ticks  (int): Number of tick marks

            contours  (bool): Controls contours being overlaid on fields

            contour_lst  (list): Specify contour values

            stress_fields  (bool): Calculate acoustic stress fields

            pdf_png  (str): File type to save, either 'png' or 'pdf'

            prefix  (str): Add a string to start of file name

            suffix  (str): Add a string to end of file name.

            modal_gains (float array): Pre-calculated gain for each acoustic mode given chosen optical fields.
    """

    # TODO: ft_emac  ,a new var to distinguish field_type string and enum
    if sim_result.is_AC():
        field_type = FieldType.AC
    else:
        try:
            # TODO:ugly that this changes from string to enum
            field_type = FieldType.from_str(field_type)
        except:
            raise ValueError("field_type must be either 'AC', 'EM_E' or 'EM_H'.")

    # Do this elsewhere by function on Mode

    # Calculate the magnetic field from the electric field
    # Move this into plot_mode_H
    fm = sim_result.fem_mesh

    if field_type == FieldType.EM_H: sim_result.make_H_fields()
 

    # assemble desired list of eigenmodes to plot
    if not ivals is None:
        ival_range = ivals
    else:
        ival_range = range(len(sim_result.nu_AC_all()))

    mode_helper = sim_result.get_mode_helper()
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
#                  'colorbar':colorbar, 'contours':contours, 'contour_lst':contour_lst, 'field_type':field_type,
#                  'prefix': prefix, 'suffix': suffix, 'pdf_png': pdf_png,
#                  'stress_fields':stress_fields, 'modal_gain':modal_gain, 'decorator': decorator,
#                  'suppress_imimre':suppress_imimre,
#              # 'n_pts_x': n_pts_x, 'n_pts_y': n_pts_y,
#               'quiver_points': quiver_points }
#

    if not prefix: prefix=numbat.NumBATApp().outprefix()

    mode_helper.set_plot_params(xlim_min=xlim_min, xlim_max=xlim_max, ylim_min=ylim_min, ylim_max=ylim_max,
                                field_type=field_type,
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

    v_modes = sim_result.get_all_modes()

    modetype = 'em'
    if field_type == FieldType.AC:
        modetype = 'acoustic'

    if len(ival_range) > 1:
        print(f'Plotting {modetype} modes m={ival_range[0]} to {ival_range[-1]}.')
    else:
        print(f'Plotting {modetype} mode m={ival_range[0]}.')

    # TODO: mode drawing and especially saving is very slow. 
    for ival in ival_range:
        v_modes[ival].plot_mode(comps, field_type)

