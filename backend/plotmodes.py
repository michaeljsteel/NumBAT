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

from collections.abc import Iterable
import copy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors

import numbat
from nbtypes import FieldTag, FieldType, SI_THz, SI_GHz, SI_um, twopi
from plottools import save_and_close_figure


def modeplot_filename_2D(field_code, plps, mode_index, label='', cut=''):
    nbapp = numbat.NumBATApp()

    pref=plps['prefix']
    suf = plps['suffix']
    comp = field_code.as_str()

    fullpref = str(nbapp.outdir_fields_path(prefix=pref))
    filestart = f'{fullpref}/{comp}_mode_{mode_index:02d}'

    if cut:
        filestart = filestart + f'_{cut}cut'
    if suf:
        filestart = filestart + f'_{suf}'

    if label: filestart += '_'+label

    return filestart + nbapp.plotfile_ext()

# TODO: allow a suffix to label multiple cut planes
def modeplot_filename_1D(field_code, plps, mode_index, cut, label=''):

    return modeplot_filename_2D(field_code, plps, mode_index, label, cut)
    # nbapp = numbat.NumBATApp()

    # pref=plps['prefix']
    # suf = plps['suffix']
    # comp = field_code.as_str()

    # fullpref = str(nbapp.outdir_fields_path(prefix=pref))
    # filestart = f'{fullpref}/{comp}_field_{mode_index:02d}{suf}_{cut}cut'

    # if label: filestart += '_'+label

    # return filestart + nbapp.plotfile_ext()

def field_type_to_intensity_code(ft):
    return {FieldType.EM_E: FieldTag('Eabs'), FieldType.EM_H: FieldTag(
        'Habs'), FieldType.AC: FieldTag('uabs')}[ft]

def field_type_to_vector_code(ft):
    return  {FieldType.EM_E: FieldTag('Et'), FieldType.EM_H: FieldTag(
        'Ht'), FieldType.AC: FieldTag('ut')}[ft]


class TidyAxes:

    def __init__(self, nax=1, props={}):

        self._nax = nax
        self._set_defaults_for_num_axes(nax)
        self._props.update(props)

    def update_property(self, k, v):
        self._props[k]=v

    def _set_defaults_for_num_axes(self, nax):
        props = {}

        if nax in (4,6):
            props['ax_label_fs'] = 10
            props['ax_label_xpad'] = ''
            props['ax_label_ypad'] = ''
            props['ax_ticklabel_fs'] = 10
            props['ax_tickwidth'] = .25
            props['ax_linewidth'] = 1
        elif nax == 2:
            props['ax_label_fs'] = 12
            props['ax_label_pad'] = 5
            props['ax_ticklabel_fs'] =  10
            props['ax_tickwidth'] = 1
            props['ax_linewidth'] = 1
        else:
            props['ax_label_fs'] =  18
            props['ax_label_pad'] = 5
            props['ax_ticklabel_fs'] =  16
            props['ax_tickwidth'] =  1
            props['ax_linewidth'] =  1

        props['aspect'] = 0.0

        props['axes_color'] = 'gray'

        props['cb_linewidth'] = 1
        props['cb_label_fs'] = 10
        props['cb_ticklabel_fs'] = 8
        props['cb_tickwidth'] = 0.25

        props['cb_edgecolor'] = 'gray'

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

            # Shape
            if pr['aspect'] >0 : ax.set_aspect(pr['aspect'])

            # Ticks
            ax.tick_params(labelsize=pr['ax_ticklabel_fs'],
                           width=pr['ax_tickwidth'])

            # Axis labels
            ax.xaxis.label.set_size(pr['ax_label_fs'])
            ax.yaxis.label.set_size(pr['ax_label_fs'])

            xpad = self._props['ax_label_xpad']
            ypad = self._props['ax_label_ypad']
            if xpad:
                xlab = ax.xaxis.get_label_text()
                ax.set_xlabel(xlab, labelpad=xpad)
            if ypad:
                ylab = ax.yaxis.get_label_text()
                ax.set_ylabel(ylab, labelpad=ypad)

            # Axes visibility
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(pr['ax_linewidth'])
                ax.spines[axis].set_color(pr['axes_color'])

    def hide_axes(self, ax):
        ax.set_axis_off()

    def apply_to_cbars(self, cbs):
        if not isinstance(cbs, Iterable):
            cbs = (cbs,)

        pr = self._props
        for cb in cbs:

            lab = cb.ax.get_ylabel()  # don't seem to be able to set label size except by setting label again
            cb.set_label(lab, size=pr['cb_label_fs'])

            cb.outline.set_linewidth(pr['cb_linewidth'])
            cb.outline.set_color(pr['cb_edgecolor'])


            cb.ax.tick_params(labelsize=pr['cb_ticklabel_fs'],
                           width=pr['cb_tickwidth'])

            # Not sure how to do these

            #if pr['cb_shrink'] >0 :
            #    cb.set_shrink(pr['cb_shrink'])
            #    cb.set_pad(pr['cb_pad'])


class Decorator(object):

    def __init__(self, multi=False):
        mp_base_fs = 18
        sp_base_fs = 24

        self._multi_props = {'figsize': (10, 8),'subplots_hspace': .2, 'subplots_wspace': .6,
                             'base_fs': mp_base_fs, 'title_fs': mp_base_fs-2, 'subtitle_fs': mp_base_fs-5,
                             'title_pad': 10, 'subtitle_pad': 10,                                     'ax_label_fs': mp_base_fs-10, 'ax_label_pad': 20, 'ax_tick_fs': mp_base_fs-10,
                             'data_label_fs': mp_base_fs-8,
                             'cbar_tick_fs': mp_base_fs-12,
                             'ax_linewidth': '.75', # used in plotgain only
                             'ax_edgecolor': 'gray',
                             'cbar_size': '5%', 'cbar_pad': '6%'}

        self._single_props = {'figsize': (10, 8),'subplots_hspace': .2, 'subplots_wspace': .2,
                              'base_fs': sp_base_fs, 'title_fs': sp_base_fs-2,
                              'subtitle_fs': sp_base_fs-3,
                              'title_pad': 20, 'subtitle_pad': 20,
                              'ax_label_fs': 20, 'ax_label_pad': 20, 'ax_tick_fs': 20,
                              'data_label_fs': sp_base_fs-8,
                              'cbar_tick_fs': 20,
                              'ax_linewidth': '.75',
                              'ax_edgecolor': 'gray',
                              'cbar_size': '5%', 'cbar_pad': '2%'}

        #self._props = self._multi_props if multi else self._single_props
        self._props = copy.deepcopy(self._multi_props) if multi else copy.deepcopy(
            self._single_props)
        self._is_single = not multi
        self._cmap_limits = {}
        self._title = ''
        self._waveguide = None

    def __getitem__(self, kw):
        return self._props[kw]

    def __setitem__(self, kw, val):
        self._props[kw]=val

    def set_waveguide(self, wg):
        self._waveguide = wg

    def add_waveguide_border(self, ax):
        plprefs = numbat.NumBATPlotPrefs()
        styles = {}

        styles['edgecolor'] = plprefs.WG_FRAME_EDGE_COLOR
        if self._is_single:
            styles['linewidth'] = plprefs.WG_FRAME_LINEWIDTH_WHOLEFIG
        else:
            styles['linewidth'] = plprefs.WG_FRAME_LINEWIDTH_SUBFIG


        #if self._waveguide is not None:
        self._waveguide.draw_mpl_frame(ax, styles)

    def _get_props(self):
        return self._props

    def set_for_single(self):
        self._is_single = True
        self._props.update(self._single_props)

    def set_for_multi(self):
        self._is_single = False
        self._props.update(self._multi_props)

    def set_cmap_limits(self, d):
        '''Specify the lower and upper contour plot limits for a field component plot.

          :param dict d: Dictionary mapping component ('x','y','z') to (zlo, zhi) tuple for contour plots.
          '''
        self._cmap_limits.update(d)

    def get_cmap_limits(self, comp):
        return self._cmap_limits.get(comp)


    def get_property(self, lab):
        ans=''
        try:
            ans = self._props[lab]
        except Exception:
            print(f'Warning: unknown fontsize label "{lab}" in Decorator::get_property()')
        return ans

    def set_property(self, label, prop):
        '''Add or override an axes property for a single plot corresponding to the given label.'''
        self._props[label] = prop

    def set_title(self, t): self._title = t

    def add_title(self, ax): # used only in plotgain , remove?
        if self._title:
            ax.set_title(self._title)

    def extra_axes_commands(self, ax):
        '''Function to make additions to a standard plot.

          This function is called after all other plotting operations are performed.
          The default implementation does nothing. By subclassing :Decorator: and implementing this function,
          users may add extra features to a plot.
          '''
        pass

class PlotParams:

    def __init__(self):
        self._d_pp = {}

    def update(self, d):
        self._d_pp.update(d)
        if self._d_pp['decorator'] is None:
            self._d_pp['decorator'] = Decorator()

    def __getitem__(self, k):
        v= self._d_pp[k]
        return v

    def __setitem__(self, k,v):
        self._d_pp[k] = v

    def get(self, k, dflt, override_None=False):
        val = self._d_pp.get(k, dflt)
        if val is None and override_None:
            val = dflt
        return val

    def __str__(self):
        return str(self._d_pp)

class PlotParams2D(PlotParams):
    def __init__(self):
        super().__init__()

        self.update({'xlim_min': None, 'xlim_max': None,
                     'ylim_min': None, 'ylim_max': None,
                      'aspect': 1.0,
                      'ticks': True, 'num_ticks': None,
                      'colorbar': True, 'contours': False, 'contour_lst': None,
                     # 'EM_AC': FieldType.EM_E,
                      'title': True, 'frame': True,
                      'hide_vector_field': False,
                      'prefix': '', 'suffix': '',
                      'decorator': Decorator(),
                      'suppress_imimre': True,
                      'quiver_points': 30
                      })

class PlotParams1D(PlotParams):

    def __init__(self):
        super().__init__()

        self.update({#'xlim_min': 0, 'xlim_max': 0, 'ylim_min': 0, 'ylim_max': 0,
                      'aspect': 1.0,
                      'ticks': True, 'num_ticks': None,
                      #'colorbar': True, 'contours': False, 'contour_lst': None,
                     # 'EM_AC': FieldType.EM_E,
                      #'hide_vector_field': False,
                      'prefix': '', 'suffix': '',
                      'decorator': Decorator(),
                      'suppress_imimre': True,
                      #'quiver_points': 30
                      })




def plot_set_ticks(ax, plps):
    if plps['ticks']:
        ax.set_xlabel('$x$ [μm]')
        ax.set_ylabel('$y$ [μm]')
    else:
        ax.set_xticks([])
        ax.set_yticks([])

#UNUSERD DELETE
# def plot_set_axes_style(ax, plps, decorator):
#     if decorator.get_property('axes.linewidth') is not None:
#         plt.setp(ax.spines.values(),
#                  linewidth=decorator.get_property('axes.linewidth'))


def plot_set_title(ax, comp_label, plps, decorator):
    if plps.get('add_title', True):
        ax.set_title(comp_label, fontsize=decorator.get_property('subtitle_fs'),
                     pad=decorator.get_property('subtitle_pad'))


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


def add_contour_plot(fig, ax, d_xy, c_field, ftag_cont, plps, decorator):

    logamp = plps['logamp']

    cmap_signed = numbat.NumBATPlotPrefs().cmap_field_signed(ftag_cont)
    cmap_unsigned = numbat.NumBATPlotPrefs().cmap_field_unsigned(ftag_cont)

    cont_signed = ftag_cont.is_signed_field()
    #cmap = cmap_signed if (cont_signed and not logamp) else cmap_unsigned
    cmap = cmap_signed if cont_signed else cmap_unsigned



    v_x, v_y, m_X, m_Y = list(d_xy.values())

    # Extract and tidy the field
    # transpose because arrays
    cont_field = (c_field).copy().T

    if ftag_cont.is_abs():
        # TODO: cleanup: plot |u| as |u|^2
        cont_field = np.abs(cont_field)**2

    logmax=5  # show up to 5 orders of magnitude contrast
    if logamp:
        #cont_field = np.log10(np.abs(cont_field))
        lcf = np.log10(np.abs(cont_field)) + logmax
        lcf = np.where(lcf>0, lcf, 0)
        cont_field = np.sign(cont_field) * lcf

    # if the data is all noise, just plot zeros
    plot_threshold = 1e-8
    if np.max(np.abs(cont_field[~np.isnan(cont_field)])) < plot_threshold:
        cont_field = np.zeros(np.shape(cont_field))

    # set imshow plot (and tick) range to match the input x and y domain
    if True or plps['ticks']:
        extents = np.array([v_x[0], v_x[-1], v_y[0], v_y[-1]])
    else:
        extents = None


    interp = None
    # interp='bilinear';

    d_kw = {'origin':'lower', 'extent':extents, 'interpolation':interp, 'cmap':cmap}

    if decorator.get_cmap_limits(ftag_cont._xyz) is not None:
        (act_zlo, act_zhi) = decorator.get_cmap_limits(ftag_cont._xyz)
        tsnorm = mplcolors.TwoSlopeNorm(
            vmin=act_zlo, vmax=act_zhi, vcenter=(act_zlo+act_zhi)/2)

        d_kw['norm'] = tsnorm
    elif logamp:
        vma = logmax
        if cont_signed:
            vmi = -logmax
        else:
            vmi=0

        d_kw['vmin'] = vmi
        d_kw['vmax'] = vma

    else:
        act_zlo = np.nanmin(cont_field)
        act_zhi = np.nanmax(cont_field)
        vma = max(abs(act_zlo), abs(act_zhi))
        vmi = -vma if cont_signed else 0.0

        d_kw['vmin'] = vmi
        d_kw['vmax'] = vma

    im_co = ax.imshow(cont_field, **d_kw)

    do_cbar = plps['colorbar']
    do_contours = plps['contours']

    cbar = None
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
        #divider = make_axes_locatable(ax)

        #cax = divider.append_axes("right", size=decorator.get_property('cbar_size'),
         #                           pad=decorator.get_property('cbar_pad'))
        #cbar = plt.colorbar(im_co, cax=cax)
        #print('pad',decorator.get_property('cbar_pad'))
        #cbpad = 0.05
        cax = ax.inset_axes([1.04, .2, 0.05, 0.6])
        cbar=fig.colorbar(im_co, cax=cax)


        if cbarticks is not None:
            cbar.set_ticks(cbarticks)
            cbarlabels = [f'{t:.2f}' for t in cbarticks]
            cbar.set_ticklabels(cbarlabels)
        cbar.ax.tick_params(labelsize=decorator.get_property('cbar_tick_fs'))

    if do_contours:
        if np.max(np.abs(cont_field[~np.isnan(cont_field)])) > plot_threshold:
            CS2 = ax.contour(m_X, m_Y, cont_field, levels=cbarticks, #colors=mycolors[::-1],
                                linewidths=(1.5,))
            if do_cbar:
                cbar.add_lines(CS2)

    return im_co, cbar

def _make_quiver_slices(v_x, v_y, plps, deftrim):
    quiver_points = plps.get('quiver_points', 20)
    n_pts_x, n_pts_y = len(v_x), len(v_y)

    # grid points to skip for each arrow
    # get y points spaced as closely as possible to x points for a square grid
    quiver_points_x = quiver_points
    quiver_points_y = round(quiver_points*(v_y[-1]-v_y[0])/(v_x[-1]-v_x[0]))

    #Scale so they look good for any aspect ratio, by adding points on the elongated axis
    asp = abs(plps['aspect'])
    if asp<1:
        quiver_points_x = round(quiver_points_x/asp)
    if asp>1:
        quiver_points_y = round(quiver_points_y*asp)

    # To keep things safe, have at least 3 points in each direction
    # If get_quiver_skip_range was safer, this wouldn't be necessary
    quiver_points_x = max(3, quiver_points_x)
    quiver_points_y = max(3, quiver_points_y)

    # We could define quiver_points vanilla as the number along the _visible_ x range accounting for trim
    #  or along the whole x range before trimming.
    # We choose the latter so trim doesn't change the spacing

    xlmi = plps.get('xlim_min', deftrim, override_None=True)
    xlma = plps.get('xlim_max', deftrim, override_None=True)
    ylmi = plps.get('ylim_min', deftrim, override_None=True)
    ylma = plps.get('ylim_max', deftrim, override_None=True)

    quiver_skip_x = int(round(n_pts_x/quiver_points_x * (1-xlmi-xlma)))
    quiver_skip_y = int(round(n_pts_y/quiver_points_y * (1-ylmi-ylma)))

    # getting a nice symmetric pattern of points to do quivers centred around the middle
    qslice_x = get_quiver_skip_range(n_pts_x, quiver_skip_x)
    qslice_y = get_quiver_skip_range(n_pts_y, quiver_skip_y)

    return qslice_x, qslice_y

def add_quiver_plot(fig, ax, d_xy, v_fields, cc, plps, decorator, do_cont):


    # give a little space around elastic profiles
    deftrim = -0.1 if cc._F=='u' else -0.05

    v_x, v_y, m_x, m_y = d_xy.values()

    # shapes:
    # m_X and other matrices are stored with y points as rows but upside down
    # mesh_grid_objects:
    # m_X =  [[x0, x1, x2, ...]
    #         [x0, x1, x2, ...]
    #         [x0, x1, x2, ...]
    # m_Y =  [[y0, y0, y0, ...]
    #         [y1, y1, y1, ...]
    #         [y2, y2, y2, ...]



    # TODO: spend a day getting truly perfect spacings of grid poitns
    #         and trying to figure out the scaling and arrow head settings once and for all

    qslice_x, qslice_y = _make_quiver_slices(v_x, v_y, plps, deftrim)

    #print('pts', n_pts_x, n_pts_y, m_x.shape, m_y.shape, v_fields['Fxr'].shape, qslice_x, qslice_y)
    v_x_q = v_x[qslice_x]
    v_y_q = v_y[qslice_y]
    m_x_q = (m_x.T)[qslice_x[:, np.newaxis], qslice_y]  # TODO: track down why m_x/y need .T but fields don't
    m_y_q = (m_y.T)[qslice_x[:, np.newaxis], qslice_y]

    # TODO: why no transpose on these fields?
    m_ReEx_q = v_fields['Fxr'][qslice_x[:, np.newaxis], qslice_y]
    m_ReEy_q = v_fields['Fyr'][qslice_x[:, np.newaxis], qslice_y]
    #m_ImEx_q = v_fields['Fxi'][qslice_x[:, np.newaxis], qslice_y]
    #m_ImEy_q = v_fields['Fyi'][qslice_x[:, np.newaxis], qslice_y]


    vecarrow_col = numbat.NumBATPlotPrefs().vector_field_arrow_color(cc)
    vecarrow_scale = numbat.NumBATPlotPrefs().vector_field_arrow_scale()
    vecarrow_lw = numbat.NumBATPlotPrefs().vector_field_arrow_linewidth()
    vecarrow_hw = numbat.NumBATPlotPrefs().vecplot_arrow_headwidth



    # Ignore all imaginary values. If there are significant imag values,
    # then instantaneous vector plots don't make much sense anyway
    d_quiv_kw = {#'linewidths': (vecarrow_lw,),
                 'width': vecarrow_lw,
                 'edgecolor': vecarrow_col,
                 'facecolor':vecarrow_col,
                 'pivot':'mid',
                 #'headlength':3,
                 'scale':vecarrow_scale,
                 'headwidth': vecarrow_hw}

    if do_cont:  # no colours in the quiver
        d_quiv_kw['color'] = 'gray'
        ax.quiver(m_x_q, m_y_q, m_ReEx_q, m_ReEy_q,  ** d_quiv_kw)
    else:
        m_arrcolour = np.sqrt(m_ReEx_q*m_ReEx_q + m_ReEy_q*m_ReEy_q)
        ax.quiver(m_x_q, m_y_q, m_ReEx_q, m_ReEy_q, m_arrcolour, ** d_quiv_kw)


    if not do_cont:
        # ax.set_aspect('equal')  # UNDO

        # this step is needed because quiver doesn't seem
        # to use its input x and y vectors to set range limits
        # clean this up so as to avoid seemingly circular calls following
        ax.set_xlim(v_x_q[0], v_x_q[-1])
        ax.set_ylim(v_y_q[0], v_y_q[-1])


def plot_contour_and_quiver(fig, ax, d_xy, v_fields, plps, ftag_scalar=None, ftag_vector=None,
                            is_single_plot=False):

    is_AC = (ftag_scalar and ftag_scalar.is_AC()) or (ftag_vector and ftag_vector.is_AC())

    do_cont = ftag_scalar is not None
    do_quiv = ftag_vector is not None

    decorator = plps['decorator']

    cbar = None
    if do_cont:
        im_co, cbar = add_contour_plot(fig, ax, d_xy, v_fields[ftag_scalar.component_as_F()],
                                       ftag_scalar, plps, decorator)

    if do_quiv:
        add_quiver_plot(fig, ax, d_xy, v_fields,
                        ftag_vector, plps, decorator, do_cont)

    # Adjustments to the visible plot domain

    deftrim = -.1 if is_AC else -0.05
    xlmi = plps.get('xlim_min', deftrim, override_None=True)
    xlma = plps.get('xlim_max', deftrim, override_None=True)
    ylmi = plps.get('ylim_min', deftrim, override_None=True)
    ylma = plps.get('ylim_max', deftrim, override_None=True)

    if xlmi != 0 or xlma != 0:
        xmin, xmax = ax.get_xlim()
        width_x = xmax-xmin
        ax.set_xlim(xmin+xlmi*width_x, xmax-xlma*width_x)

    if ylmi != 0 or ylma != 0:
        ymin, ymax = ax.get_ylim()
        width_y = ymax-ymin
        ax.set_ylim(ymin+ylmi*width_y, ymax-ylma*width_y)


    labs = []
    if do_cont: labs.append(ftag_scalar.get_tex_plot_label())
    if do_quiv: labs.append(ftag_vector.get_tex_plot_label())
    comp_label = ', '.join(labs)

    #lw = decorator.get_property('linewidth')# TODO: Unused
    ec = decorator.get_property('ax_edgecolor')

    aspect = plps.get('aspect', 1.0)

    nax = 1 if is_single_plot else 6
    tidy = TidyAxes(nax=nax,
                    props={'axes_color':ec, 'cb_edgecolor':ec,
                        'ax_label_xpad':3, 'ax_label_ypad':1,
                        'aspect':aspect })


    tidy.apply_to_axes(ax)
    if cbar:
        tidy.apply_to_cbars(cbar)

    if not plps['frame']:
        tidy.hide_axes(ax)

    plot_set_ticks(ax, plps)
    #plot_set_axes_style(ax, plps, decorator)

    if plps['title']:
        plot_set_title(ax, comp_label, plps, decorator)

    decorator.add_waveguide_border(ax)
    decorator.extra_axes_commands(ax)





def write_mode_data(ax, plps, sim_result, mode_index, field_code):  # mode data summary

    decorator = plps['decorator']
    fs = decorator.get_property('data_label_fs')

    x0 = .05
    yhi = .99
    dy = .10
    y0 = yhi-dy

    ax.set_xlim(0, 1)
    ax.set_ylim(0, yhi)
    ax.set_aspect('equal')
    ax.axis('off')

    mode = sim_result.get_all_modes()[mode_index]
    mode.analyse_mode()


    # make a compact line writer
    def _write_line(tx, ty, s, tfs=fs):
        ax.text(tx, ty, s, transform=ax.transAxes, fontsize=tfs)

    [f_x, f_y, f_t, f_z] = mode.field_fracs()
    (r0x, r0y) = mode.center_of_mass()         # In units of um
    (wx, wy, w0) = mode.second_moment_widths() # In units of um


    _write_line(x0-.05, y0, f'Mode properties: m={mode_index}', fs+2); y0 -= dy

    lines=[]

    if sim_result.is_EM():
        lines.append( r'$\omega/(2\pi)$: ' + f'{sim_result.omega_EM/(twopi*SI_THz):.5f} THz')
        lines.append( r'$k$: ' + f'{sim_result.kz_EM(mode_index)/1.e6:.5f} ' + r'μm$^{{-1}}$')

        neff = sim_result.neff(mode_index)
        if sim_result.ngroup_EM_available():
            ng = sim_result.ngroup_EM(mode_index)
            lines.append( r'$\bar{{n}}$: ' + f'{neff:.6f},' + f'$n_g$: {ng:.6f}')
        else:
            lines.append( r'$\bar{{n}}$: ' + f'{neff:.6f}')
    else:
        q_AC = sim_result.q_AC
        nu_AC = np.real(sim_result.nu_AC(mode_index))
        vp = sim_result.vp_AC(mode_index)

        lines.append(
                   f'$q$: {q_AC * SI_um:.5f} ' + r'μm$^{{-1}}$, $\lambda:$ '+ f'{twopi/q_AC/SI_um:.5f} ' + r'μm') ;y0 -= dy

        lines.append( r'$q/2\pi$: ' + f'{q_AC*SI_um/twopi:.5f} '+ r'μm$^{{-1}}$')
        lines.append( r'$\Omega/(2\pi)$: ' + f'{nu_AC/SI_GHz:.5f} GHz' )
        if sim_result.vgroup_AC_available():
            vg = sim_result.vg_AC(mode_index)
            lines.append( f'$v_p$: {vp:.2f} m/s, $v_g$: {vg:.2f} m/s')
        else:
            lines.append( f'$v_p$: {vp:.2f} m/s')
    lines.append( f'$f_x:$ {f_x:.3f}, $f_y$: {f_y:.3f}')
    lines.append( f'$f_t:$ {f_t:.3f}, $f_z$: {f_z:.3f}')
    lines.append( r'$\mathbf{{r}}_0:$ '+ f'({r0x:.3f}, {r0y:.3f}) ' +r'μm')
    lines.append( f'$(w_x, w_y):$ ({wx:.3f}, {wy:.3f}) ' + r'μm')

    if field_code.is_EM_H():
        lines.append( r'$H$ field multiplied by $Z_0=376.7\, \Omega$')

    sc = sim_result.symmetry_classification(mode_index)
    if len(sc):
        lines.append( f'Sym: {sc}')

    # TODO: reactivate
    # if sim_result.is_AC() and sim_result.Q_method != QAcMethod.NotSet:
    #     if sim_result.Q_method == QAcMethod.Intrinsic:
    #         ax.text(x0, r, 'Losses (intrinsic):',
    #                 transform=ax.transAxes, fontsize=fs)
    #         y0 -= dy
    #     else:
    #         ax.text(x0, r, 'Losses (fixed Q):',
    #                 transform=ax.transAxes, fontsize=fs)
    #         y0 -= dy

    #     ax.text(x0+.1, r, r'$\alpha$: {0:.3e} s$^{{-1}}$, {1:.2f} cm$^{{-1}}$, {2:.2f} dB/cm'.format(
    #         sim_result.alpha_t_AC(mode_index),
    #         sim_result.alpha_s_AC(mode_index)/100.,
    #         sim_result.alpha_s_AC(mode_index)/100./(np.log(10.0)/10.0)
    #     ), transform=ax.transAxes, fontsize=fs)
    #     y0 -= dy

    #     ax.text(x0+.1, r, r'$Q_m$: {0:.2f}'.format(
    #         sim_result.Qmech_AC(mode_index)), transform=ax.transAxes, fontsize=fs)
    #     y0 -= dy
    #     ax.text(x0+.1, r, r'$\Delta\Omega/(2\pi)$: {0:.4f} MHz'.format(
    #         1.e-6*sim_result.linewidth_AC(mode_index)), transform=ax.transAxes, fontsize=fs)
    #     y0 -= dy
    #     mg_pe = plps['modal_gain'].get('PE', 0)
    #     mg_mb = plps['modal_gain'].get('MB', 0)
    #     mg_tot = plps['modal_gain'].get('Tot', 0)
    #     ax.text(x0+.1, r, r'Gain (PE, MB, Tot): ({0:.4f}, {1:.4f}, {2:.4f} ) $\mathrm{{(mW)}}^{{-1}}$'.format(
    #         mg_pe, mg_mb, mg_tot), transform=ax.transAxes, fontsize=fs)
    #     y0 -= dy

    d_ext = mode.get_mode_data()
    if len(d_ext):
        lines.append( r'Extra data:')
        for (k, v) in d_ext.items():
            lines.append( f'   {k}: {v}')

    for line in lines:
        _write_line(x0, y0, line)
        y0 -= dy

def plot_all_components(field_code, d_xy, v_plots, plps, sim_result, mode_index):
    decorator = plps['decorator']
    figsz = decorator['figsize']
    ws = decorator['subplots_wspace']
    hs = decorator['subplots_hspace']

    decorator.set_waveguide(sim_result._structure.wg_geom)

    hide_minors = plps['suppress_imimre']
    rows = 2 if hide_minors else 3

    fig, axs = plt.subplots(rows, 3, figsize=figsz)
    fig.subplots_adjust(hspace=hs, wspace=ws*1.5)

    axs = axs.flat
    axi = 0

    ax = axs[axi]; axi += 1
    write_mode_data(ax, plps, sim_result, mode_index, field_code)  # mode data summary

    ft = field_code.as_field_type()

    ftag_scal = FieldTag.make_from_field_and_component(ft, 'a')
    ftag_transvec = FieldTag.make_from_field_and_component(ft, 't')

    if plps['hide_vector_field']:
        ftag_transjoint = None
    else:
        ftag_transjoint = ftag_transvec

    ax = axs[axi]; axi += 1  # the whole field (top row middle)
    plot_contour_and_quiver(fig, ax, d_xy, v_plots, plps,
                            ftag_scalar=ftag_scal, ftag_vector=ftag_transjoint)

    ax = axs[axi]; axi += 1  # the intensity (top row right)
    plot_contour_and_quiver(fig, ax, d_xy, v_plots, plps, ftag_vector=ftag_transvec)

    for flab in v_plots.keys():
        cc = FieldTag.make_from_field_and_Fcode(ft, flab)

        if not cc.is_signed_field(): continue  # already done vector and energy plots

        # don't plot what will likely be pure zero fields
        if (hide_minors and not cc.is_minor()): continue

        ax = axs[axi]; axi += 1
        plot_contour_and_quiver(fig, ax, d_xy, v_plots, plps, ftag_scalar=cc)  # the scalar plots


    fig_fname = modeplot_filename_2D(field_code, plps, mode_index)
    save_and_close_figure(fig, fig_fname)




def plot_one_component(field_code, d_xy, v_fields, plps, mode_index, cc, axis=None):
    decorator = plps['decorator']
    if axis is None:
        figsz = decorator['figsize']
        fig, ax = plt.subplots(figsize=figsz)
    else:
        ax = axis

    if cc.is_transverse():
        ftag_transvec = cc
        ft = field_code.as_field_type()
        ftag_scal = FieldTag.make_from_field_and_component(ft, 'a')
    else:
        ftag_transvec = None
        ftag_scal = cc

    plot_contour_and_quiver(fig, ax, d_xy, v_fields, plps,
                            ftag_scalar=ftag_scal, ftag_vector=ftag_transvec, is_single_plot=True)

    if axis is None:  # If user passed in the axis, they can look after saving.
        fig_fname = modeplot_filename_2D(field_code, plps, mode_index, cc._user_code)
        save_and_close_figure(fig, fig_fname)




