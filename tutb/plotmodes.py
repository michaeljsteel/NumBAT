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

from collections.abc import Iterable

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numbat
from numbattools import save_and_close_figure
from nbtypes import component_t, FieldType, SI_THz, SI_GHz, SI_um, twopi


def modeplot_filename(plps, ival, label=''):
    nbapp = numbat.NumBATApp()
    fullpref = str(nbapp.path_fields())

    comp = plps['EM_AC'].name
    suf = plps['suffix']
    filestart = f'{fullpref}/{comp}_field_{ival:02d}{suf}'

    if label: filestart += '_'+label

    return filestart + nbapp.plotfile_ext()

def field_type_to_intensity_code(ft):
    return {FieldType.EM_E: component_t('Eabs'), FieldType.EM_H: component_t(
        'Habs'), FieldType.AC: component_t('uabs')}[ft]

def field_type_to_vector_code(ft):
    return  {FieldType.EM_E: component_t('Et'), FieldType.EM_H: component_t(
        'Ht'), FieldType.AC: component_t('ut')}[ft]


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
            props['ax_linewthidth'] = 1
        else:
            props['ax_label_fs'] =  12
            props['ax_label_pad'] = 5
            props['ax_ticklabel_fs'] =  10
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

            if pr['aspect'] >0 : ax.set_aspect(pr['aspect'])

            ax.tick_params(labelsize=pr['ax_ticklabel_fs'],
                           width=pr['ax_tickwidth'])

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


            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(pr['ax_linewidth'])
                ax.spines[axis].set_color(pr['axes_color'])


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

        self._multi_props = {'figsize': (10, 8),'subplots_hspace': .2, 'subplots_wspace': .4,
                             'title_fs': mp_base_fs-2, 'subtitle_fs': mp_base_fs-5,
                             'title_pad': 10, 'subtitle_pad': 10,                                     'ax_label_fs': mp_base_fs-10, 'ax_label_pad': 20, 'ax_tick_fs': mp_base_fs-10,
                             'data_label_fs': mp_base_fs-8,
                             'cbar_tick_fs': mp_base_fs-12,
                             'linewidth': '.75', 'edgecolor': 'gray',
                             'cbar_size': '5%', 'cbar_pad': '6%'}

        self._single_props = {'figsize': (10, 8),'subplots_hspace': .2, 'subplots_wspace': .2,
                              'title_fs': sp_base_fs-2,  'subtitle_fs': sp_base_fs-3,
                              'title_pad': 20, 'subtitle_pad': 20,
                              'ax_label_fs': 20, 'ax_label_pad': 20, 'ax_tick_fs': 20,
                              'data_label_fs': sp_base_fs-8,
                              'cbar_tick_fs': 20,
                              'linewidth': '.75', 'edgecolor': 'gray',
                              'cbar_size': '5%', 'cbar_pad': '2%'}

        self._props = self._multi_props if multi else self._single_props
        #self._is_single = True
        self._cmap_limits = {}
        self._title = ''
        self._frame_drawer = None

    def set_frame_drawer(self, fd):
        self._frame_drawer = fd

    def add_frame(self, ax):
        if self._frame_drawer is not None:
            self._frame_drawer.draw_mpl_frame(ax)

    #def _fontsizes(self):
    #    if self._is_single:
    #        return self._single_fontsizes
    ##    else:
      #      return self._multi_fontsizes

    def _get_props(self):
        return self._props

    def set_for_single(self):
        self._props.update(self._single_props)

    def set_for_multi(self):
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

    # def is_single_plot(self):
    #     '''Returns True if this Decorator is for a single axes plot such as a spectrum or spatial map of a single field component.
    #        '''
    #     return self._is_single

    def set_property(self, label, prop):
        '''Add or override an axes property for a single plot corresponding to the given label.'''
        self._props[label] = prop

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



def plot_set_ticks(ax, plps):
    if plps['ticks']:
        ax.set_xlabel('$x$ [μm]')
        ax.set_ylabel('$y$ [μm]')
    else:
        ax.set_xticks([])
        ax.set_yticks([])


def plot_set_axes_style(ax, plps, decorator):
    if decorator.get_property('axes.linewidth') is not None:
        plt.setp(ax.spines.values(),
                 linewidth=decorator.get_property('axes.linewidth'))


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


def add_contour_plot(ax, d_xy, c_field, cc_cont, plps, decorator):

    cmap_signed = 'seismic'   # This should be a plot_param
    cmap_unsigned = 'OrRd'

    cont_signed = cc_cont.is_signed_field()
    cmap = cmap_signed if cont_signed else cmap_unsigned


    v_x, v_y, m_X, m_Y = list(d_xy.values())

    # Extract and tidy the field
    # transpose because arrays
    cont_field = (c_field).copy().T

    if cc_cont.is_abs():
        # TODO: cleanup: plot |u| as |u|^2
        cont_field = np.abs(cont_field)**2

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

    if decorator.get_cmap_limits(cc_cont._xyz) != None:
        (act_zlo, act_zhi) = decorator.get_cmap_limits(cc_cont._xyz)
        tsnorm = mplcolors.TwoSlopeNorm(
            vmin=act_zlo, vmax=act_zhi, vcenter=(act_zlo+act_zhi)/2)

        d_kw['norm'] = tsnorm
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
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=decorator.get_property('cbar_size'),
                                    pad=decorator.get_property('cbar_pad'))
        cbar = plt.colorbar(im_co, cax=cax)
        if not cbarticks is None:
            cbar.set_ticks(cbarticks)
            cbarlabels = ['%.2f' % t for t in cbarticks]
            cbar.set_ticklabels(cbarlabels)
        cbar.ax.tick_params(labelsize=decorator.get_property('cbar_tick_fs'))

    if do_contours:
        if np.max(np.abs(cont_field[~np.isnan(cont_field)])) > plot_threshold:
            CS2 = ax.contour(m_X, m_Y, cont_field, levels=cbarticks, #colors=mycolors[::-1],
                                linewidths=(1.5,))
            if do_cbar:
                cbar.add_lines(CS2)

    return im_co, cbar

def add_quiver_plot(ax, d_xy, v_fields, cc, plps, decorator, do_cont):

    quiver_points = plps.get('quiver_points', 20)

    # give a little space around elastic profiles
    deftrim = -.05 if cc._F=='u' else 0.0

    xlmi = plps.get('xlim_min', deftrim)
    xlma = plps.get('xlim_max', deftrim)
    ylmi = plps.get('ylim_min', deftrim)
    ylma = plps.get('ylim_max', deftrim)

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

    n_pts_x, n_pts_y = len(v_x), len(v_y)


    # grid points to skip for each arrow
    # this could probably be chosen nicer. make equally spaced on aspect=1?
    quiver_skip_x = int(round(n_pts_x/quiver_points * (1-xlmi-xlma)))
    quiver_skip_y = int(round(n_pts_y/quiver_points * (1-ylmi-ylma)))


    # getting a nice symmetric pattern of points to do quivers centred around the middle
    qrange_x = get_quiver_skip_range(n_pts_x, quiver_skip_x)
    qrange_y = get_quiver_skip_range(n_pts_y, quiver_skip_y)

    #print('pts', n_pts_x, n_pts_y, m_x.shape, m_y.shape, v_fields['Fxr'].shape, qrange_x, #qrange_y)
    v_x_q = v_x[qrange_x]
    v_y_q = v_y[qrange_y]
    m_x_q = (m_x.T)[qrange_x[:, np.newaxis], qrange_y]  # TODO: track down why m_x/y need .T but fields don't
    m_y_q = (m_y.T)[qrange_x[:, np.newaxis], qrange_y]


    # TODO: why no transpose on these fields?
    m_ReEx_q = v_fields['Fxr'][qrange_x[:, np.newaxis], qrange_y]
    m_ReEy_q = v_fields['Fyr'][qrange_x[:, np.newaxis], qrange_y]
    #m_ImEx_q = v_fields['Fxi'][qrange_x[:, np.newaxis], qrange_y]
    #m_ImEy_q = v_fields['Fyi'][qrange_x[:, np.newaxis], qrange_y]


    # Ignore all imaginary values. If there are significant imag values,
    # then instaneous vector plots don't make much sense anyway
    d_quiv_kw = {'linewidths': (0.2,), 'edgecolors': ('gray'), 'pivot':'mid', 'headlength':5}
    if do_cont:  # no colours in the quiver
        d_quiv_kw['color'] = 'gray'
        ax.quiver(m_x_q, m_y_q, m_ReEx_q, m_ReEy_q,  ** d_quiv_kw)
    else:
        m_arrcolour = np.sqrt(m_ReEx_q*m_ReEx_q + m_ReEy_q*m_ReEy_q)
        ax.quiver(m_x_q, m_y_q, m_ReEx_q, m_ReEy_q, m_arrcolour, ** d_quiv_kw)


    if not do_cont:
        ax.set_aspect('equal')
        # this step is needed because quiver doesn't seem
        # to use its input x and y vectors to set range limits
        # clean this up so as to avoid seemingly circular calls following
        ax.set_xlim(v_x_q[0], v_x_q[-1])
        ax.set_ylim(v_y_q[0], v_y_q[-1])



def plot_contour_and_quiver(ax, d_xy, v_fields, plps, cc_scalar=None, cc_vector=None):

    #v_x, v_y, m_X, m_Y = list(d_xy.values())

    do_cont = not cc_scalar is None
    do_quiv = not cc_vector is None

    decorator = plps['decorator']

    cbar = None
    if do_cont:
        im_co, cbar = add_contour_plot(ax, d_xy, v_fields[cc_scalar._f_code], cc_scalar, plps, decorator)

    if do_quiv:
        add_quiver_plot(ax, d_xy, v_fields, cc_vector, plps, decorator, do_cont)

    # Adjustments to the visible plot domain

    is_AC = (cc_scalar and cc_scalar.is_AC()) or (cc_vector and cc_vector.is_AC())
    deftrim = -.05

    xlmi = plps.get('xlim_min', deftrim)
    xlma = plps.get('xlim_max', deftrim)
    ylmi = plps.get('ylim_min', deftrim)
    ylma = plps.get('ylim_max', deftrim)

    if is_AC: # By default, give a little space around elastic profiles
        if not xlmi: xlmi = deftrim
        if not xlma: xlma = deftrim
        if not ylmi: ylmi = deftrim
        if not ylma: ylma = deftrim


    if xlmi != 0 or xlma != 0:
        xmin, xmax = ax.get_xlim()
        width_x = xmax-xmin
        ax.set_xlim(xmin+xlmi*width_x, xmax-xlma*width_x)


    if ylmi != 0 or ylma != 0:
        ymin, ymax = ax.get_ylim()
        width_y = ymax-ymin
        ax.set_ylim(ymin+ylmi*width_y, ymax-ylma*width_y)


    labs = []
    if do_cont: labs.append(cc_scalar.get_label())
    if do_quiv: labs.append(cc_vector.get_label())
    comp_label = ', '.join(labs)

    lw = decorator.get_property('linewidth')
    ec = decorator.get_property('edgecolor')

    tidy = TidyAxes(nax=6,
                    props={'axes_color':ec, 'cb_edgecolor':ec,
                        'ax_label_xpad':3, 'ax_label_ypad':1 })

    tidy.apply_to_axes(ax)
    if cbar:
        tidy.apply_to_cbars(cbar)


    plot_set_ticks(ax, plps)
    #plot_set_axes_style(ax, plps, decorator)

    plot_set_title(ax, comp_label, plps, decorator)


    decorator.add_frame(ax)

    decorator.extra_axes_commands(ax)




def write_mode_data(ax, plps, sim_result, ival):  # mode data summary

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

    mode = sim_result.get_all_modes()[ival]
    mode.analyse_mode()


    # make a compact line writer
    def _write_line(tx, ty, s, tfs=fs):
        ax.text(tx, ty, s, transform=ax.transAxes, fontsize=tfs)

    [f_x, f_y, f_t, f_z] = mode.field_fracs()
    (r0x, r0y) = mode.center_of_mass() / SI_um
    (wx, wy, w0) = mode.second_moment_widths()  / SI_um




    _write_line(x0-.05, y0, f'Mode properties: m={ival}', fs+2); y0 -= dy

    if sim_result.is_EM():
        _write_line(x0, y0, r'$\omega/(2\pi)$: ' + f'{sim_result.omega_EM/(twopi*SI_THz):.5f} THz') ; y0 -= dy
        _write_line(x0, y0, r'$k$: ' + f'{sim_result.kz_EM(ival)/1.e6:.5f} ' + r'μm$^{{-1}}$');  y0 -= dy

        neff = sim_result.neff(ival)
        if sim_result.ngroup_EM_available():
            ng = sim_result.ngroup_EM(ival)
            _write_line(x0, y0, r'$\bar{{n}}$: ' + f'{neff:.6f},' + f'$n_g$: {ng:.6f}');  y0 -= dy
        else:
            _write_line(x0, y0, r'$\bar{{n}}$: ' + f'{neff:.6f}');  y0 -= dy
    else:
        q_AC = sim_result.q_AC
        nu_AC = np.real(sim_result.nu_AC(ival))
        vp = sim_result.vp_AC(ival)

        _write_line(x0, y0,
                   f'$q$: {q_AC * SI_um:.5f} ' + r'μm$^{{-1}}$, $\lambda:$ '+ f'{twopi/q_AC/SI_um:.5f} ' + r'μm') ;y0 -= dy

        _write_line(x0, y0, r'$q/2\pi$: ' + f'{q_AC*SI_um/twopi:.5f} '+ r'μm$^{{-1}}$'); y0 -= dy
        _write_line(x0, y0, r'$\Omega/(2\pi)$: ' + f'{nu_AC/SI_GHz:.5f} GHz' ); y0 -= dy
        if sim_result.vgroup_AC_available():
            vg = sim_result.vg_AC(ival)
            _write_line(x0, y0, f'$v_p$: {vp:.2f} m/s, $v_g$: {vg:.2f} m/s'); y0 -= dy
        else:
            _write_line(x0, y0, f'$v_p$: {vp:.2f} m/s'); y0 -= dy
    _write_line(x0, y0, f'$f_x:$ {f_x:.3f}, $f_y$: {f_y:.3f}'); y0 -= dy
    _write_line(x0, y0, f'$f_t:$ {f_t:.3f}, $f_z$: {f_z:.3f}'); y0 -= dy
    _write_line(x0, y0, r'$\mathbf{{r}}_0:$ '+ f'({r0x:.2f}, {r0y:.3f}) ' +r'μm'); y0 -= dy
    _write_line(x0, y0, f'$(w_x, w_y):$ ({wx:.2f}, {wy:.2f}) ' + r'μm'); y0 -= dy

    if mode.field_type == FieldType.EM_H:
        _write_line(x0, y0, r'$H$ field multiplied by $Z_0=376.7\, \Omega$'); y0 -= dy

    sc = sim_result.symmetry_classification(ival)
    if len(sc):
        _write_line(x0, y0, f'Sym: {sc}'); y0 -= dy

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
    #         sim_result.alpha_t_AC(ival),
    #         sim_result.alpha_s_AC(ival)/100.,
    #         sim_result.alpha_s_AC(ival)/100./(np.log(10.0)/10.0)
    #     ), transform=ax.transAxes, fontsize=fs)
    #     y0 -= dy

    #     ax.text(x0+.1, r, r'$Q_m$: {0:.2f}'.format(
    #         sim_result.Qmech_AC(ival)), transform=ax.transAxes, fontsize=fs)
    #     y0 -= dy
    #     ax.text(x0+.1, r, r'$\Delta\Omega/(2\pi)$: {0:.4f} MHz'.format(
    #         1.e-6*sim_result.linewidth_AC(ival)), transform=ax.transAxes, fontsize=fs)
    #     y0 -= dy
    #     mg_pe = plps['modal_gain'].get('PE', 0)
    #     mg_mb = plps['modal_gain'].get('MB', 0)
    #     mg_tot = plps['modal_gain'].get('Tot', 0)
    #     ax.text(x0+.1, r, r'Gain (PE, MB, Tot): ({0:.4f}, {1:.4f}, {2:.4f} ) $\mathrm{{(mW)}}^{{-1}}$'.format(
    #         mg_pe, mg_mb, mg_tot), transform=ax.transAxes, fontsize=fs)
    #     y0 -= dy

    d_ext = mode.get_mode_data()
    if len(d_ext):
        _write_line(x0, y0, r'Extra data:') ; y0 -= dy
        for (k, v) in d_ext.items():
            _write_line(x0, y0, f'   {k}: {v}'); y0 -= dy


def plot_all_components(d_xy, v_plots, plps, sim_result, ival):
    decorator = plps['decorator']
    figsz = decorator.get_property('figsize')
    ws = decorator.get_property('subplots_wspace')
    hs = decorator.get_property('subplots_hspace')


    decorator.set_frame_drawer(sim_result._structure.wg_geom)


    hide_minors = plps['suppress_imimre']
    rows = 2 if hide_minors else 3

    fig, axs = plt.subplots(rows, 3, figsize=figsz, )
    fig.subplots_adjust(hspace=hs, wspace=ws)

    axs = axs.flat
    axi = 0

    ax = axs[axi]; axi += 1
    write_mode_data(ax, plps, sim_result, ival)  # mode data summary

    ft = plps['EM_AC']
    cc_scal = field_type_to_intensity_code(ft)
    cc_transvec = field_type_to_vector_code(ft)

    ax = axs[axi]; axi += 1
    plot_contour_and_quiver(ax, d_xy, v_plots, plps, cc_scalar=cc_scal, cc_vector=cc_transvec)  # the whole field

    ax = axs[axi]; axi += 1
    plot_contour_and_quiver(ax, d_xy, v_plots, plps, cc_vector=cc_transvec)  # the intensity

    for flab in v_plots.keys():
        cc = component_t.make_comp(ft, flab)

        if (hide_minors and not cc.is_dominant()) or not cc.is_signed_field():
            continue

        ax = axs[axi]; axi += 1
        plot_contour_and_quiver(ax, d_xy, v_plots, plps, cc_scalar=cc)  # the scalar plots


    fig_fname = modeplot_filename(plps, ival)
    save_and_close_figure(fig,  fig_fname)




def plot_one_component(d_xy, v_fields, plps, ival, cc, axis=None):

    #decorator = plps['decorator']

    if axis is None:
        fig, ax = plt.subplots(figsize=(12, 10))
    else:
        ax = axis

    ft = plps['EM_AC']
    cc_scal = field_type_to_intensity_code(ft)
    cc_transvec = field_type_to_vector_code(ft)


    plot_contour_and_quiver(ax, d_xy, v_fields, plps, cc_scalar=cc_scal, cc_vector=cc_transvec)

    if axis is None:  # If user passed in the axis, they can look after saving.
        fig_fname = modeplot_filename(plps, ival, cc._user_code)
        save_and_close_figure(fig, fig_fname)







###########################################################################################################
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

    print('Warning: "plt_mode_fields" is deprecated, use "plot_modes"')
    sim_result.plot_modes(ivals, n_points, quiver_points,
                     xlim_min, xlim_max, ylim_min, ylim_max,
                     field_type, num_ticks, colorbar, contours, contour_lst, stress_fields, pdf_png,
                     prefix, suffix, ticks, comps, decorator,
                     suppress_imimre, modal_gains_PE, modal_gains_MB, modal_gains)



       # # TODO: what is this for?
    # modal_gain = {}
    # if modal_gains is not None:
    #     modal_gain['Tot'] = modal_gains[ival]
    # if modal_gains_PE is not None:
    #     modal_gain['PE'] = modal_gains_PE[ival]
    # if modal_gains_MB is not None:
    #     modal_gain['MB'] = modal_gains_MB[ival]


#    plot_params={'xlim_min': xlim_min, 'xlim_max': xlim_max, 'ylim_min': ylim_min,
#                 'ylim_max': ylim_max, 'ticks': ticks, 'num_ticks':num_ticks,
#                  'colorbar':colorbar, 'contours':contours, 'contour_lst':contour_lst, 'field_type':field_type,
#                  'prefix': prefix, 'suffix': suffix, 'pdf_png': pdf_png,
#                  'stress_fields':stress_fields, 'modal_gain':modal_gain, 'decorator': decorator,
#                  'suppress_imimre':suppress_imimre,
#              # 'n_pts_x': n_pts_x, 'n_pts_y': n_pts_y,
#               'quiver_points': quiver_points }
## modal_gains_PE=modal_gains_PE, #purpose, can these go elsewhere?
    # modal_gains_MB=modal_gains_MB,
    # modal_gains=modal_gains)
