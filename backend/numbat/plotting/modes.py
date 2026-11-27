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

from numbat.reporting import OutputLevel, nbprint

import copy
import warnings

import cmasher as cmr
import matplotlib.colors as mplcolors
import matplotlib.pyplot as plt
import numpy as np

from numbat.field_selection import field_symbol, selection_from_storage_key, storage_key
from numbat.nbtypes import (
    ComplexPart,
    FieldComponent,
    FieldSelection,
    FieldType,
    SI_GHz,
    SI_nm,
    SI_THz,
    SI_um,
    twopi,
)
from numbat.plotting.axes import TidyAxes
from numbat.plotting.tools import save_and_close_figure

from .. import runtime


def modeplot_filename_2D(field_type, plps, mode_index, label='', cut=''):
    pref=plps['prefix']
    suf = plps['suffix']
    comp = FieldType.parse(field_type).name

    fullpref = str(runtime.outdir_fields_path(prefix=pref))
    filestart = f'{fullpref}/{comp}_mode_{mode_index:02d}'

    if cut:
        filestart = filestart + f'_{cut}cut'
    if suf:
        filestart = filestart + f'_{suf}'

    if label: filestart += '_'+label

    fn = filestart + runtime.plotfile_ext()
    return fn

# TODO: allow a suffix to label multiple cut planes
def modeplot_filename_1D(field_type, plps, mode_index, cut, label=''):

    return modeplot_filename_2D(field_type, plps, mode_index, label, cut)


def field_type_to_intensity_code(ft):
    return FieldSelection.intensity(ft)

def field_type_to_vector_code(ft):
    return FieldSelection.transverse(ft)


def field_display_name(field_type):
    return {
        FieldType.EM_E: "electric",
        FieldType.EM_H: "magnetic",
        FieldType.AC: "displacement",
    }[FieldType.parse(field_type)]


def domain_display_name(field_type):
    return "acoustic" if FieldType.parse(field_type) == FieldType.AC else "optical"


def field_type_label(field_type):
    return "AC" if FieldType.parse(field_type) == FieldType.AC else "EM"


def tex_label(field_sel):
    if field_sel.is_cartesian:
        primary_is_real = field_sel.component in (FieldComponent.X, FieldComponent.Y)
        is_real = (field_sel.complex_part == ComplexPart.PRIMARY) == primary_is_real
        part = "Re" if is_real else "Im"
        return f"{part}[ ${field_symbol(field_sel.field_type)}_{field_sel.component.value}$ ]"
    if field_sel.component == FieldComponent.MAGNITUDE:
        return f"$|\\vec {field_symbol(field_sel.field_type)}|$"
    if field_sel.component == FieldComponent.INTENSITY:
        return f"$|\\vec {field_symbol(field_sel.field_type)}|^2$"
    return f"$\\vec {field_symbol(field_sel.field_type)}_t$"


def line_color(field_sel):
    if field_sel.component in (FieldComponent.MAGNITUDE, FieldComponent.INTENSITY):
        return "red"
    return {
        FieldComponent.X: "blue",
        FieldComponent.Y: "green",
        FieldComponent.Z: "brown",
    }.get(field_sel.component)


def line_style(field_sel, *, show_magnitude):
    if field_sel.is_derived or not show_magnitude:
        return "solid"
    return "dotted" if field_sel.complex_part == ComplexPart.SECONDARY else "dashed"


def selection_for_component(field_type, component):
    component = FieldComponent.parse(component)
    if component in (FieldComponent.X, FieldComponent.Y, FieldComponent.Z):
        return FieldSelection.primary(field_type, component)
    return FieldSelection(field_type, component)



class Decorator(object):

    def __init__(self, multi=False):
        mp_base_fs = 18
        sp_base_fs = 24

        self._multi_props = {'figsize': (10, 8),'subplots_hspace': .2, 'subplots_wspace': 1.2,
                             'base_fs': mp_base_fs,
                             'title_fs': mp_base_fs-2,
                             'subtitle_fs': mp_base_fs-7,  # used for component labels
                             'title_pad': 10, 'subtitle_pad': 10,
                             'ax_label_fs': mp_base_fs-10, 'ax_label_pad': 20, 'ax_tick_fs': mp_base_fs-10,
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
        self._wireframe_delta_xy = (0.0, 0.0)

    def __getitem__(self, kw):
        return self._props[kw]

    def __setitem__(self, kw, val):
        self._props[kw]=val

    def set_waveguide(self, wg):
        self._waveguide = wg
        self._wireframe_delta_xy = (0.0, 0.0)

    def configure_waveguide(self, structure):
        """Set the geometry and its translation in micrometre plot coordinates."""
        self.set_waveguide(structure.wg_geom)
        translation = getattr(structure, "translation", (0.0, 0.0))
        self._wireframe_delta_xy = tuple(
            value * SI_nm / SI_um for value in translation
        )

    def add_waveguide_border(self, ax):
        plprefs = runtime.plot_prefs()
        styles = {}

        styles['edgecolor'] = plprefs.WG_FRAME_EDGE_COLOR
        if self._is_single:
            styles['linewidth'] = plprefs.WG_FRAME_LINEWIDTH_WHOLEFIG
        else:
            styles['linewidth'] = plprefs.WG_FRAME_LINEWIDTH_SUBFIG


        self._waveguide.draw_mpl_frame(ax, styles, self._wireframe_delta_xy)

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
            nbprint(f'Warning: unknown fontsize label "{lab}" in Decorator::get_property()', level=OutputLevel.ALWAYS)
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

        self.update({'trim': None,
                     'xlim_min': None, 'xlim_max': None,
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
                      'quiver_points': 30,
                      'overview_layout': 'auto',
                      'wide_mode': None,
                      'colorbar_allowance': 0.15,
                      })


def resolve_plot_trim(plps, default_trim):
    """Return left, bottom, right, and top plot-boundary fractions."""
    legacy_keys = ("xlim_min", "ylim_min", "xlim_max", "ylim_max")
    legacy_values = tuple(plps.get(key, None) for key in legacy_keys)
    if any(value is not None for value in legacy_values):
        warnings.warn(
            "xlim_min, ylim_min, xlim_max, and ylim_max are deprecated; "
            "use trim=(left, bottom, right, top) instead.",
            DeprecationWarning,
            stacklevel=2,
        )

    trim = plps.get("trim", None)
    if trim is not None:
        try:
            left, bottom, right, top = tuple(trim)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "trim must be a four-item (left, bottom, right, top) tuple."
            ) from exc
        try:
            values = tuple(float(value) for value in (left, bottom, right, top))
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "trim values must be finite numeric fractions."
            ) from exc
        if not np.all(np.isfinite(values)):
            raise ValueError("trim values must be finite numeric fractions.")
        return values

    return (
        plps.get("xlim_min", default_trim, override_None=True),
        plps.get("ylim_min", default_trim, override_None=True),
        plps.get("xlim_max", default_trim, override_None=True),
        plps.get("ylim_max", default_trim, override_None=True),
    )


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


def plot_set_title(ax, comp_label, plps, decorator, tidy):
    if plps.get('add_title', True):
        ##ax.set_title(comp_label, fontsize=decorator.get_property('subtitle_fs'),
          #           pad=decorator.get_property('subtitle_pad'))
        ax.set_title(comp_label, fontsize=tidy._props['subtitle_fs'],
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


def add_contour_plot(fig, ax, d_xy, c_field, field_sel, plps, decorator):

    logamp = plps['logamp']

    cmap_signed = runtime.plot_prefs().cmap_field_signed(field_sel)
    cmap_unsigned = runtime.plot_prefs().cmap_field_unsigned(field_sel)

    is_signed_quantity = field_sel.is_cartesian
    #cmap = cmap_signed if (is_signed_quantity and not logamp) else cmap_unsigned
    cmap = cmap_signed if is_signed_quantity else cmap_unsigned

    v_x, v_y, m_X, m_Y = list(d_xy.values())

    # Extract and tidy the field
    # transpose because arrays
    cont_field = (c_field).copy().T

    if field_sel.component == FieldComponent.MAGNITUDE:
        cont_field = np.abs(cont_field)

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

    common_amp = True and is_signed_quantity

    if decorator.get_cmap_limits(field_sel.component.value) is not None:
        (act_zlo, act_zhi) = decorator.get_cmap_limits(field_sel.component.value)
        tsnorm = mplcolors.TwoSlopeNorm(
            vmin=act_zlo, vmax=act_zhi, vcenter=(act_zlo+act_zhi)/2)

        d_kw['norm'] = tsnorm
    elif logamp:
        vma = logmax
        if is_signed_quantity:
            vmi = -logmax
        else:
            vmi=0

        d_kw['vmin'] = vmi
        d_kw['vmax'] = vma

    else:
        if common_amp:
            act_zlo = -1
            act_zhi = 1
        else:
            act_zlo = np.nanmin(cont_field)
            act_zhi = np.nanmax(cont_field)

        vma = max(abs(act_zlo), abs(act_zhi))
        vmi = -vma if is_signed_quantity else 0.0

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

    xlmi, ylmi, xlma, ylma = resolve_plot_trim(plps, deftrim)

    quiver_skip_x = int(round(n_pts_x/quiver_points_x * (1-xlmi-xlma)))
    quiver_skip_y = int(round(n_pts_y/quiver_points_y * (1-ylmi-ylma)))

    # getting a nice symmetric pattern of points to do quivers centred around the middle
    qslice_x = get_quiver_skip_range(n_pts_x, quiver_skip_x)
    qslice_y = get_quiver_skip_range(n_pts_y, quiver_skip_y)

    return qslice_x, qslice_y

def add_quiver_plot(fig, ax, d_xy, v_fields, field_sel, plps, decorator, do_cont):


    # give a little space around elastic profiles
    deftrim = -0.1 if field_sel.field_type == FieldType.AC else -0.05

    v_x, v_y, m_x, m_y = d_xy.values()

    # shapes:
    # m_X and other matrices are stored with y points as rows but upside down
    # mesh_grid_structure:
    # m_X =  [[x0, x1, x2, ...]
    #         [x0, x1, x2, ...]
    #         [x0, x1, x2, ...]
    # m_Y =  [[y0, y0, y0, ...]
    #         [y1, y1, y1, ...]
    #         [y2, y2, y2, ...]



    # TODO: spend a day getting truly perfect spacings of grid poitns
    #         and trying to figure out the scaling and arrow head settings once and for all

    qslice_x, qslice_y = _make_quiver_slices(v_x, v_y, plps, deftrim)

    v_x_q = v_x[qslice_x]
    v_y_q = v_y[qslice_y]
    m_x_q = (m_x.T)[qslice_x[:, np.newaxis], qslice_y]  # TODO: track down why m_x/y need .T but fields don't
    m_y_q = (m_y.T)[qslice_x[:, np.newaxis], qslice_y]

    # TODO: why no transpose on these fields?
    m_ReEx_q = v_fields[
        storage_key(FieldSelection.primary(field_sel.field_type, FieldComponent.X))
    ][qslice_x[:, np.newaxis], qslice_y]
    m_ReEy_q = v_fields[
        storage_key(FieldSelection.primary(field_sel.field_type, FieldComponent.Y))
    ][qslice_x[:, np.newaxis], qslice_y]


    vecarrow_col = runtime.plot_prefs().vector_field_arrow_color(field_sel)
    vecarrow_scale = runtime.plot_prefs().vecplot_arrow_scale
    vecarrow_width = runtime.plot_prefs().vecplot_arrow_width
    vecarrow_hw = runtime.plot_prefs().vecplot_arrow_headwidth



    # Ignore all imaginary values. If there are significant imag values,
    # then instantaneous vector plots don't make much sense anyway
    d_quiv_kw = {#'linewidths': (vecarrow_lw,),
                 'width': vecarrow_width,
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


def plot_contour_and_quiver(
        fig, ax, d_xy, v_fields, plps, contour_field_sel=None,
        quiver_field_sel=None, is_single_plot=False):

    is_AC = (
        contour_field_sel and contour_field_sel.field_type == FieldType.AC
    ) or (
        quiver_field_sel and quiver_field_sel.field_type == FieldType.AC
    )

    do_cont = contour_field_sel is not None
    do_quiv = quiver_field_sel is not None

    decorator = plps['decorator']

    cbar = None
    if do_cont:
        im_co, cbar = add_contour_plot(
            fig,
            ax,
            d_xy,
            v_fields[storage_key(contour_field_sel)],
            contour_field_sel,
            plps,
            decorator,
        )

    if do_quiv:
        add_quiver_plot(fig, ax, d_xy, v_fields,
                        quiver_field_sel, plps, decorator, do_cont)

    # Adjustments to the visible plot domain

    deftrim = -.1 if is_AC else -0.05
    xlmi, ylmi, xlma, ylma = resolve_plot_trim(plps, deftrim)

    if xlmi != 0 or xlma != 0:
        xmin, xmax = ax.get_xlim()
        width_x = xmax-xmin
        ax.set_xlim(xmin+xlmi*width_x, xmax-xlma*width_x)

    if ylmi != 0 or ylma != 0:
        ymin, ymax = ax.get_ylim()
        width_y = ymax-ymin
        ax.set_ylim(ymin+ylmi*width_y, ymax-ylma*width_y)


    labs = []
    if do_cont: labs.append(tex_label(contour_field_sel))
    if do_quiv: labs.append(tex_label(quiver_field_sel))
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
        plot_set_title(ax, comp_label, plps, decorator, tidy)

    decorator.add_waveguide_border(ax)
    decorator.extra_axes_commands(ax)





def _displayed_waveguide_aspect(field_type, d_xy, plps):
    v_x = d_xy["v_x"]
    v_y = d_xy["v_y"]
    default_trim = -0.1 if FieldType.parse(field_type) == FieldType.AC else -0.05
    xlim_min, ylim_min, xlim_max, ylim_max = resolve_plot_trim(
        plps, default_trim
    )
    x_fraction = 1 - xlim_min - xlim_max
    y_fraction = 1 - ylim_min - ylim_max
    waveguide_aspect = abs(
        (v_x[-1] - v_x[0]) * x_fraction / ((v_y[-1] - v_y[0]) * y_fraction)
    )
    return waveguide_aspect


def _overview_grid_shape(layout, hide_secondary_comps):
    if layout == "tall":
        return (2 if hide_secondary_comps else 3), 3
    return 3, (2 if hide_secondary_comps else 3)


def _colorbar_allowance(plps):
    if not plps["colorbar"]:
        return 0.0

    try:
        allowance = float(plps.get("colorbar_allowance", 0.15))
    except (TypeError, ValueError) as exc:
        raise ValueError("colorbar_allowance must be a non-negative number.") from exc
    if allowance < 0:
        raise ValueError("colorbar_allowance must be a non-negative number.")
    return allowance


def _overview_width_units(layout, hide_secondary_comps, plps):
    _, columns = _overview_grid_shape(layout, hide_secondary_comps)
    return (
        columns
        + (columns - 1) * plps["decorator"]["subplots_wspace"]
        + columns * _colorbar_allowance(plps)
    )


def _overview_base_figsize(layout, plps):
    figure_width, figure_height = plps["decorator"]["figsize"]
    if layout == "wide":
        figure_width *= 1.5
    return figure_width, figure_height


def _overview_cell_aspect(layout, hide_secondary_comps, plps):
    figure_width, figure_height = _overview_base_figsize(layout, plps)
    rows, _ = _overview_grid_shape(layout, hide_secondary_comps)
    return figure_width * rows / (figure_height * _overview_width_units(
        layout, hide_secondary_comps, plps
    ))


_SQUARE_TALL_ASPECT_RANGE = (0.8, 1.25)
_SQUARE_TALL_PANEL_SCALE = 1.2


def select_overview_layout(field_type, d_xy, plps, hide_secondary_comps):
    """Choose a combined-overview layout from explicit or automatic settings."""
    wide_mode = plps.get("wide_mode", None)
    if wide_mode is not None:
        return "wide" if wide_mode else "tall"

    layout = plps.get("overview_layout", "auto")
    if layout not in ("auto", "tall", "wide"):
        raise ValueError(
            "overview_layout must be one of: 'auto', 'tall', or 'wide'."
        )
    if layout != "auto" or not hide_secondary_comps:
        return "tall" if layout == "auto" else layout

    waveguide_aspect = _displayed_waveguide_aspect(field_type, d_xy, plps)
    if waveguide_aspect <= 2:
        return "tall"

    tall_score = abs(
        np.log(_overview_cell_aspect("tall", hide_secondary_comps, plps) / waveguide_aspect)
    )
    wide_score = abs(
        np.log(_overview_cell_aspect("wide", hide_secondary_comps, plps) / waveguide_aspect)
    )
    return "wide" if wide_score < tall_score else "tall"


def overview_figure_size(field_type, d_xy, plps, hide_secondary_comps, layout):
    """Return a bounded figure size that targets equal-aspect field panels."""
    base_width, base_height = _overview_base_figsize(layout, plps)
    waveguide_aspect = _displayed_waveguide_aspect(field_type, d_xy, plps)
    if waveguide_aspect <= 0:
        return base_width, base_height

    base_aspect = base_width / base_height
    rows, _ = _overview_grid_shape(layout, hide_secondary_comps)
    target_aspect = waveguide_aspect * _overview_width_units(
        layout, hide_secondary_comps, plps
    ) / rows
    bounded_aspect = np.clip(target_aspect, base_aspect / 1.5, base_aspect * 1.5)
    figure_area = base_width * base_height
    return (
        np.sqrt(figure_area * bounded_aspect),
        np.sqrt(figure_area / bounded_aspect),
    )


def configure_overview_subplots(
        fig, layout, waveguide_aspect, rows, columns, hspace, wspace):
    """Apply overview spacing while enlarging near-square tall field panels."""
    adjustments = {"hspace": hspace, "wspace": wspace}
    if (
        layout == "tall"
        and _SQUARE_TALL_ASPECT_RANGE[0]
        <= waveguide_aspect
        <= _SQUARE_TALL_ASPECT_RANGE[1]
    ):
        adjusted_wspace = (
            (columns + (columns - 1) * wspace) / _SQUARE_TALL_PANEL_SCALE
            - columns
        ) / (columns - 1)
        subplot_params = fig.subplotpars
        current_axes_height = (
            (subplot_params.top - subplot_params.bottom)
            / (rows + (rows - 1) * hspace)
        )
        adjusted_hspace = 0.05
        adjusted_height = current_axes_height * _SQUARE_TALL_PANEL_SCALE * (
            rows + (rows - 1) * adjusted_hspace
        )
        vertical_centre = (subplot_params.top + subplot_params.bottom) / 2
        adjustments.update({
            "hspace": adjusted_hspace,
            "wspace": adjusted_wspace,
            "bottom": vertical_centre - adjusted_height / 2,
            "top": vertical_centre + adjusted_height / 2,
        })

    fig.subplots_adjust(**adjustments)


def offset_mode_data_panel(fig, ax, fraction=0.1):
    """Move the mode-data panel northwest by a fraction of its physical diagonal."""
    position = ax.get_position()
    figure_width, figure_height = fig.get_size_inches()
    panel_diagonal = np.hypot(
        position.width * figure_width, position.height * figure_height
    )
    diagonal_component = fraction * panel_diagonal / np.sqrt(2)
    shift_x = diagonal_component / figure_width
    shift_y = diagonal_component / figure_height
    ax.set_position((
        max(0, position.x0 - shift_x),
        min(1 - position.height, position.y0 + shift_y),
        position.width,
        position.height,
    ))


def write_mode_data(ax, plps, sim_result, mode_index, field_type, layout):

    decorator = plps['decorator']
    fs = decorator.get_property('data_label_fs')

    wide_mode = layout == "wide"

    x0 = .05
    yhi = .99
    dy = .10 if wide_mode else .08
    y0 = yhi-dy

    ax.set_xlim(0, 1)
    ax.set_ylim(0, yhi)
    ax.set_aspect("auto")
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
    lines2=[]

    if sim_result.is_EM():
        lines.append( r'$\omega/(2\pi)$: ' + f'{sim_result.omega_EM/(twopi*SI_THz):.5f} THz')
        lines.append( r'$k$: ' + f'{sim_result.kz_EM(mode_index)/1.e6:.5f} ' + r'μm$^{{-1}}$')

        neff = sim_result.neff(mode_index)
        lines.append( r'$\bar{{n}}$: ' + f'{neff:.6f}')
        if sim_result.ngroup_EM_available():
            ng = sim_result.ngroup_EM(mode_index)
            lines.append( f'$n_g$: {ng:.6f}')
    else:
        q_AC = sim_result.q_AC
        nu_AC = np.real(sim_result.nu_AC(mode_index))
        vp = sim_result.vp_AC(mode_index)

        lines.append(f'$q$: {q_AC * SI_um:.5f} ' + r'μm$^{{-1}}$') #;y0 -= dy
        lines.append(r'$\lambda:$ '+ f'{twopi/q_AC/SI_um:.5f} ' + r'μm') #;y0 -= dy


        lines.append( r'$q/2\pi$: ' + f'{q_AC*SI_um/twopi:.5f} '+ r'μm$^{{-1}}$')
        lines.append( r'$\Omega/(2\pi)$: ' + f'{nu_AC/SI_GHz:.5f} GHz' )
        lines.append( f'$v_p$: {vp:.2f} m/s')
        if sim_result.vgroup_AC_available():
            vg = sim_result.vg_AC(mode_index)
            lines.append( f'$v_g$: {vg:.2f} m/s')

    lines2.append( f'$f_x:$ {f_x:.3f}, $f_y$: {f_y:.3f}')
    lines2.append( f'$f_t:$ {f_t:.3f}, $f_z$: {f_z:.3f}')
    lines2.append( r'$\mathbf{{r}}_0:$ '+ f'({r0x:.3f}, {r0y:.3f}) ' +r'μm')
    lines2.append( f'$(w_x, w_y):$ ({wx:.3f}, {wy:.3f}) ' + r'μm')

    if field_type == FieldType.EM_H:
        lines2.append( r'$H$ field multiplied by $Z_0=376.7\, \Omega$')

    sc = sim_result.symmetry_classification(mode_index)
    if len(sc):
        lines2.append( f'Sym: {sc}')

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
        lines2.append( r'Extra data:')
        for (k, v) in d_ext.items():
            lines2.append( f'   {k}: {v}')

    ty0 = y0
    for line in lines:
        _write_line(x0, ty0, line)
        ty0 -= dy

    for line in lines2:
        _write_line(x0, ty0, line)
        ty0 -= dy


def plot_all_components(field_type, d_xy, v_plots, plps, sim_result, mode_index):
    decorator = plps['decorator']
    hs = decorator['subplots_hspace']
    ws = decorator['subplots_wspace']

    hide_secondary_comps = plps['suppress_imimre']
    mod_is_amp = plps['mod_is_amp']

    layout = select_overview_layout(field_type, d_xy, plps, hide_secondary_comps)
    tall_layout = layout == "tall"
    figsz = overview_figure_size(field_type, d_xy, plps, hide_secondary_comps, layout)
    waveguide_aspect = _displayed_waveguide_aspect(field_type, d_xy, plps)

    if tall_layout:
        rows = 2 if hide_secondary_comps else 3
        columns = 3

        fig, axs = plt.subplots(rows, columns, figsize=figsz)

        ax_md_data = axs[0, 0]  # mode data summary
        ax_ut = axs[0, 1]  # the whole field
        ax_uv = axs[0, 2]  # the intensity
        axs_comp = axs[1:, :].flatten()  # the components

    else:
        rows = 3
        columns = 2 if hide_secondary_comps else 3

        fig, axs = plt.subplots(rows, columns, figsize=figsz)

        ax_md_data = axs[0, 0]  # mode data summary
        ax_ut = axs[1, 0]  # the whole field
        ax_uv = axs[2, 0]  # the intensity
        axs_comp = axs[:, 1:].flatten()  # the components

    configure_overview_subplots(
        fig, layout, waveguide_aspect, rows, columns, hs, ws
    )
    offset_mode_data_panel(fig, ax_md_data)
    write_mode_data(ax_md_data, plps, sim_result, mode_index, field_type, layout)

    ft = field_type

    scal_comp = 'a' if mod_is_amp else 'i'  # intensity or amplitude plot
    contour_field_sel = selection_for_component(ft, scal_comp)
    quiver_field_sel = FieldSelection.transverse(ft)

    if plps['hide_vector_field']:
        joint_quiver_field_sel = None
    else:
        joint_quiver_field_sel = quiver_field_sel

    plot_contour_and_quiver(fig, ax_ut, d_xy, v_plots, plps,
                            contour_field_sel=contour_field_sel,
                            quiver_field_sel=joint_quiver_field_sel)

    plot_contour_and_quiver(
        fig, ax_uv, d_xy, v_plots, plps, quiver_field_sel=quiver_field_sel
    )

    axi = 0
    for flab in v_plots.keys():
        field_sel = selection_from_storage_key(ft, flab)

        if not field_sel.is_cartesian:
            continue

        # don't plot what will likely be pure zero fields
        if hide_secondary_comps and field_sel.complex_part == ComplexPart.SECONDARY:
            continue

        ax = axs_comp[axi]; axi += 1
        plot_contour_and_quiver(
            fig, ax, d_xy, v_plots, plps, contour_field_sel=field_sel
        )


    fig_fname = modeplot_filename_2D(field_type, plps, mode_index)
    save_and_close_figure(fig, fig_fname)




def plot_one_component(field_type, d_xy, v_fields, plps, mode_index, field_sel, axis=None):
    decorator = plps['decorator']
    fig=None # a problem for passed in axes
    if axis is None:
        figsz = decorator['figsize']
        fig, ax = plt.subplots(figsize=figsz)
    else:
        ax = axis

    if field_sel.is_vector:
        quiver_field_sel = field_sel
        contour_field_sel = FieldSelection.magnitude(field_type)
    else:
        quiver_field_sel = None
        contour_field_sel = field_sel

    plot_contour_and_quiver(fig, ax, d_xy, v_fields, plps,
                            contour_field_sel=contour_field_sel,
                            quiver_field_sel=quiver_field_sel,
                            is_single_plot=True)

    if axis is None:  # If user passed in the axis, they can look after saving.
        component_label = f"{field_symbol(field_sel.field_type)}{field_sel.component.value}"
        fig_fname = modeplot_filename_2D(field_type, plps, mode_index, component_label)
        save_and_close_figure(fig, fig_fname)
