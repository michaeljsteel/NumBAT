
from collections.abc import Iterable
from pathlib import Path
import tomllib


import numbat
import reporting

# make the NumBATPlotPrefs call return a reference to object in plottools.py
class PlotPrefs:
    """Selection of color scales and line properties for different field plots."""

    # Add color combinations here
    # See https://matplotlib.org/stable/users/explain/colors/colormaps.html

                  #  signed, unsigned, vector arrow
    color_tup_1 = ('seismic', 'OrRd', 'dimgray')
    color_tup_2 = ('PRGn', 'GnBu', 'brown')
    color_tup_3 = ('BrBG', 'YlOrBr', 'black')
    color_tup_4 = ('PuOr', 'YlOrBr', 'dimgray')
    color_tup_5 = ('coolwarm', 'Reds', 'dimgray')
    color_tup_6 = ('PiYG', 'GnBu', 'dimgray')
    color_tup_7 = ('RdYlGn', 'GnBu', 'dimgray')


    def __init__(self, user_settings_file='', ignore_user_settings=False):

        self._plot_extension = '.png'
        self. _load_user_settings( user_settings_file, ignore_user_settings)
        self._set_defaults()

    # TODO: do this with logger
    def _load_user_settings(self, user_settings_file='', ignore_user_settings=False):
        self._user_settings = {
            'all_plots':{},
            'multi_plots':{},
            'single_plots':{}
                               }

        if ignore_user_settings:
            print('Ignoring any user setting files and using NumBAT default preferences.')
            return

        user_file=None
        if user_settings_file:
            if not Path(user_settings_file).exists():
                s=f"NumBAT can't find the user settings .toml file at {user_settings_file}."
                reporting.report_and_exit(s)

            user_file=Path(user_settings_file)
        else: # Look for standard locations
            home_dir = Path.home()
            locs = [Path('./numbat.toml'),
                    Path(home_dir / '.numbat.toml'),
                    Path(home_dir / 'numbat.toml')

                    ]
            for loc in locs:
                if loc.exists():
                    user_file = loc
                    break

        if user_file is not None:
            print(f'Loading plot preferences from {user_file}.')
            with open(user_file, 'rb') as fin:
                self._user_settings = tomllib.load(fin)

    def _set_defaults(self):

        d_all = self._user_settings['all_plots']

        # Select color combinations here
        coltup_em = self.color_tup_1
        coltup_ac = self.color_tup_7

        self.cmap_em_field_signed = d_all.get('em_colormap_signed', coltup_em[0])
        self.cmap_em_field_unsigned = d_all.get('em_colormap_unsigned', coltup_em[1])
        self.vecplot_arrow_color_em = d_all.get('em_vector_arrow_color', coltup_em[2])

        self.cmap_ac_field_signed = d_all.get('ac_colormap_signed', coltup_ac[0])
        self.cmap_ac_field_unsigned = d_all.get('ac_colormap_unsigned', coltup_ac[1])
        self.vecplot_arrow_color_ac = d_all.get('ac_vector_arrow_color', coltup_ac[2])

        # larger makes smaller arrows
        self.vecplot_arrow_scale = d_all.get('vector_arrow_scale', 10)

        # width of the shaft in fractions of plot width
        self.vecplot_arrow_linewidth = d_all.get('vector_arrow_linewidth', 0.005)

        # multiple of shaft width
        self.vecplot_arrow_headwidth = d_all.get('vector_arrow_headwidth', 2)

        #self.vecplot_arrow_head = 2
        #self.vecplot_arrow_len = 2

        # colormap for refractive index plots
        self.cmap_ref_index = d_all.get('refindex_colormap', 'cool')


        # gain plots
        self.mode_index_label_fs = d_all.get('mode_index_label_fs', 10)

        # properties of waveguide boundary frames
        # expressed in all caps since to a user mesh these should be seen as constants

        # EDGE_COLOR="gray" "dimgray" "black" "brown" "darkred"

        self.WG_FRAME_EDGE_COLOR = d_all.get('wg_frame_edge_color', 'darkred')

        # if field plot is whole figure
        self.WG_FRAME_LINEWIDTH_WHOLEFIG = d_all.get('wg_frame_linewidth_color', 0.75)
        # if field plot is a part figure
        self.WG_FRAME_LINEWIDTH_SUBFIG = d_all.get('wg_frame_linewidth_subfig', 0.25)

        self.plot_output_resolution_dpi = d_all.get('plot_output_resolution_dpi', 300)



        d_multi = self.d_multi = {}
        d_single = self.d_single = {}

        d_multi['ax_label_fs'] = 8
        d_multi['ax_ticklabel_fs'] = 8
        d_multi['ax_tickwidth'] = 0.25
        d_multi['ax_linewidth'] = 1
        d_multi['title_fs'] = 16
        d_multi['subtitle_fs'] = 11



        d_single['ax_label_fs'] = 14
        d_single['ax_ticklabel_fs'] = 14
        d_single['ax_tickwidth'] = 1
        d_single['ax_linewidth'] = 1
        d_single['title_fs'] = 22
        d_single['subtitle_fs'] = 21


        d_multi.update(self._user_settings['multi_plots'])
        d_single.update(self._user_settings['single_plots'])


    def cmap_field_signed(self, ftag):
        if ftag.is_EM():
            return self.cmap_em_field_signed
        else:
            return self.cmap_ac_field_signed

    def cmap_field_unsigned(self, ftag):
        if ftag.is_EM():
            return self.cmap_em_field_unsigned
        else:
            return self.cmap_ac_field_unsigned

    def vector_field_arrow_color(self, ftag):
        if ftag.is_EM():
            return self.vecplot_arrow_color_em
        else:
            return self.vecplot_arrow_color_ac


    # def vector_field_arrow_scale(self):
    #     return self.vecplot_arrow_scale

    # def vector_field_arrow_linewidth(self):
    #     return self.vecplot_arrow_linewidth


class TidyAxes:

    def __init__(self, nax=1, props={}):

        self._nax = nax
        self._set_defaults_for_num_axes(nax)
        self._props.update(props)

    def update_property(self, k, v):
        self._props[k]=v

    def _set_defaults_for_num_axes(self, nax):
        user_prefs = numbat.NumBATPlotPrefs()

        props = {}

        #axtags = ('ax_label_fs', 'ax_ticklabel_fs', 'ax_tickwidth', 'ax_linewidth')
        if nax in (4,6):  # includes layout for two row mode plots
            #for s in axtags:
            #    props[s] = user_prefs.d_multi[s]
            props.update(user_prefs.d_multi)


            props['ax_label_xpad'] = 3
            props['ax_label_ypad'] = 1

        #    props['ax_label_fs'] = 8
        ##    props['ax_ticklabel_fs'] = 8
         ##   props['ax_tickwidth'] = .25
          #  props['ax_linewidth'] = 1

        elif nax == 2:
            # for s in axtags:
            #     props[s] = user_prefs.d_multi[s]
            props.update(user_prefs.d_multi)


            props['ax_label_xpad'] = 3
            props['ax_label_ypad'] = 1
            props['ax_label_pad'] = 5

            ##props['ax_label_fs'] = 12
            #props['ax_ticklabel_fs'] =  10
            #props['ax_tickwidth'] = 1
            #props['ax_linewidth'] = 1
        else:
            #for s in axtags:
            #    props[s] = user_prefs.d_single[s]
            props.update(user_prefs.d_single)

            props['ax_label_pad'] = 5
            props['ax_label_xpad'] = 3
            props['ax_label_ypad'] = 1

            #props['ax_label_fs'] =  14
            #props['ax_ticklabel_fs'] =  14
            #props['ax_tickwidth'] =  1
            #props['ax_linewidth'] =  1

        props['aspect'] = 0.0

        props['axes_color'] = 'gray'

        props['cb_linewidth'] = 1
        props['cb_label_fs'] = 10
        props['cb_ticklabel_fs'] = 8
        props['cb_tickwidth'] = 0.25

        props['cb_edgecolor'] = 'gray'

        props['cb_shrink'] = 0  # possible?
        props['cb_pad'] = 0.    # possible?

        #props['mode_index_label_fs'] = user_prefs['all_plots']['mode_index_label_fs']
        props['mode_index_label_fs'] = user_prefs.mode_index_label_fs


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
