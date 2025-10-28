
from collections.abc import Iterable
from typing import Optional, Union
from matplotlib.axes import Axes
from matplotlib.colorbar import Colorbar
import numbat

class TidyAxes:
    """Utility class for applying consistent styling to matplotlib axes and colorbars.

    This class manages plot aesthetics including font sizes, tick properties,
    axis line widths, and colorbar formatting based on the number of axes in a figure.
    """

    def __init__(self, nax: int = 1, props: Optional[dict] = None) -> None:
        """Initialize TidyAxes with default or custom properties.

        Args:
            nax: Number of axes in the figure. Different defaults apply for 1, 2, 4, or 6 axes.
            props: Optional dictionary of custom properties to override defaults.
        """

        self._nax = nax
        self._set_defaults_for_num_axes(nax)
        if props:
            self._props.update(props)

    def update_property(self, k: str, v) -> None:
        """Update a single property in the properties dictionary.

        Args:
            k: Property key to update.
            v: New value for the property.
        """
        self._props[k]=v

    def _set_defaults_for_num_axes(self, nax: int) -> None:
        """Set default properties based on the number of axes in the figure.

        Different property sets are used for single plots (nax=1), dual plots (nax=2),
        and multi-panel figures (nax=4 or 6).

        Args:
            nax: Number of axes in the figure.
        """
        user_prefs = numbat.NumBATPlotPrefs()

        props = {}

        if nax in (4,6):  # includes layout for two row mode plots
            props.update(user_prefs.d_multi)
            props['ax_label_xpad'] = 3
            props['ax_label_ypad'] = 1

        elif nax == 2:
            props.update(user_prefs.d_multi)
            props['ax_label_xpad'] = 3
            props['ax_label_ypad'] = 1
            props['ax_label_pad'] = 5

        else:
            props.update(user_prefs.d_single)
            props['ax_label_pad'] = 5
            props['ax_label_xpad'] = 3
            props['ax_label_ypad'] = 1

        props['aspect'] = 0.0

        props['axes_color'] = 'gray'

        #props['cb_linewidth'] = 1
        #props['cb_label_fs'] = 10
        #props['cb_ticklabel_fs'] = 8
        #props['cb_tickwidth'] = 0.25
        #props['cb_edgecolor'] = 'gray'

        props['cb_shrink'] = 0  # possible?
        props['cb_pad'] = 0.    # possible?

        props['mode_index_label_fs'] = user_prefs.mode_index_label_fs


        self._props = props

    def apply_to_axes(self, axs: Union[Axes, Iterable[Axes]]) -> None:
        """Apply styling properties to one or more matplotlib axes.

        Sets tick parameters, axis label sizes and padding, aspect ratios,
        and spine properties (line width and color) for all provided axes.

        Args:
            axs: A single Axes object or an iterable of Axes objects to style.
        """

        if not isinstance(axs, Iterable):
            axs = (axs,)

        pr = self._props

        for ax in axs:
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

    def hide_axes(self, ax: Axes) -> None:
        """Turn off all axis visibility for the given axes.

        Args:
            ax: The Axes object to hide.
        """
        ax.set_axis_off()

    def apply_to_cbars(self, cbs: Union[Colorbar, Iterable[Colorbar]]) -> None:
        """Apply styling properties to one or more matplotlib colorbars.

        Sets colorbar label size, outline properties (line width and color),
        and tick parameters for all provided colorbars.

        Args:
            cbs: A single Colorbar object or an iterable of Colorbar objects to style.
        """
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
