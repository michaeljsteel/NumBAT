"""Rendering orchestration for sampled acoustic strain fields."""

import matplotlib.pyplot as plt
import numpy as np

from numbat.plotting import tools as plottools
from numbat.plotting.modes import resolve_plot_trim
from numbat.runtime import outdir_fields_path, plotfile_ext
from numbat.strain_selection import (
    StrainComponent,
    StrainSelection,
    storage_key,
    tex_label,
)

_OVERVIEW_SELECTIONS = tuple(
    StrainSelection.primary(component) for component in StrainComponent
)
_STRAIN_DISPLAY_SCALE = 1e-6
_STRAIN_DISPLAY_UNIT = r"$\mathrm{M}\epsilon$"


def _apply_trim_limits(ax, plot_params):
    xlim_min, ylim_min, xlim_max, ylim_max = resolve_plot_trim(
        plot_params, default_trim=-0.1
    )

    if xlim_min != 0 or xlim_max != 0:
        xmin, xmax = ax.get_xlim()
        width = xmax - xmin
        ax.set_xlim(xmin + xlim_min * width, xmax - xlim_max * width)
    if ylim_min != 0 or ylim_max != 0:
        ymin, ymax = ax.get_ylim()
        height = ymax - ymin
        ax.set_ylim(ymin + ylim_min * height, ymax - ylim_max * height)


def _strain_plot_filename(plot_params, mode_index, label=""):
    output_directory = outdir_fields_path(prefix=plot_params["prefix"])
    suffix = f"_{plot_params['suffix']}" if plot_params["suffix"] else ""
    component = f"_{label}" if label else ""
    return (
        f"{output_directory}/AC_strain_mode_{mode_index:02d}"
        f"{suffix}{component}{plotfile_ext()}"
    )


def _strain_cut_filename(plot_params, mode_index, cut):
    output_directory = outdir_fields_path(prefix=plot_params["prefix"])
    suffix = f"_{plot_params['suffix']}" if plot_params["suffix"] else ""
    return (
        f"{output_directory}/AC_strain_mode_{mode_index:02d}_{cut}cut"
        f"{suffix}{plotfile_ext()}"
    )


def _add_mode_label(fig, mode_index):
    fig.text(
        0.01,
        0.99,
        f"Mode {mode_index}\nStrain [{_STRAIN_DISPLAY_UNIT}]",
        horizontalalignment="left",
        verticalalignment="top",
    )


def _add_colorbar(fig, ax, image):
    colorbar_axis = ax.inset_axes([1.04, 0.1, 0.05, 0.8])
    return fig.colorbar(image, cax=colorbar_axis)


def _plot_component(ax, coordinates, fields, selection, plot_params):
    values = fields[storage_key(selection)] * _STRAIN_DISPLAY_SCALE
    finite_values = values[np.isfinite(values)]
    limit = np.max(np.abs(finite_values)) if len(finite_values) else 1.0
    if limit == 0:
        limit = 1.0

    v_x = coordinates["v_x"]
    v_y = coordinates["v_y"]
    image = ax.imshow(
        values.T,
        origin="lower",
        extent=(v_x[0], v_x[-1], v_y[0], v_y[-1]),
        cmap="seismic",
        vmin=-limit,
        vmax=limit,
    )
    ax.set_aspect("equal")
    _apply_trim_limits(ax, plot_params)
    if plot_params["ticks"]:
        ax.set_xlabel("$x$ [μm]")
        ax.set_ylabel("$y$ [μm]")
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    if plot_params["title"]:
        ax.set_title(tex_label(selection))
    if plot_params["frame"]:
        plot_params["decorator"].add_waveguide_border(ax)
    return image


def render_strain_2d(
    *,
    coordinates,
    fields,
    components,
    simulation_result,
    mode_index,
    plot_params,
    axis=None,
):
    """Render a six-panel strain overview or individual strain components."""
    if axis is not None and len(components) != 1:
        raise ValueError(
            "When providing an axis, specify exactly one strain component."
        )

    decorator = plot_params["decorator"]
    decorator.set_waveguide(simulation_result._structure.wg_geom)

    if not components:
        decorator.set_for_multi()
        fig, axes = plt.subplots(2, 3, figsize=decorator["figsize"])
        _add_mode_label(fig, mode_index)
        images = [
            _plot_component(ax, coordinates, fields, selection, plot_params)
            for ax, selection in zip(axes.flat, _OVERVIEW_SELECTIONS)
        ]
        if plot_params["colorbar"]:
            for ax, image in zip(axes.flat, images):
                _add_colorbar(fig, ax, image)
        fig.tight_layout()
        plottools.save_and_close_figure(
            fig, _strain_plot_filename(plot_params, mode_index)
        )
        return

    decorator.set_for_single()
    for component in components:
        selection = StrainSelection.primary(component)
        if axis is None:
            fig, component_axis = plt.subplots(figsize=decorator["figsize"])
        else:
            fig, component_axis = axis.figure, axis
        _add_mode_label(fig, mode_index)
        image = _plot_component(
            component_axis, coordinates, fields, selection, plot_params
        )
        if axis is None and plot_params["colorbar"]:
            _add_colorbar(fig, component_axis, image)
        if axis is None:
            plottools.save_and_close_figure(
                fig,
                _strain_plot_filename(
                    plot_params, mode_index, label=selection.component.value
                ),
            )


def render_strain_1d(
    *,
    x_axis,
    fields,
    components,
    x_label,
    cut,
    mode_index,
    plot_params,
):
    """Render a sampled 1D strain cut using display-primary components."""
    if not components:
        components = tuple(StrainComponent)

    fig, axis = plt.subplots()
    for component in components:
        selection = StrainSelection.primary(component)
        axis.plot(
            x_axis,
            fields[storage_key(selection)] * _STRAIN_DISPLAY_SCALE,
            label=tex_label(selection),
        )

    axis.set_xlabel(x_label)
    axis.set_ylabel(f"Strain [{_STRAIN_DISPLAY_UNIT}]")
    axis.legend()
    plottools.save_and_close_figure(
        fig, _strain_cut_filename(plot_params, mode_index, cut)
    )
