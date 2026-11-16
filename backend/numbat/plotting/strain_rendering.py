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


def _plot_component(ax, coordinates, fields, selection, plot_params):
    values = fields[storage_key(selection)]
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
        images = [
            _plot_component(ax, coordinates, fields, selection, plot_params)
            for ax, selection in zip(axes.flat, _OVERVIEW_SELECTIONS)
        ]
        if plot_params["colorbar"]:
            for ax, image in zip(axes.flat, images):
                fig.colorbar(image, ax=ax)
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
        image = _plot_component(
            component_axis, coordinates, fields, selection, plot_params
        )
        if axis is None and plot_params["colorbar"]:
            fig.colorbar(image, ax=component_axis)
        if axis is None:
            plottools.save_and_close_figure(
                fig,
                _strain_plot_filename(
                    plot_params, mode_index, label=selection.component.value
                ),
            )
