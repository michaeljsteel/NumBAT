"""High-level rendering of sampled mode fields."""

import matplotlib.pyplot as plt

from numbat.field_selection import storage_key
from numbat.nbtypes import ComplexPart, FieldComponent
from numbat.plotting import modes as mode_presentation
from numbat.plotting import tools as plottools


def selections_for_1d_plot(field_type, components, *, include_secondary):
    """Build ordered field selections for a 1D mode plot."""
    selections = [
        mode_presentation.selection_for_component(field_type, component)
        for component in components
    ]
    if not include_secondary:
        return selections

    secondary_selections = [
        selection.with_part(ComplexPart.SECONDARY)
        for selection in selections
        if selection.is_cartesian
    ]
    return selections + secondary_selections


def render_mode_1d(
    *,
    field_type,
    fields,
    x_axis,
    components,
    x_label,
    cut,
    mode_index,
    plot_params,
):
    """Render a sampled 1D mode cut and save the resulting figure."""
    logx = plot_params.get("logx", False)
    logy = plot_params.get("logy", False)

    if not components:
        components = (
            FieldComponent.X,
            FieldComponent.Y,
            FieldComponent.Z,
            FieldComponent.MAGNITUDE,
        )

    selections = selections_for_1d_plot(
        field_type,
        components,
        include_secondary=not plot_params["suppress_imimre"],
    )

    fig, axis = plt.subplots()
    show_magnitude = FieldComponent.MAGNITUDE in components
    for selection in selections:
        field_key = storage_key(selection)
        field_values = fields[field_key]
        if logy:
            field_values = abs(field_values)

        axis.plot(
            x_axis,
            field_values,
            label=mode_presentation.tex_label(selection),
            color=mode_presentation.line_color(selection),
            linestyle=mode_presentation.line_style(
                selection, show_magnitude=show_magnitude
            ),
        )

    if logx:
        axis.set_xscale("log")
    if logy:
        axis.set_yscale("log")

    axis.set_xlabel(x_label)
    axis.set_ylabel("fields")
    axis.legend()

    filename = mode_presentation.modeplot_filename_1D(
        field_type, plot_params, mode_index, cut, label=""
    )
    plottools.save_and_close_figure(fig, filename)
