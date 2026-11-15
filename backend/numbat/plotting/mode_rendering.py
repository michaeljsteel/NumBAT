"""High-level rendering of sampled mode fields."""

import matplotlib.pyplot as plt

from numbat import reporting
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


def render_mode_2d(
    *,
    field_type,
    coordinates,
    fields,
    include_overview,
    components,
    simulation_result,
    mode_index,
    plot_params,
    axis=None,
):
    """Render selected 2D mode components and optionally their overview."""
    if axis is not None and (include_overview or len(components) != 1):
        reporting.report_and_exit(
            "When providing an axis to plot on, specify exactly one individual "
            "modal component and do not request the 'e' overview component."
        )

    decorator = plot_params["decorator"]
    if include_overview:
        decorator.set_for_multi()
        mode_presentation.plot_all_components(
            field_type,
            coordinates,
            fields,
            plot_params,
            simulation_result,
            mode_index,
        )

    if components:
        decorator.set_for_single()
        for component in components:
            selection = mode_presentation.selection_for_component(
                field_type, component
            )
            mode_presentation.plot_one_component(
                field_type,
                coordinates,
                fields,
                plot_params,
                mode_index,
                selection,
                axis,
            )
