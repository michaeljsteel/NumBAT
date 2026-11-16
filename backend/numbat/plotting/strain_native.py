"""Diagnostic rendering for native element-local acoustic strain."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.tri import Triangulation


def plot_native_strain_component(samples, selection, ax=None):
    """Render one native strain component without averaging interface values."""
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    values = samples.display_values(selection)
    finite_values = values[np.isfinite(values)]
    if len(finite_values) == 0:
        raise ValueError("Cannot render a strain component with no finite values.")
    limit = np.max(np.abs(finite_values))
    if limit == 0:
        limit = 1.0

    coordinates = samples.flattened_coordinates()
    triangulation = Triangulation(
        coordinates[:, 0], coordinates[:, 1], samples.subtriangles()
    )
    image = ax.tripcolor(
        triangulation,
        values,
        shading="gouraud",
        cmap="seismic",
        vmin=-limit,
        vmax=limit,
    )
    ax.set_aspect("equal")
    ax.set_title(selection.component.value)
    colorbar = fig.colorbar(image, ax=ax)
    return fig, ax, colorbar
