"""Tests for passing structure translations to wireframe hooks."""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

sys.path.append(str(Path(__file__).resolve().parents[1] / "backend"))

from numbat.plotting.modes import Decorator


@pytest.fixture
def plot_preferences(monkeypatch):
    preferences = SimpleNamespace(
        WG_FRAME_EDGE_COLOR="red",
        WG_FRAME_LINEWIDTH_WHOLEFIG=2,
        WG_FRAME_LINEWIDTH_SUBFIG=1,
    )
    monkeypatch.setattr("numbat.plotting.modes.runtime.plot_prefs", lambda: preferences)


def points_in_data(ax, artist, points):
    return ax.transData.inverted().transform(artist.get_transform().transform(points))


@pytest.mark.parametrize("translation", [(200, -100), (-75, 120), (0, 0)])
def test_wireframe_hook_receives_delta_and_applies_it(plot_preferences, translation):
    fig = Figure()
    ax = fig.subplots()
    existing, = ax.plot([2, 3], [4, 5])
    existing_data = existing.get_xydata().copy()
    received = {}
    artists = {}

    class Geometry:
        def draw_mpl_frame(self, axis, styles, delta_xy):
            received["delta_xy"] = delta_xy
            dx, dy = delta_xy
            patch = Rectangle((dx, dy), 1, 2, edgecolor=styles["edgecolor"])
            axis.add_patch(patch)
            line, = axis.plot(
                np.array([0, 1]) + dx,
                np.array([1, 2]) + dy,
                color=styles["edgecolor"],
            )
            artists.update(patch=patch, line=line)

    decorator = Decorator()
    decorator.configure_waveguide(
        SimpleNamespace(wg_geom=Geometry(), translation=translation)
    )
    decorator.add_waveguide_border(ax)

    delta_um = np.asarray(translation) / 1000
    np.testing.assert_allclose(received["delta_xy"], delta_um)
    np.testing.assert_allclose(
        points_in_data(ax, artists["patch"], [[0, 0]]), [delta_um], atol=1e-12
    )
    np.testing.assert_allclose(
        artists["line"].get_xydata(),
        np.array([[0, 1], [1, 2]]) + delta_um,
        atol=1e-12,
    )
    np.testing.assert_array_equal(existing.get_xydata(), existing_data)
    assert artists["line"].get_color() == "red"


def test_hook_can_translate_collection_artists(plot_preferences):
    fig = Figure()
    ax = fig.subplots()

    class Geometry:
        def draw_mpl_frame(self, axis, styles, delta_xy):
            dx, dy = delta_xy
            axis.scatter([dx], [dy], color=styles["edgecolor"])

    decorator = Decorator()
    decorator.configure_waveguide(
        SimpleNamespace(wg_geom=Geometry(), translation=(250, -125))
    )
    decorator.add_waveguide_border(ax)

    np.testing.assert_allclose(
        ax.collections[0].get_offsets(), [[0.25, -0.125]], atol=1e-12
    )


def test_direct_set_waveguide_clears_previous_structure_translation(plot_preferences):
    fig = Figure()
    ax = fig.subplots()
    received = {}

    class Geometry:
        def draw_mpl_frame(self, axis, styles, delta_xy):
            received["delta_xy"] = delta_xy

    geometry = Geometry()
    decorator = Decorator()
    decorator.configure_waveguide(
        SimpleNamespace(wg_geom=geometry, translation=(200, -100))
    )
    decorator.set_waveguide(geometry)
    decorator.add_waveguide_border(ax)

    assert received["delta_xy"] == (0.0, 0.0)
