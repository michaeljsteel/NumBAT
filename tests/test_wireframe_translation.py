"""Tests for translating user-provided wireframe artists."""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.transforms import Affine2D

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
def test_new_lines_and_patches_are_shifted_in_data_coordinates(
    plot_preferences, translation
):
    fig = Figure()
    ax = fig.subplots()
    existing, = ax.plot([2, 3], [4, 5])
    existing_transform = existing.get_transform()
    artists = {}

    class Geometry:
        def draw_mpl_frame(self, axis, styles):
            local = Affine2D().rotate_deg(30).translate(0.7, -0.3)
            patch = Rectangle((0.2, 0.1), 0.4, 0.2, transform=local + axis.transData)
            axis.add_patch(patch)
            line, = axis.plot([0, 1], [1, 2], color=styles["edgecolor"])
            artists.update(patch=patch, line=line, local=local)

    decorator = Decorator()
    decorator.configure_waveguide(SimpleNamespace(wg_geom=Geometry(), translation=translation))
    decorator.add_waveguide_border(ax)

    shift = np.asarray(translation) / 1000
    expected_patch = artists["local"].transform(
        artists["patch"].get_patch_transform().transform([[0, 0], [1, 1]])
    ) + shift
    np.testing.assert_allclose(
        points_in_data(ax, artists["patch"], [[0, 0], [1, 1]]),
        expected_patch, atol=1e-12,
    )
    np.testing.assert_allclose(
        points_in_data(ax, artists["line"], [[0, 1], [1, 2]]),
        np.array([[0, 1], [1, 2]]) + shift, atol=1e-12,
    )
    assert existing.get_transform() is existing_transform
    assert artists["line"].get_color() == "red"
    assert artists["line"].get_clip_on()
    assert artists["patch"].get_clip_on()


def test_unsupported_artist_fails_explicitly_for_translated_structure(plot_preferences):
    fig = Figure()
    ax = fig.subplots()

    class Geometry:
        def draw_mpl_frame(self, axis, styles):
            axis.scatter([0], [0])

    decorator = Decorator()
    decorator.configure_waveguide(SimpleNamespace(wg_geom=Geometry(), translation=(1, 0)))
    with pytest.raises(TypeError, match="PathCollection"):
        decorator.add_waveguide_border(ax)


def test_non_data_artist_transform_fails_explicitly(plot_preferences):
    fig = Figure()
    ax = fig.subplots()

    class Geometry:
        def draw_mpl_frame(self, axis, styles):
            axis.plot([0, 1], [0, 1], transform=axis.transAxes)

    decorator = Decorator()
    decorator.configure_waveguide(SimpleNamespace(wg_geom=Geometry(), translation=(1, 0)))
    with pytest.raises(ValueError, match="data coordinates"):
        decorator.add_waveguide_border(ax)


def test_legacy_structure_without_translation_preserves_custom_artists(plot_preferences):
    fig = Figure()
    ax = fig.subplots()

    class Geometry:
        def draw_mpl_frame(self, axis, styles):
            axis.scatter([0], [0])

    decorator = Decorator()
    decorator.configure_waveguide(SimpleNamespace(wg_geom=Geometry()))
    decorator.add_waveguide_border(ax)
    assert len(ax.collections) == 1


def test_direct_set_waveguide_clears_previous_structure_translation(plot_preferences):
    fig = Figure()
    ax = fig.subplots()

    class Geometry:
        def draw_mpl_frame(self, axis, styles):
            axis.plot([0, 1], [1, 2])

    geometry = Geometry()
    decorator = Decorator()
    decorator.configure_waveguide(SimpleNamespace(wg_geom=geometry, translation=(200, -100)))
    decorator.set_waveguide(geometry)
    decorator.add_waveguide_border(ax)
    np.testing.assert_allclose(
        points_in_data(ax, ax.lines[0], [[0, 1], [1, 2]]),
        [[0, 1], [1, 2]], atol=1e-12,
    )
