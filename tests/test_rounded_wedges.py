import re
import shutil
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np
import pytest

sys.path.append("backend")

from msh.builtin_meshes import RoundedWedge3, RoundedWedge4, RoundedWedge6

WEDGE_CLASSES = [(RoundedWedge3, 3, 90.0, 60.0), (RoundedWedge4, 4, 45.0, 45.0), (RoundedWedge6, 6, 0.0, 30.0)]

BASE_PARAMS = {
    "domain_x": 2500.0,
    "domain_y": 2500.0,
    "rad_wg": 1000.0,
    "rad_tip": 10.0,
    "alpha_tip": 0.2,
    "gap_g": 20.0,
    "bowlen_l": 100.0,
    "rad_corner": 5.0,
    "lc_bkg": 0.1,
    "lc_refine_1": 1.0,
    "lc_refine_2": 2.0,
}

POINT_RE = re.compile(r"Point\((\d+)\) = \{([^,]+), ([^,]+), 0, [^}]+\};")


def make_wedge(cls, **overrides):
    geom = cls({**BASE_PARAMS, **overrides}, {})
    geom.init_geometry()
    return geom


def geo_points(geo):
    return {int(m.group(1)): (float(m.group(2)), float(m.group(3))) for m in POINT_RE.finditer(geo)}


def rotate(points, angle):
    c, s = np.cos(angle), np.sin(angle)
    return [(c * x - s * y, s * x + c * y) for x, y in points]


@pytest.mark.parametrize("cls, n, phi0_deg, alpha_limit_deg", WEDGE_CLASSES)
def test_wedges_are_centred_at_expected_angles(cls, n, phi0_deg, alpha_limit_deg):
    geo = make_wedge(cls).make_geometry("")
    pts = geo_points(geo)
    for k in range(n):
        tip_x, tip_y = pts[10 + 9 * k + 6]  # tip centre of wedge k
        angle = np.degrees(np.arctan2(tip_y, tip_x)) % 360
        expected = (phi0_deg + k * 360.0 / n) % 360
        assert np.isclose(angle, expected, atol=1e-6) or np.isclose(abs(angle - expected), 360, atol=1e-6)


@pytest.mark.parametrize("cls, n, phi0_deg, alpha_limit_deg", WEDGE_CLASSES)
def test_wedges_are_rotations_of_one_wedge(cls, n, phi0_deg, alpha_limit_deg):
    pts = geo_points(make_wedge(cls).make_geometry(""))
    wedge0 = [pts[10 + j] for j in range(9)]
    for k in range(1, n):
        expected = rotate(wedge0, 2 * np.pi * k / n)
        actual = [pts[10 + 9 * k + j] for j in range(9)]
        assert np.allclose(actual, expected, atol=1e-10)


@pytest.mark.parametrize("cls, n, phi0_deg, alpha_limit_deg", WEDGE_CLASSES)
def test_rot_deg_rotates_whole_structure(cls, n, phi0_deg, alpha_limit_deg):
    base = geo_points(make_wedge(cls).make_geometry(""))
    turned = geo_points(make_wedge(cls, rot_deg=17.0).make_geometry(""))
    for pid, (x, y) in base.items():
        if pid in (1, 2, 3, 4, 5):  # outer boundary is rotation invariant
            continue
        (rx, ry), = rotate([(x, y)], np.radians(17.0))
        assert np.allclose(turned[pid], (rx, ry), atol=1e-10)


@pytest.mark.parametrize("cls, n, phi0_deg, alpha_limit_deg", WEDGE_CLASSES)
def test_alpha_tip_limit_enforced(cls, n, phi0_deg, alpha_limit_deg):
    ok, _ = make_wedge(cls, alpha_tip=np.radians(alpha_limit_deg) * 0.9).check_dimensions()
    assert ok

    ok, msg = make_wedge(cls, alpha_tip=np.radians(alpha_limit_deg) * 1.05).check_dimensions()
    assert not ok
    assert f"{alpha_limit_deg:g} degrees" in msg


@pytest.mark.parametrize("cls, n, phi0_deg, alpha_limit_deg", WEDGE_CLASSES)
def test_material_count_is_wedge_plus_background(cls, n, phi0_deg, alpha_limit_deg):
    geo = make_wedge(cls).make_geometry("")
    assert "Physical Surface(1) = {111};" in geo
    wedge_surfaces = re.search(r"Physical Surface\(2\) = \{([^}]*)\}", geo).group(1)
    assert len(wedge_surfaces.split(",")) == n


@pytest.mark.parametrize("cls, n, phi0_deg, alpha_limit_deg", WEDGE_CLASSES)
def test_draw_mpl_frame_outlines_every_wedge(cls, n, phi0_deg, alpha_limit_deg):
    fig, ax = plt.subplots()
    try:
        make_wedge(cls).draw_mpl_frame(ax, {"edgecolor": "k", "linewidth": 1.0})
        assert len(ax.lines) == n
    finally:
        plt.close(fig)


@pytest.mark.skipif(shutil.which("gmsh") is None, reason="gmsh executable not available")
@pytest.mark.parametrize("cls, n, phi0_deg, alpha_limit_deg", WEDGE_CLASSES)
def test_generated_geometry_meshes_in_gmsh(cls, n, phi0_deg, alpha_limit_deg, tmp_path):
    geo_file = tmp_path / "wedges.geo"
    geo_file.write_text(make_wedge(cls, rot_deg=10.0).make_geometry(""))
    result = subprocess.run(
        ["gmsh", "-2", str(geo_file), "-o", str(tmp_path / "wedges.msh")],
        check=False,
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "wedges.msh").is_file()
