import json
import pickle
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

sys.path.append("backend")

from meshing import templates
from meshing.templates import _load_waveguide_templates
from msh.builtin_meshes import RoundedBowtie, RoundedWedge

PROJECT_ROOT = Path(__file__).resolve().parents[1]
BUILTIN_WAVEGUIDE_DIRECTORY = PROJECT_ROOT / "backend" / "numbat_core" / "mesh" / "builtin"
BUILTIN_TEMPLATE_DIRECTORY = BUILTIN_WAVEGUIDE_DIRECTORY / "geo"
EXTERNAL_TEMPLATE_DIRECTORY = PROJECT_ROOT / "examples" / "user_waveguides"
ROUNDED_BOWTIE_PARAMETERS = {
    "domain_x": 3000,
    "domain_y": 3000,
    "rad_tip": 10,
    "alpha_tip": 3.141592653589793 / 6,
    "gap_g": 20,
    "bowlen_l": 100,
    "rad_corner": 5,
    "lc_bkg": 0.1,
    "lc_refine_1": 1,
    "lc_refine_2": 2,
    "rad_wg": 1000,
}

def test_builtin_waveguide_index_uses_normal_import_path():
    waveguides = _load_waveguide_templates(
        BUILTIN_WAVEGUIDE_DIRECTORY / "builtin_waveguides.json"
    )

    assert waveguides
    assert all(waveguide["module"] == "msh.builtin_meshes" for waveguide in waveguides)
    assert all("wg_impl" not in waveguide for waveguide in waveguides)
    assert all(
        waveguide["wg_template_cls"].__module__.startswith("numbat_core.mesh.builtin")
        for waveguide in waveguides
    )


def test_internal_builtin_waveguide_module_path_is_supported(tmp_path):
    manifest_path = tmp_path / "numbat_waveguides.json"
    builtin_manifest = json.loads(
        (BUILTIN_WAVEGUIDE_DIRECTORY / "builtin_waveguides.json").read_text()
    )
    for waveguide in builtin_manifest["wguides"]:
        waveguide["module"] = "numbat_core.mesh.builtin"
    manifest_path.write_text(json.dumps(builtin_manifest))

    waveguides = _load_waveguide_templates(manifest_path)

    assert waveguides
    assert all(
        waveguide["wg_template_cls"].__module__.startswith("numbat_core.mesh.builtin")
        for waveguide in waveguides
    )


def test_module_based_user_waveguide_uses_neighbouring_geo_template_and_normal_import(
    tmp_path, monkeypatch
):
    index_path = tmp_path / "waveguides.json"
    index_path.write_text(
        json.dumps(
            {
                "wguides": [
                    {
                        "module": "elliptical_rod",
                        "wg_class": "EllipticalRod",
                        "inc_shape": "elliptical_rod",
                    }
                ]
            }
        )
    )
    monkeypatch.syspath_prepend(str(EXTERNAL_TEMPLATE_DIRECTORY))
    sys.modules.pop("elliptical_rod", None)

    waveguides = _load_waveguide_templates(index_path)
    geometry = waveguides[0]["wg_template_cls"](
        {
            "domain_x": 3000,
            "domain_y": 2000,
            "inc_a_x": 800,
            "inc_a_y": 400,
            "lc_bkg": 0.1,
            "lc_refine_1": 2,
            "lc_refine_2": 3,
        },
        {},
    )
    geometry.init_geometry()
    rendered_template = geometry.make_geometry(BUILTIN_TEMPLATE_DIRECTORY)

    assert type(geometry).__module__ == "elliptical_rod"
    assert geometry.gmsh_template_path(BUILTIN_TEMPLATE_DIRECTORY) == (
        EXTERNAL_TEMPLATE_DIRECTORY / "elliptical_rod_msh_template.geo"
    )
    assert "radius_x = a1 / (2 * dx_in_nm) * d;" in rendered_template
    assert "a1 = 800.000000;" in rendered_template
    assert "a1y = 400.000000;" in rendered_template

    gmsh = shutil.which("gmsh")
    if gmsh is not None:
        geometry_path = tmp_path / "elliptical_rod.geo"
        mesh_path = tmp_path / "elliptical_rod.msh"
        geometry_path.write_text(rendered_template)
        subprocess.run(
            [gmsh, "-2", str(geometry_path), "-format", "msh2", "-o", str(mesh_path)],
            check=True,
            capture_output=True,
            text=True,
        )

    serialized_geometry = pickle.dumps(geometry)
    sys.modules.pop("elliptical_rod", None)
    rebuilt_geometry = pickle.loads(serialized_geometry)

    assert type(rebuilt_geometry).__module__ == "elliptical_rod"
    assert rebuilt_geometry._d_params == geometry._d_params


def test_registered_waveguide_directory_loads_active_manifest_entries(monkeypatch):
    class Application:
        def path_mesh_templates(self):
            return BUILTIN_TEMPLATE_DIRECTORY

    monkeypatch.setattr(templates, "_registered_waveguide_directories", [])
    monkeypatch.setattr(templates, "g_waveguide_templates", [])
    monkeypatch.syspath_prepend(str(EXTERNAL_TEMPLATE_DIRECTORY))
    sys.modules.pop("elliptical_rod", None)

    templates.register_waveguide_directory(EXTERNAL_TEMPLATE_DIRECTORY)
    templates.initialise_waveguide_templates(Application())

    elliptical_rod = next(
        waveguide
        for waveguide in templates.g_waveguide_templates
        if waveguide["inc_shape"] == "elliptical_rod"
    )
    assert elliptical_rod["wg_template_cls"].__module__ == "elliptical_rod"


def test_register_waveguide_directory_rejects_incomplete_manifest_entries(tmp_path, monkeypatch):
    manifest_path = tmp_path / "numbat_waveguides.json"
    manifest_path.write_text(json.dumps({"wguides": [{"module": "elliptical_rod"}]}))
    monkeypatch.setattr(templates, "_registered_waveguide_directories", [])

    with pytest.raises(ValueError, match="missing: inc_shape, wg_class"):
        templates.register_waveguide_directory(tmp_path)


def test_waveguide_manifest_rejects_removed_wg_impl_entries(tmp_path):
    index_path = tmp_path / "waveguides.json"
    index_path.write_text(
        json.dumps(
            {
                "wguides": [
                    {
                        "wg_impl": "legacy_waveguide.py",
                        "wg_class": "LegacyWaveguide",
                        "inc_shape": "legacy_waveguide",
                    }
                ]
            }
        )
    )

    with pytest.raises(ValueError, match="missing: module"):
        _load_waveguide_templates(index_path)


def test_rounded_bowtie_is_registered_and_substitutes_its_half_angle():
    waveguides = _load_waveguide_templates(
        BUILTIN_WAVEGUIDE_DIRECTORY / "builtin_waveguides.json"
    )
    rounded_bowtie = next(
        waveguide for waveguide in waveguides if waveguide["inc_shape"] == "rounded_bowtie"
    )
    geometry = rounded_bowtie["wg_template_cls"](ROUNDED_BOWTIE_PARAMETERS, {})
    geometry.init_geometry()
    mesh_geometry = geometry.make_geometry(BUILTIN_TEMPLATE_DIRECTORY)

    assert isinstance(geometry, RoundedBowtie)
    assert "alpha_tip = 0.523599;" in mesh_geometry
    assert geometry.check_dimensions() == (True, "")


def test_rounded_bowtie_rejects_crossing_outer_corner_transitions():
    invalid_parameters = dict(ROUNDED_BOWTIE_PARAMETERS, bowlen_l=120, rad_corner=50)
    geometry = RoundedBowtie(invalid_parameters, {})

    dimensions_ok, message = geometry.check_dimensions()

    assert not dimensions_ok
    assert "outer corner transitions cross" in message


def test_rounded_bowtie_draws_mirrored_inclusion_outlines():
    geometry = RoundedBowtie(ROUNDED_BOWTIE_PARAMETERS, {})
    figure, axes = plt.subplots()
    try:
        geometry.draw_mpl_frame(axes, {"edgecolor": "black", "linewidth": 1})
        assert len(axes.lines) == 2
        right_outline = axes.lines[0].get_xydata()
        left_outline = axes.lines[1].get_xydata()
        assert right_outline.shape == left_outline.shape
        assert (right_outline[:, 0] == -left_outline[:, 0]).all()
        assert (right_outline[:, 1] == left_outline[:, 1]).all()
        assert np.isclose(right_outline[:, 0].min(), ROUNDED_BOWTIE_PARAMETERS["gap_g"] / 2000)
    finally:
        plt.close(figure)


def test_rounded_bowtie_template_creates_three_physical_surfaces(tmp_path):
    gmsh = shutil.which("gmsh")
    if gmsh is None:
        pytest.skip("Gmsh is not installed")

    geometry = RoundedBowtie(ROUNDED_BOWTIE_PARAMETERS, {})
    geometry.init_geometry()
    rendered_template = geometry.make_geometry(BUILTIN_TEMPLATE_DIRECTORY)
    geometry_path = tmp_path / "rounded_bowtie.geo"
    mesh_path = tmp_path / "rounded_bowtie.msh"
    geometry_path.write_text(rendered_template)
    subprocess.run(
        [gmsh, "-2", str(geometry_path), "-format", "msh2", "-o", str(mesh_path)],
        check=True,
        capture_output=True,
        text=True,
    )

    physical_surfaces = set()
    lines = mesh_path.read_text().splitlines()
    element_start = lines.index("$Elements") + 2
    element_end = lines.index("$EndElements")
    for element_line in lines[element_start:element_end]:
        fields = element_line.split()
        if fields[1] == "2":
            physical_surfaces.add(int(fields[3]))

    assert physical_surfaces == {1, 2, 3}


ROUNDED_WEDGE_PARAMETERS = {
    "domain_x": 3000,
    "domain_y": 3000,
    "rad_wg": 1000,
    "rad_tip": 10,
    "alpha_tip": 3.141592653589793 / 6,
    "bowlen_l": 100,
    "rad_corner": 5,
    "lc_bkg": 0.1,
    "lc_refine_1": 1,
    "lc_refine_2": 2,
}


def test_rounded_wedge_is_registered_and_substitutes_its_parameters():
    waveguides = _load_waveguide_templates(
        BUILTIN_WAVEGUIDE_DIRECTORY / "builtin_waveguides.json"
    )
    rounded_wedge = next(
        waveguide for waveguide in waveguides if waveguide["inc_shape"] == "rounded_wedge"
    )
    geometry = rounded_wedge["wg_template_cls"](ROUNDED_WEDGE_PARAMETERS, {})
    geometry.init_geometry()
    mesh_geometry = geometry.make_geometry(BUILTIN_TEMPLATE_DIRECTORY)

    assert isinstance(geometry, RoundedWedge)
    assert "alpha_tip = 0.523599;" in mesh_geometry
    assert "bowlen_l = 100.000000;" in mesh_geometry
    assert geometry.check_dimensions() == (True, "")


def test_rounded_wedge_rejects_corner_transitions_that_meet_on_axis():
    invalid_parameters = dict(ROUNDED_WEDGE_PARAMETERS, rad_corner=70)
    geometry = RoundedWedge(invalid_parameters, {})

    dimensions_ok, message = geometry.check_dimensions()

    assert not dimensions_ok
    assert "rad_corner" in message


def test_rounded_wedge_rejects_tip_overlapping_corner_transition():
    invalid_parameters = dict(ROUNDED_WEDGE_PARAMETERS, rad_tip=80, rad_corner=60)
    geometry = RoundedWedge(invalid_parameters, {})

    dimensions_ok, message = geometry.check_dimensions()

    assert not dimensions_ok
    assert "rounded tip" in message


def test_rounded_wedge_rejects_wedge_not_inside_fibre():
    invalid_parameters = dict(ROUNDED_WEDGE_PARAMETERS, rad_wg=100, bowlen_l=100)
    geometry = RoundedWedge(invalid_parameters, {})

    dimensions_ok, message = geometry.check_dimensions()

    assert not dimensions_ok
    assert "waveguide radius" in message


def test_rounded_wedge_draws_single_symmetric_outline():
    geometry = RoundedWedge(ROUNDED_WEDGE_PARAMETERS, {})
    figure, axes = plt.subplots()
    try:
        geometry.draw_mpl_frame(axes, {"edgecolor": "black", "linewidth": 1})
        assert len(axes.lines) == 1
        outline = axes.lines[0].get_xydata()
        assert np.isclose(outline[:, 0].min(), -outline[:, 0].max())
        assert np.isclose(outline[:, 1].min(), -ROUNDED_WEDGE_PARAMETERS["bowlen_l"] / 1000, atol=1e-3)
        apex_y = ROUNDED_WEDGE_PARAMETERS["bowlen_l"] - ROUNDED_WEDGE_PARAMETERS["rad_tip"] * (
            1 / np.sin(ROUNDED_WEDGE_PARAMETERS["alpha_tip"]) - 1
        )
        assert np.isclose(outline[:, 1].max(), apex_y / 1000, atol=1e-4)
    finally:
        plt.close(figure)


def test_rounded_wedge_template_creates_two_physical_surfaces(tmp_path):
    gmsh = shutil.which("gmsh")
    if gmsh is None:
        pytest.skip("Gmsh is not installed")

    geometry = RoundedWedge(ROUNDED_WEDGE_PARAMETERS, {})
    geometry.init_geometry()
    rendered_template = geometry.make_geometry(BUILTIN_TEMPLATE_DIRECTORY)
    geometry_path = tmp_path / "rounded_wedge.geo"
    mesh_path = tmp_path / "rounded_wedge.msh"
    geometry_path.write_text(rendered_template)
    subprocess.run(
        [gmsh, "-2", str(geometry_path), "-format", "msh2", "-o", str(mesh_path)],
        check=True,
        capture_output=True,
        text=True,
    )

    physical_surfaces = set()
    lines = mesh_path.read_text().splitlines()
    element_start = lines.index("$Elements") + 2
    element_end = lines.index("$EndElements")
    for element_line in lines[element_start:element_end]:
        fields = element_line.split()
        if fields[1] == "2":
            physical_surfaces.add(int(fields[3]))

    assert physical_surfaces == {1, 2}
