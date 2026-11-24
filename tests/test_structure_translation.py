"""Tests for structure-level translation during mesh generation."""

import sys
from pathlib import Path

import numpy as np
import pytest
from matplotlib.figure import Figure

sys.path.append(str(Path(__file__).resolve().parents[1] / "backend"))

import numbat
from numbat import materials, runtime
from numbat.geometry.structure import Structure
from numbat.mesh.fem import FemMesh
from numbat.mesh.gmsh_io import read_msh2_mesh
from numbat.mesh.mail import MailData
from numbat.mesh import templates
from numbat.plotting.modes import PlotParams2D
from numbat.solvers.modecalcs import ACSimulation


@pytest.fixture
def nbapp(tmp_path, monkeypatch):
    monkeypatch.setattr(numbat._SingletonMeta, "_instances", {})
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    monkeypatch.setattr(runtime, "path_scratch", lambda: scratch)
    prefs = tmp_path / "numbat.toml"
    prefs.write_text("")
    return numbat.NumBATApp(
        "translation", outdir=str(tmp_path / "out"), user_settings_file=str(prefs)
    )


def build_rectangular(app, **kwargs):
    return app.make_structure(
        "rectangular", 1500.0, 1000.0, 300.0, 280.0,
        material_bkg=materials.make_material("Vacuum"),
        material_a=materials.make_material("Si_2016_Smith"),
        lc_bkg=0.15, lc_refine_1=2.5, lc_refine_2=2.5,
        **kwargs,
    )


@pytest.mark.parametrize("translation", [
    None, (1,), (1, 2, 3), "ab", (1, "2"), (True, 2),
    (float("nan"), 0), (0, float("inf")), np.array(1.0), (10**400, 0),
])
def test_invalid_translations_are_rejected(translation):
    structure = Structure.__new__(Structure)
    with pytest.raises(ValueError, match="translation"):
        structure._clean_and_handle_args_and_parameters(
            "rectangular", 1500, 1000, lc_bkg=0.1,
            direct_call=False, translation=translation,
        )


def mesh_bounds(mesh):
    return np.array([mesh.min(axis=0), mesh.max(axis=0)])


@pytest.mark.parametrize("translation", [(200.0, -100.0), (-125.0, 230.0)])
def test_translated_rectangular_geometry_mesh_and_mail(nbapp, translation):
    original = build_rectangular(nbapp)
    with pytest.warns(UserWarning, match="point-based"):
        shifted = build_rectangular(nbapp, translation=translation)

    assert shifted.translation == translation
    assert "translation" not in shifted.all_params
    original_base = Path(original.mesh_mail_fname).with_suffix("")
    shifted_base = Path(shifted.mesh_mail_fname).with_suffix("")
    assert "Translate" not in original_base.with_suffix(".geo").read_text()
    assert "Translate" in shifted_base.with_suffix(".geo").read_text()

    normalised_shift = np.asarray(translation) / original.domain_x
    mesh_a, triangles_a, _ = read_msh2_mesh(original_base.with_suffix(".msh"))
    mesh_b, triangles_b, _ = read_msh2_mesh(shifted_base.with_suffix(".msh"))
    np.testing.assert_allclose(
        mesh_bounds(mesh_b) - mesh_bounds(mesh_a),
        np.tile(normalised_shift, (2, 1)), atol=1e-12,
    )
    assert {material for _, material in triangles_a} == {
        material for _, material in triangles_b
    }

    mail_a, mail_b = MailData(original.mesh_mail_fname), MailData(shifted.mesh_mail_fname)
    bounds_a = mesh_bounds(np.column_stack((mail_a.v_x, mail_a.v_y)))
    bounds_b = mesh_bounds(np.column_stack((mail_b.v_x, mail_b.v_y)))
    np.testing.assert_allclose(bounds_b - bounds_a, np.tile(normalised_shift, (2, 1)), atol=1e-12)
    assert set(mail_a.v_elts_mat) == set(mail_b.v_elts_mat)

    fem_a, fem_b = FemMesh(), FemMesh()
    fem_a.build_from_gmsh_mail(original)
    fem_b.build_from_gmsh_mail(shifted)
    limits_a = np.array(fem_a.get_xy_limits()).reshape(2, 2)
    limits_b = np.array(fem_b.get_xy_limits()).reshape(2, 2)
    np.testing.assert_allclose(
        limits_b - limits_a, np.repeat(np.asarray(translation)[:, None] / 1000, 2, axis=1),
        atol=1e-12,
    )


def test_explicit_zero_translation_preserves_geo(nbapp):
    default = build_rectangular(nbapp)
    zero = build_rectangular(nbapp, translation=(0, 0))
    assert default.translation == zero.translation == (0, 0)
    assert Path(default.mesh_mail_fname).with_suffix(".geo").read_text() == (
        Path(zero.mesh_mail_fname).with_suffix(".geo").read_text()
    )


@pytest.mark.parametrize(
    "shape,parameters",
    [
        ("circular", {"inc_a_x": 300.0, "inc_a_y": 280.0}),
        ("onion1", {"inc_a_x": 300.0}),
        ("slot", {
            "rib_w": 150.0, "rib_h": 100.0, "slab_w": 800.0,
            "slab_h": 100.0, "slot_w": 80.0,
        }),
        ("rounded_wedge3", {
            "rad_wg": 500.0, "rad_tip": 10.0, "alpha_tip": 0.2,
            "gap_g": 20.0, "bowlen_l": 100.0, "rad_corner": 5.0,
        }),
        ("elliptical_rod", {"inc_a_x": 300.0, "inc_a_y": 280.0}),
    ],
)
def test_translation_preserves_geometry_and_material_regions(
    nbapp, monkeypatch, shape, parameters
):
    if shape == "elliptical_rod":
        monkeypatch.syspath_prepend(
            str(Path(__file__).resolve().parents[1] / "examples" / "user_waveguides")
        )
        from elliptical_rod import EllipticalRod

        monkeypatch.setattr(
            templates, "g_waveguide_templates",
            [*templates.g_waveguide_templates, {
                "inc_shape": shape, "wg_template_cls": EllipticalRod,
            }],
        )

    silicon = materials.make_material("Si_2016_Smith")
    args = {
        "domain_x": 2000.0, "domain_y": 1600.0,
        "material_bkg": materials.make_material("Vacuum"),
        "material_a": silicon, "material_b": silicon, "material_c": silicon,
        "lc_bkg": 0.15, "lc_refine_1": 2.0, "lc_refine_2": 2.0,
        **parameters,
    }
    original = nbapp.make_structure(shape, **args)
    with pytest.warns(UserWarning, match="point-based"):
        shifted = nbapp.make_structure(shape, translation=(300.0, -150.0), **args)

    def region_bounds(wguide):
        coords, triangles, _ = read_msh2_mesh(
            Path(wguide.mesh_mail_fname).with_suffix(".msh")
        )
        result = {}
        for nodes, material in triangles:
            result.setdefault(material, set()).update(nodes)
        return {
            material: mesh_bounds(coords[sorted(nodes)])
            for material, nodes in result.items()
        }

    original_regions, shifted_regions = region_bounds(original), region_bounds(shifted)
    assert original_regions.keys() == shifted_regions.keys()
    for material in original_regions:
        np.testing.assert_allclose(
            shifted_regions[material] - original_regions[material],
            np.tile((0.15, -0.075), (2, 1)), atol=1e-10,
        )


def test_translated_mesh_reaches_em_and_both_ac_paths(nbapp):
    with pytest.warns(UserWarning, match="point-based"):
        structure = build_rectangular(nbapp, translation=(200.0, -100.0))

    em = structure.calc_EM_modes(20, wl_um=1.55, n_eff=3.38)
    em_ac = ACSimulation(structure, q_AC=1e7, simres_EM=em, shift_Hz=10e9)
    direct_ac = ACSimulation(structure, q_AC=1e7, shift_Hz=10e9)
    direct_mesh, em_mesh = direct_ac.fem_mesh, em_ac.fem_mesh
    np.testing.assert_allclose(
        direct_mesh.node_coordinates, em_mesh.node_coordinates, rtol=1e-12, atol=1e-18
    )
    np.testing.assert_allclose(
        direct_mesh.get_xy_limits(), (50e-9, 350e-9, -240e-9, 40e-9),
        rtol=1e-12, atol=1e-18,
    )
    np.testing.assert_array_equal(
        direct_mesh.element_node_connectivity, em_mesh.element_node_connectivity
    )
    np.testing.assert_array_equal(direct_mesh.element_materials, em_mesh.element_materials)
    assert direct_mesh.element_convert_table == em_mesh.element_convert_table

    em_ac.do_main_eigensolve(None)
    direct_ac.do_main_eigensolve(None)
    np.testing.assert_allclose(
        np.sort(np.real(em_ac.eigs_nu)), np.sort(np.real(direct_ac.eigs_nu)),
        rtol=1e-9, atol=1e-3,
    )

    direct_ac.make_result()
    for mode, is_strain in (
        (em.get_mode(0), False), (direct_ac.get_sim_result().get_mode(0), True)
    ):
        ax = Figure().subplots()
        params = PlotParams2D()
        params["colorbar"] = False
        params["logamp"] = False
        if is_strain:
            mode.plot_strain(comps="xx", ax=ax, n_pts=31, plot_params=params)
        else:
            mode.plot_mode(comps="x", ax=ax, n_pts=31, plot_params=params)
        patch = ax.patches[-1]
        vertices = ax.transData.inverted().transform(
            patch.get_transform().transform([[0, 0], [1, 1]])
        )
        np.testing.assert_allclose(
            vertices, [[0.05, -0.24], [0.35, 0.04]], atol=1e-12
        )
