"""Regression tests for acoustic mesh coordinate units."""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1] / "backend"))

from numbat.mesh.fem import FemMesh
from numbat.nbtypes import SI_nm


class ElasticProperties:
    typ_el_AC = {2: 1}
    n_mats_ac = 1
    v_acoustic_mats = []

    def is_elastic_material_index(self, index):
        return index == 2

    def active_material_index(self, index):
        return self.typ_el_AC[index]


def test_direct_ac_mesh_matches_em_mesh_in_metres(tmp_path):
    mesh_file = tmp_path / "mesh.mail"
    mesh_file.write_text(
        "9 2\n"
        "1 0 0 0\n"
        "2 0.2 0 0\n"
        "3 0 0.1 0\n"
        "4 0.1 0 0\n"
        "5 0.1 0.05 0\n"
        "6 0 0.05 0\n"
        "7 0.4 0 0\n"
        "8 0.3 0 0\n"
        "9 0.2 0.05 0\n"
        "1 1 2 3 4 5 6 2\n"
        "2 2 7 3 8 9 5 1\n"
    )
    waveguide = SimpleNamespace(
        mesh_mail_fname=str(mesh_file),
        domain_x=1500.0,
        optical_props=SimpleNamespace(n_mats_em=0),
        get_optical_materials=lambda: [],
        elastic_props=ElasticProperties(),
    )

    initial_mesh = FemMesh()
    initial_mesh.build_from_gmsh_mail(waveguide)
    np.testing.assert_allclose(initial_mesh.node_coordinates[:, 1], [0.3, 0.0])

    # Emulate the Fortran EM mesh output from the same normalized Mail coordinates.
    em_mesh = FemMesh()
    em_mesh.store_fortran_em_mesh_properties(
        initial_mesh.element_materials,
        np.zeros(initial_mesh.num_nodes, dtype=np.int64),
        initial_mesh.element_node_connectivity,
        initial_mesh.node_coordinates * (SI_nm / 0.001),
    )
    from_em = FemMesh()
    from_em.ac_build_from_em(waveguide, em_mesh)
    direct = FemMesh()
    direct.ac_build_from_gmsh_mail(waveguide)

    assert (direct.num_nodes, direct.num_elements) == (6, 1)
    np.testing.assert_allclose(direct.node_coordinates, from_em.node_coordinates, rtol=1e-14, atol=1e-20)
    np.testing.assert_allclose(direct.node_coordinates[:, 1], [3e-7, 0.0], rtol=1e-14, atol=1e-20)
    np.testing.assert_allclose(direct._extents, from_em._extents, rtol=1e-14, atol=1e-20)
    np.testing.assert_array_equal(direct.element_node_connectivity, from_em.element_node_connectivity)
    np.testing.assert_array_equal(direct.element_materials, from_em.element_materials)
    assert direct.element_convert_table == from_em.element_convert_table == {0: 0}
