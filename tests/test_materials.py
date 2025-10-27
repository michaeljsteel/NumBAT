import sys
sys.path.append('backend')
import copy

import numpy as np

import materials
from math import pi

mat_iso='SiO2_2023_Steel'
mat_cubic='Si_1970_Auld'

def test_isotropic_material_is_rotation_insensitive():
    mat = materials.make_material(mat_iso)

    c0 = mat.stiffness_c_IJ.value().copy()

    mat.rotate((1,2,3), np.pi*.234)
    c1 = mat.stiffness_c_IJ.value().copy()
    atol = 1e-10 * np.max(np.abs(c0))

    #with np.printoptions(precision=3):
    #    print(c0)
    #    print(c1)
    #    print(c1-c0)
    np.testing.assert_allclose(c0, c1, atol=atol)

def test_cubic_material_satisfies_symmetries():
    mat = materials.make_material(mat_cubic)

    c0 = mat.stiffness_c_IJ.value().copy()
    atol = 1e-10 * np.max(np.abs(c0))
    print(c0)

    for nvec,th in (
            ((1,0,0),pi/2),
            ((0,1,0),pi/2),
            ((0,0,1),pi/2),

            ((1,1,1),2*pi/3),
            ((1,-1,1),2*pi/3),
            ((-1,-1,1),2*pi/3),

            ((1,1,0),pi),
            ((1,0,1),pi),
            ((0,1,1),pi),

            ):

        tmat = copy.deepcopy(mat)
        tmat.rotate(nvec, th)
        c1 = tmat.stiffness_c_IJ.value().copy()

        np.testing.assert_allclose(c0, c1, atol=atol)



