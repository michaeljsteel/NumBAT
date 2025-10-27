
import sys
sys.path.append('backend')

import numpy as np


from voigt import PlainTensor2_ij, VoigtTensor3_iJ, VoigtTensor4_IJ
from voigt import make_rotation_matrix, make_rotation_Bond_matrix_M, make_rotation_Bond_matrix_N

twopi = 2 * np.pi

g_rot_vecs = ((1,0,0),(0,1,0),(0,0,1),
            (1,1,0),(1,0,1),(0,1,1), (1,1,1))

g_rot_orders = range(1,6)

# def test_plain_mat3x3_rotation():
#     mat = np.random.rand(3, 3)
#     mat_R = make_rotation_matrix((1, 0, 0), twopi / 4)  # Rotate around x-axis by 90 degrees

#     mat_rot1 = _rotate_2tensor_b(mat, mat_R)
#     mat_rot2 = _rotate_2tensor(mat, mat_R)

#     np.testing.assert_allclose(mat_rot1, mat_rot2, atol=1e-12)


def test_Bond_matrices_are_different():
    """
    Test that the bond matrices for stiffness and compliance tensors are different.
    """
    mat_R = make_rotation_matrix((1, 2, 3), twopi *.34)  # arbitrary rotation

    mat_MN_stiffness = make_rotation_Bond_matrix_M(mat_R)
    mat_MN_compliance = make_rotation_Bond_matrix_N(mat_R)

    np.testing.assert_raises(AssertionError, np.testing.assert_allclose,
                             mat_MN_stiffness, mat_MN_compliance, atol=1e-12)

def test_plain_tensor2_ij_rotations_recover_identity():

    tens = PlainTensor2_ij('testsym')
    tens.fill_random()

    for nvec in g_rot_vecs:

        for nrots in g_rot_orders:
            mat_R = make_rotation_matrix(nvec, twopi/nrots)
            tens_rot = tens.copy()

            for i in range(nrots):
                tens_rot.rotate(mat_R)
            np.testing.assert_allclose(tens.mat, tens_rot.mat, atol=1e-12)


def test_voigt_tensor3_iJ_rotations_recover_identity():

    tens = VoigtTensor3_iJ('testsym')
    tens.fill_random()

    for nvec in g_rot_vecs:
        for nrots in g_rot_orders:
            mat_R = make_rotation_matrix(nvec, twopi/nrots)
            tens_rot = tens.copy()

            for i in range(nrots):
                tens_rot.rotate(mat_R)
            np.testing.assert_allclose(tens.mat, tens_rot.mat, atol=1e-12)


def test_voigt_tensor4_IJ_rotations_recover_identity():

    for is_stiffness in (True, False):
        tens = VoigtTensor4_IJ('testsym',
                               transforms_like_stiffness=is_stiffness)
        tens.fill_random()

        for nvec in g_rot_vecs:
            for nrots in g_rot_orders:
                mat_R = make_rotation_matrix(nvec, twopi/nrots)
                tens_rot = tens.copy()

                for i in range(nrots):
                    tens_rot.rotate(mat_R)
                np.testing.assert_allclose(tens.mat, tens_rot.mat, atol=1e-12)

def test_voigt_tensor3_iJ_old_rotations_recover_identity():

    tens = VoigtTensor3_iJ('testsym')
    tens.fill_random()

    for nvec in g_rot_vecs:
        for nrots in g_rot_orders:
            mat_R = make_rotation_matrix(nvec, twopi/nrots)
            tens_rot = tens.copy()

            for i in range(nrots):
                tens_rot._old_rotate(mat_R)
            np.testing.assert_allclose(tens.mat, tens_rot.mat, atol=1e-12)

def test_voigt_tensor4_IJ_old_rotations_recover_identity():

    for is_stiffness in (True, False):
        tens = VoigtTensor4_IJ('testsym',
                               transforms_like_stiffness=is_stiffness)
        tens.fill_random()

        for nvec in g_rot_vecs:
            for nrots in g_rot_orders:
                mat_R = make_rotation_matrix(nvec, twopi/nrots)
                tens_rot = tens.copy()

                for i in range(nrots):
                    tens_rot._old_rotate(mat_R)
                np.testing.assert_allclose(tens.mat, tens_rot.mat, atol=1e-12)

def test_voigt_tensor3_iJ_old_and_new_rotations_are_equivalent_dijk():
    tens = VoigtTensor3_iJ('testsym', transforms_like_piezo_d=True)
    tens.fill_random()

    for nvec in g_rot_vecs:
        for nrots in g_rot_orders:
            mat_R = make_rotation_matrix(nvec, twopi/nrots)
            tens_rot = tens.copy()
            tens_rot2 = tens.copy()

            tens_rot.rotate(mat_R)
            tens_rot2._old_rotate(mat_R)
            np.testing.assert_allclose(tens_rot.mat, tens_rot2.mat, atol=1e-12)

def test_voigt_tensor3_iJ_old_and_new_rotations_are_equivalent_eijk():
    tens = VoigtTensor3_iJ('testsym', transforms_like_piezo_d=False)
    tens.fill_random()

    for nvec in g_rot_vecs:
        for nrots in g_rot_orders:
            mat_R = make_rotation_matrix(nvec, twopi/nrots)
            tens_rot = tens.copy()
            tens_rot2 = tens.copy()

            tens_rot.rotate(mat_R)
            tens_rot2._old_rotate(mat_R)
            np.testing.assert_allclose(tens_rot.mat, tens_rot2.mat, atol=1e-12)


def test_voigt_tensor4_IJ_stiffness_old_and_new_rotations_are_equivalent():
    tens = VoigtTensor4_IJ('testsym',
                            transforms_like_stiffness=True)
    tens.fill_random()

    for nvec in g_rot_vecs:
        for nrots in g_rot_orders:
            mat_R = make_rotation_matrix(nvec, twopi/nrots)
            tens_rot = tens.copy()
            tens_rot2 = tens.copy()

            tens_rot.rotate(mat_R)
            tens_rot2._old_rotate(mat_R)
            np.testing.assert_allclose(tens_rot.mat, tens_rot2.mat, atol=1e-12)

def test_voigt_tensor4_IJ_compliance_old_and_new_rotations_are_equivalent():

    tens = VoigtTensor4_IJ('testsym',
                            transforms_like_stiffness=False)
    tens.fill_random()

    for nvec in g_rot_vecs:
        for nrots in g_rot_orders:
            mat_R = make_rotation_matrix(nvec, twopi/nrots)
            tens_rot = tens.copy()
            tens_rot2 = tens.copy()

            tens_rot.rotate(mat_R)
            tens_rot2._old_rotate(mat_R)
            np.testing.assert_allclose(tens_rot.mat, tens_rot2.mat, atol=1e-12)


def test_rotations_preserve_stiffness_and_compliance_inverse_relations():
    """
    Test that the stiffness and compliance tensors preserve their inverse relations
    under rotations.
    """
    tens_stiffness = VoigtTensor4_IJ('testsym', transforms_like_stiffness=True)
    tens_stiffness.fill_random()
    tens_compliance = VoigtTensor4_IJ('testsym', transforms_like_stiffness=False)
    tens_compliance.set_from_matrix(np.linalg.inv(tens_stiffness.value()))

    for nvec in g_rot_vecs:
        for nrots in g_rot_orders:
            mat_R = make_rotation_matrix(nvec, twopi/nrots)

            tens_stiffness.rotate(mat_R)
            tens_compliance.rotate(mat_R)

            np.testing.assert_allclose(
                tens_stiffness.mat, np.linalg.inv(tens_compliance.mat), atol=1e-12)

def test_piezo_diJ_cubic_rotation():
    """
    Test that the piezoelectric tensor diJ transforms correctly under cubic rotations.
    """
    # Test Auld 8.37 (who uses passive rotation)
    # Start
    tens = VoigtTensor3_iJ('testsym', transforms_like_piezo_d=True)
    tens.set_elt_iJ('x4', 1.0)
    tens.set_elt_iJ('y5', 1.0)
    tens.set_elt_iJ('z6', 1.0)

    # Target
    tens2 = VoigtTensor3_iJ('testsym', transforms_like_piezo_d=True)
    tens2.set_elt_iJ('x5', 1.0)
    tens2.set_elt_iJ('y4', -1.0)
    tens2.set_elt_iJ('z1', 0.5)
    tens2.set_elt_iJ('z2', -0.5)


    mat_R = make_rotation_matrix((0,0,1), -twopi/8)  # Rotate around z-axis by -45 degrees for passive rotation
    tens.rotate(mat_R)


    np.testing.assert_allclose(tens.mat, tens2.mat, atol=1e-12)
